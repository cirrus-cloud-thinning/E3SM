#include "physics/mam/eamxx_lp05_icenuc_process_interface.hpp"

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/physics/physics_constants.hpp"
#include "share/util/eamxx_column_ops.hpp"

#include <mam4xx/conversions.hpp>
#include <mam4xx/mam4_math.hpp>
#include <mam4xx/wv_sat_methods.hpp>

#include <ekat_assert.hpp>
#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_units.hpp>

#include <cmath>

namespace scream
{

namespace {

// mam4::NucleateIce::nucleati (Liu & Penner 2005) with optional fixes for
// two cases where it returns no ice although dust would freeze
// (docs/mam4xx_lp05_issues.md in the cirrus-thinning project repo):
//  - the branch is homogeneous-only (T < regm-5, or T < -40 C with w > 1 m/s)
//    but hf gives nihf = 0 (RH below the homogeneous threshold, or T > -37 C):
//    use immersion freezing on dust at T instead;
//  - transition band (regm-5 <= T <= regm) with nihf <= niimm: mam4xx sets
//    n1 = nihf; use n1 = niimm instead, so n1 is continuous at T = regm and
//    never zero when dust froze.
// With het_fallback = false this is mam4xx's nucleati. Outputs in #/kg.
KOKKOS_INLINE_FUNCTION
void lp05_nucleati (const Real wbar, const Real tair, const Real relhum,
                    const Real cldn, const Real rhoair, const Real so4_num,
                    const Real dst3_num, const Real subgrid, const bool het_fallback,
                    Real& nuci, Real& onihf, Real& oniimm)
{
  namespace ni = mam4::nucleate_ice;
  Real nihf = 0, niimm = 0, nidep = 0, n1 = 0;
  const Real tc = tair - Real(273.15);
  constexpr Real num_threshold = 1.0e-10;
  if (so4_num >= num_threshold && dst3_num >= num_threshold && cldn > 0 &&
      tc <= Real(-35.0) &&
      relhum*mam4::wv_sat_methods::svp_water(tair)/mam4::wv_sat_methods::svp_ice(tair)
        *subgrid >= Real(1.2)) {
    Real regm = 0;
    ni::calculate_regm_nucleati(wbar, dst3_num, regm);
    const bool fast_hom = tc < Real(-40) && wbar > Real(1);
    // homogeneous freezing at T, falling back to dust if it gives nothing
    auto hom_at_tc = [&]() {
      ni::hf(tc, wbar, relhum, so4_num, subgrid, nihf);
      niimm = 0;
      n1 = nihf;
      if (het_fallback && nihf <= 0) {
        ni::hetero(tc, wbar, dst3_num, niimm, nidep);
        n1 = niimm + nidep;
      }
    };
    if (tc > regm) {
      if (fast_hom) {
        hom_at_tc();
      } else {
        ni::hetero(tc, wbar, dst3_num, niimm, nidep);
        n1 = niimm + nidep;
      }
    } else if (tc < regm - Real(5) || fast_hom) {
      hom_at_tc();
    } else {
      ni::hf(regm - Real(5), wbar, relhum, so4_num, subgrid, nihf);
      ni::hetero(regm, wbar, dst3_num, niimm, nidep);
      if (nihf <= niimm + nidep) {
        n1 = het_fallback ? niimm + nidep : nihf;
      } else {
        n1 = (niimm + nidep)*mam4::pow((niimm + nidep)/nihf, (tc - regm)/Real(5));
      }
    }
  }
  if (n1 > Real(9999.) || n1 < 0) n1 = 0;
  constexpr Real cm3_to_kg = 1.e6;  // (#/cm3)/(kg/m3)*(1e6 cm3/m3) -> #/kg
  nuci   = n1*cm3_to_kg/rhoair;
  onihf  = nihf*cm3_to_kg/rhoair;
  oniimm = niimm*cm3_to_kg/rhoair;
}

} // anonymous namespace

LP05IceNuc::LP05IceNuc (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // file: prescribed aerosol read from prescribed_aero_files
  // mam4_aci: so4_num_icenuc/dst_num_icenuc computed by mam4_aci in the same
  //           step (for checking this process against mam4_aci)
  m_aero_source = m_params.get<std::string>("aerosol_source", "file");
  EKAT_REQUIRE_MSG(m_aero_source=="file" || m_aero_source=="mam4_aci",
      "ERROR: lp05_icenuc aerosol_source must be 'file' or 'mam4_aci'.");
  EKAT_REQUIRE_MSG(m_aero_source!="file" || m_params.isParameter("prescribed_aero_files"),
      "ERROR: prescribed_aero_files is missing from lp05_icenuc parameter list.");

  // Name of the output ice nucleation number. P3 reads ni_activated; another
  // name lets this process run as a diagnostic next to mam4_aci.
  m_ni_name = m_params.get<std::string>("ni_activated_name", "ni_activated");

  m_subgrid     = m_params.get<double>("nucleate_ice_subgrid", 1.2);
  m_wsubice_min = m_params.get<double>("wsubice_min", 0.2);
  m_wsubice_max = m_params.get<double>("wsubice_max", 10.0);
  m_dust_scale  = m_params.get<double>("dust_scale_factor", 1.0);
  m_top_lev     = m_params.get<int>("top_level", 0);

  // Nucleation updraft:
  //   tke:          SHOC w variance, clipped to [wsubice_min, wsubice_max]
  //                 (same as mam4_aci; EAM icenul_wsub_scheme=1 with CLUBB)
  //   eam_gaussian: Gaussian w PDF centred on the resolved w, LP05 called once
  //                 with the PDF's positive partial mean, int_{w>0} w f(w) dw
  //                 (EAM icenul_wsub_scheme=2, subgrid_mean_updraft)
  //   pdf_integral: same PDF, LP05 called at nbins updrafts over w>0 and the
  //                 results weighted by probability, int_{w>0} n_i(w) f(w) dw
  // For both PDF methods sigma = sqrt(wsub_tke^2 + sigma_w_extra^2), with
  // sigma >= 0.001 m/s as in EAM, and wsubice_min is not applied.
  const auto method = m_params.get<std::string>("updraft_method", "tke");
  EKAT_REQUIRE_MSG(method=="tke" || method=="eam_gaussian" || method=="pdf_integral",
      "ERROR: lp05_icenuc updraft_method must be tke, eam_gaussian or pdf_integral.");
  m_updraft_method = method=="tke" ? 0 : (method=="eam_gaussian" ? 1 : 2);
  m_sigma_w_extra  = m_params.get<double>("sigma_w_extra", 0.0);
  m_pdf_num_bins   = m_params.get<int>("pdf_num_bins", 40);
  m_pdf_num_sigma  = m_params.get<double>("pdf_num_sigma", 4.0);

  // Use dust freezing where LP05 picks homogeneous freezing but it gives no
  // (or fewer) crystals; see lp05_nucleati above. Off reproduces mam4xx.
  m_het_fallback = m_params.get<bool>("het_fallback", false);
  EKAT_REQUIRE_MSG(m_sigma_w_extra>=0,
      "ERROR: lp05_icenuc sigma_w_extra must be non-negative.");
  EKAT_REQUIRE_MSG(m_pdf_num_bins>=1 && m_pdf_num_sigma>0,
      "ERROR: lp05_icenuc requires pdf_num_bins >= 1 and pdf_num_sigma > 0.");

  m_so4_file_name = m_params.get<std::string>("so4_num_var_name", "so4_num_icenuc");
  m_dst_file_name = m_params.get<std::string>("dst_num_var_name", "dst_num_icenuc");

  m_seed_num_per_liter = m_params.get<double>("seed_inp_per_liter", 0.0);
  m_seed_lat_min = m_params.get<double>("seed_lat_min", -90.0);
  m_seed_lat_max = m_params.get<double>("seed_lat_max",  90.0);
  m_seed_lon_min = m_params.get<double>("seed_lon_min",   0.0);
  m_seed_lon_max = m_params.get<double>("seed_lon_max", 360.0);
  m_seed_p_min   = m_params.get<double>("seed_p_min",     0.0);
  m_seed_p_max   = m_params.get<double>("seed_p_max",  1.1e5);

  EKAT_REQUIRE_MSG(m_subgrid>0,
      "ERROR: lp05_icenuc nucleate_ice_subgrid must be positive.");
  EKAT_REQUIRE_MSG(m_wsubice_min>0 && m_wsubice_max>=m_wsubice_min,
      "ERROR: lp05_icenuc requires 0 < wsubice_min <= wsubice_max.");
  EKAT_REQUIRE_MSG(m_dust_scale>=0,
      "ERROR: lp05_icenuc dust_scale_factor must be non-negative.");
  EKAT_REQUIRE_MSG(m_seed_num_per_liter>=0,
      "ERROR: lp05_icenuc seed_inp_per_liter must be non-negative.");
  EKAT_REQUIRE_MSG(m_seed_lat_min<=m_seed_lat_max && m_seed_p_min<=m_seed_p_max,
      "ERROR: lp05_icenuc seeding box has min > max in latitude or pressure.");

  // Model longitudes are in [0,360). Accept [-180,180) input too.
  auto wrap360 = [](const Real lon) {
    Real l = std::fmod(lon, Real(360));
    return l < 0 ? l + 360 : l;
  };
  if (not (m_seed_lon_min==0 && m_seed_lon_max==360)) {
    m_seed_lon_min = wrap360(m_seed_lon_min);
    m_seed_lon_max = wrap360(m_seed_lon_max);
  }
}

// =========================================================================================
void LP05IceNuc::create_requests ()
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  m_grid = m_grids_manager->get_grid("physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  const auto scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);
  const auto nondim = Units::nondimensional();
  constexpr auto m2 = pow(m, 2);
  constexpr auto s2 = pow(s, 2);

  // Model state
  add_field<Required>("T_mid",       scalar3d_mid, K,     grid_name);
  add_field<Required>("p_mid",       scalar3d_mid, Pa,    grid_name);
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);
  add_field<Required>("w_variance",  scalar3d_mid, m2/s2, grid_name);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);
  if (m_updraft_method!=0) {
    add_field<Required>("omega", scalar3d_mid, Pa/s, grid_name);
  }
  add_tracer<Required>("qv", m_grid, kg/kg);

  // Aerosol number used by LP05 [#/kg dry air], from file or from mam4_aci
  add_field<Computed>("lp05_so4_num", scalar3d_mid, 1/kg, grid_name);
  add_field<Computed>("lp05_dst_num", scalar3d_mid, 1/kg, grid_name);
  if (m_aero_source=="mam4_aci") {
    add_field<Required>("so4_num_icenuc", scalar3d_mid, 1/kg, grid_name);
    add_field<Required>("dst_num_icenuc", scalar3d_mid, 1/kg, grid_name);
  }

  // Output used by P3
  add_field<Computed>(m_ni_name, scalar3d_mid, 1/kg, grid_name);

  // Diagnostics
  add_field<Computed>("lp05_nihf",     scalar3d_mid, 1/kg,   grid_name);
  add_field<Computed>("lp05_niimm",    scalar3d_mid, 1/kg,   grid_name);
  add_field<Computed>("lp05_wsubice",  scalar3d_mid, m/s,    grid_name);
  add_field<Computed>("lp05_seed_num", scalar3d_mid, 1/kg,   grid_name);
  // 0: no nucleation, 1: homogeneous only, 2: heterogeneous only, 3: both
  add_field<Computed>("lp05_regime",   scalar3d_mid, nondim, grid_name);
  // T and RH wrt water as LP05 saw them (output T/qv are written after P3)
  add_field<Computed>("lp05_T",        scalar3d_mid, K,      grid_name);
  add_field<Computed>("lp05_relhum",   scalar3d_mid, nondim, grid_name);
  if (m_updraft_method!=0) {
    // resolved w = -omega/(rho g), PDF width, and probability of w > 0
    add_field<Computed>("lp05_w_resolved", scalar3d_mid, m/s,    grid_name);
    add_field<Computed>("lp05_sigma_w",    scalar3d_mid, m/s,    grid_name);
    add_field<Computed>("lp05_frac_rising", scalar3d_mid, nondim, grid_name);
  }
}

// =========================================================================================
void LP05IceNuc::initialize_impl (const RunType /* run_type */)
{
  // The LP05 kernel only needs the default configuration; the subgrid factor
  // is passed explicitly to nucleati at run time.
  mam4::AeroConfig aero_config;
  m_nucleate_ice.init(aero_config, mam4::NucleateIce::Config(m_subgrid));

  // Work arrays for the nucleation updraft
  m_dz        = view_2d("lp05_dz", m_num_cols, m_num_levs);
  m_w_sec_int = view_2d("lp05_w_sec_int", m_num_cols, m_num_levs+1);

  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<LowerBound>(get_field_out(m_ni_name), m_grid, 0);
  if (m_aero_source=="mam4_aci") {
    return;
  }

  // Prescribed aerosol input
  std::vector<Field> aero_fields = {
    get_field_out("lp05_so4_num").alias(m_so4_file_name),
    get_field_out("lp05_dst_num").alias(m_dst_file_name)
  };
  const auto files    = m_params.get<std::vector<std::string>>("prescribed_aero_files");
  EKAT_REQUIRE_MSG(files.size()>0 && files[0]!="" && files[0]!="none",
      "ERROR: lp05_icenuc needs at least one file in prescribed_aero_files.");
  const auto map_file = m_params.get<std::string>("prescribed_aero_remap_file", "");
  const auto time_interp = m_params.get<std::string>("time_interpolation_method", "linear");
  const auto pname    = m_params.get<std::string>("prescribed_aero_ps_var_name", "PS");

  m_data_interpolation = std::make_shared<DataInterpolation>(m_grid, aero_fields);
  m_data_interpolation->setup_time_database(files, time_interp);
  m_data_interpolation->create_horiz_remappers(map_file=="none" ? "" : map_file,
                                               m_iop_data_manager);
  DataInterpolation::VertRemapData vremap_data;
  vremap_data.vr_type = DataInterpolation::Dynamic3DRef;
  vremap_data.pname = pname;
  vremap_data.pmid = get_field_in("p_mid");
  m_data_interpolation->create_vert_remapper(vremap_data);
  m_data_interpolation->init_time_interpolation(start_of_step_ts(), DataInterpolation::Linear);

  add_postcondition_check<LowerBound>(get_field_out("lp05_so4_num"), m_grid, 0, true);
  add_postcondition_check<LowerBound>(get_field_out("lp05_dst_num"), m_grid, 0, true);
}

// =========================================================================================
void LP05IceNuc::run_impl (const double /* dt */)
{
  using PF = PhysicsFunctions<DefaultDevice>;
  using CO = ColumnOps<DefaultDevice, Real>;
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;

  if (m_aero_source=="file") {
    // Update prescribed aerosol to the end of this step
    m_data_interpolation->run(end_of_step_ts());
  } else {
    get_field_out("lp05_so4_num").deep_copy(get_field_in("so4_num_icenuc"));
    get_field_out("lp05_dst_num").deep_copy(get_field_in("dst_num_icenuc"));
  }

  const auto T_mid   = get_field_in("T_mid").get_view<const Real**>();
  const auto p_mid   = get_field_in("p_mid").get_view<const Real**>();
  const auto p_del   = get_field_in("pseudo_density").get_view<const Real**>();
  const auto cldfrac = get_field_in("cldfrac_tot").get_view<const Real**>();
  const auto w_var   = get_field_in("w_variance").get_view<const Real**>();
  const auto qv      = get_field_in("qv").get_view<const Real**>();
  const auto so4_num = get_field_out("lp05_so4_num").get_view<const Real**>();
  const auto dst_num = get_field_out("lp05_dst_num").get_view<const Real**>();

  const auto ni_activated = get_field_out(m_ni_name).get_view<Real**>();
  const auto nihf         = get_field_out("lp05_nihf").get_view<Real**>();
  const auto niimm        = get_field_out("lp05_niimm").get_view<Real**>();
  const auto wsubice_out  = get_field_out("lp05_wsubice").get_view<Real**>();
  const auto seed_out     = get_field_out("lp05_seed_num").get_view<Real**>();
  const auto regime       = get_field_out("lp05_regime").get_view<Real**>();
  const auto T_out        = get_field_out("lp05_T").get_view<Real**>();
  const auto relhum_out   = get_field_out("lp05_relhum").get_view<Real**>();

  const_view_2d omega;
  view_2d w_res_out, sigma_out, frac_out;
  if (m_updraft_method!=0) {
    omega     = get_field_in("omega").get_view<const Real**>();
    w_res_out = get_field_out("lp05_w_resolved").get_view<Real**>();
    sigma_out = get_field_out("lp05_sigma_w").get_view<Real**>();
    frac_out  = get_field_out("lp05_frac_rising").get_view<Real**>();
  }

  const auto lat = m_grid->get_geometry_data("lat").get_view<const Real*>();
  const auto lon = m_grid->get_geometry_data("lon").get_view<const Real*>();

  // Local copies for device capture
  const auto nucleate_ice = m_nucleate_ice;
  const auto dz           = m_dz;
  const auto w_sec_int    = m_w_sec_int;
  const Real subgrid      = m_subgrid;
  const Real wmin         = m_wsubice_min;
  const Real wmax         = m_wsubice_max;
  const Real dust_scale   = m_dust_scale;
  const int  method       = m_updraft_method;
  const Real sigma_extra  = m_sigma_w_extra;
  const int  nbins        = m_pdf_num_bins;
  const Real nsigma       = m_pdf_num_sigma;
  const bool het_fallback = m_het_fallback;
  constexpr Real gravit   = physics::Constants<Real>::gravit.value;
  const int  top_lev      = m_top_lev;
  const Real seed_cm3     = m_seed_num_per_liter*1e-3; // #/L -> #/cm3
  const Real lat_min = m_seed_lat_min, lat_max = m_seed_lat_max;
  const Real lon_min = m_seed_lon_min, lon_max = m_seed_lon_max;
  const Real p_min   = m_seed_p_min,   p_max   = m_seed_p_max;

  // Same constants as mam4::NucleateIce::compute_tendencies
  constexpr Real tmelt_m_five = 273.15 - 5.0;
  constexpr Real mincld       = 1e-4;
  constexpr Real kg_to_cm3    = 1e-6;  // (#/kg)*(kg/m3)*(m3/cm3) -> #/cm3

  const int nlev = m_num_levs;
  const auto policy = TPF::get_default_team_policy(m_num_cols, nlev);

  Kokkos::parallel_for("lp05_icenuc::run", policy,
      KOKKOS_LAMBDA(const KT::MemberType& team) {
    const int icol = team.league_rank();

    // Nucleation updraft, computed exactly as in mam4_aci
    // (mam_coupling::compute_updraft_velocities, compute_tke_at_interfaces,
    // compute_subgrid_scale_velocities): w variance is interpolated to
    // interfaces using dz, converted to TKE, and averaged back to midpoints.
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](const int k) {
      dz(icol,k) = PF::calculate_dz(p_del(icol,k), p_mid(icol,k),
                                    T_mid(icol,k), qv(icol,k));
    });
    team.team_barrier();
    const auto w_var_col = ekat::subview(w_var, icol);
    CO::compute_interface_values_linear(team, nlev, w_var_col,
                                        ekat::subview(dz, icol),
                                        w_var_col(0), w_var_col(nlev-1),
                                        ekat::subview(w_sec_int, icol));
    team.team_barrier();

    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](const int k) {
      const Real temp = T_mid(icol,k);
      const Real pmid = p_mid(icol,k);
      const Real rho  = mam4::conversions::density_of_ideal_gas(temp, pmid);

      // Subgrid (turbulent) updraft from SHOC, as in mam4_aci
      Real wsub = 0.001;  // above top_lev, as in mam4_aci
      if (k >= top_lev) {
        const Real tke_k   = (3.0/2.0)*w_sec_int(icol,k);
        const Real tke_kp1 = (3.0/2.0)*w_sec_int(icol,k+1);
        wsub = mam4::sqrt(0.5*(tke_k + tke_kp1)*(2.0/3.0));
      }

      // Gaussian w PDF for the PDF methods (EAM: w0 and wsig are 0 and 0.001
      // above top_lev)
      Real w0 = 0, sigma = 0.001;
      if (method!=0 && k >= top_lev) {
        w0 = -omega(icol,k)/(rho*gravit);
        const Real wt = mam4::min(wsub, Real(10));
        sigma = mam4::max(Real(0.001), mam4::sqrt(wt*wt + sigma_extra*sigma_extra));
      }

      Real wsubice = 0.001;
      if (method==0) {
        if (k >= top_lev) wsubice = mam4::min(wmax, mam4::max(wmin, wsub));
      } else if (method==1) {
        // EAM subgrid_mean_updraft: 50 bins over w0 +/- 3 sigma, midpoint
        // rule, positive bins only
        constexpr int nb = 50;
        const Real dw = 6*sigma/nb;
        const Real norm = dw/(sigma*mam4::sqrt(2*Real(M_PI)));
        Real ww = 0;
        for (int b=0; b<nb; ++b) {
          const Real w = w0 - 3*sigma + (b + Real(0.5))*dw;
          if (w > 0) ww += w*mam4::exp(-(w-w0)*(w-w0)/(2*sigma*sigma))*norm;
        }
        wsubice = ww > 0 ? ww : Real(0.001);
      }

      // Seeding box. Seeds act as extra heterogeneous INP; LP05 only uses
      // heterogeneous INP in the cirrus regime (T <= -35 C).
      const bool in_lon = (lon_min <= lon_max)
                          ? (lon(icol) >= lon_min && lon(icol) <= lon_max)
                          : (lon(icol) >= lon_min || lon(icol) <= lon_max); // wraps 0
      const bool in_box = lat(icol) >= lat_min && lat(icol) <= lat_max && in_lon
                          && pmid >= p_min && pmid <= p_max;
      const Real seed = in_box ? seed_cm3 : Real(0);

      Real nuci = 0, onihf = 0, oniimm = 0;
      Real frac_rising = 0;
      Real relhum_diag = 0;
      if (temp < tmelt_m_five) {
        // mam4_aci passes LP05 the dry-air water vapor mixing ratio
        const Real qv_dry = PF::calculate_drymmr_from_wetmmr(qv(icol,k), qv(icol,k));
        Real es = 0, qs = 0;
        mam4::wv_sat_methods::wv_sat_qsat_water(temp, pmid, es, qs);
        const Real relhum = qv_dry / qs;
        relhum_diag = relhum;
        const Real icldm  = mam4::max(cldfrac(icol,k), mincld);

        // Prescribed numbers are dry mixing ratios [#/kg], as in mam4_aci
        const Real so4_cm3 = mam4::max(Real(0), so4_num(icol,k))*rho*kg_to_cm3;
        const Real dst_cm3 = dust_scale*mam4::max(Real(0), dst_num(icol,k))*rho*kg_to_cm3
                             + seed;

        // LP05 at one updraft: mam4xx's kernel, or the local copy with fixes
        auto lp05 = [&](const Real w, Real& n, Real& nhf, Real& nimm) {
          if (het_fallback) {
            lp05_nucleati(w, temp, relhum, icldm, rho, so4_cm3, dst_cm3, subgrid,
                          true, n, nhf, nimm);
          } else {
            Real ndep = 0, nmey = 0;
            nucleate_ice.nucleati(w, temp, pmid, relhum, icldm, rho,
                                  so4_cm3, dst_cm3, subgrid,
                                  n, nhf, nimm, ndep, nmey);
          }
        };
        if (method!=2) {
          lp05(wsubice, nuci, onihf, oniimm);
        } else {
          // Probability-weighted LP05 over the rising part of the PDF:
          // nbins equal bins on [max(0, w0 - nsigma*sigma), w0 + nsigma*sigma],
          // each weighted by its exact Gaussian probability and evaluated at
          // its midpoint (capped at wsubice_max). Sinking air nucleates nothing.
          const Real lo = mam4::max(Real(0), w0 - nsigma*sigma);
          const Real hi = w0 + nsigma*sigma;
          const Real rs2 = 1/(mam4::sqrt(Real(2))*sigma);
          Real wsum = 0;
          if (hi > lo) {
            const Real dw = (hi - lo)/nbins;
            for (int b=0; b<nbins; ++b) {
              const Real wa = lo + b*dw;
              const Real p  = Real(0.5)*(Kokkos::erf((wa + dw - w0)*rs2)
                                         - Kokkos::erf((wa - w0)*rs2));
              if (p <= 0) continue;
              const Real w = mam4::min(wmax, wa + Real(0.5)*dw);
              Real n = 0, nhf = 0, nimm = 0;
              lp05(w, n, nhf, nimm);
              nuci   += p*n;
              onihf  += p*nhf;
              oniimm += p*nimm;
              frac_rising += p;
              wsum   += p*w;
            }
          }
          // Mean updraft of the rising air, for diagnostics only
          wsubice = frac_rising > 0 ? wsum/frac_rising : Real(0);
        }
      }
      if (method!=0) {
        if (method==1 || temp >= tmelt_m_five) {
          frac_rising = Real(0.5)*(1 + Kokkos::erf(w0/(mam4::sqrt(Real(2))*sigma)));
        }
        w_res_out(icol,k) = w0;
        sigma_out(icol,k) = sigma;
        frac_out(icol,k)  = frac_rising;
      }

      T_out(icol,k)        = temp;
      relhum_out(icol,k)   = relhum_diag;
      ni_activated(icol,k) = nuci;    // [#/kg]
      nihf(icol,k)         = onihf;   // [#/kg]
      niimm(icol,k)        = oniimm;  // [#/kg]
      wsubice_out(icol,k)  = wsubice;
      seed_out(icol,k)     = seed/(rho*kg_to_cm3);  // back to [#/kg]
      regime(icol,k)       = (onihf > 0 ? 1 : 0) + (oniimm > 0 ? 2 : 0);
    });
  });
}

} // namespace scream
