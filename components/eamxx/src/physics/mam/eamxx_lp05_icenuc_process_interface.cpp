#include "physics/mam/eamxx_lp05_icenuc_process_interface.hpp"

#include "share/algorithm/eamxx_data_interpolation.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
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
      Real wsubice = 0.001;  // above top_lev, as in mam4_aci
      if (k >= top_lev) {
        const Real tke_k   = (3.0/2.0)*w_sec_int(icol,k);
        const Real tke_kp1 = (3.0/2.0)*w_sec_int(icol,k+1);
        const Real wsub = mam4::sqrt(0.5*(tke_k + tke_kp1)*(2.0/3.0));
        wsubice = mam4::min(wmax, mam4::max(wmin, wsub));
      }

      const Real temp = T_mid(icol,k);
      const Real pmid = p_mid(icol,k);
      const Real rho  = mam4::conversions::density_of_ideal_gas(temp, pmid);

      // Seeding box. Seeds act as extra heterogeneous INP; LP05 only uses
      // heterogeneous INP in the cirrus regime (T <= -35 C).
      const bool in_lon = (lon_min <= lon_max)
                          ? (lon(icol) >= lon_min && lon(icol) <= lon_max)
                          : (lon(icol) >= lon_min || lon(icol) <= lon_max); // wraps 0
      const bool in_box = lat(icol) >= lat_min && lat(icol) <= lat_max && in_lon
                          && pmid >= p_min && pmid <= p_max;
      const Real seed = in_box ? seed_cm3 : Real(0);

      Real nuci = 0, onihf = 0, oniimm = 0, onidep = 0, onimey = 0;
      if (temp < tmelt_m_five) {
        // mam4_aci passes LP05 the dry-air water vapor mixing ratio
        const Real qv_dry = PF::calculate_drymmr_from_wetmmr(qv(icol,k), qv(icol,k));
        Real es = 0, qs = 0;
        mam4::wv_sat_methods::wv_sat_qsat_water(temp, pmid, es, qs);
        const Real relhum = qv_dry / qs;
        const Real icldm  = mam4::max(cldfrac(icol,k), mincld);

        // Prescribed numbers are dry mixing ratios [#/kg], as in mam4_aci
        const Real so4_cm3 = mam4::max(Real(0), so4_num(icol,k))*rho*kg_to_cm3;
        const Real dst_cm3 = dust_scale*mam4::max(Real(0), dst_num(icol,k))*rho*kg_to_cm3
                             + seed;

        nucleate_ice.nucleati(wsubice, temp, pmid, relhum, icldm, rho,
                              so4_cm3, dst_cm3, subgrid,
                              nuci, onihf, oniimm, onidep, onimey);
      }

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
