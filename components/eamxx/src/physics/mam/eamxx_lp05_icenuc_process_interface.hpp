#ifndef EAMXX_LP05_ICENUC_PROCESS_INTERFACE_HPP
#define EAMXX_LP05_ICENUC_PROCESS_INTERFACE_HPP

#include "share/atm_process/atmosphere_process.hpp"

#include <mam4xx/nucleate_ice.hpp>

#include <ekat_parameter_list.hpp>

#include <string>

namespace scream
{

class DataInterpolation;

/*
 * Liu & Penner (2005) cirrus ice nucleation driven by prescribed aerosol.
 *
 * This process calls the mam4xx port of LP05 (mam4::NucleateIce::nucleati)
 * without running the rest of MAM4xx. Sulfate and dust number concentrations
 * are read from file (e.g. written by a nudged MAM4xx run) with time and
 * space interpolation, so that seeded and control simulations see identical
 * background aerosol. Its main output is ni_activated, which P3 uses for
 * cirrus temperatures when use_aerosol_cirrus_ice_nucleation is on.
 *
 * Optional cirrus seeding adds a fixed number of seed ice nucleating
 * particles to the dust (heterogeneous INP) number inside a lat/lon/pressure
 * box.
 */
class LP05IceNuc : public AtmosphereProcess
{
public:
  using KT            = ekat::KokkosTypes<DefaultDevice>;
  using view_1d       = typename KT::template view_1d<Real>;
  using view_2d       = typename KT::template view_2d<Real>;
  using const_view_1d = typename KT::template view_1d<const Real>;
  using const_view_2d = typename KT::template view_2d<const Real>;

  LP05IceNuc (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::string name () const { return "lp05_icenuc"; }

  void create_requests ();

protected:

  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   () { /* Nothing to do */ }

  std::shared_ptr<const AbstractGrid> m_grid;
  int m_num_cols;
  int m_num_levs;

  std::shared_ptr<DataInterpolation> m_data_interpolation;

  mam4::NucleateIce m_nucleate_ice;

  // LP05 settings
  Real m_subgrid;       // subgrid-scale factor on RH (EAM: 1.2 for P3+SHOC)
  Real m_wsubice_min;   // floor on nucleation updraft [m/s]
  Real m_wsubice_max;   // cap on nucleation updraft [m/s]
  Real m_dust_scale;    // multiplier on prescribed dust number
  int  m_top_lev;       // levels above this use wsubice=0.001 (as mam4_aci)

  // Work arrays for the nucleation updraft
  view_2d m_dz;         // layer thickness [m]
  view_2d m_w_sec_int;  // w variance at interfaces [m2/s2]

  // Seeding settings
  Real m_seed_num_per_liter;  // seed INP number concentration [#/L]
  Real m_seed_lat_min, m_seed_lat_max;  // [deg]
  Real m_seed_lon_min, m_seed_lon_max;  // [deg], 0-360
  Real m_seed_p_min,   m_seed_p_max;    // [Pa]

  std::string m_aero_source;  // "file" or "mam4_aci"
  std::string m_ni_name;      // output field name (default ni_activated)

  // Names of the prescribed fields in the input file(s)
  std::string m_so4_file_name;
  std::string m_dst_file_name;
};

} // namespace scream

#endif // EAMXX_LP05_ICENUC_PROCESS_INTERFACE_HPP
