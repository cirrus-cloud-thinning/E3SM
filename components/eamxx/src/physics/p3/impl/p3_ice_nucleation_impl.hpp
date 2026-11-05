#ifndef P3_ICE_NUCLEATION_IMPL_HPP
#define P3_ICE_NUCLEATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_nucleation(
  const Pack& temp, const Pack& inv_rho, const Pack& ni, const Pack& ni_activated,
  const Pack& qv_supersat_i, const Scalar& inv_dt, const bool& do_predict_nc, const bool& do_prescribed_CCN,
  Pack& qv2qi_nucleat_tend, Pack& ni_nucleat_tend,
  const P3Runtime& runtime_options, 
  const Mask& context)
{
   constexpr Scalar nsmall   = C::NSMALL;
   constexpr Scalar tmelt    = C::Tmelt.value;
   constexpr Scalar T_icenuc = C::Tmelt.value-sp(15.0);
   constexpr Scalar zero     = C::ZERO;
   constexpr Scalar piov3    = C::PIOV3;
   constexpr Scalar mi0      = sp(4.0)*piov3*sp(900.0)*sp(1.e-18);

   const Scalar deposition_nucleation_exponent =
       runtime_options.deposition_nucleation_exponent;

   // Below this temperature, an aerosol-aware ice nucleation process (e.g. LP05)
   // can supply ni_activated for cirrus even when CCN are prescribed.
   constexpr Scalar T_cirrus = C::Tmelt.value-sp(35.0);

   const auto t_lt_T_icenuc = temp < T_icenuc;
   const auto qv_supersat_i_ge_005 = qv_supersat_i >= 0.05;

   // do_log: use the temperature-only deposition nucleation formula
   // (Cooper-type) instead of the aerosol-derived ni_activated.
   const auto do_log = (!do_predict_nc || do_prescribed_CCN);

   // With use_aerosol_cirrus_ice_nucleation, cirrus temperatures use ni_activated
   // even when do_log is true; warmer temperatures keep the do_log behavior.
   const Mask use_aero_here = runtime_options.use_aerosol_cirrus_ice_nucleation
                              ? Mask(temp < T_cirrus) : Mask(false);

   const auto base_mask      = t_lt_T_icenuc && qv_supersat_i_ge_005 && context;
   const auto any_if_log     = base_mask && (Mask(!do_log) || use_aero_here);
   const auto any_if_not_log = base_mask && Mask(do_log) && !use_aero_here;

   Pack dum{0.0}, N_nuc{0.0}, Q_nuc{0.0};

   // NOTE: the two branches below act on disjoint lanes. With the runtime
   // option off, at most one of them is active, which reproduces the original
   // if/else structure bit-for-bit.
   if (any_if_not_log.any()) {
     dum = sp(0.005)*exp(sp(deposition_nucleation_exponent)*(tmelt-temp))*sp(1.0e3)*inv_rho;

     dum = min(dum, sp(1.0e5)*inv_rho);

     N_nuc = max(zero, (dum-ni)*inv_dt);

     const auto n_nuc_ge_nsmall = N_nuc >= nsmall && context;

     if (n_nuc_ge_nsmall.any()) {
       Q_nuc = max(zero, (dum-ni)*mi0*inv_dt);

       qv2qi_nucleat_tend.set(any_if_not_log && n_nuc_ge_nsmall, Q_nuc);

       ni_nucleat_tend.set(any_if_not_log && n_nuc_ge_nsmall, N_nuc);
     }
   }
   if (any_if_log.any()) {
     ni_nucleat_tend.set(any_if_log, max(zero, (ni_activated-ni)*inv_dt));
     qv2qi_nucleat_tend.set(any_if_log, ni_nucleat_tend*mi0);
   }
}

} // namespace p3
} // namespace scream

#endif // P3_ICE_NUCLEATION_IMPL_HPP
