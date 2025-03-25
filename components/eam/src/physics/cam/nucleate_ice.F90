module nucleate_ice

!-------------------------------------------------------------------------------
! Purpose:
!  A parameterization of ice nucleation.
!
!  *** This module is intended to be a "portable" code layer.  Ideally it should
!  *** not contain any use association of modules that belong to the model framework.
!
!
! Method:
!  The current method is based on Liu & Penner (2005) & Liu et al. (2007)
!  It related the ice nucleation with the aerosol number, temperature and the
!  updraft velocity. It includes homogeneous freezing of sulfate & immersion
!  freezing on mineral dust (soot disabled) in cirrus clouds, and 
!  Meyers et al. (1992) deposition nucleation in mixed-phase clouds
!
!  The effect of preexisting ice crystals on ice nucleation in cirrus clouds is included, 
!  and also consider the sub-grid variability of temperature in cirrus clouds,
!  following X. Shi et al. ACP (2014).
!
!  Ice nucleation in mixed-phase clouds now uses classical nucleation theory (CNT),
!  follows Y. Wang et al. ACP (2014), Hoose et al. (2010).
!
! Authors:
!  Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
!  Xiangjun Shi & Xiaohong Liu, 01/2014.
!
!  With help from C. C. Chen and B. Eaton (2014)
!-------------------------------------------------------------------------------

use wv_saturation,  only: svp_water, svp_ice
use cam_logfile,    only: iulog
use phys_control,   only: cam_chempkg_is 

!lk+
use K2022, only: get_cirrus
!lk-

implicit none
private
save

integer, parameter :: r8 = selected_real_kind(12)

public :: nucleati_init, nucleati, nucleati_K2022

logical  :: use_preexisting_ice
logical  :: use_hetfrz_classnuc
logical  :: use_nie_nucleate
logical  :: use_dem_nucleate

real(r8) :: pi
real(r8) :: mincld

! Subgrid scale factor on relative humidity (dimensionless)
real(r8) :: subgrid

real(r8), parameter :: Shet   = 1.3_r8     ! het freezing threshold
real(r8), parameter :: rhoice = 0.5e3_r8   ! kg/m3, Wpice is not sensitive to rhoice
real(r8), parameter :: minweff= 0.001_r8   ! m/s
real(r8), parameter :: gamma1=1.0_r8 
real(r8), parameter :: gamma2=1.0_r8 
real(r8), parameter :: gamma3=2.0_r8 
real(r8), parameter :: gamma4=6.0_r8 

real(r8) :: ci

!===============================================================================
contains
!===============================================================================

subroutine nucleati_init( &
   use_preexisting_ice_in, use_hetfrz_classnuc_in, &
   use_nie_nucleate_in, use_dem_nucleate_in, &
   iulog_in, pi_in, mincld_in, subgrid_in)

   logical,  intent(in) :: use_preexisting_ice_in
   logical,  intent(in) :: use_hetfrz_classnuc_in
   logical,  intent(in) :: use_nie_nucleate_in
   logical,  intent(in) :: use_dem_nucleate_in
   integer,  intent(in) :: iulog_in
   real(r8), intent(in) :: pi_in
   real(r8), intent(in) :: mincld_in
   real(r8), intent(in) :: subgrid_in

   use_preexisting_ice = use_preexisting_ice_in
   use_hetfrz_classnuc = use_hetfrz_classnuc_in
   use_nie_nucleate    = use_nie_nucleate_in
   use_dem_nucleate    = use_dem_nucleate_in
   iulog               = iulog_in
   pi                  = pi_in
   mincld              = mincld_in
   subgrid             = subgrid_in

   ci = rhoice*pi/6._r8

end subroutine nucleati_init

!===============================================================================

subroutine nucleati(  &
   wbar, tair, pmid, relhum, cldn,      &
   qc, qi, ni_in, rhoair,               &
   so4_num, dst_num, soot_num,          &
   dst1_sfc_to_num, dst3_sfc_to_num,    &
   nuci, onihf, oniimm, onidep, onimey, &
   wpice, weff, fhom,                   &
   dst1_num,dst2_num,dst3_num,dst4_num, &
   organic_num, clim_modal_aero )

   !---------------------------------------------------------------
   ! Purpose:
   !  The parameterization of ice nucleation.
   !
   ! Method: The current method is based on Liu & Penner (2005)
   !  It related the ice nucleation with the aerosol number, temperature and the
   !  updraft velocity. It includes homogeneous freezing of sulfate, immersion
   !  freezing of soot, and Meyers et al. (1992) deposition nucleation
   !
   ! Authors: Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
   !----------------------------------------------------------------

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: pmid        ! pressure at layer midpoints (pa)
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: qi          ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
   real(r8), intent(in) :: ni_in       ! grid-mean preexisting cloud ice number conc (#/kg) 
   real(r8), intent(in) :: rhoair      ! air density (kg/m3)
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)
   real(r8), intent(in) :: dst1_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst2_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst3_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst4_num     ! dust aerosol number (#/cm^3)
   real(r8), intent(in) :: dst1_sfc_to_num
   real(r8), intent(in) :: dst3_sfc_to_num
   real(r8), intent(in) :: organic_num  !organic aerosol number (primary carbon) (#/cm^3)
   logical,  intent(in) :: clim_modal_aero !whether MAM is used or not

   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)
   real(r8), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
   real(r8), intent(out) :: fhom       ! how much fraction of cloud can reach Shom

   ! Local workspace
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, A, B, regm            ! work variable
   real(r8) :: esl, esi, deles           ! work variable

   !(09/29/2014)DeMott for mixed-phase cloud ice nucleation 
   real(r8)  :: na500, na500_1                            ! aerosol number with D>500 nm (#/cm^3)
   real(r8)  :: na500stp                                  ! aerosol number with D>500 nm (#/cm^3) at STP
   real(r8)  :: nimeystp                                  ! nucleated number from ice nucleation (meyers) at STP
   real(r8)  :: ad, bd   
   real(r8)  :: wbar1, wbar2

   ! Niemand et al. for mixed-phase cloud immersion ice nucleation (surface area based, dust)
   real(r8)  :: an
   real(r8)  :: ns_dust_imm    ! dust active site surface densities from AIDA experiments (m-2)

   ! used in SUBROUTINE Vpreice
   real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
   real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(r8) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
   real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
   real(r8) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet

   real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 
   !-------------------------------------------------------------------------------

   ! temp variables that depend on use_preexisting_ice
   wbar1 = wbar
   wbar2 = wbar

   if (use_preexisting_ice) then

      Ni_preice = ni_in*rhoair                    ! (convert from #/kg -> #/m3)
      Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density 

!!== KZ_BUGFIX 
      if (Ni_preice > 10.0_r8 .and. qi > 1.e-10_r8 ) then    ! > 0.01/L = 10/m3   
!!== KZ_BUGFIX 
         Shom = -1.5_r8   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005
         lami = (gamma4*ci*ni_in/qi)**(1._r8/3._r8)
         Ri_preice = 0.5_r8/lami                  ! radius
         Ri_preice = max(Ri_preice, 1e-8_r8)       ! >0.01micron
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shom, wpice)
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shet, wpicehet)
      else
         wpice    = 0.0_r8
         wpicehet = 0.0_r8
      endif

      weff     = max(wbar-wpice, minweff)
      wpice    = min(wpice, wbar)
      weffhet  = max(wbar-wpicehet,minweff)
      wpicehet = min(wpicehet, wbar)

      wbar1 = weff
      wbar2 = weffhet

      detaT   = wbar/0.23_r8
      RHimean = 1.0_r8
      call frachom(tair, RHimean, detaT, fhom)

   end if

   ni = 0._r8
   tc = tair - 273.15_r8
   ! initialize
   niimm = 0._r8
   nidep = 0._r8
   nihf  = 0._r8
   deles = 0._r8
   esi   = 0._r8

   if(so4_num >= 1.0e-10_r8 .and. (soot_num+dst3_num) >= 1.0e-10_r8 .and. cldn > 0._r8) then

#ifdef USE_XLIU_MOD
!++ Mod from Xiaohong is the following two line conditional.
!   It changes answers so needs climate validation.
      if ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8) then
         if ( ((tc.le.0.0_r8).and.(tc.ge.-37.0_r8).and.(qc.lt.1.e-12_r8)).or.(tc.le.-37.0_r8)) then
#else
      if((tc.le.-35.0_r8) .and. ((relhum*svp_water(tair)/svp_ice(tair)*subgrid).ge.1.2_r8)) then ! use higher RHi threshold
#endif

            A = -1.4938_r8 * log(soot_num+dst3_num) + 12.884_r8
            B = -10.41_r8  * log(soot_num+dst3_num) - 67.69_r8

            regm = A * log(wbar1) + B

            ! heterogeneous nucleation only
            if (tc .gt. regm) then

               if(tc.lt.-40._r8 .and. wbar1.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc,wbar1,relhum,so4_num,nihf)
                  niimm=0._r8
                  nidep=0._r8

                  if (use_preexisting_ice) then
                     if (nihf.gt.1e-3_r8) then ! hom occur,  add preexisting ice
                        niimm=min(dst3_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                        nihf=nihf + Ni_preice*1e-6_r8 - niimm
                     endif
                     nihf=nihf*fhom
                     n1=nihf+niimm 
                  else
                     n1=nihf
                  end if

               else

                  call hetero(tc,wbar2,soot_num+dst3_num,niimm,nidep)
                  if (use_preexisting_ice) then
                     if (niimm .gt. 1e-6_r8) then ! het freezing occur, add preexisting ice
                        niimm = niimm + Ni_preice*1e-6_r8
                        niimm = min(dst3_num, niimm)        ! niimm < dst_num 
                     end if
                  end if
                  nihf=0._r8
                  n1=niimm+nidep

               endif

            ! homogeneous nucleation only
            else if (tc.lt.regm-5._r8) then

               call hf(tc,wbar1,relhum,so4_num,nihf)
               niimm=0._r8
               nidep=0._r8

               if (use_preexisting_ice) then
                  if (nihf.gt.1e-3_r8) then !  hom occur,  add preexisting ice
                     niimm=min(dst3_num,Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                     nihf=nihf + Ni_preice*1e-6_r8 - niimm
                  endif
                  nihf=nihf*fhom
                  n1=nihf+niimm 
               else
                  n1=nihf
               end if

            ! transition between homogeneous and heterogeneous: interpolate in-between
            else

               if (tc.lt.-40._r8 .and. wbar1.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                  call hf(tc, wbar1, relhum, so4_num, nihf)
                  niimm = 0._r8
                  nidep = 0._r8

                  if (use_preexisting_ice) then
                     if (nihf .gt. 1e-3_r8) then ! hom occur,  add preexisting ice
                        niimm = min(dst3_num, Ni_preice*1e-6_r8)       ! assuming dst_num freeze firstly
                        nihf  = nihf + Ni_preice*1e-6_r8 - niimm
                     endif
                     nihf = nihf*fhom
                     n1   = nihf + niimm
                  else
                     n1 = nihf
                  end if

               else

                  call hf(regm-5._r8,wbar1,relhum,so4_num,nihf)
                  call hetero(regm,wbar2,soot_num+dst3_num,niimm,nidep)

                  if (use_preexisting_ice) then
                     nihf = nihf*fhom
                  end if

                  if (nihf .le. (niimm+nidep)) then
                     n1 = nihf
                  else
                     n1=(niimm+nidep)*((niimm+nidep)/nihf)**((tc-regm)/5._r8)
                  endif

                  if (use_preexisting_ice) then
                     if (n1 .gt. 1e-3_r8) then   ! add preexisting ice
                        n1    = n1 + Ni_preice*1e-6_r8
                        niimm = min(dst3_num, n1)  ! assuming all dst_num freezing earlier than hom  !!
                        nihf  = n1 - niimm
                     else
                        n1    = 0._r8
                        niimm = 0._r8
                        nihf  = 0._r8                         
                     endif
                  end if

               end if
            end if

            ni = n1

         end if
      end if
#ifdef USE_XLIU_MOD
   end if
#endif

   if (use_dem_nucleate) then      ! DeMott, use particles number with D>0.5 um
      !++iceMP

      if(clim_modal_aero) then
         if(cam_chempkg_is('trop_mam7') .or. cam_chempkg_is('trop_mam9')) then
           na500_1  = dst1_num*0.566_r8
         elseif(cam_chempkg_is('trop_mam3') .or. cam_chempkg_is('trop_mam4') .or. &
              cam_chempkg_is('trop_mam4_resus') .or. cam_chempkg_is('trop_mam4_resus_soag')  .or. &
              cam_chempkg_is('trop_mam4_mom') .or. cam_chempkg_is('trop_mam4_resus_mom') .or. &
              cam_chempkg_is('linoz_mam3') .or. cam_chempkg_is('linoz_mam4_resus') .or. &
              cam_chempkg_is('linoz_mam4_resus_soag') .or. cam_chempkg_is('linoz_mam4_resus_mom') .or. &
              cam_chempkg_is('linoz_mam4_resus_mom_soag') .or. &
              cam_chempkg_is('chemuci_linozv3_mam5_vbs') .or. &
              cam_chempkg_is('superfast_mam4_resus_mom_soag')) then !ASK Hailong about trop_mam4 
            na500_1 = dst1_num*0.488_r8 + dst3_num
         else
            na500_1 = dst1_num*0.488_r8 + dst2_num + dst3_num + dst4_num   ! scaled for D>0.5-1 um from 0.1-1 um
         endif
      endif

      ! prepare aerosol number and surface data for ice nucleation in mixed-phase clouds        
      na500    = ( soot_num + organic_num ) * 0.0256_r8 + na500_1  ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6

      na500stp = na500 * 101325._r8 * tair/( 273.15_r8 * pmid ) ! at STP      
      ad  = 1.25
      bd  = -0.46*tc-11.6

   !--iceMP   
   endif

   if (use_nie_nucleate) then  ! Niemand et al., surface area based
      an = -0.517*tc + 8.934
   endif

   ! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
   if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
      if (use_dem_nucleate) then          ! use DeMott et al.         
         !++iceMP
         nimeystp = 1.e-3_r8 * 3.0_r8 * (na500stp**ad) * exp(bd)               ! cm^-3 
         nimey=nimeystp*273.15_r8*pmid/(101325_r8*tair)
      else if (use_nie_nucleate) then
         ns_dust_imm = exp(an)          ! m^-2
         nimey = dst1_num * (1._r8 - exp(-1.0_r8 * dst1_sfc_to_num * ns_dust_imm))&   ! cm^-3
               + dst3_num * (1._r8 - exp(-1.0_r8 * dst3_sfc_to_num * ns_dust_imm))

      else
         !--iceMP
         esl = svp_water(tair)     ! over water in mixed clouds
         esi = svp_ice(tair)     ! over ice
         deles = (esl - esi)
         nimey=1.e-3_r8*exp(12.96_r8*deles/esi - 0.639_r8) 
      endif 
   else
      nimey=0._r8
   endif

   if (use_hetfrz_classnuc) nimey = 0._r8

   nuci=ni+nimey
   if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
      write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
      write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
      nuci=0._r8
   endif

   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair

end subroutine nucleati

!===============================================================================
!lk+
subroutine nucleati_K2022(  &
   wbar, tair, pmid, relhum, cldn,      &
   qc, qi, ni_in, rhoair,               &
   so4_num, dst_num, soot_num,          &
   nuci, onihf, oniimm, onidep, onimey, &
   wpice, weff, fhom,                   &
   call_frm_zm_in)

   ! Input Arguments
   real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: tair        ! temperature (K)
   real(r8), intent(in) :: pmid        ! pressure at layer midpoints (pa)
   real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
   real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
   real(r8), intent(in) :: qi          ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
   real(r8), intent(in) :: ni_in       ! grid-mean preexisting cloud ice number conc (#/kg) 
   real(r8), intent(in) :: rhoair      ! air density (kg/m3)
   real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
   real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
   real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)

   ! Output Arguments
   real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
   real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
   real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
   real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
   real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)
   real(r8), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
   real(r8), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
   real(r8), intent(out) :: fhom       ! how much fraction of cloud can reach Shom

   ! Optional Arguments
   logical,  intent(in), optional :: call_frm_zm_in ! true if called from ZM convection scheme

   ! Local workspace
   real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
   real(r8) :: niimm                     ! nucleated number from immersion freezing
   real(r8) :: nidep                     ! nucleated number from deposition nucleation
   real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
   real(r8) :: n1, ni                    ! nucleated number
   real(r8) :: tc, A, B                  ! work variable
   real(r8) :: esl, esi, deles           ! work variable
   real(r8) :: wbar1, wbar2

   ! used in SUBROUTINE Vpreice
   real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
   real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
   real(r8) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
   real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
   real(r8) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet

   real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 

   logical  :: call_frm_zm
!lk+  K2022
!**** Bernd Kärcher, DLR-IPA
!C**** A parameterization of cirrus cloud formation: Revisiting competing ice nucleation
!C**** J. Geophys. Res. (Atmos.) 127, https://doi.org/10.1029/2022JD036907 (2022).
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C  THE USE OF THIS MODULE IMPLIES THAT USERS AGREE NOT TO DISTRIBUTE THE MODULE OR ANY PORTION  C
!C  OF IT AND NOT TO CHANGE THE NAME OF THE MODULE, "get_cirrus".  USERS MAY MODIFY THIS MODULE  C
!C  AS NEEDED OR INCORPORATE IT INTO A (LARGER) MODEL.  ANY SPECIAL REQUESTS WITH REGARD TO THE  C
!C  ORIGINAL MODULE MAY BE DIRECTED TO bernd.kaercher@dlr.de                        Aug 25 2023  C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!c** units cgs-mb-K (output ice concentrations in #/L-air)
!c
!c** set loop=0 for single point evaluation, else for multiple evaluations in a loop
!c** set xnci=0 to prevent pre-existing cirrus ice
!c** set jinp=0 and xnci=0 for pure homogeneous freezing
!c
!c** INP types/ice activity representations:
!c
!c** inp_type = 1-10  (with analytical solutions)
!c**  1:  Heavyside step
!c**  2:  linear ramp
!c
!c** inp_type = 12-20 (requiring numerical integration)
!c**  12: hyperbolic tangent
!c
!      parameter (loop=0)
!c
      INTEGER, PARAMETER :: inpmx=3
      INTEGER, PARAMETER :: iin=7+3*inpmx, ipar=5, iout=6+inpmx
      !parameter (inpmx=3)
      integer inp_type(inpmx)
      real xinp_tot(inpmx),xinp_alpha(inpmx)
!c
!c** i/o-arrays communicating with the parameterization
!c
!      parameter (iin=7+3*inpmx, ipar=5, iout=6+inpmx)
      real zin(iin),zpar(ipar),zout(iout)
      real(r8) :: temp,press,wup,xnci,rci,alphaci,wdown,ssimx,frachomtmp,xnhom,xnhet,xnhom0,xntot
      integer  :: jinp,k,i
!c
!c** print allows screen output from the parameterization for single point evaluation
!c
      logical print!,type
!c
!c** parameters to evaluate INP spectra
!c
!      data sstar          /0.35/               ! step function
!      data sm,sp          /0.22,0.3/           ! linear ramp
!      data smid,ds        /0.4,0.02/           ! hyperbolic tangent
       real(r8):: sstar
       real(r8):: sm,sp 
       real(r8):: smid,ds

       sstar=0.35
       sm=0.22
       sp=0.3
       smid=0.4
       ds=0.02
!lk-
   !-------------------------------------------------------------------------------

   RHimean = relhum*svp_water(tair)/svp_ice(tair)*subgrid

   ! temp variables that depend on use_preexisting_ice
   wbar1 = wbar
   wbar2 = wbar

   ! If not using prexisting ice, the homogeneous freezing happens in the
   ! entire gridbox.
   fhom = 1._r8

   if (present(call_frm_zm_in)) then
     call_frm_zm = call_frm_zm_in
   else
     call_frm_zm = .false.
   end if

   if (use_preexisting_ice .and. (.not. call_frm_zm)) then

      Ni_preice = ni_in*rhoair                    ! (convert from #/kg -> #/m3)
      Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density 

      if (Ni_preice > 10.0_r8 .and. qi > 1.e-10_r8) then    ! > 0.01/L = 10/m3   
         Shom = -1.5_r8   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005
         lami = (gamma4*ci*ni_in/qi)**(1._r8/3._r8)
         Ri_preice = 0.5_r8/lami                  ! radius
         Ri_preice = max(Ri_preice, 1e-8_r8)       ! >0.01micron
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shom, wpice)
         call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shet, wpicehet)
      else
         wpice    = 0.0_r8
         wpicehet = 0.0_r8
      endif

      weff     = max(wbar-wpice, minweff)
      wpice    = min(wpice, wbar)
      weffhet  = max(wbar-wpicehet,minweff)
      wpicehet = min(wpicehet, wbar)

      wbar1 = weff
      wbar2 = weffhet

      detaT   = wbar/0.23_r8
!      if (use_incloud_nuc) then
!        call frachom(tair, 1._r8, detaT, fhom)
!      else
      call frachom(tair, RHimean, detaT, fhom)
!      end if 
!AH
   end if

   ni = 0._r8
   tc = tair - 273.15_r8

   ! initialize
   niimm = 0._r8
   nidep = 0._r8
   nihf  = 0._r8
   deles = 0._r8
   esi   = 0._r8

!lk+   K2022
      temp  = tair !  220.          K                   ! air temperature   
      press = pmid/100.   !250.      pa-> mb                       ! air pressure
      wup   = wbar*100.   !15.      m/s -> cm/s                        ! updraft speed
      wup   = max(0.01,wup)   !if (wup.lt.0.01)  pause 'updraft below 0.1 mm/s
!c
      xnci    = Ni_preice*1E-6  !0.e-3      m-3  ->cm-3                  ! cirrus total ice crystal number concentration
      rci     = Ri_preice*1E2   !30.e-4     m->cm                    ! cirrus number-mean ice crystal radius
      alphaci = 0.1                            ! deposition coefficient for cirrus ice crystals
!c
      jinp = 1!2                                 ! total number of INP types
      inp_type(1)   = 2                        ! INP type                dust
      xinp_tot(1)   = dst_num!*1E-6 !2.e-3       m-3  -> cm-3             ! total INP number concentration
      xinp_alpha(1) = 0.1                      ! deposition coefficient for INP-derived ice crystals
      !xinp_alpha(1) = 0.3                      ! deposition coefficient for INP-derived ice crystals
!      inp_type(2)   = 12
!      xinp_tot(2)   = 8.e-3
!      xinp_alpha(2) = 0.3
!c
!c** set in-array
!c


      zin(1) = temp
      zin(2) = press
      zin(3) = wup
      zin(4) = xnci
      zin(5) = rci
      zin(6) = alphaci
      zin(7) = jinp
      do k = 1, jinp
        i = 8 + 3*(k-1)
        zin(i)   = real(inp_type(k))
        zin(i+1) = xinp_tot(k)
        zin(i+2) = xinp_alpha(k)
      enddo
!c
!c** store INP parameters
!c
      zpar(1) = sp
      zpar(2) = sm
      zpar(3) = sstar
      zpar(4) = ds
      zpar(5) = smid
        print = .false.

   if ((so4_num >= 1.0e-10_r8 .or. (soot_num+dst_num) >= 1.0e-10_r8) .and. cldn > 0._r8) then
      if (RHimean.ge.1.2_r8) then
        call get_cirrus (zin,zpar,zout,iin,ipar,iout,print)


!lk-   K2022

!lk+  K2022
          wdown   = zout(1)
          wdown   = wdown * 0.01      ! cm/s  -> m/s
          ssimx   = zout(2)  +1      ! cesm 1.2     Karcher 0.2
          frachomtmp = zout(3)
          frachomtmp = frachomtmp * 0.01    ! % -> 1
          xnhom   = zout(4)              !cm-3
          xnhom0  = zout(5)               !cm-3
          xnhet   = zout(6)               !cm-3
!          write(66,916) wup,(xntot*1.e3),ssimx,wdown,xnci,
!     >                (xnhet*1.e3),max(1.e-6,(xnhom*1.e3)),
!     >                max(1.e-6,(xnhom0*1.e3)),max(1.e-6,frachom)
!AH
!  if (xnhom.eq.0._r8.and.xnhet.eq.0._r8) then
!          write (iulog,'(A,F10.4,A,F10.4,A,F10.4)') 'ALLENHU double zero',xnhom,' ',xnhet,' ',wup,' ',wdown
!          write (iulog,*) xnhom,xnhet,wup,wdown
!          write (iulog,*) xnhom,xnhet
!  end if

!AH
   xntot  = xnhet + xnhom + xnci
   nuci=xntot
   nimey=0._r8
   nidep=0._r8               ! no deposition
   niimm=xnhet               ! immersion nucleation of dust
!   nidep=xnhet
!   niimm=0._r8
   nihf=xnhom
   fhom=frachomtmp
   wpice=wdown
   weff=wup-wdown
      else
      nuci=0._r8
      onimey=0._r8
      onidep=0._r8
      oniimm=0._r8
      onihf=0._r8
      endif
endif

!lk- K2022   
   nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
   onimey = nimey*1.e+6_r8/rhoair
   onidep = nidep*1.e+6_r8/rhoair
   oniimm = niimm*1.e+6_r8/rhoair
   onihf  = nihf*1.e+6_r8/rhoair
end subroutine nucleati_K2022

!lk-
!===============================================================================

subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) B,C

!---------------------------------------------------------------------
! parameters

      A11 = 0.0263_r8
      A12 = -0.0185_r8
      A21 = 2.758_r8
      A22 = 1.3221_r8
      B11 = -0.008_r8
      B12 = -0.0468_r8
      B21 = -0.2667_r8
      B22 = -1.4588_r8

!     ice from immersion nucleation (cm^-3)

      B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
      C =  A21+B21*log(Ns)

      Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
      Nis = min(Nis,Ns)

      Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

end subroutine hetero

!===============================================================================

subroutine hf(T,ww,RH,Na,Ni)

      real(r8), intent(in)  :: T, ww, RH, Na
      real(r8), intent(out) :: Ni

      real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
      real(r8)    A2_fast,B2_fast
      real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
      real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
      real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
      real(r8)    regm
      real(r8)    A,B,C
      real(r8)    RHw

!---------------------------------------------------------------------
! parameters

      A1_fast  =0.0231_r8
      A21_fast =-1.6387_r8  !(T>-64 deg)
      A22_fast =-6.045_r8   !(T<=-64 deg)
      B1_fast  =-0.008_r8
      B21_fast =-0.042_r8   !(T>-64 deg)
      B22_fast =-0.112_r8   !(T<=-64 deg)
      C1_fast  =0.0739_r8
      C2_fast  =1.2372_r8

      A1_slow  =-0.3949_r8
      A2_slow  =1.282_r8
      B1_slow  =-0.0156_r8
      B2_slow  =0.0111_r8
      B3_slow  =0.0217_r8
      C1_slow  =0.120_r8
      C2_slow  =2.312_r8

      Ni = 0.0_r8

!----------------------------
!RHw parameters
      A = 6.0e-4_r8*log(ww)+6.6e-3_r8
      B = 6.0e-2_r8*log(ww)+1.052_r8
      C = 1.68_r8  *log(ww)+129.35_r8
      RHw=(A*T*T+B*T+C)*0.01_r8

      if((T.le.-37.0_r8) .and. ((RH*subgrid).ge.RHw)) then

        regm = 6.07_r8*log(ww)-55.0_r8

        if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
            A2_fast=A21_fast
            B2_fast=B21_fast
          else
            A2_fast=A22_fast
            B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

        else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

        endif

      end if

end subroutine hf

!===============================================================================

SUBROUTINE Vpreice(P_in, T_in, R_in, C_in, S_in, V_out)

   !  based on  Karcher et al. (2006)
   !  VERTICAL VELOCITY CALCULATED FROM DEPOSITIONAL LOSS TERM

   ! SUBROUTINE arguments
   REAL(r8), INTENT(in)  :: P_in       ! [Pa],INITIAL AIR pressure 
   REAL(r8), INTENT(in)  :: T_in       ! [K] ,INITIAL AIR temperature 
   REAL(r8), INTENT(in)  :: R_in       ! [m],INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS 
   REAL(r8), INTENT(in)  :: C_in       ! [m-3],INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY, [1/cm3]
   REAL(r8), INTENT(in)  :: S_in       ! [-],INITIAL ICE SATURATION RATIO;; if <1, use hom threshold Si 
   REAL(r8), INTENT(out) :: V_out      ! [m/s], VERTICAL VELOCITY REDUCTION (caused by preexisting ice)

   ! SUBROUTINE parameters
   REAL(r8), PARAMETER :: ALPHAc  = 0.5_r8 ! density of ice (g/cm3), !!!V is not related to ALPHAc 
   REAL(r8), PARAMETER :: FA1c    = 0.601272523_r8        
   REAL(r8), PARAMETER :: FA2c    = 0.000342181855_r8
   REAL(r8), PARAMETER :: FA3c    = 1.49236645E-12_r8        
   REAL(r8), PARAMETER :: WVP1c   = 3.6E+10_r8   
   REAL(r8), PARAMETER :: WVP2c   = 6145.0_r8
   REAL(r8), PARAMETER :: FVTHc   = 11713803.0_r8
   REAL(r8), PARAMETER :: THOUBKc = 7.24637701E+18_r8
   REAL(r8), PARAMETER :: SVOLc   = 3.23E-23_r8    ! SVOL=XMW/RHOICE
   REAL(r8), PARAMETER :: FDc     = 249.239822_r8
   REAL(r8), PARAMETER :: FPIVOLc = 3.89051704E+23_r8         
   REAL(r8) :: T,P,S,R,C
   REAL(r8) :: A1,A2,A3,B1,B2
   REAL(r8) :: T_1,PICE,FLUX,ALP4,CISAT,DLOSS,VICE

   T = T_in          ! K  , K
   P = P_in*1e-2_r8  ! Pa , hpa

   IF (S_in.LT.1.0_r8) THEN
      S = 2.349_r8 - (T/259.0_r8) ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   ELSE
      S = S_in                    ! INPUT ICE SATURATION RATIO, -,  >1
   ENDIF

   R     = R_in*1e2_r8   ! m  => cm
   C     = C_in*1e-6_r8  ! m-3 => cm-3
   T_1   = 1.0_r8/ T
   PICE  = WVP1c * EXP(-(WVP2c*T_1))
   ALP4  = 0.25_r8 * ALPHAc      
   FLUX  = ALP4 * SQRT(FVTHc*T)
   CISAT = THOUBKc * PICE * T_1   
   A1    = ( FA1c * T_1 - FA2c ) * T_1 
   A2    = 1.0_r8/ CISAT      
   A3    = FA3c * T_1 / P
   B1    = FLUX * SVOLc * CISAT * ( S-1.0_r8 ) 
   B2    = FLUX * FDc * P * T_1**1.94_r8 
   DLOSS = FPIVOLc * C * B1 * R**2 / ( 1.0_r8+ B2 * R )         
   VICE  = ( A2 + A3 * S ) * DLOSS / ( A1 * S )  ! 2006,(19)
   V_out = VICE*1e-2_r8  ! cm/s => m/s

END SUBROUTINE Vpreice

subroutine frachom(Tmean,RHimean,detaT,fhom)
   ! How much fraction of cirrus might reach Shom  
   ! base on "A cirrus cloud scheme for general circulation models",
   ! B. Karcher and U. Burkhardt 2008

   real(r8), intent(in)  :: Tmean, RHimean, detaT
   real(r8), intent(out) :: fhom

   real(r8), parameter :: seta = 6132.9_r8  ! K
   integer,  parameter :: Nbin=200          ! (Tmean - 3*detaT, Tmean + 3*detaT)

   real(r8) :: PDF_T(Nbin)    ! temperature PDF;  ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
   real(r8) :: Sbin(Nbin)     ! the fluctuations of Si that are driven by the T variations 
   real(r8) :: Sihom, deta
   integer  :: i

   Sihom = 2.349_r8-Tmean/259.0_r8   ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
   fhom  = 0.0_r8

   do i = Nbin, 1, -1

      deta     = (i - 0.5_r8 - Nbin/2)*6.0_r8/Nbin   ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
      Sbin(i)  = RHimean*exp(deta*detaT*seta/Tmean**2.0_r8)
      PDF_T(i) = exp(-deta**2.0_r8/2.0_r8)*6.0_r8/(sqrt(2.0_r8*Pi)*Nbin)
      

      if (Sbin(i).ge.Sihom) then
         fhom = fhom + PDF_T(i)
      else
         exit
      end if
   end do

   fhom=fhom/0.997_r8   ! accounting for the finite limits (-3 , 3)

end subroutine frachom

end module nucleate_ice

