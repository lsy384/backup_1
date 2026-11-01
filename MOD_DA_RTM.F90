#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_RTM
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Top-level forward observation operators for satellite brightness
!    temperature. Three schemes are provided:
!       (1) forward_CMEM  : CMEM single-layer Fresnel soil + canopy RTM
!       (2) forward_CoLM  : CoLM-Liu discrete multi-layer incoherent soil RTM
!       (3) forward_LandEM: LandEM bulk two-stream land emissivity model
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Shuyue Liu, 10/2026: Refine codes and modularize physics
!
! REFERENCES:
!   [1] Kerr, Y. H., et al. (2017). Modelling the passive microwave
!       signature from land surfaces: A review of recent results and
!       application to the L-band SMOS & SMAP soil moisture retrieval
!       algorithms. Remote Sensing of Environment, 192, 238-262.
!   [2] Weng, F., Yan, B., & Grody, N. C. (2001). A microwave land
!       emissivity model. Journal of Geophysical Research: Atmospheres,
!       106(D17), 20115-20123.
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_1DForcing
   USE MOD_DA_Const
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: nl_soil, nl_lake, N_land_classification
   USE MOD_Namelist
   USE MOD_DA_RTM_Physics
   IMPLICIT NONE
   SAVE

! Public forward operator interfaces
   PUBLIC :: forward_CMEM
   PUBLIC :: forward_CoLM
   PUBLIC :: forward_LandEM

CONTAINS

!-----------------------------------------------------------------------
   SUBROUTINE forward_CMEM( &
      patchtype, patchclass, dz_sno, &
      forc_topo, htop, &
      tref, t_soisno, tleaf, &
      wliq_soisno, wice_soisno, h2osoi, &
      snowdp, lai, sai, &
      wf_clay, wf_sand, wf_silt, BD_all, porsl, &
      sat_theta, sat_fghz, &
      tb_toa_h, tb_toa_v, &
      tb_soil_h, tb_soil_v, tb_tov_h, tb_tov_v)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Forward modeling of brightness temperature observations using the
!   CMEM scheme (single-layer Fresnel soil emission + canopy RTM).
!
! REFERENCES:
!   [1] Kerr, Y. H., et al. (2017). Modelling the passive microwave
!       signature from land surfaces: A review of recent results and
!       application to the L-band SMOS & SMAP soil moisture retrieval
!       algorithms. Remote Sensing of Environment, 192, 238-262.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      USE MOD_Vars_Global, only: nl_soil, nl_lake, maxsnl, spval, dz_soi
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                        ! land cover type
      integer, intent(in)   :: patchclass                       ! land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! topography [m]
      real(r8), intent(in)  :: htop                             ! upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! ice lens in layers [kg/m2]
      real(r8), intent(in)  :: h2osoi(nl_soil)                  ! volumetric soil water in layers [m3/m3]
      real(r8), intent(in)  :: snowdp                           ! snow depth [meter]
      real(r8), intent(in)  :: lai                              ! leaf area index
      real(r8), intent(in)  :: sai                              ! stem area index
      real(r8), intent(in)  :: wf_clay(nl_soil)                 ! gravimetric fraction of clay
      real(r8), intent(in)  :: wf_sand(nl_soil)                 ! gravimetric fraction of sand
      real(r8), intent(in)  :: wf_silt(nl_soil)                 ! gravimetric fraction of silt
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)                   ! fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out) :: tb_toa_v                         ! brightness temperature of top-of-atmosphere for V- polarization
      real(r8), intent(out) :: tb_soil_h                        ! brightness temperature of soil for H- polarization
      real(r8), intent(out) :: tb_soil_v                        ! brightness temperature of soil for V- polarization
      real(r8), intent(out) :: tb_tov_h                         ! brightness temperature of vegetation (consider snow) for H- polarization
      real(r8), intent(out) :: tb_tov_v                         ! brightness temperature of vegetation (consider snow) for V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: liq_soi(nl_soil)            ! liquid volumetric water content in layers [m3/m3]
      real(r8) :: ice_soi(nl_soil)            ! ice volumetric water content in layers [m3/m3]
      logical  :: is_low_veg                  ! flag for low vegetation (NOTE: currently unused, kept for legacy consistency)
      real(r8) :: dz_soisno(maxsnl+1:nl_soil) ! layer thickness [m]
      integer  :: lb                          ! lower bound of arrays
      real(r8) :: tau_atm                     ! atmospheric optical depth
      real(r8) :: r_r(2)                      ! rough surface reflectivity for H and V polarizations
      real(r8) :: r_sn(2)                     ! reflectivity between the snow and ground for H and V polarizations
      real(r8) :: r_snow(2)                   ! reflectivity of the snow for H and V polarizations
      real(r8) :: tb_soil(2)                  ! brightness temperature of soil for H and V polarizations
      real(r8) :: tb_tos(2)                   ! brightness temperature of snow-covered ground for H and V polarizations
      real(r8) :: tb_tov(2)                   ! brightness temperature of vegetation (consider snow) for H and V polarizations
      real(r8) :: tb_tov_noad(2)              ! brightness temperature of vegetation (no downwelling radiation) for H and V polarizations
      real(r8) :: tb_au(2)                    ! upwelling radiation (brightness temperature) of atmosphere
      real(r8) :: tb_ad(2)                    ! downwelling radiation (brightness temperature) of atmosphere
      real(r8) :: rho_snow                    ! snow density (kg/m3)
      real(r8) :: liq_snow                    ! snow liquid water content (m3/m3)
      real(r8) :: gamma_p(2)                  ! vegetation opacity for H- and V- polarization
      real(r8) :: tb_veg(2)                   ! brightness temperature of vegetation for H- and V- polarization
      real(r8) :: tb_2(2)                     ! the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_3(2)                     ! upwelling soil emission attenuated by the canopy
      real(r8) :: tb_4(2)                     ! the downwelling cosmic ray reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_toa(2)                   ! brightness temperature of top-of-atmosphere for H- and V- polarization
      real(r8) :: wf_total(nl_soil)           ! total gravimetric
      real(r8) :: BD_all_surf                 ! bulk density of soil (kg/m3) at surface
      real(r8) :: porsl_surf                  ! soil porosity at surface
      real(r8) :: t_surf                      ! soil temperature at surface (K)
      real(r8) :: t_deep                      ! soil temperature at deep layer (K)
      real(r8) :: liq_surf                    ! liquid volumetric water content at surface (m3/m3)
      real(r8) :: ice_surf                    ! ice volumetric water content at surface (m3/m3)
      real(r8) :: wf_clay_surf                ! gravimetric clay fraction at surface
      real(r8) :: wf_sand_surf                ! gravimetric sand fraction at surface
      integer  :: i

!#############################################################################
! Prepare parameters & states used in the operator
!#############################################################################
      ! get depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! calculate weighted parameters
      wf_total     = max(wf_clay + wf_sand + wf_silt, 1.0e-12_r8)
      ! surface is the (first two layers)
      wf_clay_surf = (wf_clay(1)/wf_total(1)*0.0175 + wf_clay(2)/wf_total(2)*0.0276)/0.0451
      wf_sand_surf = (wf_sand(1)/wf_total(1)*0.0175 + wf_sand(2)/wf_total(2)*0.0276)/0.0451
      BD_all_surf  = (BD_all(1)*0.0175 + BD_all(2)*0.0276)/0.0451
      porsl_surf   = (porsl(1)*0.0175 + porsl(2)*0.0276)/0.0451

      ! calculate temperature (K) at surface and deep soil layers
      t_surf = ((t_soisno(1)*(0.0175) + t_soisno(2)*(0.0451 - 0.0175))/0.0451)
      t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5)

      ! calculate liquid/ice volumetric water (first two layers)
      liq_surf = (wliq_soisno(1) + wliq_soisno(2))/(0.0451*denh2o)
      ice_surf = (wice_soisno(1) + wice_soisno(2))/(0.0451*denice)

      ! ! calculate weighted parameters
      ! wf_total     = max(wf_clay + wf_sand + wf_silt, 1.0e-12_r8)
      ! ! surface is the (first two layers)
      ! wf_clay_surf = wf_clay(1) / wf_total(1)
      ! wf_sand_surf = wf_sand(1)/wf_total(1)
      ! BD_all_surf  = BD_all(1)
      ! porsl_surf   = porsl(1)

      ! ! calculate temperature (K) at surface and deep soil layers
      ! t_surf = t_soisno(1)
      ! t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5)

      ! ! calculate liquid/ice volumetric water (first two layers)
      ! liq_surf     = wliq_soisno(1) / (dz_soisno(1)*denh2o)
      ! ice_surf     = wice_soisno(1) / (dz_soisno(1)*denice)

      ! calculate liquid/ice volumetric water profile for all soil layers
      DO i = 1, nl_soil
         liq_soi(i) = wliq_soisno(i) / (dz_soisno(i)*denh2o)
         ice_soi(i) = wice_soisno(i) / (dz_soisno(i)*denice)
      END DO

      ! calculate lower bound of snow
      lb = 0
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

!#############################################################################
! Run the forward operator
!#############################################################################
      ! check the patch type
      IF (patchtype >= 3) THEN ! ocean, lake, ice
         tb_toa = spval
      ELSE
         ! calculate parameters used in operator varied with satellite
         CALL calc_parameters(sat_theta, sat_fghz)

!#############################################################################
! atmosphere module
!#############################################################################
         CALL atm(forc_topo, tref, tau_atm, tb_au, tb_ad)

!#############################################################################
! soil module
!#############################################################################
         CALL soil_CMEM(&
            patchclass, nl_soil, dz_soi(1:nl_soil), t_soisno(1:nl_soil), &
            liq_soi(1:nl_soil), ice_soi(1:nl_soil), wf_sand(1:nl_soil), wf_clay(1:nl_soil), BD_all(1:nl_soil), porsl(1:nl_soil), &
            t_surf, t_deep, &
            liq_surf, ice_surf, &
            wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
            r_r, tb_soil)

!#############################################################################
! vegetation and snow module
!    We categorized four different cases for the calculations:
!    1) no vegetation and no snow
!    2) no vegetation with snow
!    3) vegetation without snow
!    4) vegetation with snow
!#############################################################################
         ! roughly judge low or high vegetation (only for IGBP)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         END IF

         ! ensure snow density to <= 1000 kg/m3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/dz_soisno(lb)
            liq_snow = wliq_soisno(lb)/(dz_soisno(lb)*denh2o)

            IF (liq_snow > 1.0 .or. rho_snow > 1000.0_r8) then
               rho_snow = 1000.0_r8
               liq_snow = wliq_soisno(lb)/(dz_soisno(lb)*denh2o)
            END IF
         END IF

         ! main procedures
         ! --------------------------------------------------------------------
         ! 1) no veg and no snow
         !    two components:
         !       (1) brightness temperature of soil
         !       (2) the downwelling radiation reflected by the soil
         ! --------------------------------------------------------------------
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         ! --------------------------------------------------------------------
         ! 2) no veg and has snow
         !    two components:
         !       (1) brightness temperature of snow
         !       (2) the downwelling radiation reflected by the snow
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN
            ! calculate brightness temperature of snow-covered ground
            call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         ! --------------------------------------------------------------------
         ! 3) has veg and no snow
         !    four components:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN
            ! calculate brightness temperature of vegetation
            IF (DEF_DA_RTM_veg == 0) THEN
               CALL veg_Wig07(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            ELSE IF (DEF_DA_RTM_veg == 1) THEN
               CALL veg_Jac91(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            ELSE IF (DEF_DA_RTM_veg == 2) THEN
               CALL veg_Kir79(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            END IF

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            END DO

         ! --------------------------------------------------------------------
         ! 4) has veg and has snow
         !    We need to determine the positional relationship between vegetation and snow.
         !
         !    If vegetation is higher than snow,
         !    we first calculate brightness temperature of snow (soil boundary), then calculate
         !    four components to derive brightness temperature of top of vegetation:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the snow and attenuated by the canopy layer
         !       (3) upwelling snow emission attenuated by the canopy
         !       (4) the downwelling reflected by the snow and attenuated by the canopy layer
         !
         !    If vegetation is lower than snow
         !    we first calculate brightness temperature of top of vegetation (soil boundary), then calculate
         !    four components to derive brightness temperature of top of snow:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               ! calculate brightness temperature of vegetation
               IF (DEF_DA_RTM_veg == 0) THEN
                  CALL veg_Wig07(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 1) THEN
                  CALL veg_Jac91(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 2) THEN
                  CALL veg_Kir79(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               END IF

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO

               ! calculate reflectivity between the snow and low veg (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! calculate brightness temperature of snow-covered ground
               call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

               ! calculate brightness temperature of top of snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! calculate brightness temperature of snow-covered ground
               call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

               ! calculate brightness temperature of vegetation
               IF (DEF_DA_RTM_veg == 0) THEN
                  CALL veg_Wig07(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 1) THEN
                  CALL veg_Jac91(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 2) THEN
                  CALL veg_Kir79(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               END IF

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO
            END IF
         END IF

!#############################################################################
! Calculate brightness temperature of top-of-atmosphere
!#############################################################################
         tb_toa = tb_tov*exp(-tau_atm) + tb_au
      END IF

      tb_toa_h  = tb_toa(1);  tb_toa_v  = tb_toa(2)
      tb_soil_h = tb_soil(1); tb_soil_v = tb_soil(2)
      tb_tov_h  = tb_tov(1);  tb_tov_v  = tb_tov(2)

   END SUBROUTINE forward_CMEM


!-----------------------------------------------------------------------
   SUBROUTINE forward_CoLM( &
      patchtype, patchclass, dz_sno, &
      forc_topo, htop, &
      tref, t_soisno, tleaf, &
      wliq_soisno, wice_soisno, h2osoi, &
      snowdp, lai, sai, &
      wf_clay, wf_sand, wf_silt, BD_all, porsl, &
      sat_theta, sat_fghz, &
      tb_toa_h, tb_toa_v, &
      tb_soil_h, tb_soil_v, tb_tov_h, tb_tov_v)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Forward modeling of brightness temperature observations using the
!   CoLM-Liu discrete multi-layer incoherent radiative transfer scheme.
!   The soil emission is solved by a discrete multi-layer analytical
!   incoherent RTM (soil_CoLM) instead of the single-layer Fresnel
!   scheme used in forward_CMEM. Snow, vegetation, and atmospheric
!   modules are shared with forward_CMEM.
!
! REFERENCES: Liu et al. (in preparation)
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      USE MOD_Vars_Global, only: nl_soil, nl_lake, maxsnl, spval, dz_soi
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                        ! land cover type
      integer, intent(in)   :: patchclass                       ! land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! topography [m]
      real(r8), intent(in)  :: htop                             ! upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! ice lens in layers [kg/m2]
      real(r8), intent(in)  :: h2osoi(nl_soil)                  ! volumetric soil water in layers [m3/m3]
      real(r8), intent(in)  :: snowdp                           ! snow depth [meter]
      real(r8), intent(in)  :: lai                              ! leaf area index
      real(r8), intent(in)  :: sai                              ! stem area index
      real(r8), intent(in)  :: wf_clay(nl_soil)                 ! gravimetric fraction of clay
      real(r8), intent(in)  :: wf_sand(nl_soil)                 ! gravimetric fraction of sand
      real(r8), intent(in)  :: wf_silt(nl_soil)                 ! gravimetric fraction of silt
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)                   ! fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out) :: tb_toa_v                         ! brightness temperature of top-of-atmosphere for V- polarization
      real(r8), intent(out) :: tb_soil_h                        ! brightness temperature of soil for H- polarization
      real(r8), intent(out) :: tb_soil_v                        ! brightness temperature of soil for V- polarization
      real(r8), intent(out) :: tb_tov_h                         ! brightness temperature of vegetation (consider snow) for H- polarization
      real(r8), intent(out) :: tb_tov_v                         ! brightness temperature of vegetation (consider snow) for V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: liq_soi(nl_soil), ice_soi(nl_soil)   ! liquid/ice volumetric water content in soil layers [m3/m3]
      logical  :: is_low_veg                           ! flag for low vegetation (NOTE: currently unused, kept for legacy consistency)
      real(r8) :: dz_soisno(maxsnl+1:nl_soil)          ! layer thickness [m]
      integer  :: lb                                   ! lower bound of arrays
      real(r8) :: tau_atm, r_r(2), r_sn(2), r_snow(2)  ! atmospheric optical depth and reflectivities for H/V polarizations
      real(r8) :: tb_soil(2), tb_tos(2), tb_tov(2), tb_tov_noad(2) ! brightness temperatures at soil, snow, and canopy levels
      real(r8) :: tb_au(2), tb_ad(2)                   ! up-/down-welling brightness temperature of atmosphere
      real(r8) :: rho_snow, liq_snow, gamma_p(2), tb_veg(2) ! snow and vegetation parameters
      real(r8) :: tb_2(2), tb_3(2), tb_4(2), tb_toa(2) ! multi-component brightness temperatures
      real(r8) :: wf_total(nl_soil), BD_all_surf, porsl_surf ! surface soil parameters
      real(r8) :: t_surf, t_deep, liq_surf, ice_surf        ! surface temperature and wetness
      real(r8) :: wf_clay_surf, wf_sand_surf                ! surface soil texture
      integer  :: i

!#############################################################################
! Prepare parameters & states used in the operator
!#############################################################################
      ! get depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! calculate weighted parameters
      wf_total     = max(wf_clay + wf_sand + wf_silt, 1.0e-12_r8)
      ! surface is the (first two layers)
      wf_clay_surf = wf_clay(1) / wf_total(1)
      wf_sand_surf = wf_sand(1)/wf_total(1)
      BD_all_surf  = BD_all(1)
      porsl_surf   = porsl(1)

      ! calculate temperature (K) at surface and deep soil layers
      t_surf = t_soisno(1)
      t_deep = ((t_soisno(7)*(0.8289-0.5) + t_soisno(8)*(1.0 - 0.8289))/0.5)

      ! calculate liquid/ice volumetric water (first two layers)
      liq_surf     = wliq_soisno(1) / (dz_soisno(1)*denh2o)
      ice_surf     = wice_soisno(1) / (dz_soisno(1)*denice)
      
      ! calculate liquid/ice volumetric water profile for all soil layers
      DO i = 1, nl_soil
         liq_soi(i) = wliq_soisno(i) / (dz_soisno(i)*denh2o)
         ice_soi(i) = wice_soisno(i) / (dz_soisno(i)*denice)
      END DO

      ! calculate lower bound of snow
      lb = 0
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

!#############################################################################
! Run the forward operator
!#############################################################################
      ! check the patch type
      IF (patchtype >= 3) THEN ! ocean, lake, ice
         tb_toa = spval
      ELSE
         ! calculate parameters used in operator varied with satellite
         CALL calc_parameters(sat_theta, sat_fghz)

!#############################################################################
! atmosphere module
!#############################################################################
         CALL atm(forc_topo, tref, tau_atm, tb_au, tb_ad)

!#############################################################################
! soil module
!    NOTE: soil_CoLM replaces the single-layer Fresnel scheme used in
!    forward_CMEM with a discrete multi-layer analytical incoherent
!    radiative transfer (see MOD_DA_RTM_Physics for details).
!#############################################################################
         CALL soil_CoLM(&
            patchclass, nl_soil, dz_soi(1:nl_soil), t_soisno(1:nl_soil), &
            liq_soi(1:nl_soil), ice_soi(1:nl_soil), wf_sand(1:nl_soil), wf_clay(1:nl_soil), BD_all(1:nl_soil), porsl(1:nl_soil), &
            t_surf, liq_surf, ice_surf, &
            wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
            r_r, tb_soil)

!#############################################################################
! vegetation and snow module
!    We categorized four different cases for the calculations:
!    1) no vegetation and no snow
!    2) no vegetation with snow
!    3) vegetation without snow
!    4) vegetation with snow
!#############################################################################
         ! roughly judge low or high vegetation (only for IGBP)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         END IF

         ! ensure snow density to <= 1000 kg/m3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/dz_soisno(lb)
            liq_snow = wliq_soisno(lb)/(dz_soisno(lb)*denh2o)

            IF (liq_snow > 1.0 .or. rho_snow > 1000.0_r8) then
               rho_snow = 1000.0_r8
               liq_snow = wliq_soisno(lb)/(dz_soisno(lb)*denh2o)
            END IF
         END IF

         ! main procedures
         ! --------------------------------------------------------------------
         ! 1) no veg and no snow
         !    two components:
         !       (1) brightness temperature of soil
         !       (2) the downwelling radiation reflected by the soil
         ! --------------------------------------------------------------------
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         ! --------------------------------------------------------------------
         ! 2) no veg and has snow
         !    two components:
         !       (1) brightness temperature of snow
         !       (2) the downwelling radiation reflected by the snow
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN
            ! calculate brightness temperature of snow-covered ground
            call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         ! --------------------------------------------------------------------
         ! 3) has veg and no snow
         !    four components:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN
            ! calculate brightness temperature of vegetation
            IF (DEF_DA_RTM_veg == 0) THEN
               CALL veg_Wig07(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            ELSE IF (DEF_DA_RTM_veg == 1) THEN
               CALL veg_Jac91(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            ELSE IF (DEF_DA_RTM_veg == 2) THEN
               CALL veg_Kir79(patchclass, lai, htop, 0.0_r8, tleaf, tb_veg, gamma_p)
            END IF

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            END DO

         ! --------------------------------------------------------------------
         ! 4) has veg and has snow
         !    We need to determine the positional relationship between vegetation and snow.
         !
         !    If vegetation is higher than snow,
         !    we first calculate brightness temperature of snow (soil boundary), then calculate
         !    four components to derive brightness temperature of top of vegetation:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the snow and attenuated by the canopy layer
         !       (3) upwelling snow emission attenuated by the canopy
         !       (4) the downwelling reflected by the snow and attenuated by the canopy layer
         !
         !    If vegetation is lower than snow
         !    we first calculate brightness temperature of top of vegetation (soil boundary), then calculate
         !    four components to derive brightness temperature of top of snow:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               ! calculate brightness temperature of vegetation
               IF (DEF_DA_RTM_veg == 0) THEN
                  CALL veg_Wig07(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 1) THEN
                  CALL veg_Jac91(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 2) THEN
                  CALL veg_Kir79(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               END IF

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO

               ! calculate reflectivity between the snow and low veg (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! calculate brightness temperature of snow-covered ground
               call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

               ! calculate brightness temperature of top of snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! calculate brightness temperature of snow-covered ground
               call snow_Pul97(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

               ! calculate brightness temperature of vegetation
               IF (DEF_DA_RTM_veg == 0) THEN
                  CALL veg_Wig07(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 1) THEN
                  CALL veg_Jac91(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               ELSE IF (DEF_DA_RTM_veg == 2) THEN
                  CALL veg_Kir79(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
               END IF

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO
            END IF
         END IF

!#############################################################################
! Calculate brightness temperature of top-of-atmosphere
!#############################################################################
         tb_toa = tb_tov*exp(-tau_atm) + tb_au
      END IF

      tb_toa_h  = tb_toa(1);  tb_toa_v  = tb_toa(2)
      tb_soil_h = tb_soil(1); tb_soil_v = tb_soil(2)
      tb_tov_h  = tb_tov(1);  tb_tov_v  = tb_tov(2)

   END SUBROUTINE forward_CoLM


!-----------------------------------------------------------------------
   SUBROUTINE forward_LandEM( &
      patchtype, patchclass, dz_sno, &
      forc_topo, htop, &
      tref, t_soisno, tleaf, &
      wliq_soisno, wice_soisno, h2osoi, &
      snowdp, lai, sai, &
      wf_clay, wf_sand, wf_silt, BD_all, porsl, &
      sat_theta, sat_fghz, &
      tb_toa_h, tb_toa_v, &
      tb_soil_h, tb_soil_v, tb_tov_h, tb_tov_v)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Forward operator based on the LandEM bulk two-stream microwave land
!   emissivity model (Weng et al. 2001), adapted for CoLM continuous
!   states. LandEM is a BULK model and does NOT output separate soil
!   and vegetation brightness temperatures: tb_soil_* are set to spval.
!
! REFERENCES:
!   [1] Weng, F., Yan, B., & Grody, N. C. (2001). A microwave land
!       emissivity model. Journal of Geophysical Research: Atmospheres,
!       106(D17), 20115-20123.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      USE MOD_Vars_Global, only: nl_soil, maxsnl, spval
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                        ! land cover type
      integer, intent(in)   :: patchclass                       ! land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! topography [m]
      real(r8), intent(in)  :: htop                             ! upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! ice lens in layers [kg/m2]
      real(r8), intent(in)  :: h2osoi(nl_soil)                  ! volumetric soil water in layers [m3/m3]
      real(r8), intent(in)  :: snowdp                           ! snow depth [meter]
      real(r8), intent(in)  :: lai                              ! leaf area index
      real(r8), intent(in)  :: sai                              ! stem area index
      real(r8), intent(in)  :: wf_clay(nl_soil)                 ! gravimetric fraction of clay
      real(r8), intent(in)  :: wf_sand(nl_soil)                 ! gravimetric fraction of sand
      real(r8), intent(in)  :: wf_silt(nl_soil)                 ! gravimetric fraction of silt
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)                   ! fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out) :: tb_toa_v                         ! brightness temperature of top-of-atmosphere for V- polarization
      real(r8), intent(out) :: tb_soil_h                        ! brightness temperature of soil for H- polarization
      real(r8), intent(out) :: tb_soil_v                        ! brightness temperature of soil for V- polarization
      real(r8), intent(out) :: tb_tov_h                         ! brightness temperature of vegetation (consider snow) for H- polarization
      real(r8), intent(out) :: tb_tov_v                         ! brightness temperature of vegetation (consider snow) for V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: tau_atm, tb_au(2), tb_ad(2)              ! atmospheric optical depth and up/down-welling TB
      real(r8) :: angle_deg, fghz_local                    ! LandEM inputs in degrees and GHz
      real(r8) :: smc_surf, veg_frac                       ! surface soil moisture and vegetation fraction
      real(r8) :: t_surf, t_skin                           ! surface soil and skin temperatures
      real(r8) :: snowdp_mm                                ! snow depth in mm (LandEM input)
      real(r8) :: sand_frac, clay_frac, bulk_density       ! surface soil composition
      real(r8) :: em_h, em_v                               ! bulk emissivity for H- and V- polarization
      real(r8) :: wf_total_surf                            ! total surface gravimetric fraction
      integer  :: lb, i   ! index of topmost active snow layer (dynamic search)

!-----------------------------------------------------------------------
! Non-land surface (ocean, lake, ice): set all TB outputs to spval
!-----------------------------------------------------------------------
      IF (patchtype >= 3) THEN
         tb_toa_h  = spval; tb_toa_v  = spval
         tb_soil_h = spval; tb_soil_v = spval
         tb_tov_h  = spval; tb_tov_v  = spval
         RETURN
      END IF

!#############################################################################
! Prepare surface states
!#############################################################################
      angle_deg  = sat_theta * 180.0_r8 / pi
      fghz_local = sat_fghz

      ! calculate surface soil moisture (liquid + ice, first two layers)
      smc_surf = ( (wliq_soisno(1) + wliq_soisno(2)) / denh2o &
                 + (wice_soisno(1) + wice_soisno(2)) / denice ) / 0.0451_r8
      smc_surf = max(min(smc_surf, 1.0_r8), 0.0_r8)

      ! calculate vegetation fraction from LAI + SAI
      veg_frac = 1.0_r8 - exp(-0.52_r8 * (lai + sai))
      veg_frac = max(min(veg_frac, 1.0_r8), 0.0_r8)

      ! calculate bulk surface soil temperature (first two layers)
      t_surf = (t_soisno(1)*0.0175_r8 + t_soisno(2)*0.0276_r8) / 0.0451_r8

!#############################################################################
! Determine skin temperature (t_skin)
!#############################################################################
      ! 1) Locate the topmost active snow layer (lb)
      lb = 1  ! Default to 1 (first soil layer) as a safe fallback when no snow is active
      DO i = maxsnl+1, 0
         IF (wliq_soisno(i) + wice_soisno(i) > 0.0_r8) THEN
            lb = i
            EXIT
         ENDIF
      ENDDO

      ! 2) Assign skin temperature based on the active surface type
      IF (snowdp > 0.0001_r8 .and. lb <= 0) THEN
         ! only use snow-layer temperature when the physical model has actually activated snow (lb <= 0)
         t_skin = t_soisno(lb)
      ELSE IF (lai + sai > 1.0e-6_r8) THEN
         t_skin = tleaf
      ELSE
         t_skin = t_surf
      END IF

      ! 3) Final safety guard: prevent uninitialized state / division by zero
      IF (t_skin < 150.0_r8 .or. t_skin > 350.0_r8) THEN
         t_skin = t_surf
      END IF

!#############################################################################
! Prepare surface composition parameters (first two soil layers)
!#############################################################################
      ! calculate weighted surface composition (first two soil layers)
      wf_total_surf = (wf_sand(1)+wf_clay(1)+wf_silt(1))*0.0175_r8 + &
                      (wf_sand(2)+wf_clay(2)+wf_silt(2))*0.0276_r8
      IF (wf_total_surf > 1.0e-12_r8) THEN
         sand_frac = (wf_sand(1)*0.0175_r8 + wf_sand(2)*0.0276_r8) / wf_total_surf
         clay_frac = (wf_clay(1)*0.0175_r8 + wf_clay(2)*0.0276_r8) / wf_total_surf
      ELSE
         sand_frac = 0.0_r8
         clay_frac = 0.0_r8
      END IF
      bulk_density = ((BD_all(1)*0.0175_r8 + BD_all(2)*0.0276_r8)/0.0451_r8)

!#############################################################################
! Run LandEM bulk emissivity model
!#############################################################################
      CALL calc_parameters(sat_theta, sat_fghz)
      CALL atm(forc_topo, tref, tau_atm, tb_au, tb_ad)

      IF (sat_fghz >= 80.0_r8) THEN
         ! above 80 GHz, assume a constant emissivity due to strong
         ! atmospheric opacity and limited LandEM calibration
         em_h = 0.95_r8
         em_v = 0.95_r8
      ELSE
         CALL LandEM(angle_deg, fghz_local, smc_surf, veg_frac, &
                     t_surf, t_skin, lai, &
                     sand_frac, clay_frac, bulk_density, &
                     snowdp, em_h, em_v)
      END IF

      em_h = max(min(em_h, 1.0_r8), 0.0_r8)
      em_v = max(min(em_v, 1.0_r8), 0.0_r8)

!#############################################################################
! Convert emissivity to top-of-vegetation and top-of-atmosphere TB
!#############################################################################
      tb_tov_h = em_h * t_skin + (1.0_r8 - em_h) * tb_ad(1)
      tb_tov_v = em_v * t_skin + (1.0_r8 - em_v) * tb_ad(2)

      ! LandEM is a bulk model: separate soil TB is not available
      tb_soil_h = spval
      tb_soil_v = spval

      tb_toa_h = tb_tov_h * exp(-tau_atm) + tb_au(1)
      tb_toa_v = tb_tov_v * exp(-tau_atm) + tb_au(2)

   END SUBROUTINE forward_LandEM

END MODULE MOD_DA_RTM
#endif