#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_RTM_Physics
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Physical subroutines, parameterizations, and numerical solvers for
!    microwave radiative transfer forward models (atmosphere, single-
!    layer and multi-layer soil, dielectric constants, rough boundary,
!    vegetation, and snow).
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!   Lu Li, 10/2025: Debug and clean codes
!   Shuyue Liu, 10/2026: Refined codes and decoupled multi-layer solvers
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_DA_Const
   USE MOD_Vars_Global, only: nl_soil, N_land_classification
   USE MOD_Namelist
   IMPLICIT NONE
   SAVE

! Satellite configuration and wave parameters
   real(r8), public :: fghz     ! [GHz] frequency of satellite
   real(r8), public :: theta    ! [rad] incidence angle of satellite
   real(r8), public :: f        ! [Hz] frequency
   real(r8), public :: omega    ! [rad/s] radian frequency
   real(r8), public :: lam      ! [m] wavelength
   real(r8), public :: k        ! [rad/m] wave number
   real(r8), public :: kcm      ! [rad/cm] wave number
   real(r8), public :: kr       ! [-] size parameter used in calculate single-particle albedo

! Public physics subroutines called by forward operators
   PUBLIC :: calc_parameters
   PUBLIC :: atm
   PUBLIC :: soil_CMEM
   PUBLIC :: soil_CoLM
   PUBLIC :: snow_Pul97
   PUBLIC :: veg_Wig07
   PUBLIC :: veg_Jac91
   PUBLIC :: veg_Kir79
   PUBLIC :: LandEM

CONTAINS

!-----------------------------------------------------------------------
   SUBROUTINE calc_parameters(sat_theta, sat_fghz)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the satellite-dependent wave parameters (frequency,
!   wavelength, wave number, and size parameter) stored in module
!   variables for use by all downstream physics subroutines.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_DA_Const
      IMPLICIT NONE
      real(r8), intent(in) :: sat_theta ! [rad] incidence angle of satellite
      real(r8), intent(in) :: sat_fghz  ! [GHz] frequency of satellite

      theta = sat_theta                   ! [rad] incidence angle of satellite
      fghz  = sat_fghz                    ! [GHz] frequency of satellite
      f     = fghz * 1.0e9_r8             ! [Hz] frequency
      omega = 2.0_r8 * pi * f             ! [rad/s] radian frequency
      lam   = C / f                       ! [m] wavelength
      k     = 2.0_r8 * pi / lam           ! [rad/m] wave number
      kcm   = k / 100.0_r8                ! [rad/cm] wave number
      kr    = k * (0.5_r8 * 1.0e-3_r8)    ! [-] size parameter for single-particle albedo

   END SUBROUTINE calc_parameters

!-----------------------------------------------------------------------
   SUBROUTINE atm(z, tref, tau_atm, tb_au, tb_ad)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the atmospheric opacity and up- and down-welling
!   brightness temperature of the atmosphere.
!
! REFERENCES:
!   [1] Pellarin, T., et al. (2003), Two-year global simulation of L-band brightness
!       temperature over land, IEEE Trans. Geosci. Remote Sens., 41, 2135-2139.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)  :: z            ! [m] altitude
      real(r8), intent(in)  :: tref         ! [K] 2 m air temperature
      real(r8), intent(out) :: tau_atm      ! [-] atmospheric optical depth
      real(r8), intent(out) :: tb_au(2)     ! [K] up-welling radiation (brightness temperature) of atmosphere, (1)=H,(2)=V
      real(r8), intent(out) :: tb_ad(2)     ! [K] down-welling radiation (brightness temperature) of atmosphere, (1)=H,(2)=V

      real(r8) :: t_sky = 2.7_r8            ! [K] cosmic ray radiation
      real(r8) :: t_eq                      ! [K] equivalent layer temperature
      real(r8) :: gossat                    ! [-] atmospheric transmittance

      ! calculate optical depth of atmosphere [1] eq(A1)
      tau_atm = exp(-3.9262_r8 - 0.2211_r8*z/1000.0_r8 - 0.00369_r8*tref) / cos(theta)   ! [-]
      gossat  = exp(-tau_atm)                                                              ! [-]
      ! calculate equivalent layer temperature
      t_eq    = exp(4.9274_r8 + 0.002195_r8*tref)                                         ! [K]

      ! up-welling radiation (brightness temperature) of atmosphere
      tb_au(:) = t_eq * (1.0_r8 - gossat)                                                  ! [K]
      ! down-welling radiation (brightness temperature) of atmosphere [1] eq(A2)
      tb_ad(:) = t_eq * (1.0_r8 - gossat) + t_sky * gossat                                 ! [K]

   END SUBROUTINE atm

!-----------------------------------------------------------------------
   SUBROUTINE soil_CoLM( &
      patchclass, nl_soil, dz_soi, t_soi, liq_soi, ice_soi, wf_sand, wf_clay, BD_all, porsl, &
      t_surf, liq_surf, ice_surf, &
      wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
      r_r, tb_soil)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Single-point driver of the CoLM-Liu discrete multi-layer incoherent
!   soil RTM. It builds the surface and per-layer dielectric constants,
!   computes the smooth and rough surface reflectivity, and then calls
!   tb_soil_Liu26 for the analytical two-stream solution.
!
!   NOTE: r_r is INTENT(INOUT): it is first set by rough_reflectivity
!   (Fresnel + roughness) and is then OVERWRITTEN inside
!   tb_soil_Liu26 by the effective multi-layer reflectivity.
!
! REFERENCES: Liu et al. (in preparation)
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      USE MOD_DA_Const
      IMPLICIT NONE
      integer, intent(in)     :: patchclass         ! [-] land cover class
      integer, intent(in)     :: nl_soil            ! [-] number of soil layers
      real(r8), intent(in)    :: dz_soi(nl_soil)    ! [m] layer thickness
      real(r8), intent(in)    :: t_soi(nl_soil)     ! [K] soil temperature
      real(r8), intent(in)    :: liq_soi(nl_soil)   ! [m3/m3] liquid volumetric water
      real(r8), intent(in)    :: ice_soi(nl_soil)   ! [m3/m3] ice volumetric water
      real(r8), intent(in)    :: wf_sand(nl_soil)   ! [-] gravimetric sand fraction
      real(r8), intent(in)    :: wf_clay(nl_soil)   ! [-] gravimetric clay fraction
      real(r8), intent(in)    :: BD_all(nl_soil)    ! [kg/m3] bulk density
      real(r8), intent(in)    :: porsl(nl_soil)     ! [m3/m3] porosity
      real(r8), intent(in)    :: t_surf             ! [K] surface soil temperature
      real(r8), intent(in)    :: liq_surf           ! [m3/m3] surface liquid water
      real(r8), intent(in)    :: ice_surf           ! [m3/m3] surface ice water
      real(r8), intent(in)    :: wf_sand_surf       ! [-] surface sand fraction
      real(r8), intent(in)    :: wf_clay_surf       ! [-] surface clay fraction
      real(r8), intent(in)    :: BD_all_surf        ! [kg/m3] surface bulk density
      real(r8), intent(in)    :: porsl_surf         ! [m3/m3] surface porosity
      real(r8), intent(inout) :: r_r(2)             ! [-] rough surface reflectivity, (1)=H,(2)=V; overwritten by multi-layer effective reflectivity
      real(r8), intent(out)   :: tb_soil(2)         ! [K] soil brightness temperature, (1)=H,(2)=V

      complex(r8) :: eps_prof(nl_soil)               ! [-] dielectric constant profile for each soil layer
      complex(r8) :: ew_layer                       ! [-] water dielectric constant for a layer
      complex(r8) :: eps_soil                       ! [-] surface soil dielectric constant
      complex(r8) :: ew                             ! [-] surface water dielectric constant
      real(r8)    :: ffrz_layer                     ! [-] frozen fraction per layer
      real(r8)    :: ffrz                           ! [-] surface frozen fraction
      complex(r8) :: eps_f = (5.0_r8, 0.5_r8)         ! [-] dielectric constant of frozen soil
      real(r8)    :: sal_soil = 0.0_r8                ! [psu] soil salinity
      complex(r8) :: eps_grain = (4.0_r8, 0.05_r8)    ! [-] dielectric constant of soil grain used by Liu RTM
      real(r8)    :: r_s(2)                           ! [-] smooth surface reflectivity for H and V polarizations, (1)=H,(2)=V
      integer     :: i_layer                          ! [-] layer index

!#############################################################################
! 1. Surface dielectric constant
!#############################################################################
      IF (liq_surf + ice_surf <= 0.0_r8) THEN
         ffrz = 0.0_r8
      ELSE
         ffrz = ice_surf / (liq_surf + ice_surf)
      ENDIF

      ! calculate dielectric constant of water (or pure ice if mostly frozen)
      IF (ffrz > 0.95_r8) THEN
         CALL diel_ice(t_surf, ew)
      ELSE
         CALL diel_water(2, liq_surf, t_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, ew)
      END IF

      ! calculate dielectric constant of wet soil (scheme selected by DEF_DA_RTM_diel)
      IF (DEF_DA_RTM_diel == 0) THEN
         CALL diel_soil_W80(ew, t_surf, liq_surf, wf_sand_surf, wf_clay_surf, porsl_surf, eps_soil)
         eps_soil = eps_soil*(1.0_r8 - ffrz) + eps_f*ffrz
      ELSE IF (DEF_DA_RTM_diel == 1) THEN
         CALL diel_soil_D85(ew, liq_surf, ice_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, eps_soil)
      ELSE IF (DEF_DA_RTM_diel == 2) THEN
         CALL diel_soil_M04(liq_surf, wf_clay_surf, eps_soil)
         eps_soil = eps_soil*(1.0_r8 - ffrz) + eps_f*ffrz
      ELSE IF (DEF_DA_RTM_diel == 3) THEN
         CALL diel_soil_M09(liq_surf, t_surf, wf_clay_surf, eps_soil)
         eps_soil = eps_soil*(1.0_r8 - ffrz) + eps_f*ffrz
      ENDIF

!#############################################################################
! 2. Multi-layer dielectric profile eps_prof (used by tb_soil_Liu26)
!#############################################################################
      DO i_layer = 1, nl_soil
         IF (liq_soi(i_layer) + ice_soi(i_layer) <= 0.0_r8) THEN
            ffrz_layer = 0.0_r8
         ELSE
            ffrz_layer = ice_soi(i_layer) / (liq_soi(i_layer) + ice_soi(i_layer))
         END IF

         IF (ffrz_layer > 0.95_r8) THEN
            CALL diel_ice(t_soi(i_layer), ew_layer)
         ELSE
            CALL diel_water(2, liq_soi(i_layer), t_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), BD_all(i_layer), sal_soil, ew_layer)
         END IF

         IF (DEF_DA_RTM_diel == 0) THEN
            CALL diel_soil_W80(ew_layer, t_soi(i_layer), liq_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), porsl(i_layer), eps_prof(i_layer))
            eps_prof(i_layer) = eps_prof(i_layer)*(1.0_r8 - ffrz_layer) + eps_f*ffrz_layer
         ELSE IF (DEF_DA_RTM_diel == 1) THEN
            CALL diel_soil_D85(ew_layer, liq_soi(i_layer), ice_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), BD_all(i_layer), eps_prof(i_layer))
         ELSE IF (DEF_DA_RTM_diel == 2) THEN
            CALL diel_soil_M04(liq_soi(i_layer), wf_clay(i_layer), eps_prof(i_layer))
            eps_prof(i_layer) = eps_prof(i_layer)*(1.0_r8 - ffrz_layer) + eps_f*ffrz_layer
         ELSE IF (DEF_DA_RTM_diel == 3) THEN
            CALL diel_soil_M09(liq_soi(i_layer), t_soi(i_layer), wf_clay(i_layer), eps_prof(i_layer))
            eps_prof(i_layer) = eps_prof(i_layer)*(1.0_r8 - ffrz_layer) + eps_f*ffrz_layer
         ENDIF
      END DO

!#############################################################################
! 3. Smooth and rough surface reflectivity
!#############################################################################
      CALL smooth_reflectivity(eps_soil, r_s)
      CALL rough_reflectivity(.false., patchclass, r_s, r_r)

!#############################################################################
! 4. Analytical multi-layer two-stream incoherent RTM (Liu scheme)
!    Note: r_r is overwritten here by the effective multi-layer reflectivity.
!#############################################################################
      CALL tb_soil_Liu26(nl_soil, dz_soi, t_soi, eps_prof, wf_sand, porsl, eps_grain, r_r, tb_soil)

   END SUBROUTINE soil_CoLM

!-----------------------------------------------------------------------
   SUBROUTINE tb_soil_Liu26(nl_soil, dz_soi, t_soi, eps_prof, wf_sand, porsl, eps_grain, r_r, tb_soil)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Analytically decoupled discrete multi-layer two-stream incoherent
!   RTM. For each polarization, the reflectivity of each interface is
!   computed from Fresnel laws, a linear system is assembled for the
!   up- and down-welling effective reflectivities (W1 and W2), and the
!   layer emission weights are then used to integrate the soil
!   temperature profile into a single brightness temperature.
!
! REFERENCES: Liu et al. (in preparation)
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)     :: nl_soil         ! [-] number of soil layers
      real(r8), intent(in)    :: dz_soi(:)       ! [m] soil layer thickness
      real(r8), intent(in)    :: t_soi(:)        ! [K] soil temperature
      complex(r8), intent(in) :: eps_prof(:)     ! [-] per-layer complex dielectric constant
      real(r8), intent(in)    :: wf_sand(:)      ! [-] gravimetric sand fraction
      real(r8), intent(in)    :: porsl(:)        ! [m3/m3] porosity
      complex(r8), intent(in) :: eps_grain       ! [-] soil grain dielectric constant
      real(r8), intent(inout) :: r_r(2)          ! [-] in: rough surface reflectivity; out: effective multi-layer reflectivity, (1)=H,(2)=V
      real(r8), intent(out)   :: tb_soil(2)      ! [K] soil brightness temperature, (1)=H,(2)=V

      integer  :: i                                ! [-] layer index
      integer  :: pol                              ! [-] polarization index (1=H,2=V)
      integer  :: row                              ! [-] linear-system row index
      real(r8) :: k_a                              ! [1/m] absorption coefficient
      real(r8) :: k_s                              ! [1/m] scattering coefficient
      real(r8) :: g_asym                           ! [-] asymmetry factor
      real(r8) :: omega_albedo                     ! [-] single-scattering albedo
      real(r8) :: a_j                              ! [-] similarity parameter
      real(r8) :: L_pen                            ! [m] penetration length
      real(r8) :: mu_z                             ! [-] refraction cosine
      real(r8) :: kappa_z                          ! [1/m] attenuation coefficient
      real(r8) :: r_j                              ! [m] mean grain radius
      real(r8) :: f_j                              ! [-] dry-soil volume fraction
      real(r8) :: Q_s                              ! [-] Rayleigh scattering efficiency
      real(r8) :: y_R                              ! [-] real Clausius-Mossotti factor
      real(r8) :: r_sand                           ! [m] sand grain radius
      real(r8) :: r_clay                           ! [m] clay grain radius
      real(r8) :: sin2theta                        ! [-] sin^2(theta)
      complex(r8) :: c_j                           ! [-] sqrt(eps_j - sin^2theta)
      complex(r8) :: c_jp1                         ! [-] sqrt(eps_{j+1} - sin^2theta)
      real(r8) :: R_int(nl_soil)                   ! [-] Fresnel reflectivity of internal interfaces
      real(r8) :: r_inf(nl_soil)                   ! [-] similarity-parameter reflectivity
      real(r8) :: tau_arr(nl_soil)                 ! [-] per-layer transmittance
      integer  :: n_sys                            ! [-] linear system size (= 2*nl_soil)
      real(r8) :: M(2*nl_soil, 2*nl_soil)          ! [-] coefficient matrix
      real(r8) :: M_trans(2*nl_soil, 2*nl_soil)    ! [-] transposed coefficient matrix
      real(r8) :: C_rhs(2*nl_soil)                 ! [-] right-hand sides
      real(r8) :: W1(2*nl_soil)                    ! [-] solution for unit incident up-welling
      real(r8) :: W2(2*nl_soil)                    ! [-] solution for unit incident down-welling
      real(r8) :: K_vec(2*nl_soil)                 ! [-] auxiliary vector
      real(r8) :: U_vec(nl_soil)                   ! [-] auxiliary vector
      real(r8) :: f_weights(nl_soil)               ! [-] emission weights
      real(r8) :: f_down                           ! [-] total down-welling reflectivity
      real(r8) :: R0_tilde                         ! [-] rough surface reflectivity for current polarization

!#############################################################################
! 1. Pre-compute particle-size and mixing parameters (polarization-independent)
!#############################################################################
      sin2theta = sin(theta)**2
      y_R = real((eps_grain - 1.0_r8) / (eps_grain + 2.0_r8))
      r_sand = 0.5e-3_r8                          ! sand grain radius (m)
      r_clay = 0.01e-3_r8                         ! clay grain radius (m)
      n_sys = 2 * nl_soil

!#############################################################################
! 2. Per-layer properties: refraction angle, absorption, scattering, albedo
!#############################################################################
      DO pol = 1, 2
         DO i = 1, nl_soil
            ! calculate the real part of the refraction cosine in layer i
            mu_z = sqrt((abs(eps_prof(i) - sin2theta) + real(eps_prof(i)) - sin2theta) / &
                        (abs(eps_prof(i) - sin2theta) + real(eps_prof(i)) + sin2theta))
            mu_z = max(mu_z, 0.001_r8)

            ! calculate the absorption coefficient of layer i
            k_a = 2.0_r8 * k * abs(aimag(sqrt(eps_prof(i))))
            k_a = max(k_a, 1.0e-12_r8)

            ! calculate the mean grain radius and the dry-soil volume fraction
            r_j = r_clay + (r_sand - r_clay) * wf_sand(i)
            f_j = (1.0_r8 - porsl(i)) * wf_sand(i)

            ! calculate the Rayleigh single-particle scattering efficiency (Grody & Weng 2008)
            Q_s = (8.0_r8 / 3.0_r8) * (k * r_j)**4 * (y_R**2) * &
                  ( ((1.0_r8 - f_j)**4) / ( ((1.0_r8 + 2.0_r8 * f_j)**2) * &
                  ((1.0_r8 - f_j * y_R)**1.5_r8) * ((1.0_r8 + 2.0_r8 * f_j * y_R)**0.5_r8) ) )

            ! calculate the volume scattering coefficient and the asymmetry parameter
            k_s = (3.0_r8 * f_j / (4.0_r8 * r_j)) * Q_s
            g_asym = 0.23_r8 * (k * r_j)**2

            ! calculate the single-scattering albedo and the similarity parameter
            omega_albedo = min(max(k_s / (k_s + k_a), 0.0_r8), 0.9999_r8)
            a_j = sqrt((1.0_r8 - omega_albedo) / (1.0_r8 - omega_albedo * g_asym))

            ! calculate the similarity-parameter reflectivity and per-layer transmittance
            r_inf(i) = (1.0_r8 - a_j) / (1.0_r8 + a_j)
            L_pen = a_j / k_a
            kappa_z = 1.0_r8 / (mu_z * L_pen)
            tau_arr(i) = exp(-kappa_z * dz_soi(i))
         END DO

!#############################################################################
! 3. Internal-interface Fresnel reflectivity (polarization dependent)
!#############################################################################
         DO i = 1, nl_soil - 1
            c_j = sqrt(eps_prof(i) - sin2theta)
            c_jp1 = sqrt(eps_prof(i+1) - sin2theta)
            IF (pol == 1) THEN
               R_int(i) = abs((c_j - c_jp1) / (c_j + c_jp1))**2
            ELSE
               R_int(i) = abs((c_j/eps_prof(i) - c_jp1/eps_prof(i+1)) / &
                              (c_j/eps_prof(i) + c_jp1/eps_prof(i+1)))**2
            END IF
         END DO

!#############################################################################
! 4. Assemble the linear system M * W = C for the up- and down-welling
!    effective reflectivities at each interface
!#############################################################################
         M = 0.0_r8
         R0_tilde = r_r(pol)

         M(1, 1) = r_inf(1) - R0_tilde
         M(1, 2) = 1.0_r8 - R0_tilde * r_inf(1)

         DO i = 1, nl_soil - 1
            row = 2 * i
            M(row, 2*i-1) = (1.0_r8 - R_int(i)*r_inf(i)) / tau_arr(i)
            M(row, 2*i)   = (r_inf(i) - R_int(i)) * tau_arr(i)
            M(row, 2*i+1) = -(1.0_r8 - R_int(i))
            M(row, 2*i+2) = -(1.0_r8 - R_int(i)) * r_inf(i+1)

            M(row+1, 2*i-1) = -(1.0_r8 - R_int(i)) * r_inf(i) / tau_arr(i)
            M(row+1, 2*i)   = -(1.0_r8 - R_int(i)) * tau_arr(i)
            M(row+1, 2*i+1) = r_inf(i+1) - R_int(i)
            M(row+1, 2*i+2) = 1.0_r8 - R_int(i) * r_inf(i+1)
         END DO

         M(n_sys, n_sys-1) = 1.0_r8

!#############################################################################
! 5. Solve M^T * W = C for two right-hand sides (unit incident up/down)
!#############################################################################
         M_trans = transpose(M)
         C_rhs = 0.0_r8
         C_rhs(1) = 1.0_r8
         CALL solve_linear_system(n_sys, M_trans, C_rhs, W1)

         M_trans = transpose(M)
         C_rhs = 0.0_r8
         C_rhs(2) = 1.0_r8
         CALL solve_linear_system(n_sys, M_trans, C_rhs, W2)

         DO i = 1, n_sys
            K_vec(i) = W1(i) + r_inf(1) * W2(i)
         END DO

         U_vec = 0.0_r8
         DO i = 1, nl_soil - 1
            U_vec(i) = (K_vec(2*i) - K_vec(2*i+1)) * (1.0_r8 - R_int(i))
         END DO

!#############################################################################
! 6. Layer emission weights f_weights(i) and total down-welling reflectivity
!#############################################################################
         f_down = (1.0_r8 - R0_tilde)**2 * K_vec(1) + R0_tilde

         IF (nl_soil == 1) THEN
            f_weights(1) = (1.0_r8 - R0_tilde) * (1.0_r8 - K_vec(1)*(1.0_r8 - R0_tilde))
         ELSE
            f_weights(1) = (1.0_r8 - R0_tilde) * (1.0_r8 - K_vec(1)*(1.0_r8 - R0_tilde) - U_vec(1))
            DO i = 2, nl_soil - 1
               f_weights(i) = (1.0_r8 - R0_tilde) * (U_vec(i-1) - U_vec(i))
            END DO
            f_weights(nl_soil) = (1.0_r8 - R0_tilde) * U_vec(nl_soil-1)
         END IF

         tb_soil(pol) = 0.0_r8
         DO i = 1, nl_soil
            tb_soil(pol) = tb_soil(pol) + f_weights(i) * t_soi(i)
         END DO

         r_r(pol) = f_down
      END DO

   END SUBROUTINE tb_soil_Liu26

!-----------------------------------------------------------------------
   SUBROUTINE solve_linear_system(n, A, b, x)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Solve the dense linear system A * x = b by Gaussian elimination with
!   partial pivoting followed by back substitution. A and b are
!   overwritten in place.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      integer, intent(in)     :: n                     ! [-] system size
      real(r8), intent(inout) :: A(n, n)               ! [-] coefficient matrix; units follow the equation (dimensionless in this module)
      real(r8), intent(inout) :: b(n)                  ! [-] right-hand side; units follow the equation (dimensionless in this module)
      real(r8), intent(out)   :: x(n)                  ! [-] solution; units follow the equation (dimensionless in this module)

      integer  :: i                                    ! [-] loop/index variable
      integer  :: j                                    ! [-] loop/index variable
      integer  :: k                                    ! [-] loop/index variable
      integer  :: max_idx                              ! [-] loop/index variable
      real(r8) :: temp                                 ! [-] temporary value
      real(r8) :: factor                               ! [-] elimination factor

      ! forward elimination with partial pivoting
      DO i = 1, n - 1
         max_idx = i
         DO j = i + 1, n
            IF (abs(A(j, i)) > abs(A(max_idx, i))) max_idx = j
         END DO
         IF (max_idx /= i) THEN
            DO k = i, n
               temp = A(i, k)
               A(i, k) = A(max_idx, k)
               A(max_idx, k) = temp
            END DO
            temp = b(i)
            b(i) = b(max_idx)
            b(max_idx) = temp
         END IF
         DO j = i + 1, n
            factor = A(j, i) / A(i, i)
            DO k = i, n
               A(j, k) = A(j, k) - factor * A(i, k)
            END DO
            b(j) = b(j) - factor * b(i)
         END DO
      END DO

      ! back substitution
      x(n) = b(n) / A(n, n)
      DO i = n - 1, 1, -1
         temp = b(i)
         DO j = i + 1, n
            temp = temp - A(i, j) * x(j)
         END DO
         x(i) = temp / A(i, i)
      END DO
   END SUBROUTINE solve_linear_system

!-----------------------------------------------------------------------
   SUBROUTINE soil_CMEM( &
      patchclass, nl_soil, dz_soi, t_soi, liq_soi, ice_soi, wf_sand, wf_clay, BD_all, porsl, &
      t_surf, t_deep, &
      liq_surf, ice_surf, &
      wf_sand_surf, wf_clay_surf, BD_all_surf, porsl_surf, &
      r_r, tb_soil)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Single-point driver of the CMEM single-layer soil emission scheme.
!   It builds the surface (and, for the Wilheit/Lv effective-temperature
!   schemes, per-layer) dielectric constants, computes the effective
!   temperature, applies Fresnel and rough surface reflectivity, and
!   finally computes the soil brightness temperature.
!
!   NOTE 1 (r_s source): when DEF_DA_RTM_teff == 0 (Wilheit), the
!   smooth-surface reflectivity r_s is returned directly by
!   soil_eff_temp_Wil78 as the effective reflectivity implied by the
!   multi-layer HPLD solution. For DEF_DA_RTM_teff == 1/2/3, r_s is
!   recomputed from the surface Fresnel law with eps_soil.
!
!   NOTE 2 (desert path): when the patch is classified as desert,
!   the Grody & Weng (2008) desert emissivity is used instead of the
!   t_eff*(1-r_r) formulation.
!
! REFERENCES:
!   [1] Wigneron et al. (2007), L-band Microwave Emission of the Biosphere
!       (L-MEB) model, Remote Sensing of Environment, 107, 639-655.
!   [2] Grody & Weng (2008), Microwave emission and scattering from deserts,
!       IEEE Trans. Geosci. Remote Sens., 46, 361-375.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)   :: patchclass         ! [-] land cover class
      integer, intent(in)   :: nl_soil            ! [-] number of soil layers
      real(r8), intent(in)  :: dz_soi(nl_soil)    ! [m] layer thickness
      real(r8), intent(in)  :: t_soi(nl_soil)     ! [K] soil temperature
      real(r8), intent(in)  :: liq_soi(nl_soil)   ! [m3/m3] liquid volumetric water
      real(r8), intent(in)  :: ice_soi(nl_soil)   ! [m3/m3] ice volumetric water
      real(r8), intent(in)  :: wf_sand(nl_soil)   ! [-] gravimetric sand fraction
      real(r8), intent(in)  :: wf_clay(nl_soil)   ! [-] gravimetric clay fraction
      real(r8), intent(in)  :: BD_all(nl_soil)    ! [kg/m3] bulk density
      real(r8), intent(in)  :: porsl(nl_soil)     ! [m3/m3] porosity
      real(r8), intent(in)  :: t_surf             ! [K] surface soil temperature
      real(r8), intent(in)  :: t_deep             ! [K] deep soil temperature
      real(r8), intent(in)  :: liq_surf           ! [m3/m3] surface liquid water
      real(r8), intent(in)  :: ice_surf           ! [m3/m3] surface ice water
      real(r8), intent(in)  :: wf_sand_surf       ! [-] surface sand fraction
      real(r8), intent(in)  :: wf_clay_surf       ! [-] surface clay fraction
      real(r8), intent(in)  :: BD_all_surf        ! [kg/m3] surface bulk density
      real(r8), intent(in)  :: porsl_surf         ! [m3/m3] surface porosity
      real(r8), intent(out) :: r_r(2)             ! [-] rough surface reflectivity, (1)=H,(2)=V
      real(r8), intent(out) :: tb_soil(2)         ! [K] soil brightness temperature, (1)=H,(2)=V

      complex(r8) :: eps_prof(nl_soil)              ! [-] per-layer complex dielectric constant for multi-layer schemes
      complex(r8) :: ew_layer                       ! [-] water dielectric constant for a layer
      real(r8)    :: ffrz_layer                     ! [-] frozen fraction per layer
      logical     :: is_desert_layer                ! [-] desert flag per layer
      integer     :: i_layer                        ! [-] layer index
      real(r8)    :: t_eff(2)                       ! [K] effective temperature for H and V polarizations, (1)=H,(2)=V
      complex(r8) :: eps_soil                       ! [-] surface soil dielectric constant
      complex(r8) :: ew                             ! [-] surface water dielectric constant
      real(r8)    :: r_s(2)                         ! [-] smooth surface reflectivity, (1)=H,(2)=V
      logical     :: is_desert                      ! [-] desert flag for surface
      real(r8)    :: ffrz                           ! [-] surface frozen fraction
      complex(r8) :: eps_f = (5.0, 0.5)             ! [-] frozen soil dielectric constant
      real(r8)    :: sal_soil = 0.0                 ! [psu] soil salinity

!#############################################################################
! 1. Surface dielectric constant
!#############################################################################
      ! whether this patch is desert (very dry and very sandy)
      is_desert = .false.
      IF (liq_surf < 0.02 .and. wf_sand_surf > 0.90_r8) THEN
         is_desert = .true.
      END IF

      ! calculate the frozen fraction of the surface
      IF (liq_surf + ice_surf <= 0.0d0) THEN
        ffrz = 0.0d0
      ELSE
        ffrz = ice_surf / (liq_surf + ice_surf)
      ENDIF

      ! calculate dielectric constant of soil (mixture medium)
      IF (is_desert) THEN
         ! microwave 1-10 GHz permittivity of dry sand (Matzler 1998, eq.1)
         eps_soil = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002
      ELSE
         ! calculate ice or water dielectric constant
         IF (ffrz > 0.95) THEN
            CALL diel_ice(t_surf, ew)
         ELSE
            CALL diel_water(2, liq_surf, t_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, sal_soil, ew)
         END IF

         ! calculate dielectric constant in mixed soil (scheme selected by DEF_DA_RTM_diel)
         IF (DEF_DA_RTM_diel == 0) THEN
            CALL diel_soil_W80(ew, t_surf, liq_surf, wf_sand_surf, wf_clay_surf, porsl_surf, eps_soil)
            eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz
         ELSE IF (DEF_DA_RTM_diel == 1) THEN
            CALL diel_soil_D85(ew, liq_surf, ice_surf, wf_sand_surf, wf_clay_surf, BD_all_surf, eps_soil)
         ELSE IF (DEF_DA_RTM_diel == 2) THEN
            CALL diel_soil_M04(liq_surf, wf_clay_surf, eps_soil)
            eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz
         ELSE IF (DEF_DA_RTM_diel == 3) THEN
            CALL diel_soil_M09(liq_surf, t_surf, wf_clay_surf, eps_soil)
            eps_soil = eps_soil*(1.-ffrz) + eps_f*ffrz
         ENDIF
      END IF

!#############################################################################
! 2. Multi-layer dielectric profile (only for Wilheit / Lv teff schemes)
!#############################################################################
      IF (DEF_DA_RTM_teff == 0 .or. DEF_DA_RTM_teff == 3) THEN
         DO i_layer = 1, nl_soil
            is_desert_layer = .false.
            IF (liq_soi(i_layer) < 0.02 .and. wf_sand(i_layer) > 0.90_r8) THEN
               is_desert_layer = .true.
            END IF
            IF (liq_soi(i_layer) + ice_soi(i_layer) <= 0.0d0) THEN
               ffrz_layer = 0.0d0
            ELSE
               ffrz_layer = ice_soi(i_layer) / (liq_soi(i_layer) + ice_soi(i_layer))
            END IF
            IF (is_desert_layer) THEN
               eps_prof(i_layer) = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002
            ELSE
               IF (ffrz_layer > 0.95) THEN
                  CALL diel_ice(t_soi(i_layer), ew_layer)
               ELSE
                  CALL diel_water(2, liq_soi(i_layer), t_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), BD_all(i_layer), sal_soil, ew_layer)
               END IF

               IF (DEF_DA_RTM_diel == 0) THEN
                  CALL diel_soil_W80(ew_layer, t_soi(i_layer), liq_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), porsl(i_layer), eps_prof(i_layer))
                  eps_prof(i_layer) = eps_prof(i_layer)*(1.-ffrz_layer) + eps_f*ffrz_layer
               ELSE IF (DEF_DA_RTM_diel == 1) THEN
                  CALL diel_soil_D85(ew_layer, liq_soi(i_layer), ice_soi(i_layer), wf_sand(i_layer), wf_clay(i_layer), BD_all(i_layer), eps_prof(i_layer))
               ELSE IF (DEF_DA_RTM_diel == 2) THEN
                  CALL diel_soil_M04(liq_soi(i_layer), wf_clay(i_layer), eps_prof(i_layer))
                  eps_prof(i_layer) = eps_prof(i_layer)*(1.-ffrz_layer) + eps_f*ffrz_layer
               ELSE IF (DEF_DA_RTM_diel == 3) THEN
                  CALL diel_soil_M09(liq_soi(i_layer), t_soi(i_layer), wf_clay(i_layer), eps_prof(i_layer))
                  eps_prof(i_layer) = eps_prof(i_layer)*(1.-ffrz_layer) + eps_f*ffrz_layer
               ENDIF
            END IF
         END DO
      END IF

!#############################################################################
! 3. Effective temperature of the soil
!    teff == 0 : Wilheit (1975)  - multi-layer HPLD, also returns r_s
!    teff == 1 : Wigneron (2001) - surface-moisture-dependent
!    teff == 2 : Holmes (2006)   - surface-dielectric-dependent
!    teff == 3 : Lv (2014)       - multi-layer discrete
!#############################################################################
      IF (DEF_DA_RTM_teff == 0) THEN
         CALL soil_eff_temp_Wil78(nl_soil, dz_soi, t_soi, eps_prof, t_eff, r_s)
      ELSE IF (DEF_DA_RTM_teff == 1) THEN
         CALL soil_eff_temp_Wig01(liq_surf, t_surf, t_deep, t_eff)
      ELSE IF (DEF_DA_RTM_teff == 2) THEN
         CALL soil_eff_temp_Hol06(t_surf, t_deep, eps_soil, t_eff)
      ELSE IF (DEF_DA_RTM_teff == 3) THEN
         CALL soil_eff_temp_Lv14(nl_soil, dz_soi, t_soi, eps_prof, t_eff)
      END IF

!#############################################################################
! 4. Smooth and rough surface reflectivity
!    (r_s is already returned by soil_eff_temp_Wil78 for teff == 0)
!#############################################################################
      IF (DEF_DA_RTM_teff /= 0) THEN
         CALL smooth_reflectivity(eps_soil, r_s)
      END IF

      CALL rough_reflectivity(is_desert, patchclass, r_s, r_r)

!#############################################################################
! 5. Brightness temperature of soil
!#############################################################################
      IF (is_desert) THEN
         CALL desert(t_eff, r_r, eps_soil, tb_soil)
      ELSE
         tb_soil = t_eff * (1 - r_r)
      END IF

   END SUBROUTINE soil_CMEM

!-----------------------------------------------------------------------
   SUBROUTINE soil_eff_temp_Hol06(t_surf, t_deep, eps_surf, t_eff)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil based on Holmes (2006).
!   The surface-to-deep weighting is controlled by the imaginary-to-real
!   ratio of the surface dielectric constant, mapped through a two-
!   parameter power law.
!
! REFERENCES:
!   [1] Holmes, T. R. H., et al. (2006), A new parameterization of the
!       effective temperature for L-band radiometry, Geophysical Research
!       Letters, 33, L07405.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      real(r8), intent(in)    :: t_surf           ! [K] surface soil temperature
      real(r8), intent(in)    :: t_deep           ! [K] deep soil temperature
      complex(r8), intent(in) :: eps_surf         ! [-] surface soil dielectric constant
      real(r8), intent(out)   :: t_eff(2)         ! [K] effective temperature for H and V, (1)=H,(2)=V

      real(r8) :: eps_r                           ! [-] real part of dielectric constant
      real(r8) :: eps_i                           ! [-] absolute imaginary part
      real(r8) :: C_param                         ! [-] weighting parameter
      real(r8), parameter :: eps_0 = 0.08_r8      ! [-] calibrated parameter from Holmes (2006) Table 1
      real(r8), parameter :: b_param = 0.87_r8    ! [-] calibrated parameter from Holmes (2006) Table 1

      eps_r = real(eps_surf)
      eps_i = abs(aimag(eps_surf))

      IF (eps_r > 0.0_r8) THEN
         C_param = ((eps_i / eps_r) / eps_0) ** b_param
      ELSE
         C_param = 1.0_r8
      END IF

      ! bounds constraint for numerical safety
      C_param = max(0.001_r8, C_param)

      ! calculate effective temperature (identical for H and V polarizations)
      t_eff(1) = t_deep + (t_surf - t_deep) * C_param
      t_eff(2) = t_eff(1)

   END SUBROUTINE soil_eff_temp_Hol06

!-----------------------------------------------------------------------
   SUBROUTINE soil_eff_temp_Wil78(nl_soil, dz_soi, t_soi, eps_prof, t_eff, r_s)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil based on Wilheit (1975)
!   multi-layer radiative transfer in a plane stratified dielectric.
!
!   NOTE: unlike the other teff schemes, this routine also returns the
!   smooth-surface reflectivity r_s as 1 - sumfa, i.e., the effective
!   reflectivity implied by the multi-layer HPLD/VPLD solution. The
!   caller therefore must NOT recompute r_s from the Fresnel law when
!   DEF_DA_RTM_teff == 0.
!
! REFERENCES:
!   [1] Wilheit, T. T. (1978), Radiative transfer in a plane stratified
!       dielectric, IEEE Trans. Geosci. Electron., GE-16(2), 138-143.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)     :: nl_soil                 ! [-] number of soil layers
      real(r8), intent(in)    :: dz_soi(:)               ! [m] layer thickness
      real(r8), intent(in)    :: t_soi(:)                ! [K] soil temperature
      complex(r8), intent(in) :: eps_prof(:)             ! [-] per-layer complex dielectric constant
      real(r8), intent(out)   :: t_eff(2)                ! [K] effective temperature, (1)=H,(2)=V
      real(r8), intent(out)   :: r_s(2)                  ! [-] effective smooth-surface reflectivity, (1)=H,(2)=V

      integer                 :: i                       ! [-] layer index
      complex(r8)             :: eps_layer               ! [-] layer dielectric constant
      complex(r8)             :: nref(nl_soil + 1)       ! [-] index of refraction (layer 1 = atmosphere)
      real(r8)                :: zsoil_wil(nl_soil + 1)  ! [m] layer thicknesses (layer 1 = atmosphere)
      complex(r8)             :: cp_h(nl_soil + 1)       ! [-] absorption fractions and surface reflectivity, H
      complex(r8)             :: cp_v(nl_soil + 1)       ! [-] absorption fractions and surface reflectivity, V
      real(r8)                :: sumtfa                  ! [K] weighted temperature sum
      real(r8)                :: sumfa                   ! [-] absorption fraction sum
      real(r8)                :: fa                      ! [-] layer absorption fraction

!#############################################################################
! 1. Initialize layers (layer 1 = atmosphere, layers 2..nl_soil+1 = soil)
!#############################################################################
      zsoil_wil(1) = 1.0_r8
      nref(1) = (1.0_r8, 0.0_r8)

      DO i = 1, nl_soil
         zsoil_wil(i+1) = dz_soi(i)
         eps_layer = eps_prof(i)
         nref(i+1) = sqrt(eps_layer)
      END DO

!#############################################################################
! 2. Horizontal polarization: HPLD solution -> t_eff(1) and r_s(1)
!#############################################################################
      CALL wilheit_hpld(nl_soil + 1, nref, zsoil_wil, cp_h)

      sumtfa = 0.0_r8
      sumfa  = 0.0_r8
      DO i = 1, nl_soil
         fa = real(cp_h(i+1))
         sumtfa = sumtfa + fa * t_soi(i)
         sumfa  = sumfa  + fa
      END DO

      IF (sumfa > 0.0_r8) THEN
         t_eff(1) = sumtfa / sumfa
         r_s(1)   = max(0.0_r8, min(1.0_r8, 1.0_r8 - sumfa))
      ELSE
         t_eff(1) = t_soi(1)
         r_s(1)   = 0.0_r8
      END IF

!#############################################################################
! 3. Vertical polarization: VPLD solution -> t_eff(2) and r_s(2)
!    (for nadir view, V and H are degenerate)
!#############################################################################
      IF (theta /= 0.0_r8) THEN
         CALL wilheit_vpld(nl_soil + 1, nref, zsoil_wil, cp_v)
         sumtfa = 0.0_r8
         sumfa  = 0.0_r8
         DO i = 1, nl_soil
            fa = real(cp_v(i+1))
            sumtfa = sumtfa + fa * t_soi(i)
            sumfa  = sumfa  + fa
         END DO
         IF (sumfa > 0.0_r8) THEN
            t_eff(2) = sumtfa / sumfa
            r_s(2)   = max(0.0_r8, min(1.0_r8, 1.0_r8 - sumfa))
         ELSE
            t_eff(2) = t_soi(1)
            r_s(2)   = 0.0_r8
         END IF
      ELSE
         t_eff(2) = t_eff(1)
         r_s(2)   = r_s(1)
      END IF

   END SUBROUTINE soil_eff_temp_Wil78

!-----------------------------------------------------------------------
   SUBROUTINE wilheit_hpld(n, cn, del, cp)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the surface reflectivity and the per-layer absorption
!   fractions at horizontal polarization for a plane-stratified
!   dielectric, following Wilheit (1978).
!
!   NOTE: This is a CoLM port of the ECMWF CMEM routine HPLD in
!   wilheit_sub.F90; the structure and variable names follow the
!   original implementation.
!
! Input/Output Variables:
!   n      integer  I   number of layers (= number of soil layers + 1,
!                       layer 1 = atmosphere)
!   cn     complex  I   index of refraction for atmosphere (cn(1)) and
!                       each soil layer
!   del    real     I   layer thicknesses (m) (the first and last layer
!                       are assumed to be semi-infinite; their thickness
!                       value is ignored)
!   cp     complex  O   on exit, real(cp(1)) is the surface reflectivity
!                       and real(cp(2))..real(cp(n)) are the absorption
!                       fractions of each layer
!
! REFERENCES:
!   [1] Wilheit, T. T. (1978), Radiative transfer in a plane stratified
!       dielectric, IEEE Trans. Geosci. Electron., GE-16(2), 138-143.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      integer, intent(in)     :: n                 ! [-] number of layers (= nl_soil + 1, layer 1 = atmosphere)
      complex(r8), intent(in) :: cn(n)             ! [-] index of refraction for atmosphere (cn(1)) and each soil layer
      real(r8), intent(in)    :: del(n)            ! [m] layer thicknesses; first and last layers semi-infinite
      complex(r8), intent(out):: cp(n)             ! [-] real(cp(1))=surface reflectivity; real(cp(2:n))=layer absorption fractions

      integer     :: i                             ! [-] loop/index variable
      integer     :: j                             ! [-] loop/index variable
      integer     :: ii                            ! [-] loop/index variable
      integer     :: nl_w                          ! [-] loop/index variable
      integer     :: nmax                          ! [-] loop/index variable
      integer     :: ll                            ! [-] loop/index variable
      real(r8)    :: r                             ! [-] |cp(j)|
      real(r8)    :: s_val                         ! [-] |cp(j-1)|
      real(r8)    :: arg                           ! [rad] phase
      real(r8)    :: e2                            ! [-] field intermediate
      real(r8)    :: dp                            ! [-] absorption fraction
      real(r8)    :: x                             ! [-] correction
      complex(r8) :: cs                            ! [-] sin(theta)/cn(i)
      complex(r8) :: cc                            ! [-] cos(theta) in layer
      complex(r8) :: carg                          ! [-] phase
      complex(r8) :: csj                           ! [-] sin(theta)/cn(j)
      complex(r8) :: ccj                           ! [-] cos(theta) in layer j
      complex(r8) :: csjp1                         ! [-] sin(theta)/cn(j+1)
      complex(r8) :: ccjp1                         ! [-] cos(theta) in layer j+1
      complex(r8) :: ca                            ! [-] transmission coeff
      complex(r8) :: cb                            ! [-] reflection coeff
      complex(r8) :: cx                            ! [-] normalization
      complex(r8) :: cxp                           ! [-] E+*conjg(E-)
      complex(r8) :: cep(n)                        ! [-] up-welling electric field E+
      complex(r8) :: cem(n)                        ! [-] down-welling electric field E-

!--------------------------------------------------------------
! 1. Calculate the squares of the interface propagators.
!    cp(i) here is equal to p(j)**2 in section II of the paper.
!--------------------------------------------------------------
      s_val = sin(theta)
      cp(1) = (1.0_r8, 0.0_r8)
      nl_w = n - 1
      nmax = 1

      DO i = 2, nl_w
         nmax = i + 1
         cs = cn(1) * s_val / cn(i)
         cc = sqrt((1.0_r8, 0.0_r8) - cs*cs)   ! in the paper cc is cos(theta)(j)
         arg = del(i) * k
         carg = 2.0_r8 * arg * cn(i) * cc * (0.0_r8, 1.0_r8)
         cp(i) = exp(carg) * cp(i-1)
         IF (abs(cp(i)) < 0.0001_r8) EXIT
      END DO

!---------------------------------------------------------------------
! 2. Calculate the electric fields within each layer.
!    cep(j) and cem(j) here are equal to E+(j) and E-(j) in the paper.
!---------------------------------------------------------------------
      cep(nmax) = (1.0_r8, 0.0_r8)
      cem(nmax) = (0.0_r8, 0.0_r8)

      DO ii = 2, nmax
         j = nmax - ii + 1
         csj = cn(1) * s_val / cn(j)
         ccj = sqrt((1.0_r8, 0.0_r8) - csj*csj)
         csjp1 = cn(1) * s_val / cn(j+1)
         ccjp1 = sqrt((1.0_r8, 0.0_r8) - csjp1*csjp1)
         ca = 2.0_r8 * cn(j) * ccj / (cn(j)*ccj + cn(j+1)*ccjp1)
         cb = (cn(j)*ccj - cn(j+1)*ccjp1) / ((cn(j)*ccj + cn(j+1)*ccjp1) * cp(j))
         cep(j) = cep(j+1)/ca + cb*cem(j+1)/ca
         cem(j) = cem(j+1) + (cep(j+1)-cep(j))*cp(j)
      END DO

      ! normalize so that cep(1) = 1
      cx = cep(1)
      DO j = 1, nmax
         cep(j) = cep(j) / cx
         cem(j) = cem(j) / cx
      END DO

!----------------------------------------------------------------
! 3. Calculate the absorption fraction of each layer (f(i) in the paper).
!----------------------------------------------------------------
      DO j = nmax, n
         cp(j) = (0.0_r8, 1.0_r8) / (1.E15_r8)
      END DO

      ll = nmax - 1
      DO ii = 1, ll
         j = nmax - ii + 1
         cs = sin(theta) / cn(j)
         cc = sqrt((1.0_r8, 0.0_r8) - cs*cs)
         r = abs(cp(j))
         s_val = abs(cp(j-1))
         e2 = (s_val-r) * abs(cep(j))**2 + (1.0_r8/r-1.0_r8/s_val) * abs(cem(j))**2
         dp = e2 * real(cn(j)*cc) / cos(theta)
         cxp = cep(j) * conjg(cem(j))
         x = 2.0_r8 * aimag(cn(j)*cc/cos(theta)) * &
             (aimag(cxp*cp(j-1)/abs(cp(j-1))) - aimag(cxp*cp(j)/abs(cp(j))))
         dp = dp - x
         cp(j) = cmplx(dp, 0.0_r8, kind=r8)
      END DO

      ! surface reflectivity is the squared magnitude of E-(1) times Re(n(1))
      r = abs(cem(1))**2 * real(cn(1))
      cp(1) = cmplx(r, 0.0_r8, kind=r8)

   END SUBROUTINE wilheit_hpld

!-----------------------------------------------------------------------
   SUBROUTINE wilheit_vpld(n, cn, del, cp)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the surface reflectivity and the per-layer absorption
!   fractions at vertical polarization for a plane-stratified
!   dielectric, following Wilheit (1978).
!
!   NOTE: This is a CoLM port of the ECMWF CMEM routine VPLD in
!   wilheit_sub.F90; the structure and variable names follow the
!   original implementation.
!
! Input/Output Variables:
!   n      integer  I   number of layers (= number of soil layers + 1,
!                       layer 1 = atmosphere)
!   cn     complex  I   index of refraction for atmosphere (cn(1)) and
!                       each soil layer
!   del    real     I   layer thicknesses (m) (the first and last layer
!                       are assumed to be semi-infinite; their thickness
!                       value is ignored)
!   cp     complex  O   on exit, real(cp(1)) is the surface reflectivity
!                       and real(cp(2))..real(cp(n)) are the absorption
!                       fractions of each layer
!
! REFERENCES:
!   [1] Wilheit, T. T. (1978), Radiative transfer in a plane stratified
!       dielectric, IEEE Trans. Geosci. Electron., GE-16(2), 138-143.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      integer, intent(in)     :: n                 ! [-] number of layers (= nl_soil + 1, layer 1 = atmosphere)
      complex(r8), intent(in) :: cn(n)             ! [-] index of refraction for atmosphere (cn(1)) and each soil layer
      real(r8), intent(in)    :: del(n)            ! [m] layer thicknesses; first and last layers semi-infinite
      complex(r8), intent(out):: cp(n)             ! [-] real(cp(1))=surface reflectivity; real(cp(2:n))=layer absorption fractions

      integer     :: i                             ! [-] loop/index variable
      integer     :: j                             ! [-] loop/index variable
      integer     :: ii                            ! [-] loop/index variable
      integer     :: nl_w                          ! [-] loop/index variable
      integer     :: nmax                          ! [-] loop/index variable
      integer     :: ll                            ! [-] loop/index variable
      real(r8)    :: r                             ! [-] |cp(j)|
      real(r8)    :: s_val                         ! [-] |cp(j-1)|
      real(r8)    :: arg                           ! [rad] phase
      real(r8)    :: e2                            ! [-] field intermediate
      real(r8)    :: dp                            ! [-] absorption fraction
      real(r8)    :: x                             ! [-] correction
      complex(r8) :: cs                            ! [-] sin(theta)/cn(i)
      complex(r8) :: cc                            ! [-] cos(theta) in layer
      complex(r8) :: carg                          ! [-] phase
      complex(r8) :: csj                           ! [-] sin(theta)/cn(j)
      complex(r8) :: ccj                           ! [-] cos(theta) in layer j
      complex(r8) :: csjp1                         ! [-] sin(theta)/cn(j+1)
      complex(r8) :: ccjp1                         ! [-] cos(theta) in layer j+1
      complex(r8) :: ca                            ! [-] numerator
      complex(r8) :: cb                            ! [-] numerator
      complex(r8) :: cd                            ! [-] denominator
      complex(r8) :: cr_c                          ! [-] refractive index ratio
      complex(r8) :: cx                            ! [-] normalization
      complex(r8) :: cxp                           ! [-] E+*conjg(E-)
      complex(r8) :: cep(n)                        ! [-] up-welling electric field E+
      complex(r8) :: cem(n)                        ! [-] down-welling electric field E-

!--------------------------------------------------------------
! 1. Calculate the squares of the interface propagators.
!    cp(i) here is equal to p(j)**2 in section II of the paper.
!--------------------------------------------------------------
      s_val = sin(theta)
      cp(1) = (1.0_r8, 0.0_r8)
      nl_w = n - 1
      nmax = 1

      DO i = 2, nl_w
         nmax = i + 1
         cs = cn(1) * s_val / cn(i)
         cc = sqrt((1.0_r8, 0.0_r8) - cs*cs)   ! in the paper cc is cos(theta)(j)
         arg = del(i) * k
         carg = 2.0_r8 * arg * cn(i) * cc * (0.0_r8, 1.0_r8)
         cp(i) = exp(carg) * cp(i-1)
         IF (abs(cp(i)) < 0.0001_r8) EXIT
      END DO

!---------------------------------------------------------------------
! 2. Calculate the electric fields within each layer.
!    cep(j) and cem(j) here are equal to E+(j) and E-(j) in the paper.
!---------------------------------------------------------------------
      cep(nmax) = (1.0_r8, 0.0_r8)
      cem(nmax) = (0.0_r8, 0.0_r8)

      DO ii = 2, nmax
         j = nmax - ii + 1
         csj = cn(1) * s_val / cn(j)
         ccj = sqrt((1.0_r8, 0.0_r8) - csj*csj)
         csjp1 = cn(1) * s_val / cn(j+1)
         ccjp1 = sqrt((1.0_r8, 0.0_r8) - csjp1*csjp1)
         cd = 2.0_r8 * cn(j) * ccj
         ca = cn(j)*ccjp1 + cn(j+1)*ccj
         cb = cn(j)*ccjp1 - cn(j+1)*ccj
         cep(j) = ca*cep(j+1)/cd + cb*cem(j+1)/(cd*cp(j))
         cr_c = cn(j+1) / cn(j)
         cem(j) = cr_c*cem(j+1) + (cep(j)-cep(j+1)*cr_c)*cp(j)
      END DO

      ! normalize so that cep(1) = 1
      cx = cep(1)
      DO j = 1, nmax
         cep(j) = cep(j) / cx
         cem(j) = cem(j) / cx
      END DO

!----------------------------------------------------------------
! 3. Calculate the absorption fraction of each layer (f(i) in the paper).
!----------------------------------------------------------------
      DO j = nmax, n
         cp(j) = (0.0_r8, 1.0_r8) / (1.E15_r8)
      END DO

      ll = nmax - 1
      DO ii = 1, ll
         j = nmax - ii + 1
         cs = sin(theta) / cn(j)
         cc = sqrt((1.0_r8, 0.0_r8) - cs*cs)
         r = abs(cp(j))
         s_val = abs(cp(j-1))
         e2 = (s_val-r) * abs(cep(j))**2 + (1.0_r8/r-1.0_r8/s_val) * abs(cem(j))**2
         dp = e2 * real(cn(j)*cc) / cos(theta)
         cxp = cep(j) * conjg(cem(j))
         x = 2.0_r8 * aimag(cn(j)*cc/cos(theta)) * &
             (aimag(cxp*cp(j-1)/abs(cp(j-1))) - aimag(cxp*cp(j)/abs(cp(j))))
         dp = dp - x
         cp(j) = cmplx(dp, 0.0_r8, kind=r8)
      END DO

      ! surface reflectivity is the squared magnitude of E-(1) times Re(n(1))
      r = abs(cem(1))**2 * real(cn(1))
      cp(1) = cmplx(r, 0.0_r8, kind=r8)

   END SUBROUTINE wilheit_vpld

!-----------------------------------------------------------------------
   SUBROUTINE soil_eff_temp_Wig01(wc_surf, t_surf, t_deep, t_eff)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil with the Wigneron (2001)
!   parameterization. The surface-to-deep weighting is controlled by the
!   surface volumetric water content through a two-parameter power law.
!
! REFERENCES:
!   [1] Wigneron, J.-P., Laguerre, L., Kerr, Y. H. (2001), A simple
!       parameterization of the L-band microwave emission from rough
!       agricultural soils, IEEE Trans. Geosci. Remote Sens., 39, 1697-1707.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      real(r8), intent(in)  :: wc_surf        ! [m3/m3] surface volumetric water content
      real(r8), intent(in)  :: t_surf         ! [K] surface soil temperature
      real(r8), intent(in)  :: t_deep         ! [K] deep soil temperature
      real(r8), intent(out) :: t_eff(2)       ! [K] effective temperature for H and V, (1)=H,(2)=V

      real(r8) :: C                           ! [-] weighting parameter depending on soil moisture
      real(r8) :: w0 = 0.30_r8                ! [m3/m3] reference soil moisture
      real(r8) :: bw = 0.30_r8                ! [-] exponent

      IF (wc_surf < 0.0_r8) THEN
        C = 0.001_r8
      ELSE
        C = max(0.001_r8, (wc_surf/w0)**bw)
      ENDIF
      t_eff(:) = t_deep + (t_surf - t_deep)*C

   END SUBROUTINE soil_eff_temp_Wig01

!-----------------------------------------------------------------------
   SUBROUTINE soil_eff_temp_Lv14(nl_soil, dz_soi, t_soi, eps_prof, t_eff)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil with the Lv (2014)
!   discrete multi-layer scheme. Each layer's contribution is weighted
!   by its one-way optical depth, with a cumulative attenuation product
!   accounting for the layers above.
!
! REFERENCES:
!   [1] Lv, S., Wen, J., Zeng, Y., Tian, H., & Su, Z. (2014). An improved
!       two-layer algorithm for estimating effective soil temperature in
!       microwave radiometry. Remote Sensing of Environment, 152, 356-363.
!       (based on Eq. 14 for the multi-layer discrete scheme)
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)     :: nl_soil         ! [-] number of soil layers
      real(r8), intent(in)    :: dz_soi(:)       ! [m] layer thickness
      real(r8), intent(in)    :: t_soi(:)        ! [K] soil temperature
      complex(r8), intent(in) :: eps_prof(:)     ! [-] per-layer complex dielectric constant
      real(r8), intent(out)   :: t_eff(2)        ! [K] effective temperature for H and V, (1)=H,(2)=V

      integer     :: i                           ! [-] layer index
      real(r8)    :: B_i                         ! [-] optical thickness
      real(r8)    :: eps_r                       ! [-] real part of dielectric constant
      real(r8)    :: eps_i                       ! [-] absolute imaginary part
      real(r8)    :: prod_term                   ! [-] cumulative attenuation
      complex(r8) :: eps                         ! [-] complex dielectric constant

      t_eff(1)  = 0.0_r8
      prod_term = 1.0_r8

      DO i = 1, nl_soil
         ! 1. get dielectric constant for the current layer
         eps = eps_prof(i)
         eps_r = real(eps)
         eps_i = abs(aimag(eps))

         ! 2. calculate B_i (optical thickness) for the current layer
         B_i = dz_soi(i) * (4.0_r8 * pi / lam) * (eps_i / (2.0_r8 * sqrt(eps_r)))

         ! 3. accumulate effective temperature (Lv et al. 2014, Eq. 14)
         IF (i < nl_soil) THEN
            t_eff(1) = t_eff(1) + t_soi(i) * (1.0_r8 - exp(-B_i)) * prod_term
            prod_term = prod_term * exp(-B_i)
         ELSE
            ! for the bottom layer, integrate to infinity
            t_eff(1) = t_eff(1) + t_soi(i) * prod_term
         END IF
      END DO

      ! assume identical Teff for H and V polarizations
      t_eff(2) = t_eff(1)

   END SUBROUTINE soil_eff_temp_Lv14

!-----------------------------------------------------------------------
   SUBROUTINE diel_ice(t, eps_i)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of pure ice.
!
! REFERENCES:
!   [1] Matzler, C. (2006). Thermal Microwave Radiation: Applications
!       for Remote Sensing, pp. 456-461.
!-----------------------------------------------------------------------
      USE MOD_Const_Physical
      USE MOD_Precision
      IMPLICIT NONE
      real(r8), intent(in)     :: t              ! [K] ice temperature
      complex(r8), intent(out) :: eps_i          ! [-] dielectric constant of ice

      real(r8) :: betam                          ! [-] beta parameter from Mishima
      real(r8) :: dbeta                          ! [-] delta beta
      real(r8) :: beta                           ! [-] beta parameter
      real(r8) :: t_inv                          ! [-] inverse temperature term
      real(r8) :: tk                             ! [K] temperature
      real(r8) :: t_C                            ! [C] temperature
      real(r8) :: alpha                          ! [-] alpha parameter
      real(r8) :: eps_i_r                        ! [-] real part
      real(r8) :: eps_i_i                        ! [-] imaginary part

      tk = t
      t_C = t - tfrz

      ! eq.(5.33): calculate beta parameter by Mishima et al. (1983)
      betam = (0.0207_r8/tk)*(exp(335.0_r8/tk)/((exp(335.0_r8/tk) - 1.0_r8)**2.0_r8)) + 1.16e-11_r8*(fghz**2.0_r8)
      ! eq.(5.35): calculate delta beta parameter
      dbeta = exp(-10.02_r8 + 0.0364_r8*t_C)
      ! eq.(5.34): calculate beta parameter
      beta  = betam + dbeta
      ! eq.(5.32): calculate alpha parameter
      t_inv = 300.0_r8/tk - 1.0_r8
      alpha = (0.00504_r8 + 0.0062_r8*t_inv)*exp(-22.1_r8*t_inv)

      ! eq.(5.30): real part
      eps_i_r = 3.1884_r8 + 9.1e-4_r8*t_C
      ! eq.(5.31): imaginary part
      eps_i_i = alpha/fghz + beta*fghz
      ! assemble complex dielectric constant
      eps_i   = eps_i_r - jj*eps_i_i

   END SUBROUTINE diel_ice

!-----------------------------------------------------------------------
   SUBROUTINE diel_water(type, swc, t, wf_sand, wf_clay, BD_all, sal, eps_w)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of water: pure water, sea water,
!   or soil water (including the soil-specific ionic conductivity).
!
! REFERENCES:
!   [1] Ulaby, F. T., Moore, R. K., & Fung, A. K. (1986), Microwave Remote
!       Sensing: Active and Passive, Vol. III, Artech House.
!   [2] Klein, L. A., & Swift, C. T. (1977), An improved model for the
!       dielectric constant of sea water at microwave frequencies, IEEE
!       Trans. Antennas Propag., AP-25(1), 104-111.
!-----------------------------------------------------------------------
      USE MOD_Const_Physical
      USE MOD_Precision
      IMPLICIT NONE
      integer, intent(in)      :: type    ! [-] 0: pure water, 1: sea water, 2: soil water
      real(r8), intent(in)     :: swc     ! [m3/m3] soil water content
      real(r8), intent(in)     :: t       ! [K] soil temperature
      real(r8), intent(in)     :: wf_sand ! [-] gravimetric sand fraction
      real(r8), intent(in)     :: wf_clay ! [-] gravimetric clay fraction
      real(r8), intent(in)     :: BD_all  ! [kg/m3] bulk density
      real(r8), intent(in)     :: sal     ! [psu] water salinity
      complex(r8), intent(out) :: eps_w   ! [-] dielectric constant of water

      real(r8) :: sigma                   ! [S/m] ionic conductivity
      real(r8) :: a                       ! [-] salinity correction for static permittivity
      real(r8) :: b                       ! [-] salinity correction for relaxation time
      real(r8) :: tau_w                   ! [s] relaxation time
      real(r8) :: eps_w0                  ! [-] static dielectric constant
      real(r8) :: wc                      ! [m3/m3] effective soil water content
      real(r8) :: t_C                     ! [C] temperature
      real(r8) :: BD_gcm3                 ! [g/cm3] bulk density

      t_C = t - tfrz
      BD_gcm3 = BD_all / 1000.0_r8

      ! [3] eq.16: tau(T, sal) = tau_w(T) * b(sal, T)
      ! calculate relaxation time of pure water (Stogryn)
      tau_w = 1.768e-11_r8 - 6.068e-13_r8*t_C + 1.104e-14_r8*t_C**2 - 8.111e-17_r8*t_C**3
      b     = 1.000_r8 + 2.282e-5_r8*sal*t_C - 7.638e-4_r8*sal - 7.760e-6_r8*sal**2 + 1.105e-8_r8*sal**3
      tau_w = tau_w*b

      ! [3] eq.13: eps_w0(sal, T) = eps_w0(T) * a(sal, T)
      ! static dielectric constant of pure water (Klein and Swift)
      eps_w0 = 87.134_r8 - 1.949e-1_r8*t_C - 1.276e-2_r8*t_C**2 + 2.491e-4_r8*t_C**3
      a      = 1.000_r8 + 1.613e-5_r8*sal*t_C - 3.656e-3_r8*sal + 3.210e-5_r8*sal**2 - 4.232e-7_r8*sal**3
      eps_w0 = eps_w0*a

      IF (type == 0) THEN  ! pure water
         ! [1] eq.19
         eps_w0 = 88.045_r8 - 0.4147_r8*t_C + 6.295e-4_r8*t_C**2 + 1.075e-5_r8*t_C**3
         eps_w  = eps_w_inf + (eps_w0 - eps_w_inf)/(1.0_r8 - jj*omega*tau_w)
      ELSEIF (type == 1) THEN  ! sea water
         ! calculate ionic conductivity [1] eq.27, eq.28
         sigma = sal*(0.182521_r8 - 1.46192e-3_r8*sal + 2.09324e-5_r8*sal**2 - 1.28205e-7_r8*sal**3) &
            *exp(-1.0_r8*(25.0_r8 - t_C)* &
            (2.033e-2_r8 + 1.266e-4_r8*(25.0_r8 - t_C) + 2.464e-6_r8*(25.0_r8 - t_C)**2 &
            - sal*(1.849e-5_r8 - 2.551e-7_r8*(25.0_r8 - t_C) + 2.551e-8_r8*(25.0_r8 - t_C)**2)))

         ! dielectric constant of sea water [1] eq.21
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1.0_r8 - jj*omega*tau_w) + jj*sigma/(omega*eps_0)
      ELSEIF (type == 2) THEN  ! soil water
         ! calculate soil conductivity
         sigma = -1.645_r8 + 1.939_r8*BD_gcm3 - 2.256_r8*wf_sand + 1.594_r8*wf_clay
         IF (sigma < 0.0_r8) sigma = 0.0_r8  ! negative for very sandy soils with low bulk density

         ! calculate dielectric constant of soil-water by modified Debye expression
         wc = max(0.001_r8, swc)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1.0_r8 - jj*omega*tau_w) &
            + jj*sigma/(omega*eps_0)*(rho_soil - BD_gcm3)/(rho_soil*wc)
      END IF

   END SUBROUTINE diel_water

!-----------------------------------------------------------------------
   SUBROUTINE diel_soil_W80(ew, t, wc, wf_sand, wf_clay, porsl, eps)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil with the
!   Wang & Schmugge (1980) empirical model.
!
! REFERENCES:
!   [1] Matzler, C. (1998). Microwave permittivity of dry sand, IEEE Trans.
!       Geosci. Remote Sens., 36(1), 317-319.
!   [2] Wang, J. R., & Schmugge, T. J. (1980). An empirical model for the
!       complex dielectric permittivity of soils as a function of water
!       content, IEEE Trans. Geosci. Remote Sens., GE-18(4), 288-295.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)     :: t                 ! [K] soil temperature
      real(r8), intent(in)     :: wc                ! [m3/m3] volumetric soil moisture
      real(r8), intent(in)     :: wf_sand           ! [-] gravimetric sand fraction
      real(r8), intent(in)     :: wf_clay           ! [-] gravimetric clay fraction
      real(r8), intent(in)     :: porsl             ! [m3/m3] porosity
      complex(r8), intent(in)  :: ew                ! [-] dielectric constant of water
      complex(r8), intent(out) :: eps               ! [-] dielectric constant of soil

      real(r8) :: wp                                ! [m3/m3] wilting point
      real(r8) :: wt                                ! [m3/m3] transition moisture
      real(r8) :: gamma                             ! [-] fitting parameter
      real(r8) :: ecl                               ! [-] conductivity loss
      real(r8) :: alpha                             ! [-] conductivity loss parameter
      complex(r8) :: eps_x                          ! [-] dielectric constant of bound water
      complex(r8) :: eps_a = (1.0, 0.0)             ! [-] dielectric constant of air
      complex(r8) :: eps_r = (5.5, 0.2)             ! [-] dielectric constant of rock/mineral
      complex(r8) :: eps_i = (3.2, 0.1)             ! [-] dielectric constant of ice-like bound water

      ! calculate wilting point at the soil layer
      wp    = 0.06774_r8 - 0.064_r8*wf_sand + 0.478_r8*wf_clay          ! [2](1)
      ! calculate fitting parameters
      gamma = -0.57_r8*wp + 0.481_r8                                        ! [2](8)
      ! calculate transition moisture point
      wt    = 0.49_r8*wp + 0.165_r8                                         ! [2](9)

      ! calculate dielectric constant of wet soil (when all soil freeze, eps_x = eps_i)
      IF (wc <= wt) THEN
        eps_x = eps_i + (ew - eps_i)*(wc/wt)*gamma                           ! [2](3)
        eps   = wc*eps_x + (porsl - wc)*eps_a + (1.0_r8 - porsl)*eps_r       ! [2](2)
      ELSE
        eps_x = eps_i + (ew - eps_i)*gamma                                   ! [2](5)
        eps   = wt*eps_x + (wc - wt)*ew + (porsl - wc)*eps_a + (1.0_r8 - porsl)*eps_r  ! [2](4)
      END IF

      ! add conductivity loss for imaginary part
      alpha = min(100.0_r8*wp, 26.0_r8)
      ecl   = alpha * wc**2                                                  ! [2](6)
      eps   = eps + jj*ecl                                                   ! [2](6)

   END SUBROUTINE diel_soil_W80

!-----------------------------------------------------------------------
   SUBROUTINE diel_soil_M04(wc, wf_clay, eps)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil with the Mironov
!   generalized refractive mixing dielectric model (GRMDM, 2004).
!
! REFERENCES:
!   [1] Mironov, V. L., et al. (2004), Generalized refractive mixing
!       dielectric model for moist soils, IEEE Trans. Geosci. Remote Sens.,
!       42(4), 773-785.
!   [2] Mironov, V. L., et al. (2009), Physically and mineralogically
!       based spectroscopic dielectric model for moist soils, IEEE Trans.
!       Geosci. Remote Sens., 47(7), 2059-2070.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)     :: wc                ! [m3/m3] volumetric soil moisture
      real(r8), intent(in)     :: wf_clay           ! [-] gravimetric clay fraction
      complex(r8), intent(out) :: eps               ! [-] dielectric constant of soil

      real(r8) :: znd                               ! [-] dry-soil refractive index
      real(r8) :: zkd                               ! [-] dry-soil normalized attenuation coefficient
      real(r8) :: zxmvt                             ! [m3/m3] maximum bound-water fraction
      real(r8) :: zep0b                             ! [-] bound-water static dielectric constant
      real(r8) :: ztaub                             ! [s] bound-water relaxation time
      real(r8) :: zsigmab                           ! [S/m] bound-water conductivity
      real(r8) :: zep0u                             ! [-] free-water static dielectric constant
      real(r8) :: ztauu                             ! [s] free-water relaxation time
      real(r8) :: zsigmau                           ! [S/m] free-water conductivity
      real(r8) :: zcxb                              ! [-] bound-water Debye intermediate
      real(r8) :: zepwbx                            ! [-] bound-water real permittivity
      real(r8) :: zepwby                            ! [-] bound-water loss permittivity
      real(r8) :: zcxu                              ! [-] free-water Debye intermediate
      real(r8) :: zepwux                            ! [-] free-water real permittivity
      real(r8) :: zepwuy                            ! [-] free-water loss permittivity
      real(r8) :: znb                               ! [-] bound-water refractive index
      real(r8) :: zkb                               ! [-] bound-water attenuation coefficient
      real(r8) :: znu                               ! [-] free-water refractive index
      real(r8) :: zku                               ! [-] free-water attenuation coefficient
      real(r8) :: zxmvt2                            ! [m3/m3] min(wc,zxmvt)
      real(r8) :: znm                               ! [-] soil refractive index
      real(r8) :: zkm                               ! [-] soil attenuation coefficient
      real(r8) :: zepmx                             ! [-] soil real permittivity
      real(r8) :: zepmy                             ! [-] soil loss permittivity
      integer  :: zflag                             ! [-] flag: 1 if wc >= zxmvt

      ! refractive index and normalized attenuation coefficient of dry soils
      znd   = 1.634_r8 - 0.539_r8 * wf_clay + 0.2748_r8 * wf_clay**2
      zkd   = 0.03952_r8 - 0.04038_r8 * wf_clay
      ! maximum bound water fraction
      zxmvt = 0.02863_r8 + 0.30673_r8 * wf_clay

      ! bound water parameters
      zep0b   = 79.8_r8 - 85.4_r8 * wf_clay + 32.7_r8 * wf_clay**2
      ztaub   = 1.062e-11_r8 + 3.450e-12_r8 * wf_clay
      zsigmab = 0.3112_r8 + 0.467_r8 * wf_clay

      ! unbound (free) water parameters
      zep0u   = 100.0_r8
      ztauu   = 8.5e-12_r8
      zsigmau = 0.3631_r8 + 1.217_r8 * wf_clay

      ! calculate water dielectric constant (bound and free)
      zcxb   = (zep0b - eps_w_inf) / (1.0_r8 + (2.0_r8*pi*f*ztaub)**2)
      zepwbx = eps_w_inf + zcxb
      zepwby = zcxb * (2.0_r8*pi*f*ztaub) + zsigmab / (2.0_r8*pi*eps_0*f)
      zcxu   = (zep0u - eps_w_inf) / (1.0_r8 + (2.0_r8*pi*f*ztauu)**2)
      zepwux = eps_w_inf + zcxu
      zepwuy = zcxu * (2.0_r8*pi*f*ztauu) + zsigmau / (2.0_r8*pi*eps_0*f)

      ! calculate refractive index of water (bound and free)
      znb = sqrt( sqrt( zepwbx**2 + zepwby**2) + zepwbx ) / sqrt(2.0_r8)
      zkb = sqrt( sqrt( zepwbx**2 + zepwby**2) - zepwbx ) / sqrt(2.0_r8)
      znu = sqrt( sqrt( zepwux**2 + zepwuy**2) + zepwux ) / sqrt(2.0_r8)
      zku = sqrt( sqrt( zepwux**2 + zepwuy**2) - zepwux ) / sqrt(2.0_r8)

      ! calculate soil refractive index (nm and km)
      zxmvt2 = min(wc, zxmvt)
      zflag  = 0
      IF (wc >= zxmvt) zflag = 1
      znm = znd + (znb - 1.0_r8) * zxmvt2 + (znu - 1.0_r8) * (wc - zxmvt) * zflag
      zkm = zkd + zkb * zxmvt2 + zku * (wc - zxmvt) * zflag

      ! calculate soil dielectric constant
      zepmx = znm**2 - zkm**2
      zepmy = znm * zkm * 2.0_r8
      eps   = cmplx(zepmx, zepmy, kind=r8)

   END SUBROUTINE diel_soil_M04

!-----------------------------------------------------------------------
   SUBROUTINE diel_soil_M09(wc, t, wf_clay, eps)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil with the
!   temperature- and mineralogy-dependent Mironov model (2009).
!
! REFERENCES:
!   [1] Mironov, V. L., & Fomin, S. V. (2009), Temperature and mineralogy
!       dependable model for microwave dielectric spectra of moist soils,
!       PIERS Online, 5(5), 411-415.
!   [2] Mironov, V. L., et al. (2009), Physically and mineralogically
!       based spectroscopic dielectric model for moist soils, IEEE Trans.
!       Geosci. Remote Sens., 47(7), 2059-2070.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)     :: t                 ! [K] soil temperature
      real(r8), intent(in)     :: wc                ! [m3/m3] volumetric soil moisture
      real(r8), intent(in)     :: wf_clay           ! [-] gravimetric clay fraction
      complex(r8), intent(out) :: eps               ! [-] dielectric constant of soil

      real(r8) :: nd                                ! [-] dry-soil refractive index
      real(r8) :: kd                                ! [-] dry-soil attenuation coefficient
      real(r8) :: mvt                               ! [m3/m3] maximum bound-water fraction
      real(r8) :: ts                                ! [C] reference temperature
      real(r8) :: e0b                               ! [-] bound-water static permittivity
      real(r8) :: Bb                                ! [-] bound-water temperature slope
      real(r8) :: Bsgb                              ! [-] bound-water conductivity slope
      real(r8) :: Fb                                ! [-] log factor
      real(r8) :: eb0                               ! [-] bound-water static permittivity at T
      real(r8) :: dHbR                              ! [-] bound-water enthalpy
      real(r8) :: dSbR                              ! [-] bound-water entropy
      real(r8) :: taub                              ! [s] bound-water relaxation time
      real(r8) :: sigmabt                           ! [S/m] bound-water conductivity at Ts
      real(r8) :: sigmab                            ! [S/m] bound-water conductivity at T
      real(r8) :: e0u                               ! [-] free-water static permittivity
      real(r8) :: Bu                                ! [-] free-water temperature slope
      real(r8) :: Bsgu                              ! [-] free-water conductivity slope
      real(r8) :: Fu                                ! [-] log factor
      real(r8) :: eu0                               ! [-] free-water static permittivity at T
      real(r8) :: dHuR                              ! [-] free-water enthalpy
      real(r8) :: dSuR                              ! [-] free-water entropy
      real(r8) :: tauu                              ! [s] free-water relaxation time
      real(r8) :: sigmau                            ! [S/m] free-water conductivity at T
      real(r8) :: sigmaut                           ! [S/m] free-water conductivity at Ts
      real(r8) :: cxb                               ! [-] bound-water Debye intermediate
      real(r8) :: eb_r                              ! [-] bound-water real permittivity
      real(r8) :: eb_i                              ! [-] bound-water loss permittivity
      real(r8) :: cxu                               ! [-] free-water Debye intermediate
      real(r8) :: eu_r                              ! [-] free-water real permittivity
      real(r8) :: eu_i                              ! [-] free-water loss permittivity
      real(r8) :: nb                                ! [-] bound-water refractive index
      real(r8) :: kb                                ! [-] bound-water attenuation coefficient
      real(r8) :: nu                                ! [-] free-water refractive index
      real(r8) :: ku                                ! [-] free-water attenuation coefficient
      real(r8) :: nm                                ! [-] soil refractive index
      real(r8) :: km                                ! [-] soil attenuation coefficient
      real(r8) :: eps_r                             ! [-] soil real permittivity
      real(r8) :: eps_i                             ! [-] soil loss permittivity
      real(r8) :: t_C                               ! [C] soil temperature

      t_C = t - tfrz

      ! refractive index and normalized attenuation coefficient of dry soils
      nd  = 1.634_r8 - 0.539_r8 * wf_clay + 0.2748_r8 * wf_clay**2
      kd  = 0.03952_r8 - 0.04038_r8 * wf_clay
      ! maximum bound water fraction
      mvt = 0.02863_r8 + 0.30673_r8 * wf_clay
      ! reference temperature for the parameter fit
      ts  = 20.0_r8

      ! bound water dielectric parameters and their temperature correction
      e0b  = 79.8_r8 - 85.4_r8 * wf_clay + 32.7_r8 * wf_clay**2
      Bb   = 8.67e-19_r8 - 1.26e-3_r8 * wf_clay + 1.84e-3_r8 * wf_clay**2 - 9.77e-4_r8*wf_clay**3 - 1.39e-7_r8*wf_clay**4
      Bsgb = 0.0028_r8 + 2.094e-2_r8*wf_clay - 1.229e-2_r8*wf_clay**2 - 5.03e-16_r8*wf_clay**3 + 4.163e-16_r8*wf_clay**4
      Fb   = log((e0b - 1.0_r8)/(e0b + 2.0_r8))
      eb0  = (1.0_r8 + 2.0_r8*exp(Fb-Bb*(t_C-ts))) / (1.0_r8 - exp(Fb-Bb*(t_C-ts)))

      ! bound water relaxation time and conductivity
      dHbR = 1467.0_r8 + 2697.0_r8*wf_clay - 980.0_r8*wf_clay**2 + 1.368e-4_r8*wf_clay**3 - 8.61e-5_r8*wf_clay**4
      dSbR = 0.888_r8 + 9.7_r8*wf_clay - 4.262_r8*wf_clay**2 + 6.79e-15_r8*wf_clay**3 + 4.263e-14_r8*wf_clay**4
      taub = 48e-12_r8 * exp(dHbR/t-dSbR)/t

      sigmabt = 0.3112_r8 + 0.467_r8*wf_clay
      sigmab  = sigmabt + Bsgb*(t_C-ts)

      ! unbound (free) water dielectric parameters and their temperature correction
      e0u  = 100.0_r8
      Bu   = 1.11e-4_r8 - 1.603e-5_r8*wf_clay + 1.239e-5_r8*wf_clay**2 + 8.33e-7_r8*wf_clay**3 - 1.007e-6_r8*wf_clay**4
      Bsgu = 0.00108_r8 + 0.1413_r8*wf_clay - 0.2555_r8*wf_clay**2 + 0.2147_r8*wf_clay**3 - 0.0711_r8*wf_clay**4
      Fu   = log((e0u - 1.0_r8)/(e0u + 2.0_r8))
      eu0  = (1.0_r8 + 2.0_r8*exp(Fu-Bu*(t_C-ts))) / (1.0_r8 - exp(Fu-Bu*(t_C-ts)))

      ! free water relaxation time and conductivity
      dHuR = 2231.0_r8 - 143.1_r8*wf_clay + 223.2_r8*wf_clay**2 - 142.1_r8*wf_clay**3 + 27.14_r8*wf_clay**4
      dSuR = 3.649_r8 - 0.4894_r8*wf_clay + 0.763_r8*wf_clay**2 - 0.4859_r8*wf_clay**3 + 0.0928_r8*wf_clay**4
      tauu = 48e-12_r8 * exp(dHuR/t-dSuR)/t

      sigmaut = 0.05_r8 + 1.4_r8*(1.0_r8 - (1.0_r8 - wf_clay)**4.664_r8)
      sigmau  = sigmaut + Bsgu*(t_C-ts)

      ! calculate dielectric constant of bound and free water
      cxb  = (eb0 - eps_w_inf) / (1.0_r8 + (2.0_r8*pi*f*taub)**2)
      eb_r = eps_w_inf + cxb
      eb_i = cxb*(2.0_r8*pi*f*taub) + sigmab/(2.0_r8*pi*eps_0*f)
      cxu  = (eu0 - eps_w_inf) / (1.0_r8 + (2.0_r8*pi*f*tauu)**2)
      eu_r = eps_w_inf + cxu
      eu_i = cxu*(2.0_r8*pi*f*tauu) + sigmau/(2.0_r8*pi*eps_0*f)

      ! calculate refractive index of bound and free water
      nb = sqrt(sqrt(eb_r**2 + eb_i**2) + eb_r) / sqrt(2.0_r8)
      kb = sqrt(sqrt(eb_r**2 + eb_i**2) - eb_r) / sqrt(2.0_r8)
      nu = sqrt(sqrt(eu_r**2 + eu_i**2) + eu_r) / sqrt(2.0_r8)
      ku = sqrt(sqrt(eu_r**2 + eu_i**2) - eu_r) / sqrt(2.0_r8)

      ! calculate soil refractive index (nm and km)
      IF (wc <= mvt) THEN
         nm = nd + (nb - 1.0_r8)*wc
         km = kd + kb*wc
      ELSE
         nm = nd + (nb - 1.0_r8)*mvt + (nu - 1.0_r8)*(wc - mvt)
         km = kd + kb*mvt + ku*(wc - mvt)
      ENDIF

      ! calculate soil dielectric constant
      eps_r = nm**2 - km**2
      eps_i = 2.0_r8 * nm * km
      eps   = cmplx(eps_r, eps_i, kind=r8)

   END SUBROUTINE diel_soil_M09

!-----------------------------------------------------------------------
   SUBROUTINE diel_soil_D85(ew, swc, ice_c, wf_sand, wf_clay, BD_all, eps)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet and frozen soil using
!   the four-phase (liquid water + ice + mineral + air) model of
!   Zhang & Zhao (2010), which extends Dobson (1985) to frozen soils.
!
! REFERENCES:
!   [1] Dobson, M. C., et al. (1985), Microwave dielectric behavior of wet
!       soil - Part II: Dielectric mixing models, IEEE Trans. Geosci.
!       Remote Sens., GE-23(1), 35-46.
!   [2] Zhang, L., Zhao, T., Jiang, L., & Zhao, S. (2010), Estimate of phase
!       transition water content in freeze-thaw process using microwave
!       radiometer, IEEE Trans. Geosci. Remote Sens., 48(12), 4248-4258.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      complex(r8), intent(in)  :: ew                ! [-] dielectric constant of liquid water
      real(r8), intent(in)     :: swc               ! [m3/m3] liquid soil moisture
      real(r8), intent(in)     :: ice_c             ! [m3/m3] ice volumetric content
      real(r8), intent(in)     :: wf_sand           ! [-] gravimetric sand fraction
      real(r8), intent(in)     :: wf_clay           ! [-] gravimetric clay fraction
      real(r8), intent(in)     :: BD_all            ! [kg/m3] soil bulk density
      complex(r8), intent(out) :: eps               ! [-] dielectric constant of soil

      real(r8) :: alphas = 0.65_r8                  ! [-] shape factor
      real(r8) :: beta                              ! [-] empirical exponent
      real(r8) :: eaa                               ! [-] mixing intermediate
      real(r8) :: eps_s                             ! [-] dry-soil dielectric constant
      real(r8) :: epsi                              ! [-] imaginary part
      real(r8) :: epsr                              ! [-] real part
      real(r8) :: wc                                ! [m3/m3] effective liquid water
      real(r8) :: vice                              ! [m3/m3] effective ice content
      real(r8) :: BD_gcm3                           ! [g/cm3] bulk density
      ! dielectric constant of pure ice
      real(r8), parameter :: eps_ice_r = 3.15_r8    ! [-] real part of pure ice dielectric constant
      real(r8), parameter :: eps_ice_i = 0.001_r8   ! [-] imaginary part of pure ice dielectric constant

      wc    = max(swc, 0.001_r8)
      vice  = max(ice_c, 0.0_r8)
      BD_gcm3 = BD_all / 1000.0_r8

      ! calculate the dry-soil dielectric constant
      eps_s = (1.01_r8 + 0.44_r8 * rho_soil)**2.0_r8 - 0.062_r8

      ! ---------------- real part ----------------
      beta  = 1.2748_r8 - 0.519_r8 * wf_sand - 0.152_r8 * wf_clay
      ! four-phase mixing model (Zhang & Zhao 2010):
      ! eps_m^alpha = 1 + (rho_b/rho_s)*(eps_s^alpha - 1) + m_v^beta*eps_fw^alpha - m_v + m_i*eps_i^alpha - m_i
      eaa   = 1.0_r8 + (BD_gcm3 / rho_soil) * (eps_s ** alphas - 1.0_r8) &
              + (wc ** beta) * (real(ew) ** alphas) - wc &
              + vice * (eps_ice_r ** alphas) - vice
      epsr  = eaa ** (1.0_r8/alphas)

      ! ---------------- imaginary part ----------------
      beta  = 1.33797_r8 - 0.603_r8 * wf_sand - 0.166_r8 * wf_clay
      ! the ice contribution is added to the imaginary part with its weak loss
      eaa   = (wc ** beta) * (abs(aimag(ew)) ** alphas) &
              + vice * (eps_ice_i ** alphas)
      epsi  = eaa ** (1.0_r8/alphas)

      eps   = cmplx(epsr, epsi, kind=r8)

   END SUBROUTINE diel_soil_D85

!-----------------------------------------------------------------------
   SUBROUTINE smooth_reflectivity(eps, r_s)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the smooth-surface reflectivity by the Fresnel law.
!
! REFERENCES:
!   [1] Njoku, E. G., & Kong, J. A. (1977), Theory for passive microwave
!       remote sensing of near-surface soil moisture, J. Geophys. Res.,
!       82(20), 3108-3118.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE
      complex(r8), intent(in) :: eps       ! [-] dielectric constant of the surface
      real(r8), intent(out)   :: r_s(2)    ! [-] reflectivities of flat surfaces for H and V polarizations, (1)=H,(2)=V
      complex(r8) :: g                     ! [-] parameter in Fresnel law, sqrt(eps - sin^2(theta))

      g = sqrt(eps - sin(theta)**2)
      r_s(1) = abs((cos(theta) - g)/(cos(theta) + g))**2.0_r8
      r_s(2) = abs((cos(theta)*eps - g)/(cos(theta)*eps + g))**2.0_r8

   END SUBROUTINE smooth_reflectivity

!-----------------------------------------------------------------------
   SUBROUTINE rough_reflectivity(is_desert, patchclass, r_s, r_r)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough-surface reflectivity from the smooth-surface
!   Fresnel reflectivity, using the polarization-mixing and
!   roughness-attenuation parameterizations of Kerr & Njoku (1990)
!   and Wigneron et al. (2017). Four roughness schemes are available,
!   selected by DEF_DA_RTM_rough.
!
! REFERENCES:
!   [1] Kerr, Y. H., & Njoku, E. G. (1990), A semiempirical model for
!       interpreting microwave emission from semiarid land surfaces as
!       seen from space, IEEE Trans. Geosci. Remote Sens., 28(3), 384-393.
!   [2] Wigneron, J.-P., et al. (2017), Modelling the passive microwave
!       signature from land surfaces: A review of recent results and
!       application to the L-band SMOS & SMAP soil moisture retrieval
!       algorithms, Remote Sensing of Environment, 192, 238-262.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_DA_Const
      IMPLICIT NONE
      logical, intent(in)   :: is_desert                     ! [-] flag for desert soil
      integer, intent(in)   :: patchclass                    ! [-] land cover class
      real(r8), intent(in)  :: r_s(2)                        ! [-] reflectivities of flat surfaces for H and V polarizations, (1)=H,(2)=V
      real(r8), intent(out) :: r_r(2)                        ! [-] reflectivities of rough surfaces for H and V polarizations, (1)=H,(2)=V

      real(r8) :: Q                                          ! [-] parameter for polarization mixing
      real(r8) :: hr(N_land_classification)                  ! [-] roughness parameter
      real(r8) :: nrh(N_land_classification)                 ! [-] parameter for H polarization
      real(r8) :: nrv(N_land_classification)                 ! [-] parameter for V polarization
      real(r8) :: rgh_cm                                     ! [cm] roughness height

      IF (is_desert) THEN
         r_r = r_s
      ELSE
         rgh_cm = rgh_surf*100  ! m to cm
         ! calculate parameter for polarization mixing due to surface roughness
         IF (fghz < 2.0_r8) THEN
            Q = 0.0_r8
         ELSE
            Q = 0.35_r8 * (1.0_r8 - exp(-0.6_r8 * rgh_cm**2 * fghz))    ! [1](16)
         END IF

         ! calculate rough-surface reflectivity (default settings used in [2])
         IF (DEF_DA_RTM_rough == 0) THEN
            hr(:)  = (2.0_r8 * k * rgh_surf)**2.0_r8
            nrh(:) = 0.0_r8
            nrv(:) = 0.0_r8
         ELSE IF (DEF_DA_RTM_rough == 1) THEN
            hr(:)  = hr_SMOS
            nrh(:) = 2.0_r8
            nrv(:) = 0.0_r8
         ELSE IF (DEF_DA_RTM_rough == 2) THEN
            hr(:)  = hr_SMAP
            nrh(:) = 2.0_r8
            nrv(:) = 2.0_r8
         ELSE IF (DEF_DA_RTM_rough == 3) THEN
            hr(:)  = hr_P16
            nrh(:) = -1.0_r8
            nrv(:) = -1.0_r8
         END IF

         ! rough-surface reflectivity for H and V polarizations
         r_r(1) = (Q*r_s(2) + (1.0_r8-Q)*r_s(1)) * exp(-hr(patchclass) * cos(theta)**nrh(patchclass))
         r_r(2) = (Q*r_s(1) + (1.0_r8-Q)*r_s(2)) * exp(-hr(patchclass) * cos(theta)**nrv(patchclass))
      END IF

   END SUBROUTINE rough_reflectivity

!-----------------------------------------------------------------------
   SUBROUTINE desert(t_soil, r_r, eps, tb_desert)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the desert emissivity with the Grody & Weng (2008)
!   parameterization, which accounts for the scattering by sand grains.
!
! REFERENCES:
!   [1] Grody, N. C., & Weng, F. (2008), Microwave emission and scattering
!       from deserts: Theory compared with satellite measurements, IEEE
!       Trans. Geosci. Remote Sens., 46(2), 361-375.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)    :: t_soil(2)    ! [K] desert surface temperature, (1)=H,(2)=V
      real(r8), intent(in)    :: r_r(2)       ! [-] rough surface reflectivity, (1)=H,(2)=V
      complex(r8), intent(in) :: eps          ! [-] dielectric constant of desert
      real(r8), intent(out)   :: tb_desert(2) ! [K] brightness temperature of desert, (1)=H,(2)=V

      real(r8) :: f0 = 0.7_r8                 ! [-] fractional volume of spherical particles
                                              ! (f = (4/3)*pi*r^3*n0),
                                              ! r : particle radius = 0.5 (mm)
                                              ! n0 : number of particles per unit volume
      real(r8) :: w                           ! [-] single-particle albedo
      real(r8) :: g                           ! [-] asymmetry parameter
      real(r8) :: a                           ! [-] similarity parameter
      real(r8) :: em(2)                       ! [-] desert emissivity, (1)=H,(2)=V
      real(r8) :: y_r                         ! [-] real part of y-parameter
      real(r8) :: y_i                         ! [-] imaginary part of y-parameter

      ! calculate y-parameters (eq.A15)
      y_r = (real(eps) - 1.0_r8)/(real(eps) + 2.0_r8)                 ! [1](A15)
      y_i = 3.0_r8*aimag(eps)/(real(eps) + 2.0_r8)**2                 ! [1](A15)

      ! calculate single-particle albedo (eq.A16)
      w = (1.0_r8 - f0)**4 * kr**3 * y_r**2 / &
         ((1.0_r8 - f0)**4 * kr**3 * y_r**2 + 1.5_r8*(1.0_r8 + 2.0_r8*f0)**2 * y_i)  ! [1](A16)

      ! calculate asymmetry parameter (p.374)
      g = 0.23_r8 * kr**2                                             ! [1] p.374

      ! calculate similarity parameter (eq.3b)
      a = sqrt((1.0_r8 - w)/(1.0_r8 - w*g))                           ! [1](3b)

      ! calculate desert soil emissivity (eq.A13)
      em = (1.0_r8 - r_r)*(2.0_r8*a/((1.0_r8 + a) - (1.0_r8 - a)*r_r)) ! [1](A13)

      tb_desert = t_soil * em

   END SUBROUTINE desert

!-----------------------------------------------------------------------
   SUBROUTINE veg_Wig07(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the vegetation brightness temperature and transmittance
!   using the Wigneron et al. (2007) L-MEB parameterization.
!
! REFERENCES:
!   [1] Wigneron, J.-P., et al. (2007), L-band Microwave Emission of the
!       Biosphere (L-MEB) model: Description and calibration against
!       experimental data sets over crop fields, Remote Sensing of
!       Environment, 107, 639-655.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)   :: patchclass     ! [-] land cover class
      real(r8), intent(in)  :: lai            ! [m2/m2] leaf area index
      real(r8), intent(in)  :: tleaf          ! [K] leaf temperature
      real(r8), intent(in)  :: htop           ! [m] vegetation height
      real(r8), intent(in)  :: snowdp         ! [m] snow depth
      real(r8), intent(out) :: tb_veg(2)      ! [K] vegetation brightness temperature, (1)=H,(2)=V
      real(r8), intent(out) :: gamma_p(2)     ! [-] vegetation transmittance, (1)=H,(2)=V

      real(r8) :: tau_nadir                   ! [-] nadir vegetation optical depth
      real(r8) :: tau_veg(2)                  ! [-] vegetation optical depth, (1)=H,(2)=V
      integer  :: i                           ! [-] polarization index

      ! calculate vegetation opacity (optical depth) at nadir b*VWC
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass)              ! low veg    [1](22)
      ELSE
         tau_nadir = b3(patchclass)                                   ! high veg
      END IF

      ! calculate vegetation optical depth at H and V polarizations
      tau_veg(1) = tau_nadir*(cos(theta)**2 + tth(patchclass)*sin(theta)**2)   ! [1](23)
      tau_veg(2) = tau_nadir*(cos(theta)**2 + ttv(patchclass)*sin(theta)**2)   ! [1](24)

      ! calculate brightness temperature of vegetation
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i)/cos(theta))                    ! [1](15)
         tb_veg(i)  = (1.0_r8 - w_CMEM(patchclass)) * (1.0_r8 - gamma_p(i)) * tleaf
      END DO

   END SUBROUTINE veg_Wig07

!-----------------------------------------------------------------------
   SUBROUTINE veg_Jac91(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the vegetation brightness temperature and transmittance
!   using the Jackson & Schmugge (1991) scheme. The vegetation optical
!   depth is assumed to be independent of polarization.
!
! REFERENCES:
!   [1] Jackson, T. J., & Schmugge, T. J. (1991), Vegetation effects on
!       the microwave emission of soils, Remote Sensing of Environment,
!       36(3), 203-212.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)   :: patchclass     ! [-] land cover class
      real(r8), intent(in)  :: lai            ! [m2/m2] leaf area index
      real(r8), intent(in)  :: tleaf          ! [K] leaf temperature
      real(r8), intent(in)  :: htop           ! [m] vegetation height
      real(r8), intent(in)  :: snowdp         ! [m] snow depth
      real(r8), intent(out) :: tb_veg(2)      ! [K] vegetation brightness temperature, (1)=H,(2)=V
      real(r8), intent(out) :: gamma_p(2)     ! [-] vegetation transmittance, (1)=H,(2)=V

      real(r8) :: tau_nadir                   ! [-] nadir vegetation optical depth
      integer  :: i                           ! [-] polarization index

      ! calculate nadir optical depth (equivalent to b_j * VWC in Jackson)
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass)              ! low veg
      ELSE
         tau_nadir = b3(patchclass)                                   ! high veg
      END IF

      ! calculate Jackson vegetation opacity and brightness temperature
      ! (tau is independent of polarization: tau_h = tau_v = tau_nadir)
      DO i = 1, 2
         gamma_p(i) = exp(-tau_nadir / cos(theta))
         tb_veg(i)  = (1.0_r8 - w_CMEM(patchclass)) * (1.0_r8 - gamma_p(i)) * tleaf
      END DO

   END SUBROUTINE veg_Jac91

!-----------------------------------------------------------------------
   SUBROUTINE veg_Kir79(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the vegetation opacity with the Kirdyashev et al. (1979)
!   scheme. The optical depth depends on wave number, vegetation water
!   content (VWC), and the imaginary part of the vegetation dielectric
!   constant at leaf temperature.
!
! REFERENCES:
!   [1] Kirdyashev, K. P., Chukhlantsev, A. A., & Shutko, A. M. (1979),
!       Microwave radiation of the earth's surface in the presence of a
!       vegetation cover, Radiotekhnika i Elektronika, 24, 256-264.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      integer, intent(in)   :: patchclass     ! [-] land cover class
      real(r8), intent(in)  :: lai            ! [m2/m2] leaf area index
      real(r8), intent(in)  :: tleaf          ! [K] leaf temperature
      real(r8), intent(in)  :: htop           ! [m] vegetation height
      real(r8), intent(in)  :: snowdp         ! [m] snow depth
      real(r8), intent(out) :: tb_veg(2)      ! [K] vegetation brightness temperature, (1)=H,(2)=V
      real(r8), intent(out) :: gamma_p(2)     ! [-] vegetation transmittance, (1)=H,(2)=V

      real(r8)    :: vwc                      ! [kg/m2] vegetation water content
      real(r8)    :: a_geo(2)                 ! [-] geometry structure parameter, (1)=H,(2)=V
      real(r8)    :: eps_vw_i                 ! [-] imaginary part of vegetation water dielectric constant
      real(r8)    :: tau_veg(2)               ! [-] vegetation optical depth, (1)=H,(2)=V
      complex(r8) :: eps_vw                   ! [-] vegetation water dielectric constant
      integer     :: i                        ! [-] polarization index

      ! 1. estimate VWC from LAI (generic coefficient 0.5)
      vwc = lai * 0.5_r8

      ! 2. geometry structure parameter (isotropic leaves used as default)
      a_geo(1) = 2.0_r8 / 3.0_r8
      a_geo(2) = 2.0_r8 / 3.0_r8

      ! 3. dielectric constant of pure water at leaf temperature
      CALL diel_water(0, 0.0_r8, tleaf, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, eps_vw)
      eps_vw_i = abs(aimag(eps_vw))

      ! 4. Kirdyashev optical depth
      tau_veg(1) = a_geo(1) * k * (vwc / denh2o) * eps_vw_i
      tau_veg(2) = a_geo(2) * k * (vwc / denh2o) * eps_vw_i

      ! 5. calculate gamma and tb_veg (path length correction 1/cos(theta))
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i) / cos(theta))
         tb_veg(i)  = (1.0_r8 - w_CMEM(patchclass)) * (1.0_r8 - gamma_p(i)) * tleaf
      END DO

   END SUBROUTINE veg_Kir79

!-----------------------------------------------------------------------
   SUBROUTINE snow_Pul97(t_snow, t, snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature of a snow-covered ground.
!   Snow density, liquid water content, grain size, and soil/vegetation
!   boundary reflectivity together determine the extinction coefficient
!   and the two-stream emission of the snow layer.
!
! REFERENCES:
!   [1] Matzler, C. (1987), Applications of the interaction of microwaves
!       with the natural snow cover, Remote Sensing Reviews, 2(2), 259-387.
!   [2] Anderson, E. A. (1976), A point energy and mass balance model of
!       a snow cover, NOAA Tech. Rep. NWS 19, 150 pp.
!   [3] Hallikainen, M. T., Ulaby, F. T., & Van Deventer, T. E. (1987),
!       Extinction behavior of dry snow in the 18- to 90-GHz range, IEEE
!       Trans. Geosci. Remote Sens., GE-25(6), 737-745.
!   [4] Ulaby, F. T., Moore, R. K., & Fung, A. K. (1981-1986), Microwave
!       Remote Sensing: Active and Passive, Vol. I-III, Artech House.
!   [5] Pulliainen, J., Hallikainen, M., Grandell, J., 1999. HUT snow emission model and its
!         applicability to snow water equivalent retrieval. IEEE Trans. Geosci. Remote Sens. 37,
!         1378–1390
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE
      real(r8), intent(in)  :: t_snow         ! [K] average snow temperature
      real(r8), intent(in)  :: t              ! [K] temperature at bottom of snow
      real(r8), intent(in)  :: snowdp         ! [m] snow depth
      real(r8), intent(in)  :: rho_snow       ! [kg/m3] snow density
      real(r8), intent(in)  :: liq_snow       ! [m3/m3] snow liquid water content
      real(r8), intent(in)  :: r_sn(2)        ! [-] reflectivity between snow and ground, (1)=H,(2)=V
      real(r8), intent(out) :: r_snow(2)      ! [-] reflectivity between snow and air, (1)=H,(2)=V
      real(r8), intent(out) :: tb_tos(2)      ! [K] brightness temperature of snow-covered ground, (1)=H,(2)=V

      real(r8) :: sal_snow = 0.0_r8           ! [psu] snow salinity
      real(r8) :: eps_i_r                     ! [-] real part of ice dielectric constant
      real(r8) :: eps_i_i                     ! [-] imaginary part of ice dielectric constant
      real(r8) :: eps_i_is                    ! [-] imaginary part of impure ice
      real(r8) :: eps_i_ip                    ! [-] imaginary part of pure ice
      real(r8) :: eps_ds_r                    ! [-] real part of dry-snow dielectric constant
      real(r8) :: eps_ds_i                    ! [-] imaginary part of dry-snow
      real(r8) :: eps_ws_i                    ! [-] imaginary part of wet-snow
      real(r8) :: eps_ws_r                    ! [-] real part of wet-snow
      real(r8) :: eps_w_s = 88.0_r8           ! [-] static water dielectric constant
      real(r8) :: eps_a_inf                   ! [-] infinite-frequency dielectric constant of three parts
      real(r8) :: eps_b_inf                   ! [-] infinite-frequency dielectric constant of three parts
      real(r8) :: eps_c_inf                   ! [-] infinite-frequency dielectric constant of three parts
      real(r8) :: eps_a_s                     ! [-] static dielectric constant of three parts
      real(r8) :: eps_b_s                     ! [-] static dielectric constant of three parts
      real(r8) :: eps_c_s                     ! [-] static dielectric constant of three parts
      complex(r8) :: eps_i                    ! [-] ice dielectric constant
      complex(r8) :: eps_a                    ! [-] dielectric constant of three parts
      complex(r8) :: eps_b                    ! [-] dielectric constant of three parts
      complex(r8) :: eps_c                    ! [-] dielectric constant of three parts
      complex(r8) :: eps                      ! [-] wet-snow dielectric constant
      real(r8) :: rho_ds                      ! [g/cm3] dry-snow density
      real(r8) :: rho_snow_gcm3               ! [g/cm3] snow density
      real(r8) :: rho_i = 0.916_r8            ! [g/cm3] ice density
      real(r8) :: aa = 0.005_r8               ! [-] fitting parameters
      real(r8) :: bb = 0.4975_r8              ! [-] fitting parameters
      real(r8) :: cc = 0.4975_r8              ! [-] fitting parameters
      real(r8) :: fa                          ! [GHz] relaxation frequency of three parts
      real(r8) :: fb                          ! [GHz] relaxation frequency of three parts
      real(r8) :: fc                          ! [GHz] relaxation frequency of three parts
      real(r8) :: d                           ! [mm] snow grain size
      real(r8) :: alpha                       ! [1/m] attenuation coefficient
      real(r8) :: beta                        ! [rad/m] phase coefficient
      real(r8) :: pp                          ! [1/m2] intermediate
      real(r8) :: qq                          ! [1/m2] intermediate
      real(r8) :: theta_s                     ! [rad] propagation angle in snow
      complex(r8) :: z_s                      ! [ohm] wave impedance in snow
      real(r8) :: r_sa                        ! [-] reflectivity between snow and air
      real(r8) :: tb_2                        ! [K] brightness temperature contributions from snow layer and soil
      real(r8) :: tb_3                        ! [K] brightness temperature contributions from snow layer and soil
      real(r8) :: l2_apu                      ! [-] extinction coefficient of snow (Beer's law)
      real(r8) :: q = 0.96_r8                 ! [-] parameter
      real(r8) :: ka_ws                       ! [1/m] absorption coefficient of wet snow
      real(r8) :: ka_ds                       ! [1/m] absorption coefficient of dry snow
      real(r8) :: ke                          ! [1/m] extinction coefficient of wet snow
      real(r8) :: ke_ds                       ! [1/m] extinction coefficient of dry snow
      real(r8) :: ks                          ! [1/m] scattering coefficient (all as coded; verify)
      real(r8) :: b_ds                        ! [-] loss-tangent intermediates
      real(r8) :: b_ws                        ! [-] loss-tangent intermediates
      real(r8) :: wk_h                        ! [m] as coded; verify
      integer  :: i                           ! [-] polarization index

      IF (snowdp > 0.01_r8) THEN
         rho_snow_gcm3 = rho_snow / 1000.0_r8

         ! calculate dielectric constant of ice
         CALL diel_ice(t_snow, eps_i)
         eps_i_r = real(eps_i)
         eps_i_i = aimag(eps_i)

         ! consider the effect of salinity on the dielectric constant of ice
         eps_i_is = 0.0026_r8/fghz + 0.00023_r8*(fghz**0.87_r8)   ! impure ice -5(C)
         eps_i_ip = 6.e-4_r8/fghz + 6.5e-5_r8*(fghz**1.07_r8)     ! pure ice -5(C)
         eps_i_i  = eps_i_i + (eps_i_is - eps_i_ip)*sal_snow/13.0d0

         ! dry snow dielectric constant (Polder-van Santen mixing of air and ice)
         rho_ds   = (rho_snow_gcm3 - liq_snow)/(1.0_r8 - liq_snow)
         eps_ds_r = 1.0_r8 + 1.58_r8*rho_ds/(1.0_r8 - 0.365_r8*rho_ds)
         eps_ds_i = 3.0_r8*(rho_ds/rho_i)*eps_i_i*(eps_ds_r**2)*(2.0_r8*eps_ds_r + 1.0_r8)/ &
                    ((eps_i_r + 2.0_r8*eps_ds_r)*(eps_i_r + 2.0_r8*eps_ds_r**2))

         ! wet snow dielectric constant (Matzler 1987)
         IF (liq_snow > 0.0_r8) THEN
            ! calculate relaxation frequency of three parts (eq.2.26)
            fa = f0w*(1.0_r8 + (aa*(eps_w_s - eps_w_inf)/(eps_ds_r + (aa*(eps_w_inf - eps_ds_r)))))
            fb = f0w*(1.0_r8 + (bb*(eps_w_s - eps_w_inf)/(eps_ds_r + (bb*(eps_w_inf - eps_ds_r)))))
            fc = f0w*(1.0_r8 + (cc*(eps_w_s - eps_w_inf)/(eps_ds_r + (cc*(eps_w_inf - eps_ds_r)))))

            ! infinite-frequency dielectric constant of three parts
            eps_a_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.0_r8)/(1.0_r8 + aa*((eps_w_inf/eps_ds_r) - 1.0_r8))
            eps_b_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.0_r8)/(1.0_r8 + bb*((eps_w_inf/eps_ds_r) - 1.0_r8))
            eps_c_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.0_r8)/(1.0_r8 + cc*((eps_w_inf/eps_ds_r) - 1.0_r8))

            ! static dielectric constant of three parts
            eps_a_s = (liq_snow/3.0_r8)*(eps_w_s - eps_ds_r)/(1.0_r8 + aa*((eps_w_s/eps_ds_r) - 1.0_r8))
            eps_b_s = (liq_snow/3.0_r8)*(eps_w_s - eps_ds_r)/(1.0_r8 + bb*((eps_w_s/eps_ds_r) - 1.0_r8))
            eps_c_s = (liq_snow/3.0_r8)*(eps_w_s - eps_ds_r)/(1.0_r8 + cc*((eps_w_s/eps_ds_r) - 1.0_r8))

            ! Debye equations
            eps_a = eps_a_inf + (eps_a_s - eps_a_inf)/(1.0_r8 + jj*fghz/fa)
            eps_b = eps_b_inf + (eps_b_s - eps_b_inf)/(1.0_r8 + jj*fghz/fb)
            eps_c = eps_c_inf + (eps_c_s - eps_c_inf)/(1.0_r8 + jj*fghz/fc)

            ! wet-snow dielectric constant
            eps = eps_a + eps_b + eps_c + (eps_ds_r - jj*eps_ds_i)
            eps_ws_r = real(eps)
            eps_ws_i = -1.0_r8*aimag(eps)
         ELSE
            eps_ws_r = eps_ds_r
            eps_ws_i = eps_ds_i
         END IF

         ! calculate propagation angle in snow (air -> snow)
         alpha = k*abs(aimag(sqrt(eps_ws_r - jj*eps_ws_i)))
         beta  = k*real(sqrt(eps_ws_r - jj*eps_ws_i))
         pp    = 2.0_r8*alpha*beta
         qq    = beta**2 - alpha**2 - (k*k)*(sin(theta)**2)
         theta_s = atan(k*sin(theta)/((1.0_r8/sqrt(2.0_r8))*sqrt(sqrt(pp**2 + qq**2) + qq)))

         ! calculate wave impedance in snow
         z_s = z0/sqrt(eps_ws_r - jj*eps_ws_i)

         ! calculate brightness temperature above snow for H and V polarizations
         DO i = 1, 2
            ! Fresnel reflection coefficient between snow and air
            IF (i == 1) THEN
               r_sa = abs((z_s*cos(theta) - z0*cos(theta_s))/(z_s*cos(theta) + z0*cos(theta_s)))**2
            ELSE
               r_sa = abs((z0*cos(theta) - z_s*cos(theta_s))/(z0*cos(theta) + z_s*cos(theta_s)))**2
            END IF

            ! snow grain size (mm) (Anderson 1976, eq.5.1)
            d = min(1000.0_r8*(1.6e-4_r8 + 1.1e-13_r8*(rho_snow**4)), 3.0_r8)
            ke_ds = 0.0018_r8*(fghz**2.8_r8)*(d**2)/4.3429_r8        ! [3](14)

            ! dry snow absorption coefficient
            b_ds  = (eps_ds_i/eps_ds_r)**2
            ka_ds = 2.0_r8*omega*sqrt(mu0*eps_0*eps_ds_r)*sqrt(b_ds/(2.0_r8*(sqrt(1.0_r8+b_ds) + 1.0_r8)))
            IF (ke_ds < ka_ds) ke_ds = ka_ds

            ! wet snow absorption coefficient
            b_ws  = (eps_ws_i/eps_ws_r)**2
            ka_ws = 2.0_r8*omega*sqrt(mu0*eps_0*eps_ws_r)*sqrt(b_ws/(2.0_r8*(sqrt(1.0_r8+b_ws) + 1.0_r8)))

            ! total extinction (assuming scattering is the same for dry and wet snow)
            ke = (ke_ds - ka_ds) + ka_ws
            ks = ke_ds - ka_ds

            wk_h = snowdp/rho_snow_gcm3
            IF (((ke - q*ks)*(1.0_r8/cos(theta_s))*wk_h) > (log(HUGE(1.0_r8)) - 1.0_r8)) THEN
               l2_apu = sqrt(HUGE(1.0_r8))
            ELSE
               l2_apu = min(sqrt(HUGE(1.0_r8)), exp((ke - q*ks)*(1.0_r8/cos(theta_s))*wk_h))
            END IF

            ! brightness temperature through snow-air layer ([4] pp.243, eq.4.161)
            tb_2 = (1.0_r8 + r_sn(i)/l2_apu)*(1.0_r8 - r_sa)*t_snow* &
                   (ka_ws/(ke - q*ks))*(1.0_r8 - 1.0_r8/l2_apu)/(1.0_r8 - r_sn(i)*r_sa/l2_apu**2)

            ! brightness temperature through soil-snow layer ([4] pp.243, eq.4.162)
            tb_3 = ((1.0_r8 - r_sn(i))*(1.0_r8 - r_sa)*t)/(l2_apu*(1.0_r8 - r_sn(i)*r_sa/l2_apu**2))

            ! brightness temperature of snow-covered ground
            tb_tos(i) = tb_2 + tb_3
            r_snow(i) = 1.0_r8 - (tb_2/t_snow + tb_3/t)
         END DO
      END IF

   END SUBROUTINE snow_Pul97

!=======================================================================
! LandEM: bulk two-stream microwave land emissivity model
! (Weng, Yan & Grody 2001). This is the main entry point; all
! helper subroutines used below are defined further down.
!=======================================================================
  SUBROUTINE LandEM(Angle,                 &
                    Frequency,             &
                    Soil_Moisture_Content, &
                    Vegetation_Fraction,   &
                    Soil_Temperature,      &
                    t_skin,                &
                    Lai,                   &
                    sand_frac,             &
                    clay_frac,             &
                    bulk_density,          &
                    Snow_Depth,            &
                    Emissivity_H,          &
                    Emissivity_V)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Simulate the microwave surface emissivity over land (bare soil,
!   vegetation canopy, and/or snow) with the LandEM bulk two-stream
!   model.
!
!   Two vertical branches are handled:
!     - Snow_Depth > 0.1 mm  : snow physical model (Snow_Diel +
!                              Two_Stream_Solution) or, for deep
!                              snow (> 10 mm) and freq >= 10 GHz, the
!                              AMSR-E empirical lookup (SnowEM_Default)
!     - Snow_Depth <= 0.1 mm : no-snow branch (bare soil + canopy)
!
!   The soil dielectric constant is selected by Diel_Model_Mironov
!   (default .TRUE.); if .FALSE., Dobson's model is used instead.
!   Surface roughness is selected by Roughness_Model_ChenWeng (default
!   .TRUE.); if .FALSE., the Qh model is used instead.
!
! INPUT:
!   Angle                 incidence angle (degrees)
!   Frequency             frequency (GHz)
!   Soil_Moisture_Content volumetric soil water content (m3/m3)
!   Vegetation_Fraction   vegetation fraction (0-1)
!   Soil_Temperature      soil temperature (K)
!   t_skin                skin temperature (K)
!   Lai                   leaf area index (-)
!   sand_frac             sand mass fraction (0-1)
!   clay_frac             clay mass fraction (0-1)
!   bulk_density          bulk density (kg/m3)
!   Snow_Depth            snow depth (m)
!
! OUTPUT:
!   Emissivity_H          surface emissivity at H polarization
!   Emissivity_V          surface emissivity at V polarization
!
! REFERENCES:
!   [1] Weng, F., Yan, B., & Grody, N. C. (2001). A microwave land
!       emissivity model. J. Geophys. Res. Atmos., 106(D17), 20115-20123.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8), intent(in) :: Angle                   ! [degree]
    REAL(r8), intent(in) :: Frequency               ! [GHz]
    REAL(r8), intent(in) :: Soil_Moisture_Content   ! [m3/m3]
    REAL(r8), intent(in) :: Vegetation_Fraction     ! [-]
    REAL(r8), intent(in) :: Soil_Temperature        ! [K]
    REAL(r8), intent(in) :: t_skin                  ! [K]
    REAL(r8), intent(in) :: Lai                     ! [m2/m2]
    REAL(r8), intent(in) :: sand_frac               ! [-] mass fraction 0-1
    REAL(r8), intent(in) :: clay_frac               ! [-] mass fraction 0-1
    REAL(r8), intent(in) :: bulk_density            ! [kg/m3]
    REAL(r8), intent(in) :: Snow_Depth              ! [m]
    REAL(r8), intent(out):: Emissivity_V            ! [-] surface emissivity at V polarization
    REAL(r8), intent(out):: Emissivity_H            ! [-] surface emissivity at H polarization

    REAL(r8), PARAMETER :: snow_depth_c     = 0.01_r8     ! [m] snow depth threshold for physical vs empirical model
    REAL(r8), PARAMETER :: tsoilc_undersnow = 280.0_r8    ! [K] soil temperature cap under snow
    REAL(r8), PARAMETER :: rhos             = 2650.0_r8   ! [kg/m3] density of solid soil particles

    REAL(r8) :: sand                                ! [-] mass fraction 0-1
    REAL(r8) :: clay                                ! [-] mass fraction 0-1
    REAL(r8) :: rhob                                ! [kg/m3] bulk density
    REAL(r8) :: mv                                  ! [m3/m3] soil moisture
    REAL(r8) :: veg_frac                            ! [-] vegetation fraction
    REAL(r8) :: theta_rad                           ! [rad] angles
    REAL(r8) :: theta_i                             ! [rad] angles
    REAL(r8) :: theta_t                             ! [rad] angles
    REAL(r8) :: mu                                  ! [-] cos(theta)
    REAL(r8) :: r21_h                               ! [-] reflectivities
    REAL(r8) :: r21_v                               ! [-] reflectivities
    REAL(r8) :: r23_h                               ! [-] reflectivities
    REAL(r8) :: r23_v                               ! [-] reflectivities
    REAL(r8) :: t21_v                               ! [-] transmittance
    REAL(r8) :: t21_h                               ! [-] transmittance
    REAL(r8) :: gv                                  ! [-] asymmetry factors
    REAL(r8) :: gh                                  ! [-] asymmetry factors
    REAL(r8) :: ssalb_h                             ! [-] single-scattering albedo
    REAL(r8) :: ssalb_v                             ! [-] single-scattering albedo
    REAL(r8) :: tau_h                               ! [-] optical depth
    REAL(r8) :: tau_v                               ! [-] optical depth
    REAL(r8) :: mge                                 ! [-] gravimetric water content
    REAL(r8) :: leaf_thick                          ! [m]
    REAL(r8) :: rad                                 ! [m]
    REAL(r8) :: sigma                               ! [m] roughness
    REAL(r8) :: va                                  ! [-] volume fraction
    REAL(r8) :: ep_real                             ! [-] dielectric constant parts
    REAL(r8) :: ep_imag                             ! [-] dielectric constant parts
    REAL(r8) :: t_soil                              ! [K]
    REAL(r8) :: rhoveg                              ! [kg/m3]
    REAL(r8) :: vlai                                ! [m2/m2]
    REAL(r8) :: local_snow_depth                    ! [m]
    COMPLEX(r8) :: esoil                            ! [-] dielectric constants
    COMPLEX(r8) :: eveg                             ! [-] dielectric constants
    COMPLEX(r8) :: esnow                            ! [-] dielectric constants
    COMPLEX(r8) :: eair                             ! [-] dielectric constants
    LOGICAL :: SnowEM_Physical_Model                ! [-] switch
    LOGICAL :: Diel_Model_Mironov = .TRUE.          ! [-] switch
    LOGICAL :: Roughness_Model_ChenWeng = .TRUE.    ! [-] switch

    eair = CMPLX(1.0_r8, 0.0_r8, r8)
    theta_rad = Angle*pi/180.0_r8

    ! assign local variables
    mv               = Soil_Moisture_Content
    veg_frac         = Vegetation_Fraction
    t_soil           = Soil_Temperature
    local_snow_depth = Snow_Depth

    ! surface composition parameters (already in [0,1] on input)
    sand = sand_frac
    clay = clay_frac
    rhob = bulk_density

    ! canopy parameters
    vlai       = Lai*veg_frac
    mge        = 0.30_r8
    rhoveg     = 250.0_r8
    leaf_thick = 0.12e-3_r8

    ! check soil/skin temperature consistency
    if ( (t_soil <= 100.0_r8 .OR.  t_soil >= 350.0_r8) .AND. &
         (t_skin >= 100.0_r8 .AND. t_skin <= 350.0_r8) ) t_soil = t_skin
    mv = MAX(MIN(mv,1.0_r8),0.0_r8)

!#############################################################################
! Branch 1: snow present (Snow_Depth > 0.1 mm)
!#############################################################################
    IF (local_snow_depth > 0.0001_r8) THEN
      ! physical snow model unless deep snow and high frequency
      SnowEM_Physical_Model = .TRUE.
      if (local_snow_depth > snow_depth_c) SnowEM_Physical_Model = .FALSE.
      IF (Frequency < 10.0_r8) SnowEM_Physical_Model = .TRUE.

      IF (SnowEM_Physical_Model) THEN
        ! ---- snow physical model ----
        ep_real = 3.2_r8
        ep_imag = -0.0005_r8
        sigma = 0.008_r8
        local_snow_depth = MIN(local_snow_depth,1.0_r8)
        va = 0.4_r8 + 0.4_r8*local_snow_depth
        va = MAX(MIN(va,1.0_r8),0.0_r8)
        rad = MIN((0.0005_r8 + 0.005_r8*local_snow_depth),0.001_r8)
        t_soil = MIN(t_soil,tsoilc_undersnow)

        CALL Snow_Diel(Frequency, ep_real, ep_imag, rad, va, esnow)

        ! soil dielectric constant (frozen vs. unfrozen)
        IF (t_soil > tfrz) THEN
           IF (Diel_Model_Mironov) THEN
              CALL Soil_Diel_Mironov(Frequency, t_soil, mv, clay, esoil)
           ELSE
              CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
           ENDIF
        ELSE
           CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
        ENDIF

        ! snow/air interface reflection and transmission
        theta_i = ASIN(REAL(SIN(theta_rad)*SQRT(eair)/SQRT(esnow),r8))
        CALL Smooth_Reflectance(esnow, eair, theta_i, theta_rad, r21_v, r21_h)
        CALL Transmittance(esnow, eair, theta_i, theta_rad, t21_v, t21_h)

        mu      = COS(theta_i)
        theta_t = ASIN(REAL(SIN(theta_i)*SQRT(esnow)/SQRT(esoil),r8))

        ! snow/soil interface reflection
        CALL Smooth_Reflectance(esnow, esoil, theta_i, theta_t, r23_v, r23_h)

        IF (Roughness_Model_ChenWeng) THEN
           CALL Roughness_Reflectance_ChenWeng(Frequency, sigma, theta_rad, r23_v, r23_h)
        ELSE
           CALL Roughness_Reflectance_Qh(Frequency, sigma, r23_v, r23_h)
        ENDIF

        ! snow optical properties
        CALL Snow_Optic(Frequency,rad,local_snow_depth,va,ep_real, ep_imag,gv,gh,&
                        ssalb_v,ssalb_h,tau_v,tau_h)

        ! two-stream solution -> emissivity
        CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                                 r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                                 frequency, t_soil, t_skin)
      ELSE
        ! ---- AMSR-E empirical snow lookup (deep snow, freq >= 10 GHz) ----
        CALL SnowEM_Default(Frequency,t_skin, local_snow_depth,Emissivity_V,Emissivity_H)
      END IF
    ELSE
!#############################################################################
! Branch 2: no snow (bare soil + canopy)
!#############################################################################
      veg_frac = MAX(MIN(veg_frac,1.0_r8),0.0_r8)
      mu    = COS(theta_rad)
      sigma = 0.01_r8

      ! no top half-space contribution (single canopy layer)
      r21_h = 0.0_r8
      r21_v = 0.0_r8
      t21_h = 1.0_r8
      t21_v = 1.0_r8

      ! soil dielectric constant (frozen vs. unfrozen)
      IF (t_soil >= tfrz) THEN
           IF (Diel_Model_Mironov) THEN
              CALL Soil_Diel_Mironov(Frequency, t_soil, mv, clay, esoil)
           ELSE
              CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
           ENDIF
      ELSE
          CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
      ENDIF

      ! soil surface reflection
      theta_t = ASIN(REAL(SIN(theta_rad)*SQRT(eair)/SQRT(esoil),r8))
      CALL Smooth_Reflectance(eair, esoil, theta_rad, theta_t, r23_v, r23_h)

      IF (Roughness_Model_ChenWeng) THEN
         CALL Roughness_Reflectance_ChenWeng(Frequency, sigma, theta_rad, r23_v, r23_h)
      ELSE
         CALL Roughness_Reflectance_Qh(Frequency, sigma, r23_v, r23_h)
      ENDIF

      ! canopy dielectric and optical properties
      CALL Canopy_Diel(Frequency, mge, eveg, rhoveg, t_soil)
      CALL Canopy_Optic(vlai,Frequency,theta_rad,eveg,leaf_thick,gv,gh,ssalb_v,ssalb_h,tau_v,tau_h)

      ! two-stream solution -> emissivity
      CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                               r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                               frequency, t_soil, t_skin)
    END IF

  END SUBROUTINE LandEM

!-----------------------------------------------------------------------
  SUBROUTINE Soil_Diel_Mironov(freq,t_soil,mv,clay,esm)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil with the
!   temperature- and mineralogy-dependent Mironov model (2009).
!   Same physical model as diel_soil_M09, but with the input clay
!   mass fraction in [0,1] and a different internal naming.
!
! INPUT:
!   freq      frequency (GHz)
!   t_soil    soil temperature (K)
!   mv        volumetric soil moisture (m3/m3)
!   clay      clay mass fraction (0-1)
!
! OUTPUT:
!   esm       dielectric constant of bare soil
!
! REFERENCES:
!   [1] Mironov, V. L., Kosolapova, L. G., & Fomin, S. V. (2009),
!       Physically and mineralogically based spectroscopic dielectric
!       model for moist soils, IEEE Trans. Geosci. Remote Sens., 47(7),
!       2059-2070.
!   [2] Mironov, V. L., & Fomin, S. V. (2009), Temperature and mineralogy
!       dependable model for microwave dielectric spectra of moist soils,
!       PIERS Online, 5(5), 411-415.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8) :: f                                   ! [Hz]
    REAL(r8) :: freq                                ! [GHz]
    REAL(r8) :: t_soil                              ! [K]
    REAL(r8) :: Ts                                  ! [C] reference
    REAL(r8) :: T                                   ! [C]
    REAL(r8) :: mv                                  ! [m3/m3]
    REAL(r8) :: mvt                                 ! [m3/m3]
    REAL(r8) :: clay                                ! [-] 0-1
    REAL(r8) :: Tao_wu                              ! [s]
    REAL(r8) :: Tao_wb                              ! [s]
    REAL(r8) :: Sigma_eff_wu                        ! [S/m]
    REAL(r8) :: Sigma_eff_wb                        ! [S/m]
    REAL(r8) :: E_w0u                               ! [-]
    REAL(r8) :: E_w0b                               ! [-]
    REAL(r8) :: E_wlimit                            ! [-]
    REAL(r8) :: E_0                                 ! [F/m]
    REAL(r8) :: nd                                  ! [-] refractive index
    REAL(r8) :: kd                                  ! [-] attenuation coefficient
    REAL(r8) :: nb                                  ! [-] refractive index
    REAL(r8) :: kb                                  ! [-] attenuation coefficient
    REAL(r8) :: nu                                  ! [-] refractive index
    REAL(r8) :: ku                                  ! [-] attenuation coefficient
    REAL(r8) :: nm                                  ! [-] refractive index
    REAL(r8) :: km                                  ! [-] attenuation coefficient
    REAL(r8) :: Beta_b                              ! [-] temperature-correction coefficient
    REAL(r8) :: Beta_u                              ! [-] temperature-correction coefficient
    REAL(r8) :: Delta_HbR                           ! [-] temperature-correction coefficient
    REAL(r8) :: Delta_SbR                           ! [-] temperature-correction coefficient
    REAL(r8) :: Delta_HuR                           ! [-] temperature-correction coefficient
    REAL(r8) :: Delta_SuR                           ! [-] temperature-correction coefficient
    REAL(r8) :: Beta_sigmab                         ! [-] temperature-correction coefficient
    REAL(r8) :: Beta_sigmau                         ! [-] temperature-correction coefficient
    REAL(r8) :: Fb                                  ! [-]
    REAL(r8) :: Fu                                  ! [-]
    REAL(r8) :: E_wb1                               ! [-]
    REAL(r8) :: E_wb2                               ! [-]
    REAL(r8) :: E_wu1                               ! [-]
    REAL(r8) :: E_wu2                               ! [-]
    REAL(r8) :: E_m1                                ! [-]
    REAL(r8) :: E_m2                                ! [-]
    COMPLEX(r8) :: esm                              ! [-] dielectric constant of bare soil

    f = freq*1.0e9_r8

    ! maximum combined (bound) water content
    mvt = 0.02863_r8 + 0.30673_r8*clay

    ! dry-soil refractive index and normalized attenuation coefficient
    nd = 1.634_r8 - 0.539_r8*clay + 0.2748_r8*clay**2
    kd = 0.03952_r8 - 0.04038_r8*clay

    ! bound water parameters at reference temperature Ts = 20 C
    E_w0b = 79.8_r8 - 85.4_r8*clay + 32.7_r8*clay**2
    Tao_wb = 1.062e-11_r8 + 3.450e-12_r8*clay
    Sigma_eff_wb = 0.3112_r8 + 0.467_r8*clay

    ! free water parameters at reference temperature Ts = 20 C
    E_w0u = 100.0_r8
    Tao_wu = 8.5e-12_r8

    ! reference temperature and current temperature in Celsius
    Ts = 20.0_r8
    T = t_soil - tfrz

    ! universal constants
    E_wlimit = 4.9_r8
    E_0 = 8.854e-12_r8

    ! temperature-correction coefficients for bound water
    Beta_b = 8.67e-19_r8 - 1.26e-3_r8*clay + 1.84e-3_r8*clay**2 - 9.77e-4_r8*clay**3 - 1.39e-7_r8*clay**4
    Beta_u = 1.11e-4_r8 - 1.603e-5_r8*clay + 1.239e-5_r8*clay**2 + 8.33e-7_r8*clay**3 - 1.007e-6_r8*clay**4
    Delta_HbR = 1467.0_r8 + 2697.0_r8*clay - 980.0_r8*clay**2 + 1.368e-4_r8*clay**3 - 8.61e-5_r8*clay**4
    Delta_SbR = 0.888_r8 + 9.7_r8*clay - 4.262_r8*clay**2 + 6.79e-15_r8*clay**3 + 4.263e-14_r8*clay**4
    Delta_HuR = 2231.0_r8 - 143.1_r8*clay + 223.2_r8*clay**2 - 142.1_r8*clay**3 + 27.14_r8*clay**4
    Delta_SuR = 3.649_r8 - 0.4894_r8*clay + 0.763_r8*clay**2 - 0.4859_r8*clay**3 + 0.0928_r8*clay**4
    Beta_sigmab = 0.0028_r8 + 2.094e-2_r8*clay - 1.229e-2_r8*clay**2 - 5.03e-16_r8*clay**3 + 4.163e-16_r8*clay**4
    Beta_sigmau = 0.00108_r8 + 0.1413_r8*clay - 0.2555_r8*clay**2 + 0.2147_r8*clay**3 - 0.0711_r8*clay**4
    Sigma_eff_wu = 0.05_r8 + 1.4_r8*(1.0_r8 -(1.0_r8 - clay)**4.664_r8)

    ! static dielectric constant of bound water at temperature T
    Fb = LOG((E_w0b - 1.0_r8)/(E_w0b + 2.0_r8))
    E_w0b = (1.0_r8 + 2.0_r8*EXP(Fb-Beta_b*(T - Ts)))/(1.0_r8 - EXP(Fb -Beta_b*(T -Ts)))

    ! static dielectric constant of free water at temperature T
    Fu = LOG((E_w0u - 1.0_r8)/(E_w0u + 2.0_r8))
    E_w0u = (1.0_r8 + 2.0_r8*EXP(Fb - Beta_u*(T - Ts)))/(1.0_r8 - EXP(Fb - Beta_u*(T - Ts)))

    ! relaxation time of bound water and free water at temperature T
    Tao_wb = 48e-12_r8 /( T + tfrz )*EXP( Delta_HbR / ( T + tfrz )-Delta_SbR)
    Tao_wu = 48e-12_r8 /( T + tfrz )*EXP( Delta_HuR / ( T + tfrz )-Delta_SuR)

    ! effective conductivity of bound water and free water at temperature T
    Sigma_eff_wb = Sigma_eff_wb + Beta_sigmab*(T-Ts)
    Sigma_eff_wu = Sigma_eff_wu + Beta_sigmau*(T-Ts)

    ! refractive index and normalized attenuation coefficient of bound water
    E_wb1 = E_wlimit+(E_w0b-E_wlimit)/(1.0_r8 + (2.0_r8*pi*f*Tao_wb)**2)
    E_wb2 = 2.0_r8*pi*f*Tao_wb*(E_w0b-E_wlimit)/(1.0_r8 + (2.0_r8*pi*f*Tao_wb)**2) &
            + (Sigma_eff_wb)/(2.0_r8*pi*E_0*f)
    nb = SQRT(SQRT(E_wb1**2 + E_wb2**2) + E_wb1)/SQRT(2.0_r8)
    kb = SQRT(SQRT(E_wb1**2 + E_wb2**2) - E_wb1)/SQRT(2.0_r8)

    ! refractive index and normalized attenuation coefficient of free water
    E_wu1 = E_wlimit + (E_w0u-E_wlimit)/(1.0_r8 + (2.0_r8*pi*f*Tao_wu)**2)
    E_wu2 = 2.0_r8*pi*f*Tao_wu*(E_w0u-E_wlimit)/(1.0_r8 + (2.0_r8*pi*f*Tao_wu)**2) &
            + (Sigma_eff_wu)/(2.0_r8*pi*E_0*f)
    nu = SQRT(SQRT(E_wu1**2 + E_wu2**2) + E_wu1)/SQRT(2.0_r8)
    ku = SQRT(SQRT(E_wu1**2 + E_wu2**2) - E_wu1)/SQRT(2.0_r8)

    ! refractive index and normalized attenuation coefficient of soil
    if (mv <= mvt) then
       nm = nd+(nb-1.0_r8)*mv
       km = kd+kb*mv
    else
       nm = nd+(nb-1.0_r8)*mvt+(nu-1.0_r8)*(mv-mvt)
       km = kd+kb*mvt+ku*(mv-mvt)
    endif

    ! dielectric constant of bare soil
    E_m1 = nm**2-km**2
    E_m2 = 2.0_r8*nm*km
    esm = E_m1 - CMPLX(0.0_r8,1.0_r8,r8)*E_m2

    if (AIMAG(esm) >= 0.0_r8) esm = CMPLX(REAL(esm,r8),-0.0001_r8, r8)

  END SUBROUTINE Soil_Diel_Mironov

!-----------------------------------------------------------------------
  SUBROUTINE Soil_Diel_Dobson(freq,t_soil,mv,rhob,rhos,sand,clay,esm)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet (or frozen) soil with
!   the Dobson (1985) mixing model as extended by Peplinski & Ulaby
!   (1995) to low frequencies (0.3-1.3 GHz).
!
! INPUT:
!   freq      frequency (GHz)
!   t_soil    soil temperature (K)
!   mv        volumetric soil moisture (m3/m3)
!   rhob      bulk density of soil (kg/m3)
!   rhos      density of solid soil particles (kg/m3)
!   sand      sand mass fraction (0-1)
!   clay      clay mass fraction (0-1)
!
! OUTPUT:
!   esm       dielectric constant of bare soil
!
! REFERENCES:
!   [1] Dobson, M. C., et al. (1985), Microwave dielectric behavior of wet
!       soil - Part II: Dielectric mixing models, IEEE Trans. Geosci.
!       Remote Sens., GE-23(1), 35-46.
!   [2] Peplinski, N. R., Ulaby, F. T., & Dobson, M. C. (1995), Dielectric
!       properties of soils in the 0.3-1.3-GHz range, IEEE Trans. Geosci.
!       Remote Sens., 33(3), 803-807.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8) :: f                                   ! [Hz]
    REAL(r8) :: tauw                                ! [s]
    REAL(r8) :: freq                                ! [GHz]
    REAL(r8) :: t_soil                              ! [K]
    REAL(r8) :: mv                                  ! [m3/m3]
    REAL(r8) :: rhob                                ! [kg/m3]
    REAL(r8) :: rhos                                ! [kg/m3]
    REAL(r8) :: sand                                ! [-] 0-1
    REAL(r8) :: clay                                ! [-] 0-1
    REAL(r8) :: alpha                               ! [-]
    REAL(r8) :: beta                                ! [-]
    REAL(r8) :: ess                                 ! [-]
    REAL(r8) :: rhoef                               ! [S/m]
    REAL(r8) :: t                                   ! [C]
    REAL(r8) :: eswi                                ! [-]
    REAL(r8) :: eswo                                ! [-]
    REAL(r8) :: esof                                ! [F/m]
    REAL(r8) :: vlw                                 ! [m3/m3]
    REAL(r8) :: vic                                 ! [m3/m3]
    REAL(r8) :: rhob_gcm3                           ! [g/cm3]
    REAL(r8) :: rhos_gcm3                           ! [g/cm3]
    COMPLEX(r8) :: esm                              ! [-] dielectric constant
    COMPLEX(r8) :: esw                              ! [-] dielectric constant
    COMPLEX(r8) :: es1                              ! [-] dielectric constant
    COMPLEX(r8) :: es2                              ! [-] dielectric constant
    COMPLEX(r8) :: esice                            ! [-] dielectric constant

    rhob_gcm3 = rhob / 1000.0_r8
    rhos_gcm3 = rhos / 1000.0_r8

    alpha = 0.65_r8
    beta  = 1.09_r8 - 0.11_r8*sand + 0.18_r8*clay
    ess = (1.01_r8 + 0.44_r8*rhos_gcm3)**2 - 0.062_r8

    ! effective conductivity of soil water (Dobson 1985)
    rhoef = -1.645_r8 + 1.939_r8*rhob_gcm3 - 2.25622_r8*sand + 1.594_r8*clay

    ! alternative formulation for 0.3-1.3 GHz (Peplinski & Ulaby 1995)
    if (freq < 1.4_r8) then
      rhoef = 0.0467_r8 + 0.2204_r8*rhob_gcm3 - 0.41111_r8*sand + 0.6614_r8*clay
    end if

    t = t_soil - tfrz
    f = freq*1.0e9_r8

    ! high-frequency and free-space limits of water permittivity
    eswi = 4.9_r8
    esof = 8.854e-12_r8

    ! static dielectric constant and relaxation time of water at T
    eswo = 87.134_r8 + (-1.949e-1_r8 + (-1.276e-2_r8 + 2.491e-4_r8*t)*t)*t
    tauw = 1.1109e-10_r8 + (-3.824e-12_r8 + (6.938e-14_r8 - 5.096e-16_r8*t)*t)*t

    ! soil-water dielectric constant (real and imaginary parts)
    if (mv > 0.0_r8) then
       es1 = CMPLX(eswi, -rhoef*(rhos_gcm3-rhob_gcm3)/(2.0_r8*pi*f*esof*rhos_gcm3*mv), r8)
    else
       es1 = CMPLX(eswi, 0.0_r8, r8)
    endif

    es2 = CMPLX(eswo-eswi, 0.0_r8, r8)/CMPLX(1.0_r8, f*tauw, r8)
    esw = es1 + es2

    ! unfrozen vs. frozen soil
    if (t_soil >= (tfrz + TS_FROZEN_THRESHOLD)) then
      esm = 1.0_r8 + (ess**alpha - 1.0_r8)*rhob_gcm3/rhos_gcm3 + mv**beta*esw**alpha - mv
      esm = esm**(1.0_r8/alpha)
    else
      ! frozen-soil branch: split water into unfrozen (vlw) and frozen (vic)
      esice = CMPLX(3.15_r8, 57.34_r8*(1.0_r8/f + 2.48e-14_r8*f**0.5)*EXP(0.0362_r8*t), r8)
      vic   = RATIO_FROZEN*mv
      vlw   = (1.0_r8 - RATIO_FROZEN)*mv
      esm   = 1.0_r8 + (ess**alpha - 1.0_r8)*rhob_gcm3/rhos_gcm3 + (vlw**beta*esw**alpha - vlw) &
              + vic*(esice**alpha - 1.0_r8)
      esm   = esm**(1.0_r8/alpha)
    endif

    if (AIMAG(esm) >= 0.0_r8) esm = CMPLX(REAL(esm,r8),-0.0001_r8, r8)

  END SUBROUTINE Soil_Diel_Dobson

!-----------------------------------------------------------------------
  SUBROUTINE Snow_Diel(frequency,ep_real,ep_imag,rad,frac,ep_eff)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective dielectric constant of a dense medium of
!   ice particles in air (snow), using the Polder-van Santen mixing
!   formula.
!
! INPUT:
!   frequency  frequency (GHz)
!   ep_real    real part of the dielectric constant of ice particles
!   ep_imag    imaginary part of the dielectric constant of ice particles
!   rad        particle radius (m)
!   frac       volume fraction of ice particles (0-1)
!
! OUTPUT:
!   ep_eff     effective dielectric constant of the snow medium
!
! REFERENCES:
!   [1] Matzler, C. (1987), Applications of the interaction of microwaves
!       with the natural snow cover, Remote Sensing Reviews, 2(2), 259-387.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    IMPLICIT NONE
    REAL(r8) :: ep_imag         ! [-] dielectric constant part of ice particles
    REAL(r8) :: ep_real         ! [-] dielectric constant part of ice particles
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: rad             ! [m]
    REAL(r8) :: frac            ! [-]
    REAL(r8) :: k0              ! [1/m]
    REAL(r8) :: yr              ! [-]
    REAL(r8) :: yi              ! [-]
    COMPLEX(r8) :: y            ! [-]
    COMPLEX(r8) :: ep_r         ! [-]
    COMPLEX(r8) :: ep_i         ! [-]
    COMPLEX(r8) :: ep_eff       ! [-]
    COMPLEX(r8) :: fracy        ! [-]

    ! free-space wave number at the current frequency (1/m)
    k0 = 2.0_r8*pi*frequency*1.0e9_r8/C

    ! Clausius-Mossotti (y-parameter) of the particle
    yr = (ep_real - 1.0_r8)/(ep_real + 2.0_r8)
    yi = ep_imag/(ep_real + 2.0_r8)

    y = CMPLX(yr, yi, r8)
    fracy = frac*y

    ! Polder-van Santen effective medium (real part and attenuation)
    ep_r = (1.0_r8 + 2.0_r8*fracy)/(1.0_r8 - fracy)
    ep_i = 2.0_r8*fracy*y*(k0*rad)**3*(1.0_r8-frac)**4 &
           /((1.0_r8-fracy)**2*(1.0_r8+2.0_r8*frac)**2)
    ep_eff = ep_r - CMPLX(0.0_r8,1.0_r8,r8)*ep_i

    if (AIMAG(ep_eff) >= 0.0_r8) ep_eff = CMPLX(REAL(ep_eff), -0.0001_r8, r8)

  END SUBROUTINE Snow_Diel

!-----------------------------------------------------------------------
  SUBROUTINE Canopy_Diel(frequency,mg,esv,rhoveg,t_soil)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of vegetation leaves using the
!   dual-dispersion model of Ulaby & El-Rayes (1987) at room
!   temperature, and the Mätzler (1994) low-temperature formulation
!   when the leaf temperature is below T0 + TS_LOW_THRESHOLD.
!
! INPUT:
!   frequency  frequency (GHz)
!   mg         gravimetric water content
!   rhoveg     vegetation specific density (kg/m3)
!   t_soil     soil/leaf temperature (K)
!
! OUTPUT:
!   esv        dielectric constant of the leaves
!
! REFERENCES:
!   [1] Ulaby, F. T., & El-Rayes, M. A. (1987), Microwave dielectric
!       spectrum of vegetation - Part II: Dual-dispersion model, IEEE
!       Trans. Geosci. Remote Sens., GE-25(5), 550-557.
!   [2] Matzler, C. (1994), Microwave (1-100 GHz) dielectric model of
!       leaves, IEEE Trans. Geosci. Remote Sens., 32(5), 947-949.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: t_soil          ! [K]
    REAL(r8) :: mg              ! [-]
    REAL(r8) :: md              ! [-]
    REAL(r8) :: en              ! [-]
    REAL(r8) :: vf              ! [-]
    REAL(r8) :: vb              ! [-]
    REAL(r8) :: rhoveg          ! [kg/m3]
    REAL(r8) :: rhoveg_gcm3     ! [g/cm3]
    REAL(r8) :: vmv             ! [-]
    REAL(r8) :: delta           ! [-]
    COMPLEX(r8) :: esv          ! [-]
    COMPLEX(r8) :: xx           ! [-]
    COMPLEX(r8) :: ef           ! [-]

    delta = 1.27_r8
    md = 1.0_r8-mg
    rhoveg_gcm3 = rhoveg / 1000.0_r8
    vmv = mg*rhoveg_gcm3/( 1.0_r8 - mg*(1.0_r8-rhoveg_gcm3) )
    xx = CMPLX(0.0_r8,1.0_r8,r8)

    ! dielectric constant model at room temperature
    if (t_soil >= tfrz + TS_LOW_THRESHOLD) then
      en = 1.7_r8 + (3.2_r8 + 6.5_r8*vmv)*vmv
      vf = vmv*(0.82_r8*vmv + 0.166_r8)
      vb = 31.4_r8*vmv*vmv/( 1.0_r8 + 59.5_r8*vmv*vmv)
      esv = en + vf*(4.9_r8 + 75.0_r8/(1.0_r8 + xx*frequency/18.0_r8) &
            -xx*(18.0_r8*delta/frequency)) + &
            vb*(2.9_r8 + 55.0_r8/(1.0_r8 + SQRT(xx*frequency/0.18_r8)))
    else
      ! dielectric constant model at low temperature (Matzler et al. 1994)
      ef = 4.9_r8 + 75.0_r8/(1.0_r8 + xx*frequency/18.0_r8) &
           - xx*(18.0_r8*delta/frequency)
      esv = 0.522_r8*(1.0_r8-1.32_r8*md)*ef + 0.51_r8 + 3.84_r8*md
    endif

    if (AIMAG(esv) >= 0.0_r8) esv = CMPLX(REAL(esv), -0.0001_r8, r8)

  END SUBROUTINE Canopy_Diel

!-----------------------------------------------------------------------
  SUBROUTINE Smooth_Reflectance(em1, em2, theta_i, theta_t, rv, rh)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the Fresnel reflectivity of a smooth planar interface
!   between two dielectric media.
!
! INPUT:
!   em1      dielectric constant of medium 1
!   em2      dielectric constant of medium 2
!   theta_i  incident angle (rad)
!   theta_t  transmitted angle (rad)
!
! OUTPUT:
!   rv       reflectivity at vertical polarization
!   rh       reflectivity at horizontal polarization
!
! REFERENCES: (Fresnel equations, standard textbook)
!-----------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE
    REAL(r8) :: theta_i         ! [rad] incident angle
    REAL(r8) :: theta_t         ! [rad] transmitted angle
    REAL(r8) :: rh              ! [-]
    REAL(r8) :: rv              ! [-]
    REAL(r8) :: cos_i           ! [-]
    REAL(r8) :: cos_t           ! [-]
    COMPLEX(r8) :: em1          ! [-]
    COMPLEX(r8) :: em2          ! [-]
    COMPLEX(r8) :: m1           ! [-]
    COMPLEX(r8) :: m2           ! [-]
    COMPLEX(r8) :: angle_i      ! [-]
    COMPLEX(r8) :: angle_t      ! [-]

    ! compute the refractive index ratio between medium 2 and 1
    ! using dielectric constant (n = SQRT(e))
    cos_i = COS(theta_i)
    cos_t = COS(theta_t)

    angle_i = CMPLX(cos_i, 0.0_r8, r8)
    angle_t = CMPLX(cos_t, 0.0_r8, r8)

    m1 = SQRT(em1)
    m2 = SQRT(em2)

    rv = (ABS((m1*angle_t-m2*angle_i)/(m1*angle_t+m2*angle_i)))**2
    rh = (ABS((m1*angle_i-m2*angle_t)/(m1*angle_i+m2*angle_t)))**2

  END SUBROUTINE Smooth_Reflectance

!-----------------------------------------------------------------------
  SUBROUTINE Transmittance(em1,em2,theta_i,theta_t,tv,th)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the Fresnel transmittance of a smooth planar interface
!   between two dielectric media.
!
! INPUT:
!   em1      dielectric constant of medium 1
!   em2      dielectric constant of medium 2
!   theta_i  incident angle (rad)
!   theta_t  transmitted angle (rad)
!
! OUTPUT:
!   tv       transmittance at vertical polarization
!   th       transmittance at horizontal polarization
!
! REFERENCES: (Fresnel equations, standard textbook)
!-----------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE
    REAL(r8) :: theta_i         ! [rad] incident angle
    REAL(r8) :: theta_t         ! [rad] transmitted angle
    REAL(r8) :: th              ! [-]
    REAL(r8) :: tv              ! [-]
    REAL(r8) :: rr              ! [-]
    REAL(r8) :: cos_i           ! [-]
    REAL(r8) :: cos_t           ! [-]
    COMPLEX(r8) :: em1          ! [-]
    COMPLEX(r8) :: em2          ! [-]
    COMPLEX(r8) :: m1           ! [-]
    COMPLEX(r8) :: m2           ! [-]
    COMPLEX(r8) :: angle_i      ! [-]
    COMPLEX(r8) :: angle_t      ! [-]

    cos_i = COS(theta_i)
    cos_t = COS(theta_t)

    angle_i = CMPLX(cos_i, 0.0_r8, r8)
    angle_t = CMPLX(cos_t, 0.0_r8, r8)

    m1 = SQRT(em1)
    m2 = SQRT(em2)

    rr = ABS(m2/m1)*cos_t/cos_i
    tv = rr*(ABS(2.0_r8*m1*angle_i/(m1*angle_t + m2*angle_i)))**2
    th = rr*(ABS(2.0_r8*m1*angle_i/(m1*angle_i + m2*angle_t)))**2

  END SUBROUTINE Transmittance

!-----------------------------------------------------------------------
  SUBROUTINE Roughness_Reflectance_ChenWeng(frequency,sigma,theta, rv,rh)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough-surface reflectivity using the Chen & Weng
!   (2016) model, which combines a polarization-mixing factor Q and
!   a roughness attenuation factor P.
!
! INPUT:
!   frequency  frequency (GHz)
!   sigma      standard deviation of the rough-surface height (m)
!   theta      local zenith angle (rad)
!   rv         smooth-surface reflectivity at V polarization (input)
!   rh         smooth-surface reflectivity at H polarization (input)
!
! OUTPUT:
!   rv         rough-surface reflectivity at V polarization
!   rh         rough-surface reflectivity at H polarization
!
! REFERENCES:
!   [1] Chen, M., & Weng, F. (2016), Modeling land surface roughness effect
!       on soil microwave emission in community surface emissivity model,
!       IEEE Trans. Geosci. Remote Sens., 54(3), 1716-1726.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    IMPLICIT NONE
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: sigma           ! [m]
    REAL(r8) :: theta           ! [rad]
    REAL(r8) :: f               ! [Hz]
    REAL(r8) :: lamda           ! [cm]
    REAL(r8) :: k               ! [1/cm]
    REAL(r8) :: a1              ! [-]
    REAL(r8) :: a2              ! [-]
    REAL(r8) :: a3              ! [-]
    REAL(r8) :: k1              ! [1/cm]
    REAL(r8) :: k2              ! [1/cm]
    REAL(r8) :: sigma1          ! [1/cm]
    REAL(r8) :: sigma2          ! [1/cm]
    REAL(r8) :: b1              ! [-]
    REAL(r8) :: b2              ! [-]
    REAL(r8) :: temp1           ! [-]
    REAL(r8) :: temp2           ! [-]
    REAL(r8) :: temp1_tan       ! [-]
    REAL(r8) :: temp2_tan       ! [-]
    REAL(r8) :: Q               ! [-]
    REAL(r8) :: P               ! [-]
    REAL(r8) :: rh              ! [-]
    REAL(r8) :: rv              ! [-]
    REAL(r8) :: rh_s            ! [-]
    REAL(r8) :: rv_s            ! [-]
    REAL(r8) :: sigma_cm        ! [cm]

    ! convert to CGS units
    sigma_cm = sigma * 100.0_r8
    f = frequency*1.0e9_r8
    lamda = 3.0e8_r8/f*100.0_r8
    k = 2.0_r8*pi/lamda

    ! P-factor fitting coefficients
    a1 = 1.4_r8
    a2 = 0.15_r8
    a3 = EXP(-(2.0_r8*pi*sigma_cm/5.0_r8)**sqrt(0.1_r8*cos(theta)))

    k1 = 0.0_r8
    k2 = 2.0_r8*pi/3.0_r8
    sigma1 = 2.0_r8*pi*0.3_r8
    sigma2 = 2.0_r8*pi*0.5_r8
    b1 = 0.43_r8
    b2 = 15.0_r8

    ! polarization-mixing factor Q
    Q = b1*(1.0_r8 - EXP(-b2*k*sigma_cm*(COS(theta))**2.0_r8))

    ! roughness attenuation factor P (piecewise-tanh fit)
    temp1 = (k*sigma_cm*COS(theta) - k1*sigma_cm*COS(theta))/sigma1
    temp2 = (k*sigma_cm*COS(theta) - k2*sigma_cm*COS(theta))/sigma2
    temp1_tan = (EXP(2.0_r8*temp1) - 1.0_r8)/(EXP(2.0_r8*temp1)+1.0_r8)
    temp2_tan = (EXP(2.0_r8*temp2) - 1.0_r8)/(EXP(2.0_r8*temp2)+1.0_r8)

    P = (a1+a3)/2.0_r8+(a2-a1)/2.0_r8*temp1_tan + (a3-a2)/2.0_r8*temp2_tan
    rv_s = rv
    rh_s = rh
    rv = ((1.0_r8-Q)*rv_s + Q*rh_s)*P
    rh = ((1.0_r8-Q)*rh_s + Q*rv_s)*P

  END SUBROUTINE Roughness_Reflectance_ChenWeng

!-----------------------------------------------------------------------
  SUBROUTINE Roughness_Reflectance_Qh(frequency,sigma,rv,rh)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough-surface reflectivity using the Qh
!   parameterization with Q = 0.35 (1 - exp(-0.60 f sigma^2)).
!
! INPUT:
!   frequency  frequency (GHz)
!   sigma      standard deviation of the rough-surface height (m)
!   rv         smooth-surface reflectivity at V polarization (input)
!   rh         smooth-surface reflectivity at H polarization (input)
!
! OUTPUT:
!   rv         rough-surface reflectivity at V polarization
!   rh         rough-surface reflectivity at H polarization
!
! REFERENCES:
!   [1] Wigneron, J.-P., et al. (2007), L-band Microwave Emission of the
!       Biosphere (L-MEB) model, Remote Sensing of Environment, 107,
!       639-655.
!   [2] Wang, J., & Choudhury, B. J. (1992), Passive microwave radiation
!       from soil: examples..., in Passive Microwave Remote Sensing of
!       Land-Atmosphere Interactions (eds. Choudhury et al.), VSP.
!-----------------------------------------------------------------------
    USE MOD_Precision
    IMPLICIT NONE
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: q               ! [-]
    REAL(r8) :: rh              ! [-]
    REAL(r8) :: rv              ! [-]
    REAL(r8) :: rh_s            ! [-]
    REAL(r8) :: rv_s            ! [-]
    REAL(r8) :: sigma           ! [m]
    REAL(r8) :: sigma_cm        ! [cm]

    sigma_cm = sigma * 100.0_r8
    rh_s = 0.3_r8*rh
    rv_s = 0.3_r8*rv
    q = 0.35_r8*(1.0_r8 - EXP(-0.60_r8*frequency*sigma_cm**2.0_r8))
    rh = rh_s + q*(rv_s - rh_s)
    rv = rv_s + q*(rh_s - rv_s)

  END SUBROUTINE Roughness_Reflectance_Qh

!-----------------------------------------------------------------------
  SUBROUTINE Canopy_Optic(vlai,frequency,theta,esv,d,gv,gh,&
                          ssalb_v,ssalb_h,tau_v, tau_h)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the bulk optical properties of a vegetation canopy
!   (single-scattering albedo, asymmetry factor, and optical depth)
!   from the leaf dielectric constant.
!
! INPUT:
!   vlai       vegetation-adjusted LAI (-)
!   frequency  frequency (GHz)
!   theta      incidence angle (rad)
!   esv        dielectric constant of the leaves
!   d          leaf thickness (m)
!
! OUTPUT:
!   gv,gh      asymmetry factor at V/H polarization
!   ssalb_v,h  single-scattering albedo at V/H polarization
!   tau_v,h    optical depth at V/H polarization
!
! REFERENCES:
!   [1] Weng, F., Yan, B., & Grody, N. C. (2001), A microwave land emissivity
!       model, J. Geophys. Res. Atmos., 106(D17), 20115-20123.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    IMPLICIT NONE
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: theta           ! [rad]
    REAL(r8) :: d               ! [m]
    REAL(r8) :: vlai            ! [m2/m2]
    REAL(r8) :: ssalb_v         ! [-]
    REAL(r8) :: ssalb_h         ! [-]
    REAL(r8) :: tau_v           ! [-]
    REAL(r8) :: tau_h           ! [-]
    REAL(r8) :: gv              ! [-]
    REAL(r8) :: gh              ! [-]
    REAL(r8) :: mu              ! [-]
    COMPLEX(r8) :: ix           ! [-]
    COMPLEX(r8) :: k0           ! [1/m]
    COMPLEX(r8) :: kz0          ! [1/m]
    COMPLEX(r8) :: kz1          ! [1/m]
    COMPLEX(r8) :: rhc          ! [-]
    COMPLEX(r8) :: rvc          ! [-]
    COMPLEX(r8) :: esv          ! [-]
    COMPLEX(r8) :: expval1      ! [-]
    COMPLEX(r8) :: factt        ! [-]
    COMPLEX(r8) :: factrvc      ! [-]
    COMPLEX(r8) :: factrhc      ! [-]
    REAL(r8) :: rh              ! [-]
    REAL(r8) :: rv              ! [-]
    REAL(r8) :: th              ! [-]
    REAL(r8) :: tv              ! [-]
    REAL(r8), PARAMETER :: threshold = 0.999_r8   ! [-]

    mu = COS(theta)
    ix = CMPLX(0.0_r8, 1.0_r8, r8)

    ! wave numbers in air (kz0) and inside the leaf (kz1)
    k0  = CMPLX(2.0_r8*pi*frequency*1.0e9_r8/C, 0.0_r8, r8)
    kz0 = k0*mu
    kz1 = k0*SQRT(esv - SIN(theta)**2)

    ! Fresnel reflection coefficients of a single leaf
    rhc = (kz0 - kz1)/(kz0 + kz1)
    rvc = (esv*kz0 - kz1)/(esv*kz0 + kz1)

    ! phase factor and transmission factor through the leaf
    expval1 = EXP(-2.0_r8*ix*kz1*d)
    factrvc = 1.0_r8 - rvc**2*expval1
    factrhc = 1.0_r8 - rhc**2*expval1
    factt   = 4.0_r8*kz0*kz1*EXP(ix*(kz0-kz1)*d)

    ! leaf reflectivity and transmittance
    rv = ABS(rvc*(1.0_r8 - expval1)/factrvc)**2
    rh = ABS(rhc*(1.0_r8 - expval1)/factrhc)**2

    th = ABS(factt/((kz1+kz0)**2*factrhc))**2
    tv = ABS(esv*factt/((kz1 + esv*kz0)**2*factrvc))**2

    ! asymmetry factor and canopy optical depth
    gv = 0.5_r8
    gh = 0.5_r8

    tau_v = 0.5_r8*vlai*(2.0_r8 - tv-th)
    tau_h = tau_v

    ! single-scattering albedo of the canopy
    ssalb_v = MIN((rv + rh)/(2.0_r8 -tv -th),threshold)
    ssalb_h = ssalb_v

  END SUBROUTINE Canopy_Optic

!-----------------------------------------------------------------------
  SUBROUTINE Snow_Optic(frequency,a,h,f,ep_real,ep_imag,gv,gh, &
                        ssalb_v,ssalb_h,tau_v,tau_h)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the bulk optical properties of a snow layer
!   (single-scattering albedo, asymmetry factor, and optical depth)
!   from the snow grain dielectric constant and volume fraction.
!
! INPUT:
!   frequency  frequency (GHz)
!   a          particle radius (m)
!   h          snow depth (m)
!   f          volume fraction of ice particles (0-1)
!   ep_real    real part of the particle dielectric constant
!   ep_imag    imaginary part of the particle dielectric constant
!
! OUTPUT:
!   gv,gh      asymmetry factor at V/H polarization
!   ssalb_v,h  single-scattering albedo at V/H polarization
!   tau_v,h    optical depth at V/H polarization
!
! REFERENCES:
!   [1] Weng, F., Yan, B., & Grody, N. C. (2001), A microwave land emissivity
!       model, J. Geophys. Res. Atmos., 106(D17), 20115-20123.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    IMPLICIT NONE
    REAL(r8) :: yr              ! [-]
    REAL(r8) :: yi              ! [-]
    REAL(r8) :: ep_real         ! [-]
    REAL(r8) :: ep_imag         ! [-]
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: a               ! [m]
    REAL(r8) :: h               ! [m]
    REAL(r8) :: f               ! [-]
    REAL(r8) :: ssalb_v         ! [-]
    REAL(r8) :: ssalb_h         ! [-]
    REAL(r8) :: tau_v           ! [-]
    REAL(r8) :: tau_h           ! [-]
    REAL(r8) :: gv              ! [-]
    REAL(r8) :: gh              ! [-]
    REAL(r8) :: k               ! [1/m]
    REAL(r8) :: ks1             ! [-]
    REAL(r8) :: ks2             ! [-]
    REAL(r8) :: ks3             ! [-]
    REAL(r8) :: ks              ! [1/m]
    REAL(r8) :: kr1             ! [-]
    REAL(r8) :: kr2             ! [-]
    REAL(r8) :: kr3             ! [-]
    REAL(r8) :: kr              ! [1/m]
    REAL(r8) :: ki1             ! [-]
    REAL(r8) :: ki2             ! [-]
    REAL(r8) :: ki3             ! [-]
    REAL(r8) :: ki              ! [1/m]
    REAL(r8) :: fact1           ! [-]
    REAL(r8) :: fact2           ! [-]
    REAL(r8) :: fact3           ! [-]
    REAL(r8) :: fact4           ! [-]
    REAL(r8) :: fact5           ! [-]

    k = 2.0_r8*pi*frequency*1.0e9_r8/C

    ! Clausius-Mossotti (y-parameter) of the snow grain
    yr = (ep_real - 1.0_r8)/(ep_real + 2.0_r8)
    yi = -ep_imag/(ep_real + 2.0_r8)

    ! auxiliary factors for the dense-medium scattering coefficients
    fact1 = (1.0_r8+2.0_r8*f)**2
    fact2 = 1.0_r8-f*yr
    fact3 = (1.0_r8-f)**4
    fact4 = f*(k*a)**3
    fact5 = 1.0_r8 + 2.0_r8*f*yr

    ! scattering coefficient
    ks1 = k*SQRT(fact2/fact5)
    ks2 = fact4*fact3/fact1
    ks3 = (yr/fact2)**2
    ks  = ks1*ks2*ks3

    ! real part of the effective propagation coefficient
    kr1 = fact5/fact2
    kr2 = 2.0_r8*ks2
    kr3 = 2.0_r8*yi*yr/(fact2**3)
    kr  = k*SQRT(kr1 + kr2*kr3)

    ! imaginary part of the effective propagation coefficient
    ki1 = 3.0_r8*f*yi/fact2**2
    ki2 = kr2
    ki3 = ks3
    ki  = k**2/(2.0_r8*kr)*(ki1 + ki2*ki3)

    ! asymmetry factor and optical depth
    gv = 0.5_r8
    gh = 0.5_r8

    ssalb_v = MIN(ks/ki, 0.999_r8)
    ssalb_h = ssalb_v
    tau_v = 2.0_r8*ki*h
    tau_h = tau_v

  END SUBROUTINE Snow_Optic

!-----------------------------------------------------------------------
  SUBROUTINE Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
        r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,esv,esh,frequency,t_soil,t_skin)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Two-stream solution for the emissivity of a single scattering layer
!   bounded by two Fresnel interfaces. The formulation follows the
!   total-up-welling-radiance form of Weng et al. (2001).
!
! INPUT:
!   mu                   cos(incidence angle)
!   gv,gh                asymmetry factor at V/H polarization
!   ssalb_v,ssalb_h      single-scattering albedo at V/H polarization
!   tau_v,tau_h          optical depth at V/H polarization
!   r21_v,r21_h          reflectivity of the upper interface
!   r23_v,r23_h          reflectivity of the lower interface
!   t21_v,t21_h          transmittance of the upper interface
!   frequency            frequency (GHz)
!   t_soil               soil temperature (K)
!   t_skin               skin temperature (K)
!
! OUTPUT:
!   esv,esh              emissivity at V/H polarization
!
! REFERENCES:
!   [1] Weng, F., Yan, B., & Grody, N. C. (2001), A microwave land emissivity
!       model, J. Geophys. Res. Atmos., 106(D17), 20115-20123.
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_Const_Physical
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8) :: mu              ! [-] cos(theta)
    REAL(r8) :: gv              ! [-] asymmetry factor V
    REAL(r8) :: gh              ! [-] asymmetry factor H
    REAL(r8) :: ssalb_h         ! [-] single-scattering albedo H
    REAL(r8) :: ssalb_v         ! [-] single-scattering albedo V
    REAL(r8) :: tau_h           ! [-] optical depth H
    REAL(r8) :: tau_v           ! [-] optical depth V
    REAL(r8) :: r21_h           ! [-] reflectivity
    REAL(r8) :: r21_v           ! [-] reflectivity
    REAL(r8) :: r23_h           ! [-] reflectivity
    REAL(r8) :: r23_v           ! [-] reflectivity
    REAL(r8) :: t21_v           ! [-] transmittance
    REAL(r8) :: t21_h           ! [-] transmittance
    REAL(r8) :: esv             ! [-] emissivity V
    REAL(r8) :: esh             ! [-] emissivity H
    REAL(r8) :: alfa_v          ! [-]
    REAL(r8) :: alfa_h          ! [-]
    REAL(r8) :: kk_h            ! [-]
    REAL(r8) :: kk_v            ! [-]
    REAL(r8) :: gamma_h         ! [-]
    REAL(r8) :: gamma_v         ! [-]
    REAL(r8) :: beta_v          ! [-]
    REAL(r8) :: beta_h          ! [-]
    REAL(r8) :: fact1           ! [-]
    REAL(r8) :: fact2           ! [-]
    REAL(r8) :: frequency       ! [GHz]
    REAL(r8) :: t_soil          ! [K]
    REAL(r8) :: t_skin          ! [K]
    REAL(r8) :: gsect0          ! [-]
    REAL(r8) :: gsect1_h        ! [-]
    REAL(r8) :: gsect1_v        ! [-]
    REAL(r8) :: gsect2_h        ! [-]
    REAL(r8) :: gsect2_v        ! [-]

    ! similarity parameters and extinction coefficients (H)
    alfa_h  = SQRT((1.0_r8 - ssalb_h)/(1.0_r8 - gh*ssalb_h))
    kk_h    = SQRT((1.0_r8 - ssalb_h)*(1.0_r8 -  gh*ssalb_h))/mu
    beta_h  = (1.0_r8 - alfa_h)/(1.0_r8 + alfa_h)
    gamma_h = (beta_h -r23_h)/(1.0_r8-beta_h*r23_h)

    ! similarity parameters and extinction coefficients (V)
    alfa_v  = SQRT((1.0_r8-ssalb_v)/(1.0_r8 - gv*ssalb_v))
    kk_v    = SQRT((1.0_r8-ssalb_v)*(1.0_r8 - gv*ssalb_v))/mu
    beta_v  = (1.0_r8 - alfa_v)/(1.0_r8 + alfa_v)
    gamma_v = (beta_v - r23_v)/(1.0_r8 - beta_v*r23_v)

    ! layer effective reflection factors
    fact1 = gamma_h*EXP(-2.0_r8*kk_h*tau_h)
    fact2 = gamma_v*EXP(-2.0_r8*kk_v*tau_v)

    ! Planck-ratio between layer and soil temperature
    gsect0   = (EXP(C_2*frequency/t_skin) - 1.0_r8) / (EXP(C_2*frequency/t_soil) - 1.0_r8)
    gsect1_h = (1.0_r8 - r23_h)*(gsect0 - 1.0_r8)
    gsect2_h = ((1.0_r8 - beta_h*beta_h)/(1.0_r8 - beta_h*r23_h))*EXP(-kk_h*tau_h)

    gsect1_v = (1.0_r8 - r23_v)*(gsect0 - 1.0_r8)
    gsect2_v = ((1.0_r8 - beta_v*beta_v)/(1.0_r8 - beta_v*r23_v))*EXP(-kk_v*tau_v)

    ! total up-welling emissivity at H and V polarization
    esh = t21_h*((1.0_r8 - beta_h)*(1.0_r8 + fact1) + gsect1_h*gsect2_h) &
          /(1.0_r8 - beta_h*r21_h - (beta_h - r21_h)*fact1)
    esv = t21_v*((1.0_r8 - beta_v)*(1.0_r8 + fact2) + gsect1_v*gsect2_v) &
          /(1.0_r8 - beta_v*r21_v - (beta_v - r21_v)*fact2)

    ! physical bounds (emissivity must be in [0, 1])
    if (esh < EMISSH_DEFAULT) esh = EMISSH_DEFAULT
    if (esv < EMISSV_DEFAULT) esv = EMISSV_DEFAULT
    if (esh > 1.0_r8) esh = 1.0_r8
    if (esv > 1.0_r8) esv = 1.0_r8

  END SUBROUTINE Two_Stream_Solution

!-----------------------------------------------------------------------
  SUBROUTINE SnowEM_Default(frequency,ts, depth,Emissivity_V,Emissivity_H)
!-----------------------------------------------------------------------
! DESCRIPTION:
!   Preliminary estimate of snow emissivity based on the surface
!   temperature and snow depth, using the tabulated AMSR-E snow
!   emissivity spectra (in ARMS_Common_Parameters).
!
! INPUT:
!   frequency  frequency (GHz)
!   ts         surface/skin temperature (K)
!   depth      snow depth (m)
!
! OUTPUT:
!   Emissivity_H  snow emissivity at H polarization
!   Emissivity_V  snow emissivity at V polarization
!
! REFERENCES:
!   AMSR-E ATBD (tabulated coefficients in ARMS_Common_Parameters:
!   FREQUENCY_AMSRE, GRASS_AFTER_SNOW_*_AMSRE, POWDER_SNOW_*_AMSRE,
!   WET_SNOW_*_AMSRE, DEEP_SNOW_*_AMSRE).
!-----------------------------------------------------------------------
    USE MOD_Precision
    USE MOD_DA_Const
    IMPLICIT NONE
    REAL(r8) :: frequency                             ! [GHz]
    REAL(r8) :: ts                                    ! [K]
    REAL(r8) :: depth                                 ! [m]
    REAL(r8) :: Emissivity_V                          ! [-]
    REAL(r8) :: Emissivity_H                          ! [-]
    INTEGER , PARAMETER :: new = N_AMSRE              ! [-]
    INTEGER , PARAMETER :: NFRESH_SHALLOW_SNOW = 1    ! [-]
    INTEGER , PARAMETER :: NPOWDER_SNOW        = 2    ! [-]
    INTEGER , PARAMETER :: NWET_SNOW           = 3    ! [-]
    INTEGER , PARAMETER :: NDEEP_SNOW          = 4    ! [-]
    REAL(r8), PARAMETER :: twet    = 270.0_r8         ! [K] wet-snow temperature threshold
    REAL(r8), PARAMETER :: tcrust  = 235.0_r8         ! [K] crust-snow temperature threshold
    REAL(r8), PARAMETER :: depth_s =  0.05_r8         ! [m] shallow-snow depth threshold
    REAL(r8), PARAMETER :: depth_c =  0.10_r8         ! [m] deep-snow depth threshold
    INTEGER :: ich                                    ! [-] channel index
    INTEGER :: basic_snow_type                        ! [-] snow type index
    REAL(r8), DIMENSION(new) :: ev                    ! [-] emissivity
    REAL(r8), DIMENSION(new) :: eh                    ! [-] emissivity
    REAL(r8), DIMENSION(new) :: freq                  ! [GHz]
    REAL(r8) :: df                                    ! [GHz]
    REAL(r8) :: df0                                   ! [GHz]

    freq = FREQUENCY_AMSRE(1:new)

!#############################################################################
! 1. Classify the snow type from temperature and depth
!#############################################################################
    basic_snow_type = NFRESH_SHALLOW_SNOW
    if (ts >= twet .and. depth <= depth_s) then
      basic_snow_type = NWET_SNOW
    else
      if (depth <= depth_s) then
        basic_snow_type = NFRESH_SHALLOW_SNOW
      else
        basic_snow_type = NPOWDER_SNOW
      endif
    endif
    if (ts <= tcrust .and. depth >= depth_c) basic_snow_type = NDEEP_SNOW

!#############################################################################
! 2. Assign the emissivity spectrum (AMSR-E lookup)
!#############################################################################
    SELECT CASE (basic_snow_type)
      CASE (NFRESH_SHALLOW_SNOW)
        ev = GRASS_AFTER_SNOW_EV_AMSRE(1:new)
        eh = GRASS_AFTER_SNOW_EH_AMSRE(1:new)
      CASE (NPOWDER_SNOW)
        ev = POWDER_SNOW_EV_AMSRE(1:new)
        eh = POWDER_SNOW_EH_AMSRE(1:new)
      CASE (NWET_SNOW)
        ev = WET_SNOW_EV_AMSRE(1:new)
        eh = WET_SNOW_EH_AMSRE(1:new)
      CASE (NDEEP_SNOW)
        ev = DEEP_SNOW_EV_AMSRE(1:new)
        eh = DEEP_SNOW_EH_AMSRE(1:new)
    END SELECT

!#############################################################################
! 3. Extrapolate or interpolate to the requested frequency
!#############################################################################
    ! below the lowest tabulated channel: use the first value
    if (frequency <= freq(1)) then
      Emissivity_H = eh(1)
      Emissivity_V = ev(1)
      return
    end if
    ! above the highest tabulated channel: use the last value
    if (frequency >= freq(new)) then
      Emissivity_H = eh(new)
      Emissivity_V = ev(new)
      return
    end if

    ! linear interpolation between AMSR-E channels
    Channel_loop: do ich = 2, new
      if (frequency <= freq(ich)) then
        df  = frequency - freq(ich-1)
        df0 = freq(ich) - freq(ich-1)
        Emissivity_H = eh(ich-1) + (df*(eh(ich) - eh(ich-1))/df0)
        Emissivity_V = ev(ich-1) + (df*(ev(ich) - ev(ich-1))/df0)
        exit Channel_loop
      end if
    end do Channel_loop

  END SUBROUTINE SnowEM_Default

!-----------------------------------------------------------------------
END MODULE MOD_DA_RTM_Physics
#endif