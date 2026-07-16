!##############################################################################
! Description:
! Module to compute the surface optical properties for Land��SNOW and ICE surfaces at
! microwave frequencies required for determining the surface
! contribution to the radiative transfer.
!
! Current Code Owner: ARMS Team
!
! History:
! Version          Date           Comment             
! 1.0              2020/12/30     Original
! 1.1              2021/08/04     Add solver HRTS         
! 1.2              2022/04/06     Optimize code
! 2.0              2023/04/04     Improve relevant algorithms
!
! Main Subroutine/Function:       
!     MW_Land_SfcOptics
!     MW_Land_SfcOptics_TL
!     MW_Land_SfcOptics_AD
!     MW_Snow_SfcOptics
!     MW_Snow_SfcOptics_TL
!     MW_Snow_SfcOptics_AD
!     MW_Ice_SfcOptics
!     MW_Ice_SfcOptics_TL
!     MW_Ice_SfcOptics_AD
!
!############################################################################## 

MODULE ARMS_SfcEM_MWLand

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE ARMS_Common_Basic  
  USE ARMS_Common_Parameters,     ONLY: ZERO, ONE,TWO,THREE,FOUR,&
                                        PI_AUS,TWOPI,&
                                        T0,C_2, &
                                        FREQUENCY_AMSRE, &
                                        GRASS_AFTER_SNOW_EV_AMSRE, &
                                        GRASS_AFTER_SNOW_EH_AMSRE, &
                                        POWDER_SNOW_EV_AMSRE, &
                                        POWDER_SNOW_EH_AMSRE, &
                                        WET_SNOW_EV_AMSRE, &
                                        WET_SNOW_EH_AMSRE, &
                                        DEEP_SNOW_EV_AMSRE, &
                                        DEEP_SNOW_EH_AMSRE
  USE ARMS_SpcCoeff,               ONLY: SC
  USE ARMS_Surface_Setting,        ONLY: Surface_ARMS
  USE ARMS_GeometryInfo_Setting,   ONLY: GeometryInfo_ARMS, &
                                        ARMS_GeometryInfo_GetValue
  USE ARMS_SfcOptics_Setting,      ONLY: SfcOptics_ARMS
  ! Disable implicit typing
  IMPLICIT NONE

  ! Everything private by default
  PRIVATE
  PUBLIC :: MW_Land_SfcOptics
  PUBLIC :: MW_Land_SfcOptics_TL
  PUBLIC :: MW_Land_SfcOptics_AD
  ! Science routines
  PUBLIC :: MW_Snow_SfcOptics
  PUBLIC :: MW_Snow_SfcOptics_TL
  PUBLIC :: MW_Snow_SfcOptics_AD
  PUBLIC :: MW_Ice_SfcOptics
  PUBLIC :: MW_Ice_SfcOptics_TL
  PUBLIC :: MW_Ice_SfcOptics_AD
  ! Subroutine and Function
  PUBLIC :: LandEM
  PUBLIC :: Soil_Diel_Mironov
  PUBLIC :: Soil_Diel_Dobson
  PUBLIC :: Snow_Diel
  PUBLIC :: Smooth_Reflectance
  PUBLIC :: Transmittance
  PUBLIC :: Roughness_Reflectance_ChenWeng
  PUBLIC :: Roughness_Reflectance_Qh
  PUBLIC :: Canopy_Diel
  PUBLIC :: Canopy_Optic
  PUBLIC :: Snow_Optic
  PUBLIC :: SnowEM_Default
  PUBLIC :: Two_Stream_Solution
  
  ! -----------------
  ! Module parameters
  ! -----------------

  ! Message length
  INTEGER, PARAMETER :: ML = 256
  ! Valid type indices for the microwave land emissivity model
  ! ...The soil types
  INTEGER, PARAMETER :: N_VALID_SOIL_TYPES = 8
  INTEGER, PARAMETER :: INVALID_SOIL    =  0
  INTEGER, PARAMETER :: COARSE          =  1
  INTEGER, PARAMETER :: MEDIUM          =  2
  INTEGER, PARAMETER :: FINE            =  3
  INTEGER, PARAMETER :: COARSE_MEDIUM   =  4
  INTEGER, PARAMETER :: COARSE_FINE     =  5
  INTEGER, PARAMETER :: MEDIUM_FINE     =  6
  INTEGER, PARAMETER :: COARSE_MED_FINE =  7
  INTEGER, PARAMETER :: ORGANIC         =  8
  ! ...The vegetation types
  INTEGER, PARAMETER :: N_VALID_VEGETATION_TYPES       = 12
  INTEGER, PARAMETER :: INVALID_VEGETATION             =  0
  INTEGER, PARAMETER :: BROADLEAF_EVERGREEN_TREES      =  1
  INTEGER, PARAMETER :: BROADLEAF_DECIDUOUS_TREES      =  2
  INTEGER, PARAMETER :: BROADLEAF_NEEDLELEAF_TREES     =  3
  INTEGER, PARAMETER :: NEEDLELEAF_EVERGREEN_TREES     =  4
  INTEGER, PARAMETER :: NEEDLELEAF_DECIDUOUS_TREES     =  5
  INTEGER, PARAMETER :: BROADLEAF_TREES_GROUNDCOVER    =  6
  INTEGER, PARAMETER :: GROUNDCOVER                    =  7
  INTEGER, PARAMETER :: GROADLEAF_SHRUBS_GROUNDCOVER   =  8
  INTEGER, PARAMETER :: BROADLEAF_SHRUBS_BARE_SOIL     =  9
  INTEGER, PARAMETER :: DWARF_TREES_SHRUBS_GROUNDCOVER = 10
  INTEGER, PARAMETER :: BARE_SOIL                      = 11
  INTEGER, PARAMETER :: CULTIVATIONS                   = 12
  !LANDEM
  REAL(AUS), PARAMETER :: POINT1 = 0.1_AUS
  REAL(AUS), PARAMETER :: POINT5 = 0.5_AUS
  REAL(AUS), PARAMETER :: EMISSH_DEFAULT = 0.25_AUS
  REAL(AUS), PARAMETER :: EMISSV_DEFAULT = 0.30_AUS
  REAL(AUS), PARAMETER :: ONE_TENTH = POINT1
  REAL(AUS), PARAMETER :: HALF      = POINT5

  LOGICAL, PARAMETER :: dbg = .TRUE.
  
  !LOW TEMPERATURE and FROZEN SOIL
  REAL(AUS), PARAMETER :: TS_LOW_THRESHOLD  =  5.0_AUS    !��
  REAL(AUS), PARAMETER :: TS_FROZEN_THRESHOLD  =  0.0_AUS    !��
  REAL(AUS), PARAMETER :: RATIO_FROZEN  =  0.3_AUS    
CONTAINS


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Land_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a land surface.
!
! INPUT:
!       Surface:        Structure containing the surface emissivity data.
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!
! OUTPUT:
!       SfcOptics:      Structure containing the surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Land_SfcOptics( &
    Surface     , &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   ) &  ! Output
  RESULT ( err_stat )
    ! Arguments
    TYPE(Surface_ARMS),      INTENT(IN)     :: Surface
    INTEGER,                 INTENT(IN)     :: SensorIndex
    INTEGER,                 INTENT(IN)     :: ChannelIndex
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Land_SfcOptics'
    REAL(AUS),     PARAMETER :: FREQUENCY_CUTOFF   = 80.0_AUS  ! GHz
    REAL(AUS),     PARAMETER :: DEFAULT_EMISSIVITY = 0.95_AUS
    REAL(AUS),     PARAMETER :: EPS = 1.0E-8_AUS
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER :: i


    ! Set up
    err_stat = SUCCESS
    ! ...Check the soil type...
    IF ( Surface%Soil_Type < 1 .OR. &
         Surface%Soil_Type > N_VALID_SOIL_TYPES ) THEN
      SfcOptics%Emissivity   = ZERO
      SfcOptics%Reflectivity = ZERO
      err_stat = FAILURE
      msg = 'Invalid soil type index specified'
      CALL Msg_Print( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF
    ! ...and the vegetation type
    IF ( Surface%Vegetation_Type < 1 .OR. &
         Surface%Vegetation_Type > N_VALID_VEGETATION_TYPES ) THEN
      SfcOptics%Emissivity   = ZERO
      SfcOptics%Reflectivity = ZERO
      err_stat = FAILURE
      msg = 'Invalid vegetation type index specified'
      CALL Msg_Print( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF

    ! Compute the surface optical parameters
    IF ( SC(SensorIndex)%Frequency(ChannelIndex) < FREQUENCY_CUTOFF ) THEN
      ! Frequency is low enough for the model
      DO i = 1, SfcOptics%n_Angles
        CALL        LandEM(SfcOptics%Angle(i),            & ! Input, Degree
                           SC(SensorIndex)%Frequency(ChannelIndex),   & ! Input, GHz
                           Surface%Soil_Moisture_Content, & ! Input, g.cm**-3
                           Surface%Vegetation_Fraction,   & ! Input
                           Surface%Soil_Temperature,      & ! Input, K
                           Surface%Land_Temperature,      & ! Input, K
                           Surface%Lai,                   & ! Input, Leaf Area Index
                           Surface%Soil_Type,             & ! Input, Soil Type (1 -  9)
                           Surface%Vegetation_Type,       & ! Input, Vegetation Type (1 - 13)
                           ZERO,                          & ! Input, Snow depth, mm
                           SfcOptics%Emissivity(i,2),     & ! Output, H component
                           SfcOptics%Emissivity(i,1)      ) ! Output, V component
         IF ( SfcOptics%Emissivity(i,2) <= ZERO+EPS .OR. &
             SfcOptics%Emissivity(i,2) >= ONE -EPS .OR. &
             SfcOptics%Emissivity(i,1) <= ZERO+EPS .OR. &
             SfcOptics%Emissivity(i,1) >= ONE -EPS ) THEN
          SfcOptics%Emissivity   = ZERO
          SfcOptics%Reflectivity = ZERO
          err_stat = FAILURE
          WRITE(msg,'(A,1X,F7.2,A,1X,F8.3,A,1X,2F12.6)') &
            'LandEM emissivity out of (0,1): angle=', SfcOptics%Angle(i), &
            ', freq=', SC(SensorIndex)%Frequency(ChannelIndex), ' GHz,(Eh,Ev)=', &
            SfcOptics%Emissivity(i,2), SfcOptics%Emissivity(i,1)
          CALL Msg_Print( ROUTINE_NAME, msg, err_stat ); RETURN
        END IF

        ! 2) Polarization ordering check: Ev >= Eh (strict with EPS)
        IF ( SfcOptics%Emissivity(i,1) + EPS < SfcOptics%Emissivity(i,2) ) THEN
          SfcOptics%Emissivity   = ZERO
          SfcOptics%Reflectivity = ZERO
          err_stat = FAILURE
          WRITE(msg,'(A,1X,F7.2,A,1X,F8.3,A,1X,2F12.6)') &
            'LandEM physics check failed (Ev < Eh): angle=', SfcOptics%Angle(i),&
            ', freq=', SC(SensorIndex)%Frequency(ChannelIndex), ' GHz,(Eh,Ev)=', &
            SfcOptics%Emissivity(i,2), SfcOptics%Emissivity(i,1)
          CALL Msg_Print( ROUTINE_NAME, msg, err_stat ); RETURN
        END IF

        ! Assume specular surface
        SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
        SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
      END DO
    ELSE
      ! Frequency is too high for model. Use default.
      DO i = 1, SfcOptics%n_Angles
        SfcOptics%Emissivity(i,1:2)         = DEFAULT_EMISSIVITY
        SfcOptics%Reflectivity(i,1:2,i,1:2) = ONE-DEFAULT_EMISSIVITY
      END DO
    END IF
!         print*, '----------------------'
!        print*, 'SfcOptics%Angle(i)=', SfcOptics%Angle(1)
!        print*, 'Surface%Soil_Moisture_Content=', Surface%Soil_Moisture_Content
!        print*, 'SfcOptics%Emissivity(H,1)=',SfcOptics%Emissivity(1,2)
!        print*, 'SfcOptics%Emissivity(V,1)=',SfcOptics%Emissivity(1,1)
!         print*, '----------------------'
  END FUNCTION MW_Land_SfcOptics


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Land_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a land surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_TL:   Structure containing the tangent-linear surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Land_SfcOptics_TL( &
    SfcOptics_TL) &  ! TL  Output
  RESULT ( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS), INTENT(IN OUT) :: SfcOptics_TL
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Land_SfcOptics_TL'
    ! Local variables

    ! Set up
    err_stat = SUCCESS

    ! Compute the tangent-linear surface optical parameters
    ! ***No TL models yet, so default TL output is zero***
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO

  END FUNCTION MW_Land_SfcOptics_TL


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Land_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a land surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_AD:   Structure containing the adjoint surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Land_SfcOptics_AD( &
    SfcOptics_AD) &  ! AD  Input
  RESULT( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics_AD
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Land_SfcOptics_AD'
    ! Local variables


    ! Set up
    err_stat = SUCCESS


    ! Compute the adjoint surface optical parameters
    ! ***No AD models yet, so there is no impact on AD result***
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO

  END FUNCTION MW_Land_SfcOptics_AD
  
!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Snow_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! INPUT:
!       Surface:        Structure containing the surface emissivity data.
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!
! OUTPUT:
!       SfcOptics:      Structure containing the surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------


  FUNCTION MW_Snow_SfcOptics( &
    Surface     , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   ) &  ! Output
  RESULT( Error_Status )
    ! Arguments
    TYPE(Surface_ARMS),      INTENT(IN)     :: Surface
    TYPE(GeometryInfo_ARMS), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*),  PARAMETER :: ROUTINE_NAME = 'MW_Snow_SfcOptics'
    REAL(AUS), PARAMETER :: FREQUENCY_THRESHOLD            =  80.0_AUS  ! GHz
    REAL(AUS), PARAMETER :: DEFAULT_EMISSIVITY             =   0.90_AUS
    REAL(AUS), PARAMETER :: NOT_USED(4)                    = -99.9_AUS
    ! Local variables
    INTEGER :: i
    REAL(AUS) :: Sensor_Zenith_Angle
    REAL(AUS) :: Alpha


    ! Set up
    Error_Status = SUCCESS
    CALL ARMS_GeometryInfo_GetValue( GeometryInfo, Sensor_Zenith_Angle = Sensor_Zenith_Angle )


    ! Compute the surface emissivities
        IF ( SC(SensorIndex)%Frequency(ChannelIndex) < FREQUENCY_THRESHOLD ) THEN
          DO i = 1, SfcOptics%n_Angles
            CALL LandEM( SfcOptics%Angle(i),                      & ! Input, Degree
                                SC(SensorIndex)%Frequency(ChannelIndex), & ! Input, GHz
                                NOT_USED(1),                             & ! Input, Soil_Moisture_Content, g.cm**-3
                                NOT_USED(1),                             & ! Input, Vegetation_Fraction
                                Surface%Snow_Temperature,                & ! Input, K
                                Surface%Snow_Temperature,                & ! Input, K
                                Surface%Lai,                             & ! Input, Leaf Area Index
                                Surface%Soil_Type,                       & ! Input, Soil Type (1 -  9)
                                Surface%Vegetation_Type,                 & ! Input, Vegetation Type (1 - 13)
                                Surface%Snow_Depth,                      & ! Input, mm
                                SfcOptics%Emissivity(i,2),               & ! Output, H component
                                SfcOptics%Emissivity(i,1)                ) ! Output, V component
          END DO
        ELSE
          SfcOptics%Emissivity(1:SfcOptics%n_Angles,1:2) = DEFAULT_EMISSIVITY
        END IF


    ! Compute the surface reflectivities,
    ! assuming a specular surface
    SfcOptics%Reflectivity = ZERO
    DO i = 1, SfcOptics%n_Angles
      SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
      SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
    END DO

  END FUNCTION MW_Snow_SfcOptics


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Snow_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a snow surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_TL:   Structure containing the tangent-linear surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Snow_SfcOptics_TL( &
    SfcOptics_TL) &  ! TL  Output
  RESULT ( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS), INTENT(IN OUT) :: SfcOptics_TL
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Snow_SfcOptics_TL'
    ! Local variables


    ! Set up
    err_stat = SUCCESS


    ! Compute the tangent-linear surface optical parameters
    ! ***No TL models yet, so default TL output is zero***
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO

  END FUNCTION MW_Snow_SfcOptics_TL


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Snow_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a snow surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_AD:   Structure containing the adjoint surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Snow_SfcOptics_AD( &
    SfcOptics_AD) &  ! AD  Input
  RESULT( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics_AD
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Snow_SfcOptics_AD'
    ! Local variables


    ! Set up
    err_stat = SUCCESS


    ! Compute the adjoint surface optical parameters
    ! ***No AD models yet, so there is no impact on AD result***
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO

  END FUNCTION MW_Snow_SfcOptics_AD

!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Ice_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over an ice surface.
!
!       This function is a wrapper for third party code.
!
! INPUT:
!       Surface:        Structure containing the surface emissivity data.
!
!       GeometryInfo:   Structure containing the view geometry data.
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!
! OUTPUT:
!       SfcOptics:      Structure containing the surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------


  FUNCTION MW_Ice_SfcOptics( &
    Surface     , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   ) &  ! Output
  RESULT( Error_Status )
    ! Arguments
    TYPE(Surface_ARMS),      INTENT(IN)     :: Surface
    TYPE(GeometryInfo_ARMS), INTENT(IN)     :: GeometryInfo
    INTEGER,                 INTENT(IN)     :: SensorIndex
    INTEGER,                 INTENT(IN)     :: ChannelIndex
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Ice_SfcOptics'
    REAL(AUS),     PARAMETER :: DEFAULT_EMISSIVITY = 0.92_AUS
    REAL(AUS),     PARAMETER :: NOT_USED(4) = -99.9_AUS
    ! Local variables
    INTEGER :: i
    REAL(AUS) :: Sensor_Zenith_Angle


    ! Set up
    Error_Status = SUCCESS
    CALL ARMS_GeometryInfo_GetValue( GeometryInfo, Sensor_Zenith_Angle = Sensor_Zenith_Angle )


    ! Compute the surface emissivities

      ! Default physical model
        DO i = 1, SfcOptics%n_Angles
!          CALL NESDIS_SIce_Phy_EM( SC(SensorIndex)%Frequency(ChannelIndex), &  ! Input, GHz
!                                   SfcOptics%Angle(i),                      &  ! Input, Degree
!                                   Surface%Ice_Temperature,                 &  ! Input, K
!                                   Surface_Dummy%Salinity,                  &  ! Input
!                                   SfcOptics%Emissivity(i,2),               &  ! Output, H component
!                                   SfcOptics%Emissivity(i,1)                )  ! Output, V component
          SfcOptics%Emissivity(i,1:2) = DEFAULT_EMISSIVITY
        END DO



    ! Compute the surface reflectivities,
    ! assuming a specular surface
    SfcOptics%Reflectivity = ZERO
    DO i = 1, SfcOptics%n_Angles
      SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
      SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
    END DO

  END FUNCTION MW_Ice_SfcOptics


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Ice_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over an ice surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_TL:   Structure containing the tangent-linear surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Ice_SfcOptics_TL( &
    SfcOptics_TL) &  ! TL  Output
  RESULT ( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS), INTENT(IN OUT) :: SfcOptics_TL
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Ice_SfcOptics_TL'
    ! Local variables


    ! Set up
    err_stat = SUCCESS


    ! Compute the tangent-linear surface optical parameters
    ! ***No TL models yet, so default TL output is zero***
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO

  END FUNCTION MW_Ice_SfcOptics_TL


!----------------------------------------------------------------------------------
!
! FUNCTION NAME:
!       MW_Ice_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over an ice surface.
!
! INPUT:      
!
! OUTPUT:
!       SfcOptics_AD:   Structure containing the adjoint surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!
! FUNCTION RESULT:
!       err_stat:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred!
!
!----------------------------------------------------------------------------------

  FUNCTION MW_Ice_SfcOptics_AD( &
    SfcOptics_AD) &  ! AD  Input
  RESULT( err_stat )
    ! Arguments
    TYPE(SfcOptics_ARMS),    INTENT(IN OUT) :: SfcOptics_AD
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MW_Ice_SfcOptics_AD'
    ! Local variables


    ! Set up
    err_stat = SUCCESS


    ! Compute the adjoint surface optical parameters
    ! ***No AD models yet, so there is no impact on AD result***
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO

  END FUNCTION MW_Ice_SfcOptics_AD

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       LandEM
!
! PURPOSE:
!       Subroutine to simulate microwave emissivity over land conditions.
!
! REFERENCES:
!       Weng, F. et al, 2001: "A microwave land emissivity model".
!
!
! INPUT:
!         Angle:                   The angle values in degree.
!
!         Frequency                Frequency User defines
!
!         Soil_Moisture_Content:   The volumetric water content of the soil (0:1).
!
!         Vegetation_Fraction:     The vegetation fraction of the surface (0:1).
!
!         Soil_Temperature:        The soil temperature(K).
!
!         Land_Temperature:        The land surface temperature(K).
!
!         Snow_Depth:              The snow depth(mm).
!
! OUTPUT:
!         Emissivity_H:            The surface emissivity at a horizontal
!                                  polarization.
!
!         Emissivity_V:            The surface emissivity at a vertical polarization.
!
!
! INTERNAL:
!       theta       -  local zenith angle in radian
!       rhob        -  bulk volume density of the soil (1.18-1.12)
!       rhos        -  density of the solids (2.65 g.cm**3 for solid soil material)
!       sand        -  sand fraction (sand + clay = 1.0)
!       clay        -  clay fraction (0-1.0)
!       lai         -  leaf area index (eg. lai = 4.0 for corn leaves)
!       sigma       -  surface roughness formed between medium 0.5 and 3.0,
!                      expressed as the standard deviation of roughtness height (cm)
!       leaf_thick  --  leaf thickness (mm)
!       rad         -  radius of dense medium scatterers (mm)
!       va          -  fraction volume of dense medium scatterers(0.0 - 1.0)
!       ep          -  dielectric constant of ice or sand particles, complex value
!                               (e.g, 3.0+i0.0)
!------------------------------------------------------------------------------------------------------------

  SUBROUTINE        LandEM(Angle,                 &   ! Input
                           Frequency,             &   ! Input
                           Soil_Moisture_Content, &   ! Input
                           Vegetation_Fraction,   &   ! Input
                           Soil_Temperature,      &   ! Input
                           t_skin,                &   ! Input
                           Lai,                   &   ! Input
                           Soil_Type,             &   ! Input
                           Vegetation_Type,       &   ! Input
                           Snow_Depth,            &   ! Input
                           Emissivity_H,          &   ! Output
                           Emissivity_V)              ! Output
    ! Arguments
    REAL(AUS), intent(in) :: Angle
    REAL(AUS), intent(in) :: Frequency
    REAL(AUS), intent(in) :: Soil_Moisture_Content
    REAL(AUS), intent(in) :: Vegetation_Fraction
    REAL(AUS), intent(in) :: Soil_Temperature
    REAL(AUS), intent(in) :: t_skin
    REAL(AUS), intent(in) :: Lai
    INTEGER,  intent(in) :: Soil_Type
    INTEGER,  intent(in) :: Vegetation_Type
    REAL(AUS), intent(in) :: Snow_Depth
    REAL(AUS), intent(out):: Emissivity_V,Emissivity_H
    ! Local parameters
    REAL(AUS), PARAMETER :: snow_depth_c     = 10.0_AUS    ! unit:mm
    REAL(AUS), PARAMETER :: tsoilc_undersnow = 280.0_AUS   ! unit:k
    REAL(AUS), PARAMETER :: rhos = 2.65_AUS                ! unit: g/cm**3
    REAL(AUS)            :: sand, clay, rhob
    REAL(AUS), PARAMETER, dimension(0:9) :: frac_sand = (/ 0.80_AUS,     &
                          0.92_AUS, 0.10_AUS, 0.20_AUS, 0.51_AUS, 0.50_AUS, &
                          0.35_AUS, 0.60_AUS, 0.42_AUS,  0.92_AUS  /)
    REAL(AUS), PARAMETER, dimension(0:9) :: frac_clay = (/ 0.20_AUS,     &
                          0.06_AUS, 0.34_AUS, 0.63_AUS, 0.14_AUS, 0.43_AUS, &
                          0.34_AUS, 0.28_AUS, 0.085_AUS, 0.06_AUS /)
    REAL(AUS), PARAMETER, dimension(0:9) :: rhob_soil = (/ 1.48_AUS,     &
                          1.68_AUS, 1.27_AUS, 1.21_AUS, 1.48_AUS, 1.31_AUS, &
                          1.32_AUS, 1.40_AUS, 1.54_AUS, 1.68_AUS /)
! Specific Density
    REAL(AUS), PARAMETER, dimension(0:13) :: veg_rho  = (/ 0.33_AUS,     &
                          0.40_AUS, 0.40_AUS, 0.40_AUS, 0.40_AUS, 0.40_AUS, &
                          0.25_AUS, 0.25_AUS, 0.40_AUS, 0.40_AUS, 0.40_AUS, &
                          0.40_AUS, 0.33_AUS, 0.33_AUS            /)
! MGE
    REAL(AUS), PARAMETER, dimension(0:13) :: veg_mge  = (/ 0.50_AUS,     &
                          0.45_AUS, 0.45_AUS, 0.45_AUS, 0.40_AUS, 0.40_AUS, &
                          0.30_AUS, 0.35_AUS, 0.30_AUS, 0.30_AUS, 0.40_AUS, &
                          0.30_AUS, 0.50_AUS, 0.40_AUS            /)
! LAI
    REAL(AUS), PARAMETER, dimension(0:13) :: lai_min  = (/ 0.52_AUS,     &
                          3.08_AUS, 1.85_AUS, 2.80_AUS, 5.00_AUS, 1.00_AUS, &
                          0.50_AUS, 0.52_AUS, 0.60_AUS, 0.50_AUS, 0.60_AUS, &
                          0.10_AUS, 1.56_AUS, 0.01_AUS            /)
    REAL(AUS), PARAMETER, dimension(0:13) :: lai_max  = (/ 2.90_AUS,     &
                          6.48_AUS, 3.31_AUS, 5.50_AUS, 6.40_AUS, 5.16_AUS, &
                          3.66_AUS, 2.90_AUS, 2.60_AUS, 3.66_AUS, 2.60_AUS, &
                          0.75_AUS, 5.68_AUS, 0.01_AUS            /)
! Leaf_thickness
    REAL(AUS), PARAMETER, dimension(0:13) :: leaf_th  = (/ 0.07_AUS,     &
                          0.18_AUS, 0.18_AUS, 0.18_AUS, 0.18_AUS, 0.18_AUS, &
                          0.12_AUS, 0.12_AUS, 0.12_AUS, 0.12_AUS, 0.12_AUS, &
                          0.12_AUS, 0.15_AUS, 0.12_AUS            /)
    ! Local variables
    REAL(AUS) :: mv,veg_frac,theta,theta_i,theta_t,mu,r21_h,r21_v,r23_h,r23_v,  &
                t21_v,t21_h,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v,mge, &
                    leaf_thick,rad,sigma,va,ep_real,ep_imag
    REAL(AUS) :: t_soil
    REAL(AUS) :: rhoveg, vlai
    REAL(AUS) :: local_snow_depth
    COMPLEX(AUS) :: esoil, eveg, esnow, eair
    LOGICAL :: SnowEM_Physical_Model
    LOGICAL :: Diel_Model_Mironov
    LOGICAL :: Roughness_Model_ChenWeng

    eair = CMPLX(ONE,-ZERO,AUS)
    theta = Angle*PI_AUS/180.0_AUS

    ! By default use the assign local variable
    mv               = Soil_Moisture_Content
    veg_frac         = Vegetation_Fraction
    t_soil           = Soil_Temperature
    sand = frac_sand(Soil_Type)
    clay = frac_clay(Soil_Type )
    rhob = rhob_soil(Soil_Type )
    
    local_snow_depth = Snow_Depth
    
    ! By default use the Mironov soil dielectric constant model,otherwise use the Dobson model
    Diel_Model_Mironov = .TRUE.
    
    !By default use the Chen-Weng surface roughness reflectance model,otherwise use the Qh model
    Roughness_Model_ChenWeng = .TRUE.

    ! ---- Method trace: default configuration ----
    IF (dbg) THEN
    WRITE(*,'(A)') '--- LandEM default model configuration ---'
    IF (Diel_Model_Mironov) THEN
     WRITE(*,'(A)') 'Default soil dielectric model: Mironov'
    ELSE
     WRITE(*,'(A)') 'Default soil dielectric model: Dobson'
    END IF

    IF (Roughness_Model_ChenWeng) THEN
    WRITE(*,'(A)') 'Default surface roughness model: C-W'
    ELSE
    WRITE(*,'(A)') 'Default surface roughness model: Qh'
    END IF
    WRITE(*,'(A)') '------------------------------------------'
    END IF
    
    ! Check soil/skin temperature
    if ( (t_soil <= 100.0_AUS .OR.  t_soil >= 350.0_AUS) .AND. &
         (t_skin >= 100.0_AUS .AND. t_skin <= 350.0_AUS) ) t_soil = t_skin

    ! Check soil moisture content range
    mv = MAX(MIN(mv,ONE),ZERO)

    ! Surface type based on snow depth
    IF (local_snow_depth > POINT1) THEN
      
      ! By default use the physical model for snow (Snow_Depth > 0.1mm and <=10mm )
        
      SnowEM_Physical_Model = .TRUE.
      if (local_snow_depth > snow_depth_c) SnowEM_Physical_Model = .FALSE.

      ! Compute the snow emissivity 
      IF ( SnowEM_Physical_Model ) THEN

        ep_real = 3.2_AUS   
        ep_imag = -0.0005_AUS
        sigma = 0.8_AUS 

        ! For deep snow, the performance of the model is poor
        local_snow_depth = MIN(local_snow_depth,1000.0_AUS)

        ! The fraction volume of dense medium
        ! scatterers must be in the range (0-1)
        va = 0.4_AUS + 0.0004_AUS*local_snow_depth
        va = MAX(MIN(va,ONE),ZERO)

        ! Limit for snow grain size
        rad = MIN((POINT5 + 0.005_AUS*local_snow_depth),ONE)

        ! Limit for soil temperature
        t_soil = MIN(t_soil,tsoilc_undersnow)

        CALL Snow_Diel(Frequency, ep_real, ep_imag, rad, va, esnow)
                
        IF ( t_soil > ( T0+ TS_FROZEN_THRESHOLD) ) THEN
           ! Use the default Mironov soil model
           IF ( Diel_Model_Mironov ) THEN
              CALL Soil_Diel_Mironov(Frequency, t_soil, mv, clay, esoil)
           ELSE
              CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
           ENDIF
        ELSE
           ! Use the Dobson frozen soil model
           CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
              
        ENDIF
        
        theta_i = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(esnow),AUS))

        CALL Smooth_Reflectance(esnow, eair, theta_i,  theta, r21_v, r21_h)
        CALL Transmittance(esnow, eair, theta_i, theta, t21_v, t21_h)

        mu      = COS(theta_i)
        theta_t = ASIN(REAL(SIN(theta_i)*SQRT(esnow)/SQRT(esoil),AUS))

        CALL Smooth_Reflectance(esnow, esoil, theta_i, theta_t, r23_v, r23_h)
           
        ! Use the default Chen-Weng model
        IF ( Roughness_Model_ChenWeng ) THEN
           CALL Roughness_Reflectance_ChenWeng (Frequency, sigma, theta, r23_v, r23_h)
        ELSE
           CALL Roughness_Reflectance_Qh(Frequency, sigma, r23_v, r23_h)
        ENDIF
        
        CALL Snow_Optic(Frequency,rad,local_snow_depth,va,ep_real, ep_imag,gv,gh,&
                        ssalb_v,ssalb_h,tau_v,tau_h)
                        
        CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                                 r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                               frequency, t_soil, t_skin)
      ELSE
        ! Use the empirical method (Snow_Depth > 10mm)
        CALL SnowEM_Default(Frequency,t_skin, local_snow_depth,Emissivity_V,Emissivity_H)
      END IF

    ELSE

      ! No snow, so we're going to compute vegetation canopy emissivities....
      
      ! Limit for vegetation fraction
      veg_frac = MAX(MIN(veg_frac,ONE),ZERO)

!     lai = THREE*veg_frac + POINT5
!     mge = POINT5*veg_frac
!     leaf_thick = 0.07_AUS
      mu  = COS(theta)
      sigma = ONE
      
      vlai = Lai*veg_frac
      mge = veg_mge(Vegetation_Type)
      rhoveg = veg_rho(Vegetation_Type)
      leaf_thick = leaf_th(Vegetation_Type)
     
      r21_h    = ZERO
      r21_v    = ZERO
      t21_h    = ONE
      t21_v    = ONE

      IF ( t_soil >= T0+TS_FROZEN_THRESHOLD ) THEN
           ! Use the default Mironov soil model
           IF ( Diel_Model_Mironov ) THEN
              CALL Soil_Diel_Mironov(Frequency, t_soil, mv, clay, esoil)
           ELSE
              CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)
           ENDIF
      ELSE
           ! Use the Dobson frozen soil model
          CALL Soil_Diel_Dobson(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)   
      ENDIF

      theta_t = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(esoil),AUS))
      CALL Smooth_Reflectance(eair, esoil, theta, theta_t, r23_v, r23_h)

      ! Use the default Chen-Weng model
      IF ( Roughness_Model_ChenWeng ) THEN
         CALL Roughness_Reflectance_ChenWeng (Frequency, sigma, theta, r23_v, r23_h)
      ELSE
         CALL Roughness_Reflectance_Qh(Frequency, sigma, r23_v, r23_h)
      ENDIF
      
      CALL Canopy_Diel(Frequency, mge, eveg, rhoveg, t_soil)
      
      CALL Canopy_Optic(vlai,Frequency,theta,eveg,leaf_thick,gv,gh,ssalb_v,ssalb_h,tau_v,tau_h)
      
      CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                               r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                               frequency, t_soil, t_skin)
                      
    END IF

  END SUBROUTINE LandEM



!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       SnowEM_Default
!
! PURPOSE:
!       Preliminary estimate of snow emissivity using surface temperature and snow depth
! INPUT:
!         ts:                   Surface temperature.
!
!         Frequency             Frequency (GHz)
!
! OUTPUT:
!         Emissivity_H:            The snow emissivity at a horizontal polarization.
!
!         Emissivity_V:            The snow emissivity at a vertical polarization.
!
!----------------------------------------------------------------------------------

SUBROUTINE SnowEM_Default(frequency,ts, depth,Emissivity_V,Emissivity_H)

  ! Arguments
  REAL(AUS) :: frequency,ts, depth,Emissivity_V,Emissivity_H
  ! Local parameters  
  INTEGER , PARAMETER :: new = 7     !AMSR-E seven Frequency
  INTEGER , PARAMETER :: NFRESH_SHALLOW_SNOW = 1
  INTEGER , PARAMETER :: NPOWDER_SNOW        = 2
  INTEGER , PARAMETER :: NWET_SNOW           = 3
  INTEGER , PARAMETER :: NDEEP_SNOW          = 4
  REAL(AUS), PARAMETER :: twet    = 270.0_AUS   ! unit:K
  REAL(AUS), PARAMETER :: tcrust  = 235.0_AUS   ! unit:K
  REAL(AUS), PARAMETER :: depth_s =  50.0_AUS   ! unit:mm
  REAL(AUS), PARAMETER :: depth_c = 100.0_AUS   ! unit:mm
  ! Local variables
  INTEGER :: ich,basic_snow_type
  REAL(AUS), DIMENSION(new) :: ev, eh, freq
  REAL(AUS) :: df, df0


  freq = FREQUENCY_AMSRE(1:new)
  
  ! Determine the snow type based on temperatures
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

  ! Assign the emissivity spectrum
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

  ! Handle possible extrapolation
  if (frequency <= freq(1)) then
    Emissivity_H = eh(1)
    Emissivity_V = ev(1)
    return
  end if
  if (frequency >= freq(new)) then
    Emissivity_H = eh(new)
    Emissivity_V = ev(new)
    return
  end if

  ! Interpolate emissivity at a certain frequency
  Channel_loop: do ich=2,new
    if (frequency <= freq(ich)) then
      df  = frequency-freq(ich-1)
      df0 = freq(ich)-freq(ich-1)
      Emissivity_H = eh(ich-1) + (df*(eh(ich)-eh(ich-1))/df0)
      Emissivity_V = ev(ich-1) + (df*(ev(ich)-ev(ich-1))/df0)
      exit Channel_loop
    end if
  end do Channel_loop

end subroutine SnowEM_Default

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Canopy_Optic
!
! PURPOSE:
!       Compute optic parameters for canopy
! INPUT:
!      lai:           Leaf area index
!      frequency:     Frequency (GHz)
!      theta:         Incident angle
!      esv:           Leaf dielectric constant
!      d:             Leaf thickness (mm)
!
! OUTPUT:
!      gv:             Asymmetry factor for v pol
!      gh:             Asymmetry factor for h pol
!      ssalb_v:        Single scattering albedo at v. polarization
!      ssalb_h:        Single scattering albedo at h. polarization
!      tau_v:          Optical depth at v. polarization
!      tau_h:          Optical depth at h. polarization
! INTERNAL:
!       k0:            wave number (1/mm)
!       rhc:           Horizontal polarization reflectance factor
!       rvc:           Vertical polarization reflectance factor
!       factt:         Transmission factor
!       rv:            Vertical polarization reflectivity
!       rh:            Horizontal polarization reflectivity
!       th:            Horizontal polarization transmittance
!       tv:            Vertical polarization transmittance
!----------------------------------------------------------------------------------


subroutine Canopy_Optic(vlai,frequency,theta,esv,d,gv,gh,&
                        ssalb_v,ssalb_h,tau_v, tau_h)

  REAL(AUS) :: frequency,theta,d,vlai,ssalb_v,ssalb_h,tau_v,tau_h,gv, gh, mu
  COMPLEX(AUS) :: ix,k0,kz0,kz1,rhc,rvc,esv,expval1,factt,factrvc,factrhc
  REAL(AUS) :: rh,rv,th,tv
  REAL(AUS), PARAMETER :: threshold = 0.999_AUS

  mu = COS(theta)
  ix = CMPLX(ZERO, ONE, AUS)

  k0  = CMPLX(TWOPI*frequency/300.0_AUS, ZERO, AUS)    
  kz0 = k0*mu
  kz1 = k0*SQRT(esv - SIN(theta)**2)

  rhc = (kz0 - kz1)/(kz0 + kz1)            
  rvc = (esv*kz0 - kz1)/(esv*kz0 + kz1)    

  expval1 = EXP(-TWO*ix*kz1*d)
  factrvc = ONE - rvc**2*expval1
  factrhc = ONE - rhc**2*expval1
  factt   = FOUR*kz0*kz1*EXP(ix*(kz0-kz1)*d)   

  rv = ABS(rvc*(ONE - expval1)/factrvc)**2     
  rh = ABS(rhc*(ONE - expval1)/factrhc)**2     

  th = ABS(factt/((kz1+kz0)**2*factrhc))**2    
  tv = ABS(esv*factt/((kz1 + esv*kz0)**2*factrvc))**2  

  gv = POINT5
  gh = POINT5

  tau_v = POINT5*vlai*(TWO - tv-th)
  tau_h = tau_v

  ssalb_v = MIN((rv + rh)/(TWO -tv -th),threshold)
  ssalb_h = ssalb_v

end subroutine Canopy_Optic

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Snow_Optic
!
! PURPOSE:
!       Compute optic parameters for snow
! INPUT:
!      theta:          Local zenith angle (degree)
!      frequency:      Frequency (GHz)
!      ep_real:        Real part of dielectric constant of particles
!      ep_imag:        Imaginary part of dielectric constant of particles
!      a:              Particle radiu (mm)
!      h:              Snow depth(mm)
!      f:              Fraction volume of snow (0.0 - 1.0)
!
! OUTPUT:
!       ssalb:         Single scattering albedo
!       tau:           Optical depth
!       g:             Asymmetry factor
! INTERNAL:
!       k :            wave number (/mm)
!       ks:            scattering coeffcient (/mm)
!       ka:            absorption coeffient (/mm)
!       kr:            the real part of the effective propagation coefficient K
!       ki:            the imaginary part of the effective propagation coefficient K
!       kp:            eigenvalue of two-stream approximation
!       y:             yr+iyi
!----------------------------------------------------------------------------------

subroutine Snow_Optic(frequency,a,h,f,ep_real,ep_imag,gv,gh, ssalb_v,ssalb_h,tau_v,tau_h)

  REAL(AUS) :: yr,yi,ep_real,ep_imag
  REAL(AUS) :: frequency,a,h,f,ssalb_v,ssalb_h,tau_v,tau_h,gv,gh,k
  REAL(AUS) :: ks1,ks2,ks3,ks,kr1,kr2,kr3,kr,ki1,ki2,ki3,ki
  REAL(AUS) :: fact1,fact2,fact3,fact4,fact5

  k = TWOPI/(300._AUS/frequency)

  yr = (ep_real - ONE)/(ep_real + TWO)
  yi = -ep_imag/(ep_real + TWO)

  fact1 = (ONE+TWO*f)**2
  fact2 = ONE-f*yr
  fact3 = (ONE-f)**4
  fact4 = f*(k*a)**3
  fact5 = ONE + TWO*f*yr

  ks1 = k*SQRT(fact2/fact5)
  ks2 = fact4*fact3/fact1
  ks3 = (yr/fact2)**2
  ks = ks1*ks2*ks3

  kr1 = fact5/fact2
  kr2 = TWO*ks2
  kr3 = TWO*yi*yr/(fact2**3)
  kr = k*SQRT(kr1 + kr2*kr3)

  ki1 = THREE*f*yi/fact2**2
  ki2 = kr2
  ki3 = ks3
  ki  = k**2/(TWO*kr)*(ki1 + ki2*ki3)
  
  ! ka=ki-ks

  gv = POINT5
  gh = POINT5

  ssalb_v = MIN(ks/ki, 0.999_AUS)
  ssalb_h = ssalb_v
  tau_v = TWO*ki*h
  tau_h = tau_v

end subroutine Snow_Optic

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Soil_Diel_Mironov 
! REFERENCES:
!      Mironov V L, Kosolapova L G, Fomin S V.2009.Physically and Mineralogically Based Spectroscopic Dielectric Model for Moist Soils[J].
!                 IEEE Transactions on Geoscience & Remote Sensing, 47(7):2059-2070.
!    Mironov V L ,  Fomin S V. 2009. Temperature and Mineralogy Dependable Model for Microwave Dielectric Spectra of Moist Soils[J].
!              Piers Online, 5(5):411-415.
! PURPOSE:
!       Calculate the dielectric properties of soil based on Mironov model (2009)
! INPUT:
!      frequency:      Frequency (GHz)
!      t_soil:         Soil temperature (K)
!      mv:            Volumetric moisture content (0-1.0, demensionless or cm**3/cm**3)
!      clay:           Clay fraction (0-1.0)
!
! OUTPUT:
!       esm:           Dielectric constant for bare soil
! INTERNAL:
!      E_0:            The permittivity of free space (F/m)
!      E_wlimit        The permittivity at the high frequency limit

!      E_w0u:            Static dieletric constant of free water
!      E_w0b:            Static dieletric constant of bound water
!      Tao_wu:           Relaxation time of free water
!      Tao_wb:           Relaxation time of bound water
!      Sigma_eff_wu      Effective conductivity of free water
!      Sigma_eff_wb      Effective conductivity of bound water

!      mvt             Maximum combined water content (demensionless or cm3/cm3)
!      Ts              Reference temperature (�棩

!      nd              Refractive index of dry soil
!      kd              Normalized attenuation coefficient of dry soil
!      nb              Refractive index of bound water
!      kb              Normalized attenuation coefficient of bound water
!      nu              Refractive index of free water
!      ku              Normalized attenuation coefficient of free water
!      nm              Refractive index of soil
!      km              Normalized attenuation coefficient of soil


!----------------------------------------------------------------------------------

subroutine Soil_Diel_Mironov(freq,t_soil,mv,clay,esm)

  REAL(AUS) :: f, freq, t_soil, Ts, T, mv, mvt, clay, clay_per
  REAL(AUS) :: Tao_wu, Tao_wb, Sigma_eff_wu, Sigma_eff_wb, E_w0u, E_w0b, E_wlimit, E_0
  REAL(AUS) :: nd, kd, nb, kb, nu, ku, nm, km
  REAL(AUS) :: Beta_b, Beta_u, Delta_HbR, Delta_SbR, Delta_HuR, Delta_SuR, Beta_sigmab, Beta_sigmau
  REAL(AUS) :: Fb, Fu, E_wb1, E_wb2, E_wu1, E_wu2, E_m1, E_m2
  COMPLEX(AUS) :: esm

  clay_per = clay*100_AUS
  f = freq*1.0e9_AUS
  
  mvt = 0.02863_AUS + 0.30673e-2*clay_per
  
  ! Estimation of dry soil parameters 
  nd = 1.634_AUS - 0.539e-2*clay_per + 0.2748e-4*clay_per**2
  kd = 0.03952_AUS - 0.04038e-2*clay_per
    
  ! Bound water parameter estimation (Ts is 20 ��)
  E_w0b = 79.8_AUS - 85.4e-2*clay_per + 32.7e-4*clay_per**2
  Tao_wb = 1.062e-11 + 3.450e-12*1.0e-2*clay_per
  Sigma_eff_wb = 0.3112_AUS + 0.467e-2*clay_per
  
  ! Free water parameters (Ts is 20 ��)
  E_w0u = 100_AUS
  Tao_wu = 8.5e-12_AUS
  
  Ts = 20.0_AUS  
  T = t_soil - T0
  
  E_wlimit = 4.9_AUS  
  E_0 = 8.854e-12_AUS
      
  !------Correction parameters when the reference temperature Ts is 20 ��-----------
  Beta_b = 8.67e-19 - 0.00126e-2*clay_per + 0.00184e-4*clay_per**2 - 9.77e-10*clay_per**3 - 1.39e-15*clay_per**4   
  Beta_u = 1.11e-4 - 1.603e-7*clay_per + 1.239e-9*clay_per**2 + 8.33e-13*clay_per**3 - 1.007e-14*clay_per**4
  Delta_HbR = 1467_AUS + 2697e-2*clay_per - 980e-4*clay_per**2 + 1.368e-10*clay_per**3 - 8.61e-13*clay_per**4
  Delta_SbR = 0.888_AUS + 9.7e-2*clay_per - 4.262e-4*clay_per**2 + 6.79e-21*clay_per**3 + 4.263e-22*clay_per**4
  Delta_HuR = 2231_AUS - 143.1e-2*clay_per + 2232e-4*clay_per**2 - 142.1e-6*clay_per**3 + 27.14e-8*clay_per**4
  Delta_SuR = 3.649_AUS - 0.4894e-2*clay_per + 0.763e-4*clay_per**2 - 0.4859e-6*clay_per**3 + 0.0928e-8*clay_per**4
  Beta_sigmab = 0.0028_AUS + 0.02094e-2*clay_per - 0.01229e-4*clay_per**2 - 5.03e-22*clay_per**3 + 4.163e-24*clay_per**4  
  Beta_sigmau = 0.00108_AUS + 0.1413e-2*clay_per - 0.2555e-4*clay_per**2 + 0.2147e-6*clay_per**3 - 0.0711e-8*clay_per**4  
  Sigma_eff_wu = 0.05_AUS + 1.4*(1 -(1 -1.0e-2*clay_per)**4.664) 
      
  !------Static dielectric constants of bound water and free water at arbitrary temperature T ------
  Fb = LOG((E_w0b - ONE)/(E_w0b + TWO))
  E_w0b = (ONE + TWO*EXP(Fb-Beta_b*(T - Ts)))/(ONE - EXP(Fb -Beta_b*(T -Ts)))
 
  Fu = LOG((E_w0u - ONE)/(E_w0u + TWO))
  E_w0u = (ONE + TWO*EXP(Fb - Beta_u*(T - Ts)))/(ONE - EXP(Fb - Beta_u*(T - Ts)))
 
  !------Polarization relaxation time of bound water and free water at arbitrary temperature T------
  Tao_wb = 48e-12 /( T + T0 )*EXP( Delta_HbR / ( T + T0 )-Delta_SbR)
  Tao_wu = 48e-12 /( T + T0 )*EXP( Delta_HuR / ( T + T0 )-Delta_SuR)
 
  !------Effective conductivity of bound water and free water at any temperature T------
  Sigma_eff_wb = Sigma_eff_wb + Beta_sigmab*(T-Ts)
  Sigma_eff_wu = Sigma_eff_wu + Beta_sigmau*(T-Ts)
 
  !------Estimation of refractive index and normalized attenuation coefficient of bound water------
  E_wb1 = E_wlimit+(E_w0b-E_wlimit)/(ONE + (TWOPI*f*Tao_wb)**2)
  E_wb2 = TWOPI*f*Tao_wb*(E_w0b-E_wlimit)/(ONE + (TWOPI*f*Tao_wb)**2) + (Sigma_eff_wb)/(TWOPI*E_0*f)
  nb = SQRT(SQRT(E_wb1**2 + E_wb2**2) + E_wb1)/SQRT(TWO)
  kb = SQRT(SQRT(E_wb1**2 + E_wb2**2) - E_wb1)/SQRT(TWO)
 
  !------Estimation of refractive index and normalized attenuation coefficient of free water------
  E_wu1 = E_wlimit + (E_w0u-E_wlimit)/(ONE + (TWOPI*f*Tao_wu)**2)
  E_wu2 = TWOPI*f*Tao_wu*(E_w0u-E_wlimit)/(ONE + (TWOPI*f*Tao_wu)**2) + (Sigma_eff_wu)/(TWOPI*E_0*f)
  nu = SQRT(SQRT(E_wu1**2 + E_wu2**2) + E_wu1)/SQRT(TWO)
  ku = SQRT(SQRT(E_wu1**2 + E_wu2**2) - E_wu1)/SQRT(TWO)
 
  !------Estimation of refractive index and normalized attenuation coefficient of soil------
  if ( mv <= mvt ) then 
     nm = nd+(nb-1)*mv
     km = kd+kb*mv
  else
     nm = nd+(nb-1)*mvt+(nu-1)*(mv-mvt)
     km = kd+kb*mvt+ku*(mv-mvt)
  endif
 
  !------ Estimation of dielectric constant for bare soil------
  E_m1 = nm**2-km**2
  E_m2 = TWO*nm*km
  esm = E_m1 - CMPLX(ZERO,ONE)*E_m2
 
  if(AIMAG(esm) >= ZERO) esm = CMPLX(REAL(esm,AUS),-0.0001_AUS, AUS)
 
end subroutine Soil_Diel_Mironov


!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Soil_Diel_Dobson
! REFERENCES:
!         Dobson M C . Microwave Dielectric Behavior of Wet Soil-Part II: Dielectric Mixing Models[J].
!                       Geoscience & Remote Sensing IEEE Transactions on, 1985, GE-23(1):35-46.
!         Peplinski N R , Ulaby F T . Dielectric properties of soils in the 0.3-1.3-GHz range[J]. 
!                     IEEE Trans.on Geosci. & Remote Sensing, 1995, 33(3):803-807.
! PURPOSE:
!       Calculate the dielectric properties of soil based on Dobson model
! INPUT:
!      theta:          Local zenith angle (degree)
!      frequency:      Frequency (GHz)
!      t_soil:         Soil temperature
!      mv:            Volumetric moisture content (demensionless)
!      rhob:           Bulk volume density of the soil (1.18-1.12)
!      rhos:           Density of the solids (2.65 g.cm**3 for
!                      Solid soil material)
!      sand:           Sand fraction (sand + clay = 1.0)
!      clay:           Clay fraction (0-1.0)
!
! OUTPUT:
!       esm:           Dielectric constant for bare soil
! INTERNAL:
!      esof:           The permittivity of free space
!      eswo:           Static dieletric constant of water
!      eswi            The permittivity at the high frequency limit
!      tauw:           Relaxation time of water
!      s   :           Salinity
!      alpha:          Shape factor
!      beta            Implicitly considering the empirical coefficient of bound water
!      ess             Relative dielectric constant of soil parent material (dry soil)
!      rhoef           Effective conductivity of water
!      vlw             Volume content of unfrozen water (0-1.0)
!      vic             Volume content of frozen water (0-1.0)

!----------------------------------------------------------------------------------
subroutine Soil_Diel_Dobson(freq,t_soil,mv,rhob,rhos,sand,clay,esm)


  REAL(AUS) :: f,tauw,freq,t_soil,mv,rhob,rhos,sand,clay
  REAL(AUS) :: alpha,beta,ess,rhoef,t,eswi,eswo
  REAL(AUS) :: esof,vlw,vic
  COMPLEX(AUS) :: esm,esw,es1,es2,esice

  alpha = 0.65_AUS
  beta  = 1.09_AUS - 0.11_AUS*sand + 0.18_AUS*clay
  ess = (1.01_AUS + 0.44_AUS*rhos)**2 - 0.062_AUS
  
  !rhoef = -1.645_AUS + 1.939_AUS*rhob - 0.020213_AUS*sand + 0.01594_AUS*clay
  rhoef = -1.645_AUS + 1.939_AUS*rhob - 0.0225622_AUS*sand*100 + 0.01594_AUS*clay*100    
    
  ! Effective conductivity of water for 0.3-1.3 GHz
  if (freq < 1.4) then
    rhoef = 0.0467 + 0.2204*rhob - 0.41111*sand + 0.6614*clay
  end if
  
  t = t_soil - T0
  f = freq*1.0e9_AUS

  !eswi = 5.5_AUS
  eswi = 4.9_AUS  ! Correction of the original 'eswi' value
  esof = 8.854e-12_AUS

  eswo = 87.134_AUS + (-1.949e-1_AUS + (-1.276e-2_AUS + 2.491e-4_AUS*t)*t)*t
  tauw = 1.1109e-10_AUS + (-3.824e-12_AUS + (6.938e-14_AUS - 5.096e-16_AUS*t)*t)*t

  if (mv > ZERO) then
     es1 = CMPLX(eswi, -rhoef*(rhos-rhob)/(TWOPI*f*esof*rhos*mv), AUS)
  else
     es1 = CMPLX(eswi, ZERO, AUS)
  endif

  es2 = CMPLX(eswo-eswi, ZERO, AUS)/CMPLX(ONE, f*tauw, AUS)
  esw = es1 + es2
  
  if ( t_soil >= ( T0+TS_FROZEN_THRESHOLD) ) then
    esm = ONE + (ess**alpha - ONE)*rhob/rhos + mv**beta*esw**alpha - mv
    esm = esm**(ONE/alpha)
    
  else            ! Estimation of the dielectric constant of frozen soil
  
    esice=CMPLX(3.15, 57.34*(1.0/f + 2.48e-14*f**0.5)*EXP(0.0362*t), AUS)
    vic = RATIO_FROZEN*mv
    vlw = ( ONE-RATIO_FROZEN )*mv
    esm = ONE + (ess**alpha - ONE)*rhob/rhos + (vlw**beta*esw**alpha - vlw)+vic*(esice**alpha - ONE)
    esm = esm**(ONE/alpha)
  endif
   
  if(AIMAG(esm) >= ZERO) esm = CMPLX(REAL(esm,AUS),-0.0001_AUS, AUS)

end subroutine Soil_Diel_Dobson

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Snow_Diel
!
! PURPOSE:
!       Calculate the dielectric properties of snow
! INPUT:
!      frequency:      Frequency (GHz)
!      ep_real:        Real part of dielectric constant of particle
!      ep_imag:        Imaginary part of dielectric constant of particle
!      rad:            Particle radiu (mm)
!      frac:           Fraction volume of snow (0.0 - 1.0)
!
! OUTPUT:
!      ep_eff:         Dielectric constant for the dense medium
! INTERNAL:
!        k0:           Wave number (1/mm) 
!----------------------------------------------------------------------------------

subroutine Snow_Diel(frequency,ep_real,ep_imag,rad,frac,ep_eff)

  REAL(AUS) :: ep_imag,ep_real
  REAL(AUS) :: frequency,rad,frac,k0,yr,yi
  COMPLEX(AUS) :: y,ep_r,ep_i,ep_eff,fracy

  k0 = TWOPI/(300.0_AUS/frequency)

  yr = (ep_real - ONE)/(ep_real + TWO)
  yi = ep_imag/(ep_real + TWO)

  y = CMPLX(yr, yi, AUS)
  fracy=frac*y

  ep_r = (ONE + TWO*fracy)/(ONE - fracy)
  ep_i = TWO*fracy*y*(k0*rad)**3*(ONE-frac)**4/((ONE-fracy)**2*(ONE+TWO*frac)**2)
  ep_eff = ep_r - CMPLX(ZERO,ONE,AUS)*ep_i

  if (AIMAG(ep_eff) >= ZERO) ep_eff = CMPLX(REAL(ep_eff), -0.0001_AUS, AUS)

end subroutine Snow_Diel

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Canopy_Diel
! REFERENCES:
!     Ulaby and el-rayer, 1987: microwave dielectric spectrum of vegetation part ii,
!           dual-dispersion model, ieee trans geosci. remote sensing, 25, 550-557
!      Matzler C.1994. Microwave (1-100 GHz) dielectric model of leaves.
!            IEEE Transactions on Geoscience & Remote Sensing, 32(5):9447
! PURPOSE:
!       Calculate the dielectric properties of vegetation canopy
! INPUT:
!      frequency:      Frequency (GHz)
!      mg:             Gravimetric water content
!
! OUTPUT:
!      esv:            Dielectric constant for the leaves

! INTERNAL:
!     delta:           Ionic conductivity of aqueous solution
!     vmv:             Volume water content of vegetation
!     en               Dielectric constant for nondispersive residual component
!     vf               Volume water content of free water  
!     vb               Volume components of bound water and binding substances             
!----------------------------------------------------------------------------------

subroutine Canopy_Diel(frequency,mg,esv,rhoveg,t_soil)

  REAL(AUS) :: frequency, t_soil, mg, md, en, vf, vb
  REAL(AUS) :: rhoveg, vmv, delta
  COMPLEX(AUS) :: esv, xx, ef

  delta = 1.27_AUS
  md = ONE-mg
  vmv = mg*rhoveg/( ONE - mg*(ONE-rhoveg) )
  xx = CMPLX(ZERO,ONE,AUS)
  
  !-----------------Dielectric constant model estimation at room temperature----
  if (t_soil >= T0 + TS_LOW_THRESHOLD) then
    en = 1.7_AUS + (3.2_AUS + 6.5_AUS*vmv)*vmv
    vf = vmv*(0.82_AUS*vmv + 0.166_AUS)
    vb = 31.4_AUS*vmv*vmv/( ONE + 59.5_AUS*vmv*vmv)
    esv = en + vf*(4.9_AUS + 75.0_AUS/(ONE + xx*frequency/18.0_AUS)-xx*(18.0_AUS*delta/frequency)) + &
         vb*(2.9_AUS + 55.0_AUS/(ONE + SQRT(xx*frequency/0.18_AUS)))
  else
  !-----------------Dielectric constant model estimation at low temperature��-5��< =TS <=5��,Matzler et al,1994��----       
    ef = 4.9_AUS + 75.0_AUS/(ONE + xx*frequency/18.0_AUS)- xx*(18.0_AUS*delta/frequency) 
    esv = 0.522_AUS*(ONE-1.32_AUS*md)*ef + 0.51_AUS + 3.84_AUS*md
  endif
  
  if (AIMAG(esv) >= ZERO) esv = CMPLX(REAL(esv), -0.0001_AUS, AUS)

end subroutine Canopy_Diel

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Smooth_Reflectance
! PURPOSE:
!       Compute the smooth surface reflectivity using Fresnel equation
! INPUT:
!      theta_i:        Incident angle (degree)
!      theta_t:        Transmitted angle (degree)
!      em1:            Dielectric constant of the medium 1
!      em2:            Dielectric constant of the medium 2
!
! OUTPUT:
!      rv:             Reflectivity at vertical polarization
!      rh:             Reflectivity at horizontal polarization
!
! INTERNAL:
!     m1:             Refractive index of the medium 1
!     m2:             Refractive index of the medium 2
!----------------------------------------------------------------------------------

subroutine Smooth_Reflectance(em1, em2, theta_i, theta_t, rv, rh)

  REAL(AUS) :: theta_i, theta_t
  REAL(AUS) :: rh, rv,cos_i,cos_t
  COMPLEX(AUS) :: em1, em2, m1, m2, angle_i, angle_t

  ! compute the refractive index ratio between medium 2 and 1
  ! using dielectric constant (n = SQRT(e))
  cos_i = COS(theta_i)
  cos_t = COS(theta_t)

  angle_i = CMPLX(cos_i, ZERO, AUS)
  angle_t = CMPLX(cos_t, ZERO, AUS)

  m1 = SQRT(em1)
  m2 = SQRT(em2)

  rv = (ABS((m1*angle_t-m2*angle_i)/(m1*angle_t+m2*angle_i)))**2
  rh = (ABS((m1*angle_i-m2*angle_t)/(m1*angle_i+m2*angle_t)))**2

end subroutine Smooth_Reflectance

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Transmittance
! PURPOSE:
!       Compute Transmittance  using Fresnel equation
! INPUT:
!      theta:          Local zenith angle (degree)
!      theta_i:        Incident angle (degree)
!      theta_t:        Transmitted angle (degree)
!      em1:            Dielectric constant of the medium 1
!      em2:            Dielectric constant of the medium 2
!
! OUTPUT:
!      tv:             Transmisivity at vertical polarization
!      th:             Transmisivity at horizontal polarization
! INTERNAL:
!     m1:             Refractive index of the medium 1
!     m2:             Refractive index of the medium 2
!
!----------------------------------------------------------------------------------

subroutine Transmittance(em1,em2,theta_i,theta_t,tv,th)

  REAL(AUS) :: theta_i, theta_t
  REAL(AUS) :: th, tv, rr, cos_i,cos_t
  COMPLEX(AUS) :: em1, em2, m1, m2, angle_i, angle_t

  ! compute the refractive index ratio between medium 2 and 1
  ! using dielectric constant (n = SQRT(e))
  cos_i = COS(theta_i)
  cos_t = COS(theta_t)

  angle_i = CMPLX(cos_i, ZERO, AUS)
  angle_t = CMPLX(cos_t, ZERO, AUS)

  m1 = SQRT(em1)
  m2 = SQRT(em2)

  rr = ABS(m2/m1)*cos_t/cos_i
  tv = rr*(ABS(TWO*m1*angle_i/(m1*angle_t + m2*angle_i)))**2
  th = rr*(ABS(TWO*m1*angle_i/(m1*angle_i + m2*angle_t)))**2

end subroutine Transmittance

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Roughness_Reflectance_ChenWeng
! PURPOSE:
!       Compute surface relectivity using Chen-Weng model
! REFERENCES:
!    Chen M ,  Weng F . 2016.Modeling Land Surface Roughness Effect on Soil Microwave Emission in Community Surface Emissivity Model[J].
!             IEEE Transactions on Geoscience & Remote Sensing, 54(3):1716-1726
! INPUT:
!      frequency:      Frequency (GHz)
!      theta:          Local zenith angle (degree)
!      sigma:          Standard deviation of rough surface height
!                      smooth surface:0.38, medium: 1.10, rough:2.15 cm
!      rv:             Reflectivity at vertical polarization for smooth surface
!      rh:             Reflectivity at horizontal polarization for smooth surface
!
! OUTPUT:
!      rv:             Reflectivity at vertical polarization for roughness surface
!      rh:             Reflectivity at horizontal polarization for roughness surface

!INTERNAl:
!       f:             frequency (cm) 
!     lamda:           wave length (cm) 
!       k:             wave number (1/cm)
!       Q:             Polarization coupling factor, which describes the effect of roughness on orthogonal polarization
!       P:             Roughness attenuation factor, which describes the impact of surface roughness on reflection attenuation
!    rv_s:             Reflectivity at vertical polarization for smooth surface
!    rh_s:             Reflectivity at horizontal polarization for smooth surface     
!----------------------------------------------------------------------------------

subroutine Roughness_Reflectance_ChenWeng(frequency,sigma,theta, rv,rh)

  REAL(AUS) :: frequency,  sigma,theta
  REAL(AUS) :: f, lamda, k, a1, a2, a3, k1, k2, sigma1, sigma2, b1, b2, temp1, temp2, temp1_tan, temp2_tan
  REAL(AUS) :: Q, P, rh, rv, rh_s, rv_s

  f = frequency*1.0e9_AUS  
  lamda = 3.0e8_AUS/f*100_AUS    
  k = TWOPI/lamda  
  
  a1 = 1.4_AUS
  a2 = 0.15_AUS
  a3 = EXP(-(TWOPI*sigma/5.0_AUS)**sqrt(POINT1*cos(theta)))
  
  k1 = ZERO            
  k2 = TWOPI/3.0_AUS    
  sigma1 = TWOPI*0.3_AUS   
  sigma2 = TWOPI*0.5_AUS
  b1 = 0.43_AUS
  b2 = 15.0_AUS
  
  Q = b1*(ONE - EXP(-b2*k*sigma*(COS(theta))**2.0))  
  temp1 = (k*sigma*COS(theta) - k1*sigma*COS(theta))/sigma1
  temp2 = (k*sigma*COS(theta) - k2*sigma*COS(theta))/sigma2
  temp1_tan = (EXP(TWO*temp1) - ONE)/(EXP(TWO*temp1)+ONE)
  temp2_tan = (EXP(TWO*temp2) - ONE)/(EXP(TWO*temp2)+ONE)
  
  P = (a1+a3)/TWO+(a2-a1)/2*temp1_tan + (a3-a2)/2*temp2_tan 
  rv_s = rv
  rh_s = rh
  rv = ((1-Q)*rv_s + Q*rh_s)*P
  rh = ((1-Q)*rh_s + Q*rv_s)*P
  
end subroutine Roughness_Reflectance_ChenWeng

!--------------------------------------------------------------------------------
!
!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Roughness_Reflectance_Qh
! PURPOSE:
!       Compute surface relectivity
! INPUT:
!      frequency:      Frequency (GHz)
!      theta:          Local zenith angle (degree)
!      sigma:          Standard deviation of rough surface height
!                      smooth surface:0.38, medium: 1.10, rough:2.15 cm
!
! OUTPUT:
!      rv:             Reflectivity at vertical polarization
!      rh:             Reflectivity at horizontal polarization
! REFERENCES:
!    Wigneron, J P, Kerr Y, Waldteufel P��2007. L-band microwave emission of the biosphere (L-MEB) model:
!          Description and calibration against experimental data sets over crop fields. Remote Sensing of Environment, 107, 639�C655
!    Wang, j. and b. j. choudhury, 1992: passive microwave radiation from soil: examples...
!    passive microwave remote sensing of .. ed. b. j. choudhury, etal vsp.
!    also wang and choudhury (1982)
!
!----------------------------------------------------------------------------------

subroutine Roughness_Reflectance_Qh(frequency,sigma,rv,rh)

  REAL(AUS) :: frequency
  REAL(AUS) :: q, rh, rv, rh_s, rv_s, sigma

  rh_s = 0.3_AUS*rh
  rv_s = 0.3_AUS*rv
  q = 0.35_AUS*(ONE - EXP(-0.60_AUS*frequency*sigma**TWO))
  rh = rh_s + q*(rv_s - rh_s)
  rv = rv_s + q*(rh_s - rv_s)

end subroutine Roughness_Reflectance_Qh

!--------------------------------------------------------------------------------
!
! SUBROUTINE NAME:
!       Two_Stream_Solution
! PURPOSE:
!       Two stream solution 
!       Updated with the more accurate formula of total upwelling radiance emanating from the surface.
! INPUT:
!      b:              Scattering layer temperature (k)         (gdas)   (not used here)
!      mu:             cos(theta)
!      gv:             Asymmetry factor for v pol
!      gh:             Asymmetry factor for h pol
!      ssalb_v:        Single scattering albedo at v. polarization
!      ssalb_h:        Single scattering albedo at h. polarization
!      tau_v:          Optical depth at v. polarization
!      tau_h:          Optical depth at h. polarization
!      r12_v:          Reflectivity at vertical polarization   (not used here)
!      r12_h:          Reflectivity at horizontal polarization (not used here)
!      r21_v:          Reflectivity at vertical polarization
!      r21_h:          Reflectivity at horizontal polarization
!      r23_v:          Reflectivity at vertical polarization
!      r23_h:          Reflectivity at horizontal polarization
!      t21_v:          Transmisivity at vertical polarization
!      t21_h:          Transmisivity at horizontal polarization
!      t12_v:          Transmisivity at vertical polarization   (not used here)
!      t12_h:          Transmisivity at horizontal polarization (not used here)
!      Frequency:      Frequency
!      t_soil:         Soil temperature (K)
!      t_skin:         Land surface temperature (K)
!
! OUTPUT:
!      esv:             Emissivity at vertical polarization
!      esh:             Emissivity at horizontal polarization
! REFERENCES:
!    Weng, F., B. Yan, and N. Grody, 2001: "A microwave land emissivity model",
!     J. Geophys. Res., 106, 20, 115-20, 123
!
!----------------------------------------------------------------------------------

subroutine Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
      r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,esv,esh,frequency,t_soil,t_skin)


  REAL(AUS) :: mu, gv, gh, ssalb_h, ssalb_v, tau_h,tau_v,                 &
              r21_h, r21_v, r23_h, r23_v, t21_v, t21_h, esv, esh
  REAL(AUS) :: alfa_v, alfa_h, kk_h, kk_v, gamma_h, gamma_v, beta_v, beta_h
  REAL(AUS) :: fact1,fact2
  REAL(AUS) :: frequency, t_soil, t_skin
  REAL(AUS) :: gsect0, gsect1_h, gsect1_v, gsect2_h, gsect2_v

  alfa_h  = SQRT((ONE - ssalb_h)/(ONE - gh*ssalb_h))
  kk_h    = SQRT((ONE - ssalb_h)*(ONE -  gh*ssalb_h))/mu
  beta_h  = (ONE - alfa_h)/(ONE + alfa_h)
  gamma_h = (beta_h -r23_h)/(ONE-beta_h*r23_h)

  alfa_v  = SQRT((ONE-ssalb_v)/(ONE - gv*ssalb_v))
  kk_v    = SQRT((ONE-ssalb_v)*(ONE - gv*ssalb_v))/mu
  beta_v  = (ONE - alfa_v)/(ONE + alfa_v)
  gamma_v = (beta_v - r23_v)/(ONE - beta_v*r23_v)

  fact1=gamma_h*EXP(-TWO*kk_h*tau_h)
  fact2=gamma_v*EXP(-TWO*kk_v*tau_v)

  gsect0  =(EXP(C_2*frequency/t_skin) -ONE)/(EXP(C_2*frequency/t_soil) -ONE)

  gsect1_h=(ONE-r23_h)*(gsect0-ONE)
  gsect2_h=((ONE-beta_h*beta_h)/(ONE-beta_h*r23_h))*EXP(-kk_h*tau_h)

  gsect1_v=(ONE-r23_v)*(gsect0-ONE)
  gsect2_v=((ONE-beta_v*beta_v)/(ONE-beta_v*r23_v))*EXP(-kk_h*tau_v)

  esh  = t21_h*((ONE - beta_h)*(ONE + fact1)+gsect1_h*gsect2_h) /(ONE-beta_h*r21_h-(beta_h-r21_h)*fact1)
  esv  = t21_v*((ONE - beta_v)*(ONE + fact2)+gsect1_v*gsect2_v) /(ONE-beta_v*r21_v-(beta_v-r21_v)*fact2)

  if (esh < EMISSH_DEFAULT) esh = EMISSH_DEFAULT
  if (esv < EMISSV_DEFAULT) esv = EMISSV_DEFAULT

  if (esh > ONE) esh = ONE
  if (esv > ONE) esv = ONE

end subroutine Two_Stream_Solution

END MODULE ARMS_SfcEM_MWLand
