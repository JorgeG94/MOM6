!> Initial conditions for the 2D Rossby front test
module Rossby_front_2d_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use regrid_consts, only : coordinateMode, DEFAULT_COORDINATE_MODE
use regrid_consts, only : REGRIDDING_LAYER, REGRIDDING_ZSTAR
use regrid_consts, only : REGRIDDING_RHO, REGRIDDING_SIGMA

use MOM_datatypes, only : wp

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mdl = "Rossby_front_2d_initialization" !< This module's name.
! This include declares and sets the variable "version".
#include "version_variable.h"

public Rossby_front_initialize_thickness
public Rossby_front_initialize_temperature_salinity
public Rossby_front_initialize_velocity

! Parameters defining the initial conditions of this test case
real(wp), parameter :: frontFractionalWidth = 0.5_wp !< Width of front as fraction of domain [nondim]
real(wp), parameter :: HMLmin = 0.25_wp !< Shallowest ML as fractional depth of ocean [nondim]
real(wp), parameter :: HMLmax = 0.75_wp !< Deepest ML as fractional depth of ocean [nondim]

contains

!> Initialization of thicknesses in 2D Rossby front test
subroutine Rossby_front_initialize_thickness(h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G           !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV          !< Vertical grid structure
  type(unit_scale_type),   intent(in)  :: US          !< A dimensional unit scaling type
  real(wp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: h           !< The thickness that is being initialized [H ~> m or kg m-2]
  type(param_file_type),   intent(in)  :: param_file  !< A structure indicating the open file
                                                      !! to parse for model parameter values.
  logical,                 intent(in)  :: just_read   !< If true, this call will only read
                                                      !! parameters without changing h.

  ! Local variables
  real(wp)    :: Tz         ! Vertical temperature gradient [C H-1 ~> degC m2 kg-1]
  real(wp)    :: Dml        ! Mixed layer depth [H ~> m or kg m-2]
  real(wp)    :: eta        ! An interface height depth [H ~> m or kg m-2]
  real(wp)    :: stretch    ! A nondimensional stretching factor [nondim]
  real(wp)    :: h0         ! The stretched thickness per layer [H ~> m or kg m-2]
  real(wp)    :: T_range    ! Range of temperatures over the vertical [C ~> degC]
  real(wp)    :: dRho_dT    ! The partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real(wp)    :: max_depth  ! Maximum depth of the model bathymetry [H ~> m or kg m-2]
  character(len=40) :: verticalCoordinate
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.just_read) &
    call MOM_mesg("Rossby_front_2d_initialization.F90, Rossby_front_initialize_thickness: setting thickness")

  if (.not.just_read) call log_version(param_file, mdl, version, "")
  ! Read parameters needed to set thickness
  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units="degC", default=0.0_wp, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "DRHO_DT", dRho_dT, &
                 units="kg m-3 degC-1", default=-0.2_wp, scale=US%kg_m3_to_R*US%C_to_degC, do_not_log=.true.)
  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9_wp, scale=GV%m_to_H, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0_wp) call MOM_error(FATAL, &
      "Rossby_front_initialize_thickness, Rossby_front_initialize_thickness: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")

  Tz = T_range / max_depth

  if (GV%Boussinesq) then
    select case ( coordinateMode(verticalCoordinate) )

      case (REGRIDDING_LAYER, REGRIDDING_RHO)
        ! This code is identical to the REGRIDDING_ZSTAR case but probably should not be.
        do j = G%jsc,G%jec ; do i = G%isc,G%iec
          Dml = Hml( G, G%geoLatT(i,j), max_depth )
          eta = -( -dRho_dT / GV%Rho0 ) * Tz * 0.5_wp * ( Dml * Dml )
          stretch = ( ( max_depth + eta ) / max_depth )
          h0 = ( max_depth / real(nz, wp) ) * stretch
          do k = 1, nz
            h(i,j,k) = h0
          enddo
        enddo ; enddo

      case (REGRIDDING_ZSTAR, REGRIDDING_SIGMA)
        do j = G%jsc,G%jec ; do i = G%isc,G%iec
          Dml = Hml( G, G%geoLatT(i,j), max_depth )
          ! The free surface height is set so that the bottom pressure gradient is 0.
          eta = -( -dRho_dT / GV%Rho0 ) * Tz * 0.5_wp * ( Dml * Dml )
          stretch = ( ( max_depth + eta ) / max_depth )
          h0 = ( max_depth / real(nz, wp) ) * stretch
          do k = 1, nz
            h(i,j,k) = h0
          enddo
        enddo ; enddo

      case default
        call MOM_error(FATAL,"Rossby_front_initialize: "// &
        "Unrecognized i.c. setup - set REGRIDDING_COORDINATE_MODE")

    end select
  else
    ! In non-Boussinesq mode with a flat bottom, the only requirement for no bottom pressure
    ! gradient and no abyssal flow is that all columns have the same mass.
    h0 = max_depth / real(nz, wp)
    do k=1,nz ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
      h(i,j,k) = h0
    enddo ; enddo ; enddo
  endif

end subroutine Rossby_front_initialize_thickness


!> Initialization of temperature and salinity in the Rossby front test
subroutine Rossby_front_initialize_temperature_salinity(T, S, h, G, GV, US, &
                   param_file, just_read)
  type(ocean_grid_type),                     intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),                   intent(in)  :: GV !< The ocean's vertical grid structure.
  real(wp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: T  !< Potential temperature [C ~> degC]
  real(wp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: S  !< Salinity [S ~> ppt]
  real(wp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(unit_scale_type),                     intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),                     intent(in)  :: param_file   !< Parameter file handle
  logical,                                   intent(in)  :: just_read !< If true, this call will
                                                      !! only read parameters without changing T & S.
  ! Local variables
  real(wp)      :: T_ref        ! Reference temperature within the surface layer [C ~> degC]
  real(wp)      :: S_ref        ! Reference salinity within the surface layer [S ~> ppt]
  real(wp)      :: T_range      ! Range of temperatures over the vertical [C ~> degC]
  real(wp)      :: zc           ! Position of the middle of the cell [H ~> m or kg m-2]
  real(wp)      :: zi           ! Bottom interface position relative to the sea surface [H ~> m or kg m-2]
  real(wp)      :: dTdz         ! Vertical temperature gradient [C H-1 ~> degC m-1 or degC m2 kg-1]
  real(wp)      :: Dml          ! Mixed layer depth [H ~> m or kg m-2]
  real(wp)      :: max_depth    ! Maximum depth of the model bathymetry [H ~> m or kg m-2]
  character(len=40) :: verticalCoordinate
  integer   :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call get_param(param_file, mdl,"REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 default=35.0_wp, units="ppt", scale=US%ppt_to_S, do_not_log=just_read)
  call get_param(param_file, mdl, "T_REF", T_ref, 'Reference temperature', &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units="degC", default=0.0_wp, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9_wp, scale=GV%m_to_H, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0_wp) call MOM_error(FATAL, &
      "Rossby_front_initialize_thickness, Rossby_front_initialize_temperature_salinity: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")

  T(:,:,:) = 0.0_wp
  S(:,:,:) = S_ref
  dTdz = T_range / max_depth

  ! This sets the temperature to the value at the base of the specified mixed layer
  ! depth from a horizontally uniform constant thermal stratification.
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
    zi = 0._wp
    Dml = Hml(G, G%geoLatT(i,j), max_depth)
    do k = 1, nz
      zi = zi - h(i,j,k)           ! Bottom interface position
      zc = zi - 0.5_wp*h(i,j,k)       ! Position of middle of cell
      T(i,j,k) = T_ref + dTdz * min( zc, -Dml ) ! Linear temperature profile below the mixed layer
    enddo
  enddo ; enddo

end subroutine Rossby_front_initialize_temperature_salinity


!> Initialization of u and v in the Rossby front test
subroutine Rossby_front_initialize_velocity(u, v, h, G, GV, US, param_file, just_read)
  type(ocean_grid_type),      intent(in)  :: G  !< Grid structure
  type(verticalGrid_type),    intent(in)  :: GV !< Vertical grid structure
  real(wp), dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                              intent(out) :: u  !< i-component of velocity [L T-1 ~> m s-1]
  real(wp), dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                              intent(out) :: v  !< j-component of velocity [L T-1 ~> m s-1]
  real(wp), dimension(SZI_(G),SZJ_(G), SZK_(GV)), &
                              intent(in)  :: h  !< Thickness [H ~> m or kg m-2]
  type(unit_scale_type),      intent(in)  :: US !< A dimensional unit scaling type
  type(param_file_type),      intent(in)  :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  logical,                    intent(in)  :: just_read !< If present and true, this call will only
                                                !! read parameters without setting u & v.

  real(wp)    :: T_range      ! Range of temperatures over the vertical [C ~> degC]
  real(wp)    :: T_ref        ! Reference temperature within the surface layer [C ~> degC]
  real(wp)    :: S_ref        ! Reference salinity within the surface layer [S ~> ppt]
  real(wp)    :: dUdT         ! Factor to convert dT/dy into dU/dz, g*alpha/f with rescaling
                          ! [L2 H-1 T-1 C-1 ~> m s-1 degC-1 or m4 kg-1 s-1 degC-1]
  real(wp)    :: Rho_T0_S0    ! The density at T=0, S=0 [R ~> kg m-3]
  real(wp)    :: dRho_dT      ! The partial derivative of density with temperature [R C-1 ~> kg m-3 degC-1]
  real(wp)    :: dRho_dS      ! The partial derivative of density with salinity [R S-1 ~> kg m-3 ppt-1]
  real(wp)    :: dSpV_dT      ! The partial derivative of specific volume with temperature [R-1 C-1 ~> m3 kg-1 degC-1]
  real(wp)    :: T_here       ! The temperature in the middle of a layer [C ~> degC]
  real(wp)    :: dTdz         ! Vertical temperature gradient [C H-1 ~> degC m-1 or degC m2 kg-1]
  real(wp)    :: Dml          ! Mixed layer depth [H ~> m or kg m-2]
  real(wp)    :: zi, zc, zm   ! Depths in thickness units [H ~> m or kg m-2].
  real(wp)    :: f            ! The local Coriolis parameter [T-1 ~> s-1]
  real(wp)    :: I_f          ! The Adcroft reciprocal of the local Coriolis parameter [T ~> s]
  real(wp)    :: Ty           ! The meridional temperature gradient [C L-1 ~> degC m-1]
  real(wp)    :: hAtU         ! Interpolated layer thickness in height units [H ~> m or kg m-2].
  real(wp)    :: u_int        ! The zonal velocity at an interface [L T-1 ~> m s-1]
  real(wp)    :: max_depth    ! Maximum depth of the model bathymetry [H ~> m or kg m-2]
  integer :: i, j, k, is, ie, js, je, nz
  character(len=40) :: verticalCoordinate

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call get_param(param_file, mdl, "REGRIDDING_COORDINATE_MODE", verticalCoordinate, &
                 default=DEFAULT_COORDINATE_MODE, do_not_log=just_read)
  call get_param(param_file, mdl, "T_RANGE", T_range, 'Initial temperature range', &
                 units="degC", default=0.0_wp, scale=US%degC_to_C, do_not_log=just_read)
  call get_param(param_file, mdl, "S_REF", S_ref, 'Reference salinity', &
                 default=35.0_wp, units="ppt", scale=US%ppt_to_S, do_not_log=.true.)
  call get_param(param_file, mdl, "T_REF", T_ref, 'Reference temperature', &
                 units="degC", scale=US%degC_to_C, fail_if_missing=.not.just_read, do_not_log=.true.)
  call get_param(param_file, mdl, "RHO_T0_S0", Rho_T0_S0, &
                 units="kg m-3", default=1000.0_wp, scale=US%kg_m3_to_R, do_not_log=.true.)
  call get_param(param_file, mdl, "DRHO_DT", dRho_dT, &
                 units="kg m-3 degC-1", default=-0.2_wp, scale=US%kg_m3_to_R*US%C_to_degC, do_not_log=.true.)
  call get_param(param_file, mdl, "DRHO_DS", dRho_dS, &
                 units="kg m-3 ppt-1", default=0.8_wp, scale=US%kg_m3_to_R*US%S_to_ppt, do_not_log=.true.)
  call get_param(param_file, mdl, "MAXIMUM_DEPTH", max_depth, &
                 units="m", default=-1.e9_wp, scale=GV%m_to_H, do_not_log=.true.)

  if (just_read) return ! All run-time parameters have been read, so return.

  if (max_depth <= 0.0_wp) call MOM_error(FATAL, &
      "Rossby_front_initialize_thickness, Rossby_front_initialize_velocity: "//&
      "This module requires a positive value of MAXIMUM_DEPTH.")
  if (G%grid_unit_to_L <= 0._wp) call MOM_error(FATAL, 'Rossby_front_2d_initialization.F90: '// &
          "dTdy() is only set to work with Cartesian axis units.")

  v(:,:,:) = 0.0_wp
  u(:,:,:) = 0.0_wp

  if (GV%Boussinesq) then
    do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
      f = 0.5_wp* (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1) )
      dUdT = 0.0_wp ; if (abs(f) > 0.0_wp) &
        dUdT = ( GV%H_to_Z*GV%g_Earth*dRho_dT ) / ( f * GV%Rho0 )
      Dml = Hml( G, G%geoLatCu(I,j), max_depth )
      Ty = dTdy( G, T_range, G%geoLatCu(I,j), US )
      zi = 0._wp
      do k = 1, nz
        hAtU = 0.5_wp * (h(i,j,k) + h(i+1,j,k))
        zi = zi - hAtU             ! Bottom interface position
        zc = zi - 0.5_wp*hAtU         ! Position of middle of cell
        zm = max( zc + Dml, 0._wp )    ! Height above bottom of mixed layer
        u(I,j,k) = dUdT * Ty * zm   ! Thermal wind starting at base of ML
      enddo
    enddo ; enddo
  else
    ! With an equation of state that is linear in density, the nonlinearies in
    ! specific volume require that temperature be calculated for each layer.

    dTdz = T_range / max_depth

    do j = G%jsc,G%jec ; do I = G%isc-1,G%iec+1
      f = 0.5_wp* (G%CoriolisBu(I,J) + G%CoriolisBu(I,J-1) )
      I_f = 0.0_wp ; if (abs(f) > 0.0_wp) I_f = 1.0_wp / f
      Dml = Hml( G, G%geoLatCu(I,j), max_depth )
      Ty = dTdy( G, T_range, G%geoLatCu(I,j), US )
      zi = -max_depth
      u_int = 0.0_wp ! The velocity at an interface
      ! Work upward in non-Boussinesq mode
      do k = nz, 1, -1
        hAtU = 0.5_wp * (h(i,j,k) + h(i+1,j,k))
        zc = zi + 0.5_wp*hAtU         ! Position of middle of cell
        T_here = T_ref + dTdz * min(zc, -Dml) ! Linear temperature profile below the mixed layer
        dSpV_dT = -dRho_dT / (Rho_T0_S0 + (dRho_dS * S_ref + dRho_dT * T_here) )**2
        dUdT = -( GV%H_to_RZ * GV%g_Earth * dSpV_dT ) * I_f

        ! There is thermal wind shear only within the mixed layer.
        u(I,j,k) = u_int + dUdT * Ty * min(max((zi + Dml) + 0.5_wp*hAtU, 0.0_wp), 0.5_wp*hAtU)
        u_int = u_int + dUdT * Ty * min(max((zi + Dml) + hAtU, 0.0_wp), hAtU)

        zi = zi + hAtU             ! Update the layer top interface position
      enddo
    enddo ; enddo
  endif
end subroutine Rossby_front_initialize_velocity

!> Pseudo coordinate across domain used by Hml() and dTdy()
!! returns a coordinate from -PI/2 .. PI/2 squashed towards the
!! center of the domain [radians].
real(wp) function yPseudo( G, lat )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real(wp),                  intent(in) :: lat !< Latitude in arbitrary units, often [km]
  ! Local
  real(wp) :: PI   ! The ratio of the circumference of a circle to its diameter [nondim]

  PI = 4.0_wp * atan(1.0_wp)
  yPseudo = ( ( lat - G%south_lat ) / G%len_lat ) - 0.5_wp ! -1/2 .. 1/.2
  yPseudo = PI * max(-0.5_wp, min(0.5_wp, yPseudo / frontFractionalWidth))
end function yPseudo


!> Analytic prescription of mixed layer depth in 2d Rossby front test,
!! in the same units as max_depth (usually [Z ~> m] or [H ~> m or kg m-2])
real(wp) function Hml( G, lat, max_depth )
  type(ocean_grid_type), intent(in) :: G   !< Grid structure
  real(wp),                  intent(in) :: lat !< Latitude in arbitrary units, often [km]
  real(wp),                  intent(in) :: max_depth !< The maximum depth of the ocean [Z ~> m] or [H ~> m or kg m-2]
  ! Local
  real(wp) :: dHML, HMLmean ! The range and mean of the mixed layer depths [Z ~> m] or [H ~> m or kg m-2]

  dHML = 0.5_wp * ( HMLmax - HMLmin ) * max_depth
  HMLmean = 0.5_wp * ( HMLmin + HMLmax ) * max_depth
  Hml = HMLmean + dHML * sin( yPseudo(G, lat) )
end function Hml


!> Analytic prescription of mixed layer temperature gradient in [C L-1 ~> degC m-1] in 2d Rossby front test
real(wp) function dTdy( G, dT, lat, US )
  type(ocean_grid_type), intent(in) :: G     !< Grid structure
  real(wp),                  intent(in) :: dT    !< Top to bottom temperature difference [C ~> degC]
  real(wp),                  intent(in) :: lat   !< Latitude in the same units as geoLat, often [km]
  type(unit_scale_type), intent(in) :: US    !< A dimensional unit scaling type
  ! Local
  real(wp) :: PI   ! The ratio of the circumference of a circle to its diameter [nondim]
  real(wp) :: dHML ! The range of the mixed layer depths [Z ~> m]
  real(wp) :: dHdy ! The mixed layer depth gradient [Z L-1 ~> m m-1]

  PI = 4.0_wp * atan(1.0_wp)
  dHML = 0.5_wp * ( HMLmax - HMLmin ) * G%max_depth
  dHdy = dHML * ( PI / ( frontFractionalWidth * G%len_lat * G%grid_unit_to_L ) ) * cos( yPseudo(G, lat) )
  dTdy = -( dT / G%max_depth ) * dHdy

end function dTdy


!> \namespace rossby_front_2d_initialization
!!
!! \section section_Rossby_front_2d Description of the 2d Rossby front initial conditions
!!
!! Consistent with a linear equation of state, the system has a constant stratification
!! below the mixed layer, stratified in temperature only. Isotherms are flat below the
!! mixed layer and vertical within. Salinity is constant. The mixed layer has a half sine
!! form so that there are no mixed layer or temperature gradients at the side walls.
!!
!! Below the mixed layer the potential temperature, \f$\theta(z)\f$, is given by
!! \f[ \theta(z) = \theta_0 - \Delta \theta \left( z + h_{ML} \right) \f]
!! where \f$ \theta_0 \f$ and \f$ \Delta \theta \f$ are external model parameters.
!!
!! The depth of the mixed layer, \f$H_{ML}\f$ is
!! \f[ h_{ML}(y) = h_{min} + \left( h_{max} - h_{min} \right) \cos{\pi y/L} \f].
!! The temperature in mixed layer is given by the reference temperature at \f$z=h_{ML}\f$
!! so that
!! \f{eqnarray} \theta(y,z) =
!!     \theta_0 - \Delta \theta \left( z + h_{ML} \right) & \forall & z < h_{ML}(y) T.B.D.
!! \f}

end module Rossby_front_2d_initialization
