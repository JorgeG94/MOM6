!> Initialization for the "Neverworld" configuration
module Neverworld_initialization

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_sponge, only : sponge_CS, set_up_sponge_field, initialize_sponge
use MOM_dyn_horgrid, only : dyn_horgrid_type
use MOM_error_handler, only : MOM_mesg, MOM_error, FATAL, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_get_input, only : directories
use MOM_grid, only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

use random_numbers_mod, only: initializeRandomNumberStream, getRandomNumbers, randomNumberStream

use MOM_datatypes, only : wp

implicit none ; private

#include <MOM_memory.h>

public Neverworld_initialize_topography
public Neverworld_initialize_thickness

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> This subroutine sets up the Neverworld test case topography.
subroutine Neverworld_initialize_topography(D, G, param_file, max_depth)
  type(dyn_horgrid_type),  intent(in)  :: G !< The dynamic horizontal grid type
  real(wp), dimension(G%isd:G%ied,G%jsd:G%jed), &
                           intent(out) :: D !< Ocean bottom depth in the units of depth_max [A]
  type(param_file_type),   intent(in)  :: param_file !< Parameter file structure
  real(wp),                    intent(in)  :: max_depth !< Maximum ocean depth in arbitrary units [A]

  ! Local variables
  real(wp) :: PI                   ! 3.1415926... calculated as 4*atan(1) [nondim]
  real(wp) :: x, y ! Lateral positions normalized by the domain size [nondim]
  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "Neverworld_initialize_topography" ! This subroutine's name.
  real(wp) :: nl_top_amp       ! Amplitude of large-scale topographic features as a fraction of the maximum depth [nondim]
  real(wp) :: nl_roughness_amp ! Amplitude of topographic roughness as a fraction of the maximum depth [nondim]
  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_topography: setting topography", 5)

  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "NL_ROUGHNESS_AMP", nl_roughness_amp, &
                 "Amplitude of wavy signal in bathymetry.", units="nondim", default=0.05_wp)
  call get_param(param_file, mdl, "NL_CONTINENT_AMP", nl_top_amp, &
                 "Scale factor for topography - 0.0 for no continents.", units="nondim", default=1.0_wp)

  PI = 4.0_wp*atan(1.0_wp)

!  Calculate the depth of the bottom.
  do j=js,je ; do i=is,ie
    x = (G%geoLonT(i,j)-G%west_lon) / G%len_lon
    y = (G%geoLatT(i,j)-G%south_lat) / G%len_lat
!  This sets topography that has a reentrant channel to the south.
    D(i,j) = 1.0_wp - 1.1_wp * spike(y-1,0.12_wp) - 1.1_wp * spike(y,0.12_wp) - & !< The great northern wall and Antarctica
              nl_top_amp*( &
                (1.2_wp * spike(x,0.2_wp) + 1.2_wp * spike(x-1.0_wp,0.2_wp)) * spike(MIN(0.0_wp,y-0.3_wp),0.2_wp) & !< South America
              +  1.2_wp * spike(x-0.5_wp,0.2_wp) * spike(MIN(0.0_wp,y-0.55_wp),0.2_wp)       & !< Africa
              +  1.2_wp * (spike(x,0.12_wp)  + spike(x-1,0.12_wp)) * spike(MAX(0.0_wp,y-0.06_wp),0.12_wp)    & !< Antarctic Peninsula
              +  0.1_wp * (cosbell(x,0.1_wp) + cosbell(x-1,0.1_wp))                 & !< Drake Passage ridge
              +  0.5_wp * cosbell(x-0.16_wp,0.05_wp) * (cosbell(y-0.18_wp,0.13_wp)**0.4_wp)  & !< Scotia Arc East
              +  0.4_wp * (cosbell(x-0.09_wp,0.08_wp)**0.4_wp) * cosbell(y-0.26_wp,0.05_wp)  & !< Scotia Arc North
              +  0.4_wp * (cosbell(x-0.08_wp,0.08_wp)**0.4_wp) * cosbell(y-0.1_wp,0.05_wp))   & !< Scotia Arc South
              -  nl_roughness_amp * cos(14*PI*x) * sin(14*PI*y)            & !< roughness
              -  nl_roughness_amp * cos(20*PI*x) * cos(20*PI*y)              !< roughness
    if (D(i,j) < 0.0_wp) D(i,j) = 0.0_wp
    D(i,j) = D(i,j) * max_depth
  enddo ; enddo

end subroutine Neverworld_initialize_topography

!> Returns the value of a cosine-bell function evaluated at x/L [nondim]
real(wp) function cosbell(x, L)
  real(wp) , intent(in) :: x       !< Position in arbitrary units [A]
  real(wp) , intent(in) :: L       !< Width in arbitrary units [A]
  real(wp)              :: PI      !< 3.1415926... calculated as 4*atan(1) [nondim]

  PI      = 4.0_wp*atan(1.0_wp)
  cosbell = 0.5_wp * (1 + cos(PI*MIN(ABS(x/L),1.0_wp)))
end function cosbell

!> Returns the value of a sin-spike function evaluated at x/L [nondim]
real(wp) function spike(x, L)

  real(wp) , intent(in) :: x       !< Position in arbitrary units [A]
  real(wp) , intent(in) :: L       !< Width in arbitrary units [A]
  real(wp)              :: PI      !< 3.1415926... calculated as 4*atan(1) [nondim]

  PI    = 4.0_wp*atan(1.0_wp)
  spike = (1 - sin(PI*MIN(ABS(x/L),0.5_wp)))
end function spike

!> Returns the value of a triangular function centered at x=x0 with value 1
!! and linearly decreasing to 0 at x=x0+/-L, and 0 otherwise [nondim].
!! If clip is present the top of the cone is cut off at "clip", which
!! effectively defaults to 1.
real(wp) function cone(x, x0, L, clip)
  real(wp),           intent(in) :: x    !< Coordinate in arbitrary units [A]
  real(wp),           intent(in) :: x0   !< position of peak in arbitrary units [A]
  real(wp),           intent(in) :: L    !< half-width of base of cone in arbitrary units [A]
  real(wp), optional, intent(in) :: clip !< clipping height of cone [nondim]

  cone = max( 0._wp, 1._wp - abs(x - x0) / L )
  if (present(clip)) cone = min(clip, cone)
end function cone

!> Returns an s-curve s(x) s.t. s(x0)<=0, s(x0+L)>=1 and cubic in between [nondim].
real(wp) function scurve(x, x0, L)
  real(wp), intent(in) :: x       !< Coordinate in arbitrary units [A]
  real(wp), intent(in) :: x0      !< position of peak in arbitrary units [A]
  real(wp), intent(in) :: L       !< half-width of base of cone in arbitrary units [A]
  real(wp) :: s ! A rescaled position [nondim]

  s = max( 0._wp, min( 1._wp,( x - x0 ) / L ) )
  scurve = ( 3._wp - 2._wp*s ) * ( s * s )
end function scurve

! None of the following 7 functions appear to be used.

!> Returns a "coastal" profile [nondim].
real(wp) function cstprof(x, x0, L, lf, bf, sf, sh)
  real(wp), intent(in) :: x       !< Coordinate in arbitrary units [A]
  real(wp), intent(in) :: x0      !< position of peak in arbitrary units [A]
  real(wp), intent(in) :: L       !< width of profile in arbitrary units [A]
  real(wp), intent(in) :: lf      !< fraction of width that is "land" [nondim]
  real(wp), intent(in) :: bf      !< fraction of width that is "beach" [nondim]
  real(wp), intent(in) :: sf      !< fraction of width that is "continental slope" [nondim]
  real(wp), intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real(wp) :: s ! A rescaled position [nondim]

  s = max( 0._wp, min( 1._wp,( x - x0 ) / L ) )
  cstprof = sh * scurve(s-lf,0._wp,bf) + (1._wp-sh) * scurve(s - (1._wp-sf),0._wp,sf)
end function cstprof

!> Distance between points x,y and a line segment (x0,y0) and (x0,y1) in arbitrary units [A].
real(wp) function dist_line_fixed_x(x, y, x0, y0, y1)
  real(wp), intent(in) :: x       !< X-coordinate in arbitrary units [A]
  real(wp), intent(in) :: y       !< Y-coordinate in arbitrary units [A]
  real(wp), intent(in) :: x0      !< x-position of line segment in arbitrary units [A]
  real(wp), intent(in) :: y0      !< y-position of line segment end in arbitrary units [A]
  real(wp), intent(in) :: y1      !< y-position of line segment end in arbitrary units [A]
  real(wp) :: dx, yr, dy ! Relative positions in arbitrary units [A]

  dx = x - x0
  yr = min( max(y0,y1), max( min(y0,y1), y ) ) ! bound y by y0,y1
  dy = y - yr ! =0 within y0<y<y1, =y0-y for y<y0, =y-y1 for y>y1
  dist_line_fixed_x = sqrt( (dx*dx) + (dy*dy) )
end function dist_line_fixed_x

!> Distance between points x,y and a line segment (x0,y0) and (x1,y0) in arbitrary units [A].
real(wp) function dist_line_fixed_y(x, y, x0, x1, y0)
  real(wp), intent(in) :: x       !< X-coordinate in arbitrary units [A]
  real(wp), intent(in) :: y       !< Y-coordinate in arbitrary units [A]
  real(wp), intent(in) :: x0      !< x-position of line segment end in arbitrary units [A]
  real(wp), intent(in) :: x1      !< x-position of line segment end in arbitrary units [A]
  real(wp), intent(in) :: y0      !< y-position of line segment in arbitrary units [A]

  dist_line_fixed_y = dist_line_fixed_x(y, x, y0, x0, x1)
end function dist_line_fixed_y

!> A "coast profile" applied in an N-S line from lon0,lat0 to lon0,lat1 [nondim].
real(wp) function NS_coast(lon, lat, lon0, lat0, lat1, dlon, sh)
  real(wp), intent(in) :: lon     !< Longitude [degrees_E]
  real(wp), intent(in) :: lat     !< Latitude [degrees_N]
  real(wp), intent(in) :: lon0    !< Longitude of coast [degrees_E]
  real(wp), intent(in) :: lat0    !< Latitude of coast end [degrees_N]
  real(wp), intent(in) :: lat1    !< Latitude of coast end [degrees_N]
  real(wp), intent(in) :: dlon    !< "Radius" of coast profile [degrees]
  real(wp), intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real(wp) :: r  ! A relative position [nondim]

  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_coast = cstprof(r, 0._wp, dlon, 0.125_wp, 0.125_wp, 0.5_wp, sh)
end function NS_coast

!> A "coast profile" applied in an E-W line from lon0,lat0 to lon1,lat0 [nondim].
real(wp) function EW_coast(lon, lat, lon0, lon1, lat0, dlat, sh)
  real(wp), intent(in) :: lon     !< Longitude [degrees_E]
  real(wp), intent(in) :: lat     !< Latitude [degrees_N]
  real(wp), intent(in) :: lon0    !< Longitude of coast end [degrees_E]
  real(wp), intent(in) :: lon1    !< Longitude of coast end [degrees_E]
  real(wp), intent(in) :: lat0    !< Latitude of coast [degrees_N]
  real(wp), intent(in) :: dlat    !< "Radius" of coast profile [degrees]
  real(wp), intent(in) :: sh      !< depth of shelf as fraction of full depth [nondim]
  real(wp) :: r  ! A relative position [nondim]

  r = dist_line_fixed_y( lon, lat, lon0, lon1, lat0 )
  EW_coast = cstprof(r, 0._wp, dlat, 0.125_wp, 0.125_wp, 0.5_wp, sh)
end function EW_coast

!> A NS ridge [nondim]
real(wp) function NS_ridge(lon, lat, lon0, lat0, lat1, dlon, rh)
  real(wp), intent(in) :: lon     !< Longitude [degrees_E]
  real(wp), intent(in) :: lat     !< Latitude [degrees_N]
  real(wp), intent(in) :: lon0    !< Longitude of ridge center [degrees_E]
  real(wp), intent(in) :: lat0    !< Latitude of ridge end [degrees_N]
  real(wp), intent(in) :: lat1    !< Latitude of ridge end [degrees_N]
  real(wp), intent(in) :: dlon    !< "Radius" of ridge profile [degrees]
  real(wp), intent(in) :: rh      !< depth of ridge as fraction of full depth [nondim]
  real(wp) :: r ! A distance from a point [degrees]

  r = dist_line_fixed_x( lon, lat, lon0, lat0, lat1 )
  NS_ridge = 1._wp - rh * cone(r, 0._wp, dlon)
end function NS_ridge


!> A circular ridge [nondim]
real(wp) function circ_ridge(lon, lat, lon0, lat0, ring_radius, ring_thickness, ridge_height)
  real(wp), intent(in) :: lon            !< Longitude [degrees_E]
  real(wp), intent(in) :: lat            !< Latitude [degrees_N]
  real(wp), intent(in) :: lon0           !< Longitude of center of ring [degrees_E]
  real(wp), intent(in) :: lat0           !< Latitude of center of ring [degrees_N]
  real(wp), intent(in) :: ring_radius    !< Radius of ring [degrees]
  real(wp), intent(in) :: ring_thickness !< Radial thickness of ring [degrees]
  real(wp), intent(in) :: ridge_height   !< Ridge height as fraction of full depth [nondim]
  real(wp) :: r ! A relative position [degrees]
  real(wp) :: frac_ht ! The fractional height of the topography [nondim]

  r = sqrt( ((lon - lon0)**2) + ((lat - lat0)**2) ) ! Pseudo-distance from a point
  r = abs( r - ring_radius) ! Pseudo-distance from a circle
  frac_ht = cone(r, 0._wp, ring_thickness, ridge_height) ! 0 .. frac_ridge_height
  circ_ridge = 1._wp - frac_ht ! Fractional depths (1-frac_ridge_height) .. 1
end function circ_ridge

!> This subroutine initializes layer thicknesses for the Neverworld test case,
!! by finding the depths of interfaces in a specified latitude-dependent
!! temperature profile with an exponentially decaying thermocline on top of a
!! linear stratification.
subroutine Neverworld_initialize_thickness(h, depth_tot, G, GV, US, param_file, P_ref)
  type(ocean_grid_type),   intent(in) :: G                    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in) :: GV                   !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in) :: US                   !< A dimensional unit scaling type
  real(wp), dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(out) :: h !< The thickness that is being
                                                              !! initialized [Z ~> m]
  real(wp), dimension(SZI_(G),SZJ_(G)), &
                           intent(in) :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in) :: param_file           !< A structure indicating the open
                                                              !! file to parse for model
                                                              !! parameter values.
  real(wp),                    intent(in) :: P_Ref                !< The coordinate-density
                                                              !! reference pressure [R L2 T-2 ~> Pa].
  ! Local variables
  real(wp) :: e0(SZK_(GV)+1)    ! The resting interface heights, in depth units [Z ~> m],
                            ! usually negative because it is positive upward.
  real(wp), dimension(SZK_(GV)) :: h_profile ! Vector of initial thickness profile [Z ~> m].
  real(wp) :: e_interface ! Current interface position [Z ~> m].
  real(wp) :: x, y    ! horizontal coordinates for computation of the initial perturbation normalized
                  ! by the domain sizes [nondim]
  real(wp) :: r1, r2  ! radial coordinates for computation of initial perturbation, normalized
                  ! by the domain sizes [nondim]
  real(wp) :: pert_amp ! Amplitude of perturbations as a fraction of layer thicknesses [nondim]
  real(wp) :: h_noise ! Amplitude of noise to scale h by [nondim]
  real(wp) :: noise   ! Fractional noise in the layer thicknesses [nondim]
  type(randomNumberStream) :: rns ! Random numbers for stochastic tidal parameterization
  character(len=40)  :: mdl = "Neverworld_initialize_thickness" ! This subroutine's name.
  integer :: i, j, k, is, ie, js, je, nz

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  call MOM_mesg("  Neverworld_initialization.F90, Neverworld_initialize_thickness: setting thickness", 5)
  call get_param(param_file, mdl, "INIT_THICKNESS_PROFILE", h_profile, &
                 "Profile of initial layer thicknesses.", units="m", scale=US%m_to_Z, &
                 fail_if_missing=.true.)
  call get_param(param_file, mdl, "NL_THICKNESS_PERT_AMP", pert_amp, &
                 "Amplitude of finite scale perturbations as fraction of depth.", &
                 units="nondim", default=0._wp)
  call get_param(param_file, mdl, "NL_THICKNESS_NOISE_AMP", h_noise, &
                 "Amplitude of noise to scale layer by.", units="nondim", default=0._wp)

  ! e0 is the notional position of interfaces
  e0(1) = 0._wp ! The surface
  do k=1,nz
    e0(k+1) = e0(k) - h_profile(k)
  enddo

  do j=js,je ; do i=is,ie
    e_interface = -depth_tot(i,j)
    do k=nz,2,-1
      h(i,j,k) = e0(k) - e_interface ! Nominal thickness
      x = (G%geoLonT(i,j)-G%west_lon)/G%len_lon
      y = (G%geoLatT(i,j)-G%south_lat)/G%len_lat
      r1 = sqrt(((x-0.7_wp)**2) + ((y-0.2_wp)**2))
      r2 = sqrt(((x-0.3_wp)**2) + ((y-0.25_wp)**2))
      h(i,j,k) = h(i,j,k) + pert_amp * (e0(k) - e0(nz+1)) * &
                            (spike(r1,0.15_wp)-spike(r2,0.15_wp)) ! Prescribed perturbation
      if (h_noise /= 0._wp) then
        rns = initializeRandomNumberStream( int( 4096*(x + (y+1._wp)) ) )
        call getRandomNumbers(rns, noise) ! x will be in (0,1)
        noise = h_noise * 2._wp * ( noise - 0.5_wp ) ! range -h_noise to h_noise
        h(i,j,k) = ( 1._wp + noise ) * h(i,j,k)
      endif
      h(i,j,k) = max( GV%Angstrom_Z, h(i,j,k) ) ! Limit to non-negative
      e_interface = e_interface + h(i,j,k) ! Actual position of upper interface
    enddo
    h(i,j,1) = e0(1) - e_interface ! Nominal thickness
    h(i,j,1) = max( GV%Angstrom_Z, h(i,j,1) ) ! Limit to non-negative
  enddo ; enddo

end subroutine Neverworld_initialize_thickness

end module Neverworld_initialization
