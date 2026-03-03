!> Piecewise Parabolic Method 1D reconstruction with h4 interpolation for edges
!!
!! This implementation of PPM follows White and Adcroft 2008 \cite white2008, with cells
!! resorting to PCM for extrema including first and last cells in column.
!! This scheme differs from Colella and Woodward, 1984 \cite colella1984, in the method
!! of first estimating the fourth-order accurate edge values.
!! This uses numerical expressions refactored at the beginning of 2019.
!! The first and last cells are always limited to PCM.
module Recon1d_PPM_H4_2019

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

use MOM_datatypes, only : wp

implicit none ; private

public PPM_H4_2019, testing

!> PPM reconstruction following White and Adcroft, 2008
!!
!! The source for the methods ultimately used by this class are:
!! - init()                    *locally defined
!! - reconstruct()             *locally defined
!! - average()                 *locally defined
!! - f()                       *locally defined
!! - dfdx()                    *locally defined
!! - check_reconstruction()    *locally defined
!! - unit_tests()              *locally defined
!! - destroy()                 *locally defined
!! - remap_to_sub_grid()    -> recon1d_type.remap_to_sub_grid()
!! - init_parent()          -> init()
!! - reconstruct_parent()   -> reconstruct()
type, extends (Recon1d) :: PPM_H4_2019

  real(wp), allocatable :: ul(:) !< Left edge value [A]
  real(wp), allocatable :: ur(:) !< Right edge value [A]

contains
  !> Implementation of the PPM_H4_2019 initialization
  procedure :: init => init
  !> Implementation of the PPM_H4_2019 reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PPM_H4_2019 average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PPM_H4_2019 reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PPM_H4_2019 reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PPM_H4_2019
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PPM_H4_2019 reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PPM_H4_2019 reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PPM_H4_2019

contains

!> Initialize a 1D PPM_H4_2019 reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PPM_H4_2019),     intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real(wp), optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H]
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

  this%n = n

  allocate( this%u_mean(n) )
  allocate( this%ul(n) )
  allocate( this%ur(n) )

  this%h_neglect = tiny( this%u_mean(1) )
  if (present(h_neglect)) this%h_neglect = h_neglect
  this%check = .false.
  if (present(check)) this%check = check

end subroutine init

!> Calculate a 1D PPM_H4_2019 reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PPM_H4_2019), intent(inout) :: this !< This reconstruction
  real(wp),          intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real(wp),          intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  real(wp) :: slp ! The PLM slopes (difference across cell) [A]
  real(wp) :: sigma_l, sigma_c, sigma_r ! Left, central and right slope estimates as
                                    ! differences across the cell [A]
  real(wp) :: u_min, u_max ! Minimum and maximum value across cell [A]
  real(wp) :: u_l, u_r, u_c ! Left, right, and center values [A]
  real(wp) :: h_l, h_c, h_r ! Thickness of left, center and right cells [H]
  real(wp) :: h_c0 ! Thickness of center with h_neglect added [H]
  real(wp) :: h0, h1, h2, h3        ! temporary thicknesses [H]
  real(wp) :: h_min                 ! A minimal cell width [H]
  real(wp) :: f1                    ! An auxiliary variable [H]
  real(wp) :: f2                    ! An auxiliary variable [A H]
  real(wp) :: f3                    ! An auxiliary variable [H-1]
  real(wp) :: et1, et2, et3         ! terms the expression for edge values [A H]
  real(wp) :: I_h12                 ! The inverse of the sum of the two central thicknesses [H-1]
  real(wp) :: I_h012, I_h123        ! Inverses of sums of three successive thicknesses [H-1]
  real(wp) :: I_den_et2, I_den_et3  ! Inverses of denominators in edge value terms [H-2]
  real(wp) :: dx                    ! Difference of successive values of x [H]
  real(wp) :: f                     ! value of polynomial at x in arbitrary units [A]
  real(wp) :: edge_l, edge_r        ! Edge values (left and right) [A]
  real(wp) :: expr1, expr2          ! Temporary expressions [A2]
  real(wp) :: slope_x_h             ! retained PLM slope times  half grid step [A]
  real(wp) :: u0_avg                ! avg value at given edge [A]
  real(wp), parameter :: hMinFrac = 1.e-5_wp  !< A minimum fraction for min(h)/sum(h) [nondim]
  real(wp) :: edge_values(this%n,2) ! Edge values [A]
  real(wp) :: ppoly_coef(this%n,3)  ! Polynomial coefficients [A]
  real(wp) :: dz(4)                 ! A temporary array of limited layer thicknesses [H]
  real(wp) :: u_tmp(4)              ! A temporary array of cell average properties [A]
  real(wp) :: A(4,4)     ! Differences in successive positions raised to various powers,
                     ! in units that vary with the second (j) index as [H^j]
  real(wp) :: B(4)       ! The right hand side of the system to solve for C [A H]
  real(wp) :: C(4)       ! The coefficients of a fit polynomial in units that vary
                     ! with the index (j) as [A H^(j-1)]
  integer :: k, n, km1, kp1

  n = this%n

  ! Loop on interior cells
  do K = 3, n-1

    h0 = h(k-2)
    h1 = h(k-1)
    h2 = h(k)
    h3 = h(k+1)

    ! Avoid singularities when consecutive pairs of h vanish
    if (h0+h1==0.0_wp .or. h1+h2==0.0_wp .or. h2+h3==0.0_wp) then
      h_min = hMinFrac*max( this%h_neglect, (h0+h1)+(h2+h3) )
      h0 = max( h_min, h0 )
      h1 = max( h_min, h1 )
      h2 = max( h_min, h2 )
      h3 = max( h_min, h3 )
    endif

    I_h12 = 1.0_wp / (h1+h2)
    I_den_et2 = 1.0_wp / ( ((h0+h1)+h2)*(h0+h1) ) ; I_h012 = (h0+h1) * I_den_et2
    I_den_et3 = 1.0_wp / ( (h1+(h2+h3))*(h2+h3) ) ; I_h123 = (h2+h3) * I_den_et3

    et1 = ( 1.0_wp + (h1 * I_h012 + (h0+h1) * I_h123) ) * I_h12 * (h2*(h2+h3)) * u(k-1) + &
          ( 1.0_wp + (h2 * I_h123 + (h2+h3) * I_h012) ) * I_h12 * (h1*(h0+h1)) * u(k)
    et2 = ( h1 * (h2*(h2+h3)) * I_den_et2 ) * (u(k-1)-u(k-2))
    et3 = ( h2 * (h1*(h0+h1)) * I_den_et3 ) * (u(k) - u(k+1))
    edge_values(k,1) = (et1 + (et2 + et3)) / ((h0 + h1) + (h2 + h3))
    edge_values(k-1,2) = edge_values(k,1)

  enddo ! end loop on interior cells

  ! Determine first two edge values
  do k=1,4 ; dz(k) = max(this%h_neglect, h(k) ) ; u_tmp(k) = u(k) ; enddo
  call end_value_h4(dz, u_tmp, C)

  ! Set the edge values of the first cell
  edge_values(1,1) = C(1)
  edge_values(1,2) = C(1) + dz(1) * ( C(2) + dz(1) * ( C(3) + dz(1) * C(4) ) )
  edge_values(2,1) = edge_values(1,2)

  ! Determine two edge values of the last cell
  do k=1,4 ; dz(k) = max(this%h_neglect, h(n+1-k) ) ; u_tmp(k) = u(n+1-k) ; enddo
  call end_value_h4(dz, u_tmp, C)

  ! Set the last and second to last edge values
  edge_values(n,2) = C(1)
  edge_values(n,1) = C(1) + dz(1) * ( C(2) + dz(1) * ( C(3) + dz(1) * C(4) ) )
  edge_values(n-1,2) = edge_values(n,1)

  ! Loop on cells to bound edge value
  do k = 1, n

    ! For the sake of bounding boundary edge values, the left neighbor of the left boundary cell
    ! is assumed to be the same as the left boundary cell and the right neighbor of the right
    ! boundary cell is assumed to be the same as the right boundary cell. This effectively makes
    ! boundary cells look like extrema.
    km1 = max(1,k-1) ; kp1 = min(k+1,N)

    slope_x_h = 0.0_wp
    sigma_l = ( u(k) - u(km1) )
    if ( (h(km1) + h(kp1)) + 2.0_wp*h(k) > 0._wp ) then
      sigma_c = ( u(kp1) - u(km1) ) * ( h(k) / ((h(km1) + h(kp1)) + 2.0_wp*h(k)) )
    else
      sigma_c = 0._wp
    endif
    sigma_r = ( u(kp1) - u(k) )

    ! The limiter is used in the local coordinate system to each cell, so for convenience store
    ! the slope times a half grid spacing.  (See White and Adcroft JCP 2008 Eqs 19 and 20)
    if ( (sigma_l * sigma_r) > 0.0_wp ) &
      slope_x_h = sign( min(abs(sigma_l),abs(sigma_c),abs(sigma_r)), sigma_c )

    ! Limit the edge values
    if ( (u(km1)-edge_values(k,1)) * (edge_values(k,1)-u(k)) < 0.0_wp ) then
      edge_values(k,1) = u(k) - sign( min( abs(slope_x_h), abs(edge_values(k,1)-u(k)) ), slope_x_h )
    endif

    if ( (u(kp1)-edge_values(k,2)) * (edge_values(k,2)-u(k)) < 0.0_wp ) then
      edge_values(k,2) = u(k) + sign( min( abs(slope_x_h), abs(edge_values(k,2)-u(k)) ), slope_x_h )
    endif

    ! Finally bound by neighboring cell means in case of roundoff
    edge_values(k,1) = max( min( edge_values(k,1), max(u(km1), u(k)) ), min(u(km1), u(k)) )
    edge_values(k,2) = max( min( edge_values(k,2), max(u(kp1), u(k)) ), min(u(kp1), u(k)) )

  enddo ! loop on interior edges

  do k = 1, n-1
    if ( (edge_values(k+1,1) - edge_values(k,2)) * (u(k+1) - u(k)) < 0.0_wp ) then
      u0_avg = 0.5_wp * ( edge_values(k,2) + edge_values(k+1,1) )
      u0_avg = max( min( u0_avg, max(u(k), u(k+1)) ), min(u(k), u(k+1)) )
      edge_values(k,2) = u0_avg
      edge_values(k+1,1) = u0_avg
    endif
  enddo ! end loop on interior edges

  ! Loop on interior cells to apply the standard
  ! PPM limiter (Colella & Woodward, JCP 84)
  do k = 2,N-1

    ! Get cell averages
    u_l = u(k-1)
    u_c = u(k)
    u_r = u(k+1)

    edge_l = edge_values(k,1)
    edge_r = edge_values(k,2)

    if ( (u_r - u_c)*(u_c - u_l) <= 0.0_wp) then
      ! Flatten extremum
      edge_l = u_c
      edge_r = u_c
    else
      expr1 = 3.0_wp * (edge_r - edge_l) * ( (u_c - edge_l) + (u_c - edge_r))
      expr2 = (edge_r - edge_l) * (edge_r - edge_l)
      if ( expr1 > expr2 ) then
        ! Place extremum at right edge of cell by adjusting left edge value
        edge_l = u_c + 2.0_wp * ( u_c - edge_r )
        edge_l = max( min( edge_l, max(u_l, u_c) ), min(u_l, u_c) ) ! In case of round off
      elseif ( expr1 < -expr2 ) then
        ! Place extremum at left edge of cell by adjusting right edge value
        edge_r = u_c + 2.0_wp * ( u_c - edge_l )
        edge_r = max( min( edge_r, max(u_r, u_c) ), min(u_r, u_c) ) ! In case of round off
      endif
    endif
    ! This checks that the difference in edge values is representable
    ! and avoids overshoot problems due to round off.
    !### The 1.e-60 needs to have units of [A], so this dimensionally inconsistent.
    if ( abs( edge_r - edge_l )<max(1.e-60_wp,epsilon(u_c)*abs(u_c)) ) then
      edge_l = u_c
      edge_r = u_c
    endif

    edge_values(k,1) = edge_l
    edge_values(k,2) = edge_r

  enddo ! end loop on interior cells

  ! PCM within boundary cells
  edge_values(1,:) = u(1)
  edge_values(N,:) = u(N)

  ! Store reconstruction
  do k = 1, n
    this%u_mean(k) = u(k)
    this%ul(k) = edge_values(k,1)
    this%ur(k) = edge_values(k,2)
  enddo

end subroutine reconstruct

!> Determine a one-sided 4th order polynomial fit of u to the data points for the purposes of specifying
!! edge values, as described in the appendix of White and Adcroft JCP 2008.
subroutine end_value_h4(dz, u, Csys)
  real(wp), intent(in)  :: dz(4)    !< The thicknesses of 4 layers, starting at the edge [H].
                                !! The values of dz must be positive.
  real(wp), intent(in)  :: u(4)     !< The average properties of 4 layers, starting at the edge [A]
  real(wp), intent(out) :: Csys(4)  !< The four coefficients of a 4th order polynomial fit
                                !! of u as a function of z [A H-(n-1)]

  ! Local variables
  real(wp) :: Wt(3,4)         ! The weights of successive u differences in the 4 closed form expressions.
                          ! The units of Wt vary with the second index as [H-(n-1)].
  real(wp) :: h1, h2, h3, h4  ! Copies of the layer thicknesses [H]
  real(wp) :: h12, h23, h34   ! Sums of two successive thicknesses [H]
  real(wp) :: h123, h234      ! Sums of three successive thicknesses [H]
  real(wp) :: h1234           ! Sums of all four thicknesses [H]
  ! real :: I_h1          ! The inverse of the a thickness [H-1]
  real(wp) :: I_h12, I_h23, I_h34 ! The inverses of sums of two thicknesses [H-1]
  real(wp) :: I_h123, I_h234  ! The inverse of the sum of three thicknesses [H-1]
  real(wp) :: I_h1234         ! The inverse of the sum of all four thicknesses [H-1]
  real(wp) :: I_denom         ! The inverse of the denominator some expressions [H-3]
  real(wp) :: I_denB3         ! The inverse of the product of three sums of thicknesses [H-3]
  real(wp) :: min_frac = 1.0e-6_wp  ! The square of min_frac should be much larger than roundoff [nondim]
  real(wp), parameter :: C1_3 = 1.0_wp / 3.0_wp   ! A rational parameter [nondim]

 ! if ((dz(1) == dz(2)) .and. (dz(1) == dz(3)) .and. (dz(1) == dz(4))) then
 !   ! There are simple closed-form expressions in this case
 !   I_h1 = 0.0 ; if (dz(1) > 0.0) I_h1 = 1.0 / dz(1)
 !   Csys(1) = u(1) + (-13.0 * (u(2)-u(1)) + 10.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25*C1_3)
 !   Csys(2) = (35.0 * (u(2)-u(1)) - 34.0 * (u(3)-u(2)) + 11.0 * (u(4)-u(3))) * (0.25*C1_3 * I_h1)
 !   Csys(3) = (-5.0 * (u(2)-u(1)) + 8.0 * (u(3)-u(2)) - 3.0 * (u(4)-u(3))) * (0.25 * I_h1**2)
 !   Csys(4) = ((u(2)-u(1)) - 2.0 * (u(3)-u(2)) + (u(4)-u(3))) * (0.5*C1_3)
 ! else

  ! Express the coefficients as sums of the differences between properties of successive layers.

  h1 = dz(1) ; h2 = dz(2) ; h3 = dz(3) ; h4 = dz(4)
  ! Some of the weights used below are proportional to (h1/(h2+h3))**2 or (h1/(h2+h3))*(h2/(h3+h4))
  ! so h2 and h3 should be adjusted to ensure that these ratios are not so large that property
  ! differences at the level of roundoff are amplified to be of order 1.
  if ((h2+h3) < min_frac*h1) h3 = min_frac*h1 - h2
  if ((h3+h4) < min_frac*h1) h4 = min_frac*h1 - h3

  h12 = h1+h2 ; h23 = h2+h3 ; h34 = h3+h4
  h123 = h12 + h3 ; h234 = h2 + h34 ; h1234 = h12 + h34
  ! Find 3 reciprocals with a single division for efficiency.
  I_denB3 = 1.0_wp / (h123 * h12 * h23)
  I_h12 = (h123 * h23) * I_denB3
  I_h23 = (h12 * h123) * I_denB3
  I_h123 = (h12 * h23) * I_denB3
  I_denom = 1.0_wp / ( h1234 * (h234 * h34) )
  I_h34 = (h1234 * h234) * I_denom
  I_h234 = (h1234 * h34) * I_denom
  I_h1234 = (h234 * h34) * I_denom

  ! Calculation coefficients in the four equations

  ! The expressions for Csys(3) and Csys(4) come from reducing the 4x4 matrix problem into the following 2x2
  ! matrix problem, then manipulating the analytic solution to avoid any subtraction and simplifying.
  !  (C1_3 * h123 * h23) * Csys(3) + (0.25 * h123 * h23 * (h3 + 2.0*h2 + 3.0*h1)) * Csys(4) =
  !            (u(3)-u(1)) - (u(2)-u(1)) * (h12 + h23) * I_h12
  !  (C1_3 * ((h23 + h34) * h1234 + h23 * h3)) * Csys(3) +
  !  (0.25 * ((h1234 + h123 + h12 + h1) * h23 * h3 + (h1234 + h12 + h1) * (h23 + h34) * h1234)) * Csys(4) =
  !            (u(4)-u(1)) - (u(2)-u(1)) * (h123 + h234) * I_h12
  ! The final expressions for Csys(1) and Csys(2) were derived by algebraically manipulating the following expressions:
  !  Csys(1) = (C1_3 * h1 * h12 * Csys(3) + 0.25 * h1 * h12 * (2.0*h1+h2) * Csys(4)) - &
  !            (h1*I_h12)*(u(2)-u(1)) + u(1)
  !  Csys(2) = (-2.0*C1_3 * (2.0*h1+h2) * Csys(3) - 0.5 * (h1**2 + h12 * (2.0*h1+h2)) * Csys(4)) + &
  !            2.0*I_h12 * (u(2)-u(1))
  ! These expressions are typically evaluated at x=0 and x=h1, so it is important that these are well behaved
  ! for these values, suggesting that h1/h23 and h1/h34 should not be allowed to be too large.

  Wt(1,1) = -h1 * (I_h1234 + I_h123 + I_h12)                              ! > -3
  Wt(2,1) =  h1 * h12 * ( I_h234 * I_h1234 + I_h23 * (I_h234 + I_h123) )  ! < (h1/h234) + (h1/h23)*(2+(h1/h234))
  Wt(3,1) = -h1 * h12 * h123 * I_denom                                    ! > -(h1/h34)*(1+(h1/h234))

  Wt(1,2) =  2.0_wp * (I_h12*(1.0_wp + (h1+h12) * (I_h1234 + I_h123)) + h1 * I_h1234*I_h123) ! < 10/h12
  Wt(2,2) = -2.0_wp * ((h1 * h12 * I_h1234) *       (I_h23 * (I_h234 + I_h123)) + &       ! > -(10+6*(h1/h234))/h23
                    (h1+h12) * ( I_h1234*I_h234 + I_h23 * (I_h234 + I_h123) ) )
  Wt(3,2) =  2.0_wp * ((h1+h12) * h123 + h1*h12 ) * I_denom                               ! < (2+(6*h1/h234)) / h34

  Wt(1,3) = -3.0_wp * I_h12 * I_h123* ( 1.0_wp + I_h1234 * ((h1+h12)+h123) )                 ! > -12 / (h12*h123)
  Wt(2,3) =  3.0_wp * I_h23 * ( I_h123 + I_h1234 * ((h1+h12)+h123) * (I_h123 + I_h234) )  ! < 12 / (h23^2)
  Wt(3,3) = -3.0_wp * ((h1+h12)+h123) * I_denom                                           ! > -9 / (h234*h23)

  Wt(1,4) =  4.0_wp * I_h1234 * I_h123 * I_h12                          ! Wt*h1^3 < 4
  Wt(2,4) = -4.0_wp * I_h1234 * (I_h23 * (I_h123 + I_h234))             ! Wt*h1^3 > -4* (h1/h23)*(1+h1/h234)
  Wt(3,4) =  4.0_wp * I_denom  ! = 4.0*I_h1234 * I_h234 * I_h34         ! Wt*h1^3 < 4 * (h1/h234)*(h1/h34)

  Csys(1) = ((u(1) + Wt(1,1) * (u(2)-u(1))) + Wt(2,1) * (u(3)-u(2))) + Wt(3,1) * (u(4)-u(3))
  Csys(2) = (Wt(1,2) * (u(2)-u(1)) + Wt(2,2) * (u(3)-u(2))) + Wt(3,2) * (u(4)-u(3))
  Csys(3) = (Wt(1,3) * (u(2)-u(1)) + Wt(2,3) * (u(3)-u(2))) + Wt(3,3) * (u(4)-u(3))
  Csys(4) = (Wt(1,4) * (u(2)-u(1)) + Wt(2,4) * (u(3)-u(2))) + Wt(3,4) * (u(4)-u(3))

  ! endif ! End of non-uniform layer thickness branch.

end subroutine end_value_h4

!> Value of PPM_H4_2019 reconstruction at a point in cell k [A]
real(wp) function f(this, k, x)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,            intent(in) :: k    !< Cell number
  real(wp),               intent(in) :: x    !< Non-dimensional position within element [nondim]
  real(wp) :: xc ! Bounded version of x [nondim]
  real(wp) :: du ! Difference across cell [A]
  real(wp) :: a6 ! Collela and Woordward curvature parameter [A]
  real(wp) :: u_a, u_b ! Two estimate of f [A]
  real(wp) :: lmx ! 1 - x [nondim]
  real(wp) :: wb ! Weight based on x [nondim]

  du = this%ur(k) - this%ul(k)
  a6 = 3.0_wp * ( ( this%u_mean(k) - this%ul(k) ) + ( this%u_mean(k) - this%ur(k) ) )
  xc = max( 0._wp, min( 1._wp, x ) )
  lmx = 1.0_wp - xc

  ! This expression for u_a can overshoot u_r but is good for x<<1
  u_a = this%ul(k) + xc * ( du + a6 * lmx )
  ! This expression for u_b can overshoot u_l but is good for 1-x<<1
  u_b = this%ur(k) + lmx * ( - du + a6 * xc )

  ! Since u_a and u_b are both side-bounded, using weights=0 or 1 will preserve uniformity
  wb = 0.5_wp + sign(0.5_wp, xc - 0.5_wp ) ! = 1 @ x=0, = 0 @ x=1
  f = ( ( 1._wp - wb ) * u_a ) + ( wb * u_b )

end function f

!> Derivative of PPM_H4_2019 reconstruction at a point in cell k [A]
real(wp) function dfdx(this, k, x)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,            intent(in) :: k    !< Cell number
  real(wp),               intent(in) :: x    !< Non-dimensional position within element [nondim]
  real(wp) :: xc ! Bounded version of x [nondim]
  real(wp) :: du ! Difference across cell [A]
  real(wp) :: a6 ! Collela and Woordward curvature parameter [A]

  du = this%ur(k) - this%ul(k)
  a6 = 3.0_wp * ( ( this%u_mean(k) - this%ul(k) ) + ( this%u_mean(k) - this%ur(k) ) )
  xc = max( 0._wp, min( 1._wp, x ) )

  dfdx = du + a6 * ( 2.0_wp * xc - 1.0_wp )

end function dfdx

!> Average between xa and xb for cell k of a 1D PPM reconstruction [A]
real(wp) function average(this, k, xa, xb)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  integer,       intent(in) :: k    !< Cell number
  real(wp),          intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real(wp),          intent(in) :: xb   !< End of averaging interval on element (0 to 1)
  real(wp) :: xapxb                 ! A sum of fracional positions [nondim]
  real(wp) :: mx, Ya, Yb, my        ! Various fractional positions [nondim]
  real(wp) :: u_a, u_b ! Values at xa and xb [A]
  real(wp) :: xa2pxb2,  xa2b2ab, Ya2b2ab  ! Sums of squared fractional positions [nondim]
  real(wp) :: a_L, a_R, u_c, a_c    ! Values of the polynomial at various locations [A]

  mx = 0.5_wp * ( xa + xb )
  a_L = this%ul(k)
  a_R = this%ur(k)
  u_c = this%u_mean(k)
  a_c = 0.5_wp * ( ( u_c - a_L ) + ( u_c - a_R ) ) ! a_6 / 6
  if (mx<0.5_wp) then
    ! This integration of the PPM reconstruction is expressed in distances from the left edge
    xa2b2ab = (xa*xa+xb*xb)+xa*xb
    average = a_L + ( ( a_R - a_L ) * mx &
                    + a_c * ( 3._wp * ( xb + xa ) - 2._wp*xa2b2ab ) )
  else
    ! This integration of the PPM reconstruction is expressed in distances from the right edge
    Ya = 1._wp - xa
    Yb = 1._wp - xb
    my = 0.5_wp * ( Ya + Yb )
    Ya2b2ab = (Ya*Ya+Yb*Yb)+Ya*Yb
    average = a_R  + ( ( a_L - a_R ) * my &
                     + a_c * ( 3._wp * ( Yb + Ya ) - 2._wp*Ya2b2ab ) )
  endif

end function average

!> Deallocate the PPM_H4_2019 reconstruction
subroutine destroy(this)
  class(PPM_H4_2019), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean, this%ul, this%ur )

end subroutine destroy

!> Checks the PPM_H4_2019 reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PPM_H4_2019), intent(in) :: this !< This reconstruction
  real(wp),          intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real(wp),          intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  ! Simply checks the internal copy of "u" is exactly equal to "u"
  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0._wp ) check_reconstruction = .true.
  enddo

  ! If (u - ul) has the opposite sign from (ur - u), then this cell has an interior extremum
  do k = 1, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ur(k) - this%u_mean(k) ) < 0._wp ) check_reconstruction = .true.
  enddo

  ! Check bounding of right edges, w.r.t. the cell means
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%u_mean(k+1) - this%ur(k) ) < 0._wp ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. the cell means
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%u_mean(k-1) ) < 0._wp ) check_reconstruction = .true.
  enddo

  ! Check bounding of right edges, w.r.t. this cell mean and the next cell left edge
  do K = 1, this%n-1
    if ( ( this%ur(k) - this%u_mean(k) ) * ( this%ul(k+1) - this%ur(k) ) < 0._wp ) check_reconstruction = .true.
  enddo

  ! Check bounding of left edges, w.r.t. this cell mean and the previous cell right edge
  do K = 2, this%n
    if ( ( this%u_mean(k) - this%ul(k) ) * ( this%ul(k) - this%ur(k-1) ) < 0._wp ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PPM_H4_2019 reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PPM_H4_2019), intent(inout) :: this    !< This reconstruction
  logical,       intent(in)    :: verbose !< True, if verbose
  integer,       intent(in)    :: stdout  !< I/O channel for stdout
  integer,       intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real(wp), allocatable :: ul(:), ur(:), um(:) ! test values [A]
  real(wp), allocatable :: ull(:), urr(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  if (verbose) write(stdout,'(a)') 'PPM_H4_2019:unit_tests testing with linear fn'

  call this%init(5)
  call test%test( this%n /= 5, 'Setting number of levels')
  allocate( um(5), ul(5), ur(5), ull(5), urr(5) )

  ! Straight line, f(x) = x , or  f(K) = 2*K
  call this%reconstruct( (/2._wp,2._wp,2._wp,2._wp,2._wp/), (/1._wp,3._wp,5._wp,7._wp,9._wp/) )
  call test%real_arr(5, this%u_mean, (/1._wp,3._wp,5._wp,7._wp,9._wp/), 'Setting cell values')
  call test%real_arr(5, this%ul, (/1._wp,2._wp,4._wp,6._wp,9._wp/), 'Left edge values', robits=2)
  call test%real_arr(5, this%ur, (/1._wp,4._wp,6._wp,8._wp,9._wp/), 'Right edge values')
  do k = 1, 5
    um(k) = this%u_mean(k)
  enddo
  call test%real_arr(5, um, (/1._wp,3._wp,5._wp,7._wp,9._wp/), 'Return cell mean')

  do k = 1, 5
    ul(k) = this%f(k, 0._wp)
    um(k) = this%f(k, 0.5_wp)
    ur(k) = this%f(k, 1._wp)
  enddo
  call test%real_arr(5, ul, this%ul, 'Evaluation on left edge')
  call test%real_arr(5, um, (/1._wp,3._wp,5._wp,7._wp,9._wp/), 'Evaluation in center')
  call test%real_arr(5, ur, this%ur, 'Evaluation on right edge')

  do k = 1, 5
    ul(k) = this%dfdx(k, 0._wp)
    um(k) = this%dfdx(k, 0.5_wp)
    ur(k) = this%dfdx(k, 1._wp)
  enddo
  call test%real_arr(5, ul, (/0._wp,2._wp,2._wp,2._wp,0._wp/), 'dfdx on left edge', robits=3)
  call test%real_arr(5, um, (/0._wp,2._wp,2._wp,2._wp,0._wp/), 'dfdx in center', robits=2)
  call test%real_arr(5, ur, (/0._wp,2._wp,2._wp,2._wp,0._wp/), 'dfdx on right edge', robits=6)

  do k = 1, 5
    um(k) = this%average(k, 0.5_wp, 0.75_wp) ! Average from x=0.25 to 0.75 in each cell
  enddo
  call test%real_arr(5, um, (/1._wp,3.25_wp,5.25_wp,7.25_wp,9._wp/), 'Return interval average')

  if (verbose) write(stdout,'(a)') 'PPM_H4_2019:unit_tests testing with parabola'

  ! x = 3 i   i=0 at origin
  ! f(x) = x^2 / 3   = 3 i^2
  ! f[i] = [ ( 3 i )^3 - ( 3 i - 3 )^3 ]    i=1,2,3,4,5
  ! means:   1, 7, 19, 37, 61
  ! edges:  0, 3, 12, 27, 48, 75
  call this%reconstruct( (/3._wp,3._wp,3._wp,3._wp,3._wp/), (/1._wp,7._wp,19._wp,37._wp,61._wp/) )
  do k = 1, 5
    ul(k) = this%f(k, 0._wp)
    um(k) = this%f(k, 0.5_wp)
    ur(k) = this%f(k, 1._wp)
  enddo
  call test%real_arr(5, ul, (/1._wp,3._wp,12._wp,27._wp,61._wp/), 'Return left edge', robits=2)
  call test%real_arr(5, ur, (/1._wp,12._wp,27._wp,48._wp,61._wp/), 'Return right edge', robits=1)

  call this%destroy()
  deallocate( um, ul, ur, ull, urr )

  unit_tests = test%summarize('PPM_H4_2019:unit_tests')

end function unit_tests

!> \namespace recon1d_ppm_c4_2019
!!

end module Recon1d_PPM_H4_2019
