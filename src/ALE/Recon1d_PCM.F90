!> 1D reconstructions using the Piecewise Constant Method (PCM)
module Recon1d_PCM

! This file is part of MOM6. See LICENSE.md for the license.

use Recon1d_type, only : Recon1d, testing

use MOM_datatypes, only : wp

implicit none ; private

public PCM

!> PCM (piecewise constant) reconstruction
!!
!! The source for the methods ultimately used by this class are:
!!   init()                    *locally defined
!!   reconstruct()             *locally defined
!!   average()                 *locally defined
!!   f()                       *locally defined
!!   dfdx()                    *locally defined
!!   check_reconstruction()    *locally defined
!!   unit_tests()              *locally defined
!!   destroy()                 *locally defined
!!   remap_to_sub_grid()    -> Recon1d%remap_to_sub_grid()
!!   init_parent()          -> init()
!!   reconstruct_parent()   -> parent()
type, extends (Recon1d) :: PCM

contains
  !> Implementation of the PCM initialization
  procedure :: init => init
  !> Implementation of the PCM reconstruction
  procedure :: reconstruct => reconstruct
  !> Implementation of the PCM average over an interval [A]
  procedure :: average => average
  !> Implementation of evaluating the PCM reconstruction at a point [A]
  procedure :: f => f
  !> Implementation of the derivative of the PCM reconstruction at a point [A]
  procedure :: dfdx => dfdx
  !> Implementation of deallocation for PCM
  procedure :: destroy => destroy
  !> Implementation of check reconstruction for the PCM reconstruction
  procedure :: check_reconstruction => check_reconstruction
  !> Implementation of unit tests for the PCM reconstruction
  procedure :: unit_tests => unit_tests

  !> Duplicate interface to init()
  procedure :: init_parent => init
  !> Duplicate interface to reconstruct()
  procedure :: reconstruct_parent => reconstruct

end type PCM

contains

!> Initialize a 1D PCM reconstruction for n cells
subroutine init(this, n, h_neglect, check)
  class(PCM),        intent(out) :: this      !< This reconstruction
  integer,           intent(in)  :: n         !< Number of cells in this column
  real(wp), optional,    intent(in)  :: h_neglect !< A negligibly small width used in cell reconstructions [H].
                                              !! Not used by PCM.
  logical, optional, intent(in)  :: check     !< If true, enable some consistency checking

  if (present(h_neglect)) this%n = n ! no-op to avoid compiler warning about unused dummy argument
  if (present(check)) this%check = check

  this%n = n

  allocate( this%u_mean(n) )

end subroutine init

!> Calculate a 1D PCM reconstructions based on h(:) and u(:)
subroutine reconstruct(this, h, u)
  class(PCM), intent(inout) :: this !< This reconstruction
  real(wp),       intent(in)    :: h(*) !< Grid spacing (thickness) [typically H]
  real(wp),       intent(in)    :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  this%u_mean(1) = h(1) ! no-op to avoid compiler warning about unused dummy argument

  do k = 1, this%n
    this%u_mean(k) = u(k)
  enddo

end subroutine reconstruct

!> Value of PCM reconstruction at a point in cell k [A]
real(wp) function f(this, k, x)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real(wp),       intent(in) :: x    !< Non-dimensional position within element [nondim]

  f = this%u_mean(k)

end function f

!> Derivative of PCM reconstruction at a point in cell k [A]
real(wp) function dfdx(this, k, x)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real(wp),       intent(in) :: x    !< Non-dimensional position within element [nondim]

  dfdx = 0._wp

end function dfdx

!> Average between xa and xb for cell k of a 1D PCM reconstruction [A]
real(wp) function average(this, k, xa, xb)
  class(PCM), intent(in) :: this !< This reconstruction
  integer,    intent(in) :: k    !< Cell number
  real(wp),       intent(in) :: xa   !< Start of averaging interval on element (0 to 1)
  real(wp),       intent(in) :: xb   !< End of averaging interval on element (0 to 1)

  average = xb + xa ! no-op to avoid compiler warnings about unused dummy argument
  average = this%u_mean(k)

end function average

!> Deallocate the PCM reconstruction
subroutine destroy(this)
  class(PCM), intent(inout) :: this !< This reconstruction

  deallocate( this%u_mean )

end subroutine destroy

!> Checks the PCM reconstruction for consistency
logical function check_reconstruction(this, h, u)
  class(PCM), intent(in) :: this !< This reconstruction
  real(wp),       intent(in) :: h(*) !< Grid spacing (thickness) [typically H]
  real(wp),       intent(in) :: u(*) !< Cell mean values [A]
  ! Local variables
  integer :: k

  check_reconstruction = .false.

  do k = 1, this%n
    if ( abs( this%u_mean(k) - u(k) ) > 0._wp ) check_reconstruction = .true.
  enddo

end function check_reconstruction

!> Runs PCM reconstruction unit tests and returns True for any fails, False otherwise
logical function unit_tests(this, verbose, stdout, stderr)
  class(PCM), intent(inout) :: this    !< This reconstruction
  logical,    intent(in)    :: verbose !< True, if verbose
  integer,    intent(in)    :: stdout  !< I/O channel for stdout
  integer,    intent(in)    :: stderr  !< I/O channel for stderr
  ! Local variables
  real(wp), allocatable :: ul(:), ur(:), um(:) ! test values [A]
  type(testing) :: test ! convenience functions
  integer :: k

  call test%set( stdout=stdout ) ! Sets the stdout channel in test
  call test%set( stderr=stderr ) ! Sets the stderr channel in test
  call test%set( verbose=verbose ) ! Sets the verbosity flag in test

  call this%init(3)
  call test%test( this%n /= 3, 'Setting number of levels')
  allocate( um(3), ul(3), ur(3) )

  call this%reconstruct( (/2._wp,2._wp,2._wp/), (/1._wp,3._wp,5._wp/) )
  call test%real_arr(3, this%u_mean, (/1._wp,3._wp,5._wp/), 'Setting cell values')

  do k = 1, 3
    ul(k) = this%f(k, 0._wp)
    um(k) = this%f(k, 0.5_wp)
    ur(k) = this%f(k, 1._wp)
  enddo
  call test%real_arr(3, ul, (/1._wp,3._wp,5._wp/), 'Evaluation on left edge')
  call test%real_arr(3, um, (/1._wp,3._wp,5._wp/), 'Evaluation in center')
  call test%real_arr(3, ur, (/1._wp,3._wp,5._wp/), 'Evaluation on right edge')

  do k = 1, 3
    ul(k) = this%dfdx(k, 0._wp)
    um(k) = this%dfdx(k, 0.5_wp)
    ur(k) = this%dfdx(k, 1._wp)
  enddo
  call test%real_arr(3, ul, (/0._wp,0._wp,0._wp/), 'dfdx on left edge')
  call test%real_arr(3, um, (/0._wp,0._wp,0._wp/), 'dfdx in center')
  call test%real_arr(3, ur, (/0._wp,0._wp,0._wp/), 'dfdx on right edge')

  do k = 1, 3
    um(k) = this%average(k, 0.5_wp, 0.75_wp)
  enddo
  call test%real_arr(3, um, (/1._wp,3._wp,5._wp/), 'Return interval average')

  unit_tests = test%summarize('PCM:unit_tests')

end function unit_tests

!> \namespace recon1d_pcm
!!

end module Recon1d_PCM
