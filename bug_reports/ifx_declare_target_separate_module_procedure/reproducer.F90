! ifx error reproducer:
!   `!$omp declare target` in both the interface of a separate module procedure
!   (in the parent module) and its body (in a submodule).
!
! Compiler: ifx (IFX) 2026.1.0 20260617 (Aurora)
!
! Fails at compile time:
!   ifx -qopenmp -fopenmp-targets=spir64 reproducer.F90
!   -> error #8682: Multiple !$OMP DECLARE TARGET declarations.
!
! Compiles with -DWORKAROUND (directive only in the interface).  Also compiles with
! the directive only in the body.
!
! OpenMP's Fortran restrictions for declare target say a directive in an interface
! block must match one in the definition of the procedure, so the failing form
! should be the conforming one.
module parent
  implicit none ; private
  public :: f

  interface
    !> Doubles x
    elemental module function f(x) result(a)
      real, intent(in) :: x !< Input [arbitrary]
      real :: a             !< 2*x [arbitrary]
      !$omp declare target
    end function f
  end interface
end module parent

submodule (parent) child
  implicit none
contains
  module procedure f
#ifndef WORKAROUND
    !$omp declare target    ! <-- error #8682
#endif
    a = 2.0 * x
  end procedure f
end submodule child

program reproducer
  use parent, only : f
  implicit none
  real :: v(8)
  integer :: i

  v(:) = 1.0
  !$omp target teams distribute parallel do map(tofrom: v)
  do i = 1, 8
    v(i) = f(v(i))
  enddo

  if (all(v == 2.0)) then
    print '(a)', 'PASS: f evaluated correctly in the target region'
  else
    print '(a)', 'FAIL: wrong values'
    error stop 1
  endif
end program reproducer
