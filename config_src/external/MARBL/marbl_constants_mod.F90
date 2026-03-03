!> A non-functioning template of the MARBL constants module
module marbl_constants_mod

  use MOM_datatypes, only : wp
  implicit none
  private

  !> Molecular weight of iron
  real(wp), public, parameter :: molw_Fe = 55.845_wp

end module marbl_constants_mod

