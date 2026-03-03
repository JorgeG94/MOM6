!> Defines portable kind parameters for integer and real types
!!
!! This module provides consistent type definitions across MOM6:
!!
!! - \b int32, \b int64 : 32-bit and 64-bit integers from iso_fortran_env
!! - \b real32, \b real64 : IEEE single and double from iso_fortran_env
!! - \b sp : Single precision real kind (approximately 6 decimal digits)
!! - \b wp : Working precision real kind (approximately 15 decimal digits)
!!
!! Using explicit kind parameters allows the code to be compiled without
!! the -r8 compiler flag while maintaining double precision throughout.
module MOM_datatypes

! This file is part of MOM6. See LICENSE.md for the license.

use, intrinsic :: iso_fortran_env, only : int32, int64, real32, real64

implicit none ; private

! Integer kinds re-exported from iso_fortran_env
public :: int32, int64

! Real kinds re-exported from iso_fortran_env
public :: real32, real64

integer, parameter, public :: sp = selected_real_kind(6, 37)
  !< Single precision real kind (approximately 6 decimal digits, range 1e37)
integer, parameter, public :: wp = selected_real_kind(15, 307)
  !< Working precision real kind (approximately 15 decimal digits, range 1e307)

end module MOM_datatypes
