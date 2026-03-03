!> Provides a few physical constants
module MOM_constants

! This file is part of MOM6. See LICENSE.md for the license.

use constants_mod, only : FMS_HLV => HLV
use constants_mod, only : FMS_HLF => HLF

use MOM_datatypes, only : wp

implicit none ; private

real(wp), public, parameter :: CELSIUS_KELVIN_OFFSET = 273.15_wp
  !< The constant offset for converting temperatures in Kelvin to Celsius [K]
real(wp), public, parameter :: HLV = real(FMS_HLV, kind=kind(1.0_wp))
  !< Latent heat of vaporization [J kg-1]
real(wp), public, parameter :: HLF = real(FMS_HLF, kind=kind(1.0_wp))
  !< Latent heat of fusion [J kg-1]

end module MOM_constants
