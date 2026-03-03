module MOM_MEKE_types

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_datatypes, only : wp
implicit none ; private

!> This type is used to exchange information related to the MEKE calculations.
type, public :: MEKE_type
  ! Variables
  real(wp), allocatable :: MEKE(:,:)    !< Vertically averaged eddy kinetic energy [L2 T-2 ~> m2 s-2].
  real(wp), allocatable :: GM_src(:,:)  !< MEKE source due to thickness mixing (GM) [R Z L2 T-3 ~> W m-2].
  real(wp), allocatable :: mom_src(:,:) !< MEKE source from lateral friction in the
                                    !! momentum equations [R Z L2 T-3 ~> W m-2].
  real(wp), allocatable :: mom_src_bh(:,:) !< MEKE source from the biharmonic part of the lateral friction in the
                                    !! momentum equations [R Z L2 T-3 ~> W m-2].
  real(wp), allocatable :: GME_snk(:,:) !< MEKE sink from GME backscatter in the momentum equations [R Z L2 T-3 ~> W m-2].
  real(wp), allocatable :: Kh(:,:)      !< The MEKE-derived lateral mixing coefficient [L2 T-1 ~> m2 s-1].
  real(wp), allocatable :: Kh_diff(:,:) !< Uses the non-MEKE-derived thickness diffusion coefficient to diffuse
                                    !! MEKE [L2 T-1 ~> m2 s-1].
  real(wp), allocatable :: Rd_dx_h(:,:) !< The deformation radius compared with the grid spacing [nondim].
                                    !! Rd_dx_h is copied from VarMix_CS.
  real(wp), allocatable :: Ku(:,:)      !< The MEKE-derived lateral viscosity coefficient
                                    !! [L2 T-1 ~> m2 s-1]. This viscosity can be negative when representing
                                    !! backscatter from unresolved eddies (see Jansen and Held, 2014).
  real(wp), allocatable :: Au(:,:)      !< The MEKE-derived lateral biharmonic viscosity
                                    !! coefficient [L4 T-1 ~> m4 s-1].
  real(wp), allocatable :: Le(:,:)      !< Eddy length scale [L ~> m]

  ! Parameters
  real(wp) :: KhTh_fac = 1.0_wp !< Multiplier to map Kh(MEKE) to KhTh [nondim]
  real(wp) :: KhTr_fac = 1.0_wp !< Multiplier to map Kh(MEKE) to KhTr [nondim].
  real(wp) :: backscatter_Ro_pow = 0.0_wp !< Power in Rossby number function for backscatter [nondim].
  real(wp) :: backscatter_Ro_c = 0.0_wp !< Coefficient in Rossby number function for backscatter [nondim].

end type MEKE_type

end module MOM_MEKE_types
