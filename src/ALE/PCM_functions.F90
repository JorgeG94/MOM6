! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Piecewise constant reconstruction functions
module PCM_functions

implicit none ; private

public PCM_reconstruction, PCM_reconstruction_columns

contains

!> Reconstruction by constant polynomials within each cell. There is nothing to
!! do but this routine is provided to ensure a homogeneous interface
!! throughout the regridding toolbox.
!!
!! It is assumed that the dimension of 'u' is equal to the number of cells
!! defining 'grid' and 'ppoly'. No consistency check is performed.
pure subroutine PCM_reconstruction( N, u, edge_values, ppoly_coef )
  !$omp declare target
  integer,              intent(in)    :: N !< Number of cells
  real, dimension(N),   intent(in)    :: u !< cell averages in arbitrary units [A]
  real, dimension(N,2), intent(inout) :: edge_values !< Edge value of polynomial,
                                           !! with the same units as u [A].
  real, dimension(N,1), intent(inout) :: ppoly_coef !< Coefficients of polynomial,
                                           !! with the same units as u [A].

  ! Local variables
  integer :: k

  ! The coefficients of the piecewise constant polynomial are simply
  ! the cell averages.
  ppoly_coef(:,1) = u(:)

  ! The edge values are equal to the cell average
  do k = 1,N
    edge_values(k,:) = u(k)
  enddo

end subroutine PCM_reconstruction

!> Reconstruction by constant polynomials in a block of columns, each with N cells, using
!! PCM_reconstruction.  The input and output arrays must be current on the device.
subroutine PCM_reconstruction_columns( N, nib, njb, iie, jje, active, u, edge_values, ppoly_coef )
  integer, intent(in) :: N   !< Number of cells in each column
  integer, intent(in) :: nib !< The i-size of the block arrays
  integer, intent(in) :: njb !< The j-size of the block arrays
  integer, intent(in) :: iie !< The number of valid columns in the i-direction
  integer, intent(in) :: jje !< The number of valid columns in the j-direction
  logical, dimension(nib,njb),  intent(in)    :: active !< True for the columns to reconstruct
  real, dimension(N,nib,njb),   intent(in)    :: u !< cell averages in arbitrary units [A]
  real, dimension(N,2,nib,njb), intent(inout) :: edge_values !< Edge value of polynomial,
                                           !! with the same units as u [A].
  real, dimension(N,1,nib,njb), intent(inout) :: ppoly_coef !< Coefficients of polynomial,
                                           !! with the same units as u [A].
  ! Local variables
  integer :: ii, jj

  do concurrent (jj=1:jje, ii=1:iie, active(ii,jj))
    call PCM_reconstruction( N, u(:,ii,jj), edge_values(:,:,ii,jj), ppoly_coef(:,:,ii,jj) )
  enddo

end subroutine PCM_reconstruction_columns

!> \namespace PCM_functions
!!
!! Date of creation: 2008.06.06
!! L. White
!!
!! This module contains routines that handle one-dimensional finite volume
!! reconstruction using the piecewise constant method (PCM).

end module PCM_functions
