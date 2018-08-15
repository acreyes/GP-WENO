subroutine sim_GPFinit()

#include "definition.h"

  use grid_data, only: gr_dx
  use sim_data, only: sim_intFlux
  use gp_data
  use linalg
  use GP

  implicit none

  real, dimension(2*gp_radius+1,2*gp_radius+1) :: C, L
  real, dimension(2*gp_radius+1) :: T, Z
  real, dimension(2*gp_radius+1)   :: stencil, u

  integer :: R, N, i, j
  logical :: switch
  
  R = gp_radius
  N = 2*R+1
  C = 0.
  T = 0.
  u = 1.

  switch = .true. !true for the Dongwook way

  !make the stencil
  !stencil is centered on x_iph
  if (sim_intFlux) then
     !samples are equally spaced
     do i = -R, R, 1
        stencil(i+R+1) = REAL(i)
     end do
  else
     !first sample is at interface
     !rest are at cell centers and equally spaced
     stencil(1) = 0.5
     do i = 1, 2*R
        stencil(i+1) = REAL(i-R+1)
     end do
  end if

  if (switch) then
     !we will do this by reconstructing the primitive function of the fluxes as the numerical flux
     !this is like regular FDM but using a Riemann solver instead of flux splitting for upwinding
     do i = 1, N
        do j = 1, N
           C(i,j) = intg_kernel(stencil(i), stencil(j), gp_eldel)
        end do
        T(i) = int_egrand(stencil(i), 0., gp_eldel)
     end do

  else
     !Adam's way (little harder so do it later)
  end if

  call chol(C,N,L)
  call solve_Axb(C, Z, T, L, N)
  allocate(gp_Z(2*R+1)); gp_Z = Z

end subroutine sim_GPFinit
