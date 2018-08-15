subroutine gp_Fluxinit()

#include "definition.h"

  use grid_data, only: gr_dx
  use sim_data, only: sim_intFlux, sim_DongwookFlux
  use gp_data
  use linalg
  use GP

  implicit none

  real, dimension(3,3) :: C2, L2
  real, dimension(5,5) :: C4, L4
  real, dimension(3)   :: T2, Z2
  real, dimension(5)   :: T4, Z4
  real, dimension(3)   :: stencil2, u2
  real, dimension(5)   :: stencil4, u4

  real :: eldel
  integer :: R, N, i, j, N2, N4
  logical :: switch
  
  R  = 2
  N2 = 3
  N4 = 5
  C2 = 0.; C4 = 0.
  T2 = 0.; T4 = 0.
  u2 = 1.; u4 = 0.
  eldel = gp_eldel


  !make the stencil
  !stencil is centered on x_iph
  if (sim_intFlux) then
     !samples are equally spaced
     do i = -R, R, 1
        stencil4(i+R+1) = REAL(i)
     end do
  else
     !first sample is at interface
     !rest are at cell centers and equally spaced
     stencil4(1) = 0.5
     do i = 1, 2*R
        stencil4(i+1) = REAL(i-R+1)
     end do
  end if
  stencil2 = stencil4(2:4)

  if (sim_DongwookFlux) then
     !we will do this by reconstructing the primitive function of the fluxes as the numerical flux
     !this is like regular FDM but using a Riemann solver instead of flux splitting for upwinding
!!$     do i = 1, N
!!$        do j = 1, N
!!$           C(i,j) = intg_kernel(stencil(i), stencil(j), eldel)
!!$        end do
!!$        T(1,i) = int_egrand(stencil(i), 0., eldel)
!!$     end do
!!$
!!$     call chol(C,N,L)
!!$     call solve_Axb(C, Z(1,:), T(1,:), L, N)
!!$     allocate(gpF_Z(1,2*R+1)); gpF_Z(1,:) = Z(1,:)
     !this hasn't been working and I do not understand what needs to change for it
     
  else
     !Adam's way using the flux derivatives to get the numerical flux just like Chen/DelZanna
    do i = 1, N2
        do j = 1, N2
            C2(i,j) = SE(stencil2(i), stencil2(j), eldel)
         end do
        T2(i) = d2_SE(0., stencil2(i), eldel)
     end do
     call chol(C2,N2,L2)
     call solve_Axb(C2, Z2, T2, L2, N2)
     allocate(gpF_Z2(N2)); gpF_Z2(:) = Z2(:)
     do i = 1, N4
        do j = 1, N4
            C4(i,j) = SE(stencil4(i), stencil4(j), eldel)
         end do
        T4(i) = d4_SE(0., stencil4(i), eldel)
     end do
     call chol(C4,N4,L4)
     call solve_Axb(C4, Z4, T4, L4, N4)
     allocate(gpF_Z4(N4)); gpF_Z4(:) = Z4(:)
  end if

     
  
end subroutine gp_Fluxinit
