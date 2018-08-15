subroutine soln_RK4(dt, m, Vm, Flux)
  !performs the mth step of the RK4 algorithm
  !Total Flux = (k1 + 2k2 + 2k3 + k4)/6
#include "definition.h"
  
  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: m
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)), intent(INOUT) :: Vm
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM), intent(INOUT) :: Flux

  real :: A, F, dtx, dty, dtz
  integer :: i, j, k
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) :: Um

  if (m == 1 .OR. m == 2) then
     A = 0.5
  else 
     A = 1.
  end if

  if (m == 1 .OR. m == 4) then
     F = 1./6.
  else
     F = 1./3.
  end if

  dtx = A*dt/gr_dx
  dty = A*dt/gr_dy
  dtz = A*dt/gr_dz
  
  call soln_reconstruct(dt, Vm)
  call soln_getFlux(Vm)

  if (m .NE. 4) then
     !update cons variables to mth step only if not 4th step
     do i = gr_ibeg(XDIM), gr_iend(XDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
           do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
              Um(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - &
                   dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
                   dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
                   dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))

              !convert to primitive vars
              call cons2prim(Um(DENS_VAR:ENER_VAR,i,j,k),Vm(DENS_VAR:GAME_VAR,i,j,k))
           end do
        end do
     end do


     !! get updated primitive vars from the updated conservative vars
!!$     do i = gr_ibeg(XDIM), gr_iend(XDIM)
!!$        do j = gr_ibeg(YDIM), gr_iend(YDIM)
!!$           ! Eos is automatically callled inside cons2prim
!!$           call cons2prim(Um(DENS_VAR:ENER_VAR,i,j),Vm(DENS_VAR:GAME_VAR,i,j))
!!$        end do
!!$     end do
     
  end if

  !finally add to total flux
  Flux(:,:,:,:,:) = Flux(:,:,:,:,:) + F*gr_flux(:,:,:,:,:)
  
end subroutine soln_RK4
