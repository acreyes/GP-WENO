subroutine soln_RK3(dt, m, Vm, Flux)
  !performs the mth step of the RK3 algorithm
  !Total Flux = (k1 + k2 + 4k3)/6

#include "definition.h"
  
  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: m
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)), intent(INOUT) :: Vm
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM), intent(INOUT) :: Flux

  real :: F, dtx, dty, dtz
  integer :: i,j,k
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) :: Um

  !F is the factor that multiplies the Km flux
  if (m == 3) then
     F = 2./3.
  else
     F = 1./6.
  end if
  
  call soln_reconstruct(dt, Vm)
  call soln_getFlux(Vm)

  dtx = dt/gr_dx
  dty = dt/gr_dy
  dtz = dt/gr_dz

  if (m == 1) then
     do i = gr_ibeg(XDIM), gr_iend(XDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
           do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
              !getting the U1 state
              !U1 = U0 + k1
              Um(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - &
                   dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
                   dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
                   dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))
              call cons2prim(Um(DENS_VAR:ENER_VAR,i,j,k),Vm(DENS_VAR:GAME_VAR,i,j,k))
           end do
        end do
     end do
  elseif (m == 2) then
     do i = gr_ibeg(XDIM), gr_iend(XDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
           do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
              !getting the U2 state
              !U2 = U0 + 1/4(k1 +k2)
              !k1 = 6.*Flux
              Um(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - 0.25*( &
                   dtx*( &
                   6.*(Flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) -   Flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) + &
                   (gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) -gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) &
                   ) + dty*( &
                   6.*(Flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) -   Flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) +&
                   (gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) -gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM))  &
                   ) + dtz*(&
                   6.*(Flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) -   Flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM)) + &
                   (gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) -gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))&
                   ))
              call cons2prim(Um(DENS_VAR:ENER_VAR,i,j,k),Vm(DENS_VAR:GAME_VAR,i,j,k))
           end do
        end do
     end do
  end if

  !get updated prim vars if needed for m+1 step
  !really this should be done in the above loops?
!!$  if (m .NE. 3) then
!!$     do i = gr_ibeg(XDIM), gr_iend(XDIM)
!!$        do j = gr_ibeg(YDIM), gr_iend(YDIM)
!!$           ! Eos is automatically callled inside cons2prim
!!$           call cons2prim(Um(DENS_VAR:ENER_VAR,i,j),Vm(DENS_VAR:GAME_VAR,i,j))
!!$        end do
!!$     end do
!!$  end if

  !finally add km flux to total flux with weights
  Flux(:,:,:,:,:) = Flux(:,:,:,:,:) + F*gr_flux(:,:,:,:,:)

end subroutine soln_RK3
