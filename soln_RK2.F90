subroutine soln_RK2(dt, m, Vm, Flux)
  !performs the jth step of the RK3 algorithm
  !Total Flux = (k1 + k2 + 2k3)/6

#include "definition.h"

  use grid_data
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: m
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)), intent(INOUT) :: Vm
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM), intent(INOUT) :: Flux !(var,i,j,ndim)

  real :: dtx,dty,dtz
  integer :: i,j,k
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) :: Uk
  call soln_reconstruct(dt, Vm)
  call soln_getFlux(Vm)
  dtx = dt/gr_dx
  dty = dt/gr_dy
  dtz = dt/gr_dz
  Uk = 0.
  if (m == 1) then
     do i = gr_ibeg(XDIM), gr_iend(XDIM)
        do j = gr_ibeg(YDIM), gr_iend(YDIM)
           do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
              Uk(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k)                                - &
                   dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
                   dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
                   dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))
              !print *, "x", gr_flux(:,i,j+1,2) - gr_flux(:,i,j,2), j
              call cons2prim(Uk(DENS_VAR:ENER_VAR,i,j,k),Vm(DENS_VAR:GAME_VAR,i,j,k))
           end do
        end do
     end do
     
  end if

  
  Flux(:,:,:,:,:) = Flux(:,:,:,:,:) + 0.5*gr_flux(:,:,:,:,:)
end subroutine soln_RK2
