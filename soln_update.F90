subroutine soln_update(dt)

#include "definition.h"
  
  use grid_data
  use primconsflux, only : cons2prim

  implicit none
  real, intent(IN) :: dt
  integer :: i,j,k
  real :: dtx, dty, dtz

  dtx = dt/gr_dx
  dty = dt/gr_dy
  dtz = dt/gr_dz

  !! update conservative vars
  do i = gr_ibeg(XDIM), gr_iend(XDIM)
     do j = gr_ibeg(YDIM), gr_iend(YDIM)
        do k = gr_ibeg(ZDIM), gr_iend(ZDIM)
#ifdef BDRY_VAR
           if (gr_V(BDRY_VAR,i,j,k) == -1.0) then
              gr_U(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - &
                   dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
                   dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
                   dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))
              call cons2prim(gr_U(DENS_VAR:ENER_VAR,i,j,k),gr_V(DENS_VAR:GAME_VAR,i,j,k))
           end if
#else           
           gr_U(DENS_VAR:ENER_VAR,i,j,k) = gr_U(DENS_VAR:ENER_VAR,i,j,k) - &
                dtx*(gr_flux(DENS_VAR:ENER_VAR,i+1,j,k,XDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,XDIM)) - &
                dty*(gr_flux(DENS_VAR:ENER_VAR,i,j+1,k,YDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,YDIM)) - &
                dtz*(gr_flux(DENS_VAR:ENER_VAR,i,j,k+1,ZDIM) - gr_flux(DENS_VAR:ENER_VAR,i,j,k,ZDIM))
           call cons2prim(gr_U(DENS_VAR:ENER_VAR,i,j,k),gr_V(DENS_VAR:GAME_VAR,i,j,k))
#endif           
        end do
     end do
  end do


!!$  !! get updated primitive vars from the updated conservative vars
!!$  do i = gr_ibeg(XDIM), gr_iend(XDIM)
!!$     do j = gr_ibeg(YDIM), gr_iend(YDIM)
!!$        ! Eos is automatically callled inside cons2prim
!!$        call cons2prim(gr_U(DENS_VAR:ENER_VAR,i,j),gr_V(DENS_VAR:GAME_VAR,i,j))
!!$     end do
!!$  end do
  

  return
end subroutine soln_update
