subroutine soln_ReconEvolveAvg(dt)

  
#include "definition.h"  

  use grid_data
  use sim_data
  use primconsflux
  use bc

  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) :: Vj
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM) :: Flux !(var,i,j,k,ndim)
  integer :: i,j,k

  
  if (sim_Torder == 1) then
     !this is just forward euler
     !should also be called w/ sim_RK set to true
     call soln_reconstruct(dt, gr_V)
     call soln_getFlux(gr_V)
  elseif (sim_Torder == 2) then
     !second order in time
     !do RK2
     Flux = 0.
     !initialize Vj as V0
     Vj(:,:,:,:) = gr_V(:,:,:,:)
     !do RK2 steps
     do j = 1,2
        call soln_RK2(dt, j, Vj, Flux)
        call bc_apply(Vj)
     end do
     gr_flux(DENS_VAR:ENER_VAR,:,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:,:)
     
  elseif (sim_Torder == 3) then
     !RK3 in time
     Flux = 0.
     !initialize Vj as V0

     Vj(:,:,:,:) = gr_V(:,:,:,:)
     !do RK3 steps

     do j = 1,3
        call soln_RK3(dt, j, Vj, Flux)
        select case (j)
        case(1)
           sim_bcT = sim_time + dt
        case(2)
           sim_bcT = sim_time + 0.5*dt
        case(3)
           sim_bcT = sim_time + dt
        end select
        call bc_apply(Vj)
     end do

     gr_flux(DENS_VAR:ENER_VAR,:,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:,:)

!!$     
  elseif (sim_Torder == 4) then
     !RK4 in time
     Flux = 0.
     !initialize Vj as V0
     Vj(:,:,:,:) = gr_V(:,:,:,:)
     !do rk4 steps
     do j = 1, 4
        !each step adds the flux at the kj'th state with the appropriate weight
        !to the total flux
        call soln_RK4(dt, j, Vj, Flux)
        call bc_apply(Vj)
     end do
     gr_flux(DENS_VAR:ENER_VAR,:,:,:,:) = Flux(DENS_VAR:ENER_VAR,:,:,:,:)
     
  end if
     

  return
end subroutine soln_ReconEvolveAvg
