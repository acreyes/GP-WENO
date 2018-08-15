subroutine hll(vL,vR,Flux, dir)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons


  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux
  integer, intent(IN)                   :: dir

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR
  real :: sL,sR,aL,aR
  integer :: VELC_VAR
  VELC_VAR = 1
  if (dir == XDIM) then
     VELC_VAR = VELX_VAR
  elseif (dir == YDIM) then
     VELC_VAR = VELY_VAR
  elseif (dir == ZDIM) then
     VELC_VAR = VELZ_VAR
  end if

  call prim2flux(vL,FL,dir)
  call prim2flux(vR,FR,dir)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)

  
  ! left and right sound speed a
  aL = sqrt(vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR))
  aR = sqrt(vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR))

  ! fastest left and right going velocities
  sL = min(vL(VELC_VAR) - aL,vR(VELC_VAR) - aR)
  sR = max(vL(VELC_VAR) + aL,vR(VELC_VAR) + aR)

  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)
  elseif ( (sL < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = (    sR*FL(DENS_VAR:ENER_VAR) &
                                   -sL*FR(DENS_VAR:ENER_VAR) &
                               +sR*sL*(uR(DENS_VAR:ENER_VAR) &
                                      -uL(DENS_VAR:ENER_VAR)))/(sR-sL)
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif

  return
end subroutine hll
