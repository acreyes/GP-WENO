subroutine hllc(vL,vR,Flux,dir)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons,cons2flux


  implicit none
  integer, intent(IN) :: dir
  real, dimension(NUMB_VAR), intent(INOUT) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR) :: FL,FR,uL,uR,Uhll,UstarL,UstarR
  real :: sL,sR,aL,aR,uStar,pStar,pTot
  real :: numerL, denomL, numerR, denomR
  real :: dStarL, dStarR

  integer :: VELC_VAR
  VELC_VAR = 1
  if (.false.) then
     uL = vL(DENS_VAR:ENER_VAR)
     uR = vR(DENS_VAR:ENER_VAR)
!!$     call cons2flux(uL,FL,dir,vL)
!!$     call cons2flux(uR,FR,dir,vR)
  else
     call prim2flux(vL,FL,dir)
     call prim2flux(vR,FR,dir)
     call prim2cons(vL,uL)
     call prim2cons(vR,uR)
  end if

  if (dir == XDIM) then
     VELC_VAR = VELX_VAR
  elseif (dir == YDIM) then
     VELC_VAR = VELY_VAR
  elseif (dir == ZDIM) then
     VELC_VAR = VELZ_VAR
  end if

  
  ! left and right sound speed a
  aL = sqrt(vL(GAMC_VAR)*vL(PRES_VAR)/vL(DENS_VAR))
  aR = sqrt(vR(GAMC_VAR)*vR(PRES_VAR)/vR(DENS_VAR))

  ! fastest left and right going velocities
  sL = min(vL(VELC_VAR) - aL,vR(VELC_VAR) - aR)
  sR = max(vL(VELC_VAR) + aL,vR(VELC_VAR) + aR)

  ! Get HLL states for later use
  if (sL > 0.) then
     Uhll(DENS_VAR:ENER_VAR) = uL(DENS_VAR:ENER_VAR)
  elseif ((sL <= 0.) .and. (sR >= 0.)) then
     Uhll(DENS_VAR:ENER_VAR) = &
          ( sR*uR(DENS_VAR:ENER_VAR) &
           -sL*uL(DENS_VAR:ENER_VAR) &
             - FR(DENS_VAR:ENER_VAR) &
             + FL(DENS_VAR:ENER_VAR)&
           )/(sR - sL)
  else
     Uhll(DENS_VAR:ENER_VAR) = uR(DENS_VAR:ENER_VAR)
  endif

  ! Get uStar
  uStar = vR(DENS_VAR)*vR(VELC_VAR)*(sR-vR(VELC_VAR)) &
         -vL(DENS_VAR)*vL(VELC_VAR)*(sL-vL(VELC_VAR)) &
         +vL(PRES_VAR)-vR(PRES_VAR)
  uStar = uStar/( vR(DENS_VAR)*(sR-vR(VELC_VAR)) &
                 -vL(DENS_VAR)*(sL-vL(VELC_VAR)))

  ! Convenient parameters
  numerL = sL-vL(VELC_VAR)
  denomL = sL-uStar
  numerR = sR-vR(VELC_VAR)
  denomR = sR-uStar

  ! Get pStar
  pStar = vL(DENS_VAR)*numerL*(uStar-vL(VELC_VAR)) + vL(PRES_VAR)


  ! density
  dStarL = uL(DENS_VAR)*numerL/denomL
  dStarR = uR(DENS_VAR)*numerR/denomR
  select case (dir)
     case(XDIM)
        ! left and right star regions
        UstarL(DENS_VAR) = dStarL
        UstarL(MOMX_VAR) = dStarL*uStar
        UstarL(MOMY_VAR) = uL(MOMY_VAR)*numerL/denomL
        UstarL(MOMZ_VAR) = uL(MOMZ_VAR)*numerL/denomL
        UstarL(ENER_VAR) = uL(ENER_VAR)*numerL/denomL+&
             (pStar*uStar - vL(PRES_VAR)*vL(VELC_VAR))/denomL

        UstarR(DENS_VAR) = dStarR
        UstarR(MOMX_VAR) = dStarR*uStar
        UstarR(MOMY_VAR) = uR(MOMY_VAR)*numerR/denomR
        UstarR(MOMZ_VAR) = uR(MOMZ_VAR)*numerR/denomR
        UstarR(ENER_VAR) = uR(ENER_VAR)*numerR/denomR+&
             (pStar*uStar - vR(PRES_VAR)*vR(VELC_VAR))/denomR

     case(YDIM)
        ! left and right star regions
        UstarL(DENS_VAR) = dStarL
        UstarL(MOMY_VAR) = dStarL*uStar
        UstarL(MOMX_VAR) = uL(MOMX_VAR)*numerL/denomL
        UstarL(MOMZ_VAR) = uL(MOMZ_VAR)*numerL/denomL
        UstarL(ENER_VAR) = uL(ENER_VAR)*numerL/denomL+&
             (pStar*uStar - vL(PRES_VAR)*vL(VELC_VAR))/denomL

        UstarR(DENS_VAR) = dStarR
        UstarR(MOMY_VAR) = dStarR*uStar
        UstarR(MOMX_VAR) = uR(MOMX_VAR)*numerR/denomR
        UstarR(MOMZ_VAR) = uR(MOMZ_VAR)*numerR/denomR
        UstarR(ENER_VAR) = uR(ENER_VAR)*numerR/denomR+&
             (pStar*uStar - vR(PRES_VAR)*vR(VELC_VAR))/denomR
     case(ZDIM)
        ! left and right star regions
        UstarL(DENS_VAR) = dStarL
        UstarL(MOMZ_VAR) = dStarL*uStar
        UstarL(MOMX_VAR) = uL(MOMX_VAR)*numerL/denomL
        UstarL(MOMY_VAR) = uL(MOMY_VAR)*numerL/denomL
        UstarL(ENER_VAR) = uL(ENER_VAR)*numerL/denomL+&
             (pStar*uStar - vL(PRES_VAR)*vL(VELC_VAR))/denomL

        UstarR(DENS_VAR) = dStarR
        UstarR(MOMZ_VAR) = dStarR*uStar
        UstarR(MOMX_VAR) = uR(MOMX_VAR)*numerR/denomR
        UstarR(MOMY_VAR) = uR(MOMY_VAR)*numerR/denomR
        UstarR(ENER_VAR) = uR(ENER_VAR)*numerR/denomR+&
             (pStar*uStar - vR(PRES_VAR)*vR(VELC_VAR))/denomR
     end select




  ! numerical flux
  if (sL >= 0.) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR)

  elseif ( (sL < 0.) .and. (uStar >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FL(DENS_VAR:ENER_VAR) &
          + sL*(UstarL(DENS_VAR:ENER_VAR) - uL(DENS_VAR:ENER_VAR))
  elseif ( (uStar < 0.) .and. (sR >= 0.) ) then
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR) &
          + sR*(UstarR(DENS_VAR:ENER_VAR) - uR(DENS_VAR:ENER_VAR))
  else
     Flux(DENS_VAR:ENER_VAR) = FR(DENS_VAR:ENER_VAR)
  endif

  return
end subroutine hllc
