module primconsflux

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_gamma, sim_smallPres
  use eos, only : eos_cell
  
contains

  subroutine prim2cons(V,U)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: U

    real :: ekin, eint

    U(DENS_VAR) = V(DENS_VAR)
    U(MOMX_VAR) = V(DENS_VAR)*V(VELX_VAR)
    U(MOMY_VAR) = V(DENS_VAR)*V(VELY_VAR)
    U(MOMZ_VAR) = V(DENS_VAR)*V(VELZ_VAR)
    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2+V(VELY_VAR)**2+V(VELZ_VAR)**2)
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    U(ENER_VAR) = ekin + eint
     
  end subroutine prim2cons


  subroutine cons2prim(U,V)
    implicit none
    real, dimension(NSYS_VAR), intent(IN)  :: U
    real, dimension(NUMB_VAR), intent(OUT) :: V
    real :: eint, ekin, pres
    
    V(DENS_VAR) = U(DENS_VAR)
    V(VELX_VAR) = U(MOMX_VAR)/U(DENS_VAR)
    V(VELY_VAR) = U(MOMY_VAR)/U(DENS_VAR)
    V(VELZ_VAR) = U(MOMZ_VAR)/U(DENS_VAR)
    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2 + V(VELY_VAR)**2 + V(VELZ_VAR)**2)
    eint = max(U(ENER_VAR) - ekin, sim_smallPres) !eint=rho*e
    eint = eint/U(DENS_VAR)
    ! get pressure by calling eos
    call eos_cell(U(DENS_VAR),eint,sim_gamma,pres)
    V(PRES_VAR) = pres
    V(EINT_VAR) = eint*U(DENS_VAR)
    V(GAMC_VAR) = sim_gamma
    V(GAME_VAR) = sim_gamma
    
  end subroutine cons2prim

  subroutine prim2flux(V,Flux,dir)
    implicit none
    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NSYS_VAR), intent(OUT) :: Flux
    integer, intent(IN) :: dir

    real :: ekin,eint,ener
    integer :: VEL1_VAR, VEL2_VAR, VEL3_VAR

    if (dir == XDIM) then
       VEL1_VAR = VELX_VAR
       VEL2_VAR = VELY_VAR
       VEL3_VAR = VELZ_VAR
    elseif (dir == YDIM) then
       VEL1_VAR = VELY_VAR
       VEL2_VAR = VELZ_VAR
       VEL3_VAR = VELX_VAR
    elseif (dir == ZDIM) then
       VEL1_VAR = VELZ_VAR
       VEL2_VAR = VELX_VAR
       VEL3_VAR = VELY_VAR
    end if
    
    
    ekin = 0.5*V(DENS_VAR)*(V(VELX_VAR)**2 + V(VELY_VAR)**2 + V(VELZ_VAR)**2)
    eint = V(PRES_VAR)/(V(GAME_VAR)-1.)
    ener = ekin + eint
    
    Flux(DENS_VAR) = V(DENS_VAR)*V(VEL1_VAR)
    Flux(VEL1_VAR) = Flux(DENS_VAR)*V(VEL1_VAR) + V(PRES_VAR)
    Flux(VEL2_VAR) = Flux(DENS_VAR)*V(VEL2_VAR)
    Flux(VEL3_VAR) = Flux(DENS_VAR)*V(VEL3_var)
    Flux(ENER_VAR) = V(VEL1_VAR)*(ener + V(PRES_VAR))
    
  end subroutine prim2flux

  subroutine cons2flux
    implicit none
  end subroutine cons2flux
  
end module primconsflux
