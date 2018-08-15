subroutine soln_gpWENO(dt, Nx, V, vL, vR, dir)
  
#include "definition.h"

  use gp_data
  use eigensystem
  use WENO
  use sim_data, only: sim_WENeps, sim_mval, sim_charLimiting

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: Nx, dir 
  real, dimension(Nx, NUMB_VAR), intent(IN   ) :: V
  real, dimension(NUMB_VAR),   intent(INOUT) :: vL, vR

  integer :: var, s, R
  real    :: w_norm, gamm_HI, gamm_LO

  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig

  real, dimension(Nx) :: stencil
  real, dimension(gp_radius+1) :: smth_ind
  real, dimension(gp_radius+1,2) :: nonLin_w, ENO_intp
  real, dimension(4) :: lin_w, nonLin
  real, dimension(NSYS_VAR) :: tempL, tempR

  R = gp_radius
  vL = V(R+1,:)
  vR = V(R+1,:)

  gamm_HI = 0.85
  gamm_LO = 0.85

  

  if (sim_charLimiting) then
     call eigenvectors(V(R+1,:), .false., reig, leig, dir)
  end if

  !build interpolation stencil
  do var = 1, NSYS_VAR
     do s = 1, Nx
        if (sim_charLimiting) then
           stencil(s) = dot_product(leig(:,var), V(s,DENS_VAR:ENER_VAR))
        else
           stencil(s) = V(s,var)
        end if
     end do

     !get smoothness indicators
     call gp_betas(stencil, R, smth_ind)
     !smth_ind = 1.

     !compute non-linear weights
     do s = 1, R+1
        nonLin_w(s,1) = gp_w(1,s)/(sim_WENeps+smth_ind(s))**sim_mval
        nonLin_w(s,2) = gp_w(2,s)/(sim_WENeps+smth_ind(s))**sim_mval
     end do
     !normalize weights
     w_norm = SUM(nonLin_w(:,1))
     nonLin_w(:,1) = nonLin_w(:,1)/w_norm
     w_norm = SUM(nonLin_w(:,2))
     nonLin_w(:,2) = nonLin_w(:,2)/w_norm

     !calculate ENO GP interpoplations
     do s = 1, R+1
        ENO_intp(s,1) = dot_product(gp_zk(1,1:R+1,s), stencil(s:s+R))
        ENO_intp(s,2) = dot_product(gp_zk(2,1:R+1,s), stencil(s:s+R))
     end do

     !take convex combination of ENO states
     tempL(var) = dot_product(nonLin_w(:,1), ENO_intp(:,1))
     tempR(var) = dot_product(nonLin_w(:,2), ENO_intp(:,2))
     
  end do

  do var = DENS_VAR, PRES_VAR
     if (sim_charLimiting) then
        vL(var) = dot_product(tempL(:), reig(var,:))
        vR(var) = dot_product(tempR(:), reig(var,:))
     else
        vL(var) = tempL(var)
        vR(var) = tempR(var)
     end if
  end do
  
  return
end subroutine soln_gpWENO
