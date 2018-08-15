subroutine soln_WENO(dt, Nx, V, vL, vR, dir)

#include "definition.h"

  use eigensystem
  use WENO
  use sim_data, only: sim_WENeps, sim_mval, sim_charLimiting

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: Nx, dir 
  real, dimension(Nx, NUMB_VAR), intent(IN   ) :: V
  real, dimension(NUMB_VAR),   intent(INOUT) :: vL, vR

  integer :: var, s, r
  real    :: w_norm

  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  
  real, dimension(3  ) :: lin_w, smth_ind
  real, dimension(3,3) :: ENO_coeff
  real, dimension(3,2) :: ENO_intp, nonLin_w
  real, dimension(Nx, NSYS_VAR) :: stencil
  real, dimension(NSYS_VAR) :: tempL, tempR

  lin_w(:)       = (/1./16., 10./16., 5./16 /)
  ENO_coeff(:,1) = (/ 3., -10.,  15. /)
  ENO_coeff(:,2) = (/-1.,   6.,  3.  /)
  ENO_coeff(:,3) = (/ 3.,   6., -1.  /)

  ENO_coeff = ENO_coeff/8.

  vL = V(3,:)
  vR = V(3,:)

  if (sim_charLimiting) then
     call eigenvectors(V(3,:), .false., reig, leig, dir)
  end if

  do var = 1, NSYS_VAR
     do s = 1, 5
        if (sim_charLimiting) then
           stencil(s,var) = dot_product(leig(:,var), V(s,DENS_VAR:ENER_VAR))
        else
           stencil(s,var) = V(s,var)
        end if
     end do
     !get smoothness indicators
     call betas(stencil(:,var), 2, smth_ind)
     !smth_ind = 1.
     !compute non-linear weights
     do s = 1, 3
        r = 4-s
        !WENO-JS
        nonLin_w(s,1) = lin_w(r)/(sim_WENeps + smth_ind(s))**sim_mval
        nonLin_w(s,2) = lin_w(s)/(sim_WENeps + smth_ind(s))**sim_mval
        !WENO-Z
!!$        nonLin_w(s,1) = lin_w(s)*(1.+&
!!$             abs(smth_ind(3)-smth_ind(1))/(sim_WENeps + smth_ind(r)))**sim_mval
!!$        nonLin_w(s,2) = lin_w(s)*(1.+&
!!$             abs(smth_ind(3)-smth_ind(1))/(sim_WENeps + smth_ind(s)))**sim_mval
     end do
     w_norm = SUM(nonLin_w(:,1))
     nonLin_w(:,1) = nonLin_w(:,1)/w_norm
     w_norm = SUM(nonLin_w(:,2))
     nonLin_w(:,2) = nonLin_w(:,2)/w_norm

     !calculate ENO interpolations
     do s = 1, 3
        r = 4 - s
        ENO_intp(s,1) = dot_product(stencil(s+2:s:-1,var), ENO_coeff(:,r))
        !ENO_intp(s,1) = dot_product(stencil(r+2:r:-1,var), ENO_coeff(:,s))
        ENO_intp(s,2) = dot_product(stencil(s:s+2   ,var), ENO_coeff(:,s))
     end do
     
     !take convex combination of ENO statets
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
end subroutine soln_WENO
