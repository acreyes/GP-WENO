subroutine soln_GP(dt, V, dir, nFlux)
#include "definition.h"

  use grid_data
  use sim_data
  use eigensystem
  use primconsflux

  implicit none
  real, intent(IN) :: dt
  integer, intent(IN) :: dir
  real, dimension(NUMB_VAR,gr_imax(dir), 2*gr_Tcells+1), intent(IN) :: V
  real, dimension(NSYS_VAR, gr_imax(dir)), intent(INOUT) :: nFlux

  integer :: i, var, s, x, y
  logical :: conservative

  real ::  f0L, f0R
  real, dimension(NSYS_VAR) :: vecL, vecR, Us, Uss, Fs, Fss
  real, dimension(NSYS_VAR, gr_gp_StencilPts, 2) :: flux
  real, dimension(NUMB_VAR) :: Viph

  real, dimension(NUMB_WAVE) :: lambda, maxalphas
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  real, dimension(gr_GP_stencilPts) :: un


  conservative = .true.
  
  MaxAlphas = 0.
  do i = gr_ibeg(dir), gr_iend(dir)
     !find the max wavespeeds throughout the domain
     call eigenvalues(V(:,i,gr_Tcells+1), lambda, dir)
     do var = SHOCKLEFT,SHOCKRGHT
        MaxAlphas(var) = MAX(MaxAlphas(var),abs(lambda(var)))
     end do
  end do

  
  do i = gr_ibeg(dir)-1,gr_iend(dir)

     call eigenvalues(V(DENS_VAR:GAME_VAR,i,gr_Tcells+1),lambda,dir)
     call prim2cons(V(:,i,gr_Tcells+1), Us)
     call prim2cons(V(:,i+1,gr_Tcells+1), Uss)
     call cons2prim(0.5*(Us + Uss), Viph)
     
     
     call left_eigenvectors (Viph(DENS_VAR:GAME_VAR),conservative,leig,dir)
     call right_eigenvectors(Viph(DENS_VAR:GAME_VAR),conservative,reig,dir)



     !do flux-splitting
     do s = 1,gr_gp_StencilPts
        x = i + INT(gr_gp_stencil(s,1))
        y = gr_Tcells + 1 + INT(gr_gp_stencil(s,2))
        ![s] stencil
        call prim2cons(V(DENS_VAR:NUMB_VAR,x,y), Us)
        call prim2flux(V(:,x,y), Fs, dir)
        ![s'] stencil
        call prim2cons(V(:,x+1,y), Uss)
        call prim2flux(V(:,x+1,y), Fss, dir)

        do var = 1, NUMB_WAVE
           flux(var, s, 2) = 0.5*dot_product(leig(:,var), Fs(:)  + gr_maxAlphas(var,dir)*Us(:) )
           flux(var, s, 1) = 0.5*dot_product(leig(:,var), Fss(:)  - gr_maxAlphas(var,dir)*Uss(:) )
        end do
     end do


     !if (dir == YDIM) stop
     !done w/ flux splitting

     do var = 1,NUMB_WAVE
     

        !now we begin the GP reconstruction
        un = 1.


        !calculate the most probable mean f0 (see eq 33)
        f0L = dot_product(gr_GPv, flux(var, :, 1))/dot_product(gr_GPv, un)
        f0R = dot_product(gr_GPv, flux(var, :, 2))/dot_product(gr_GPv, un)

        !calculate the left and right states
        vecL(var) = f0L + dot_product(gr_GPZ(1, :),flux(var, :, 1)) - f0L*dot_product(gr_GPZ(1, :),un)
        vecR(var) = f0R + dot_product(gr_GPZ(2, :),flux(var, :, 2)) - f0R*dot_product(gr_GPZ(2, :),un)
     end do !var

     !now all that is left is to recombine the split fluxes into the interface flux
     do var = DENS_VAR,PRES_VAR
        nflux(var,i+1) = dot_product( vecL(:) + vecR(:), reig(var,:))
     end do

  end do !i

end subroutine soln_GP
