module char_tracing

#include "definition.h"

  use grid_data
  use sim_data

contains

subroutine soln_PPMtracing(dt, V, C0, C1, C2, vL, vR, lambda, leig0, reig0)
  !subroutine to do characteristic tracing for a particular cell for PPM and put the edge states into vL and vR

  implicit none
  real, intent(IN)                                 :: dt
  real, dimension(NSYS_VAR          ), intent(IN ) :: C0, C1, C2
  real, dimension(NSYS_VAR          ), intent(OUT) :: vL, vR
  real, dimension(NUMB_VAR          ), intent(IN ) :: V
  real, dimension(NUMB_WAVE         ), intent(IN ) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE), intent(IN) :: reig0, leig0

  integer :: kWaveNum, i, j
  real :: lambdaDtDx, delC1, delC2
  real, dimension(NSYS_VAR) :: sigL, sigR, vecL, vecR

  !set initial sum to 0
  sigL(DENS_VAR:ENER_VAR) = 0.
  sigR(DENS_VAR:ENER_VAR) = 0.

  do kWaveNum = 1, NUMB_WAVE
     lambdaDtDx = lambda(kWaveNum)*dt/gr_dx
     delC1 = dot_product(leig0(DENS_VAR:PRES_VAR, kWaveNum), C1(DENS_VAR:PRES_VAR))
     delC2 = dot_product(leig0(DENS_VAR:PRES_VAR, kWaveNum), C2(DENS_VAR:PRES_VAR))

     if (sim_riemann == 'roe') then
        if (lambda(kWaveNum) > 0.) then
           vecR(DENS_VAR:PRES_VAR) = &
                0.5*( 1. -    lambdaDtDx                        )*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC1 + &
                .25*( 1. - 2.*lambdaDtDx + 4./3.*(lambdaDtDx**2))*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC2
           sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)
        elseif (lambda(kWaveNum) < 0.) then
           vecL(DENS_VAR:PRES_VAR) = &
                0.5*(-1. -    lambdaDtDx                        )*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC1 + &
                .25*( 1. + 2.*lambdaDtDx + 4./3.*(lambdaDtDx**2))*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC2
           sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
        end if
     else
        vecR(DENS_VAR:PRES_VAR) = &
             0.5*( 1. -    lambdaDtDx                        )*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC1 + &
             .25*( 1. - 2.*lambdaDtDx + 4./3.*(lambdaDtDx**2))*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC2
        sigR(DENS_VAR:PRES_VAR) = sigR(DENS_VAR:PRES_VAR) + vecR(DENS_VAR:PRES_VAR)
        
        vecL(DENS_VAR:PRES_VAR) = &
             0.5*(-1. -    lambdaDtDx                        )*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC1 + &
             .25*( 1. + 2.*lambdaDtDx + 4./3.*(lambdaDtDx**2))*reig0(DENS_VAR:PRES_VAR, kWaveNum)*delC2
        sigL(DENS_VAR:PRES_VAR) = sigL(DENS_VAR:PRES_VAR) + vecL(DENS_VAR:PRES_VAR)
     end if
  enddo
  
  vL(DENS_VAR:PRES_VAR) = C0(DENS_VAR:PRES_VAR) + sigL(DENS_VAR:PRES_VAR) 
  vR(DENS_VAR:PRES_VAR) = C0(DENS_VAR:PRES_VAR) + sigR(DENS_VAR:PRES_VAR) 


end subroutine soln_PPMtracing

end module char_tracing
