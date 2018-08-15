module soln_PPM_system
 
#include "definition.h"

  use grid_data
  use sim_data
  use slopeLimiter
  use eigensystem

  subroutine del_W_i(i, delW)

    use grid_data
    use sim_data
    use slopeLimiter
    use eigensystem
  
    implicit none

    integer, intent(IN) :: i
    real, dimension(NUMB_WAVE), intent(OUT) :: delW
    
    real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
    integer :: kWaveNum
    real, dimension(NUMB_VAR)  :: delV,delL,delR
    integer :: nVar
    logical :: conservative

    conservative = .false.

    call left_eigenvectors (gr_V(DENS_VAR:GAME_VAR,i),conservative,leig)
    call right_eigenvectors(gr_V(DENS_VAR:GAME_VAR,i),conservative,reig)
    
    ! primitive limiting
    if (.not. sim_charLimiting) then
       do kWaveNum = 1, NUMB_WAVE
          ! slope limiting
          ! deltas in primitive vars
          delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
          delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
          do nVar = DENS_VAR,PRES_VAR
             if (sim_limiter == 'minmod') then
                call minmod(delL(nVar),delR(nVar),delV(nVar))
             elseif (sim_limiter == 'vanLeer') then
                call vanLeer(delL(nVar),delR(nVar),delV(nVar))
             elseif (sim_limiter == 'mc') then
                call mc(delL(nVar),delR(nVar),delV(nVar))
             endif
          enddo
          ! project primitive delta to characteristic vars
          delW(kWaveNum) = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum),delV(DENS_VAR:PRES_VAR))
       enddo
    elseif (sim_charLimiting) then
       do kWaveNum = 1, NUMB_WAVE
          delL(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i  )-gr_V(DENS_VAR:PRES_VAR,i-1)
          delR(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR,i+1)-gr_V(DENS_VAR:PRES_VAR,i  )
          ! project onto characteristic vars
          pL = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum), delL(DENS_VAR:PRES_VAR))
          pR = dot_product(leig(DENS_VAR:PRES_VAR,kWaveNum), delR(DENS_VAR:PRES_VAR))
          if (sim_limiter == 'minmod')then
             call minmod(pL, pR, delW(kWaveNum))
          elseif (sim_limiter == 'vanLeer') then
             call vanLeer(pL, pR, delW(kWaveNum))
          elseif (sim_limiter == 'mc') then
             call mc(pL, pR, delW(kWaveNum))
          endif
       enddo

    endif
    return
  end subroutine del_W_i

  subroutine mono_check(vL, vR, vi)
    implicit none
    real, dimension(NSYS_VAR), intent(inout) :: vL, vR, vi

    real :: RHS, delv2
    integer :: var

    do var = PRES_VAR, DENS_VAR
       delv2 = (vL(var) - vR(var))**2
       RHS   = 6.*(vR(var) - vL(var))*(vi(var) - .5*(vR(var) + vL(var)) )

       if ((vR(var) - vi(var))*(vi(var) - vL(var)) <= 0) then
          vR(var) = vi(var)
          vL(var) = vi(var)
       elseif (-1.*delv2 > RHS) then
          vR(var) = 3.*vi(var) - 2.*vL(var)
       elseif (delv2 < RHS) then
          vL(var) = 3.*vi(var) - 2.*vR(var)
       end if
    end do
  end subroutine mono_check

  subroutine coefficients(C0, C1, C2)
    implicit none
    real, dimension(NSYS_VAR), intent(OUT) :: C0, C1, C2

    C2(DENS_VAR:PRES_VAR) = 6.*( .5*(vR(DENS_VAR:PRES_VAR) + vL(DENS_VAR:PRES_VAR)) &
         - gr_V(DENS_VAR:PRES_VAR))/(gr_dx**2)
    C1(DENS_VAR:PRES_VAR) = (vR(DENS_VAR:PRES_VAR) - vL(DENS_VAR:PRES_VAR)/gr_dx
    C0(DENS_VAR:PRES_VAR) = gr_V(DENS_VAR:PRES_VAR) - C2(DENS_VAR:PRES_VAR)*(gr_dx**2)/12.
    return
  end subroutine coefficients
    

    


end module SOLN_PPM_SYSTEM
