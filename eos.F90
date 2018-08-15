module eos

#include "definition.h"
  
  use grid_data
  use sim_data, only : sim_smallPres

contains

  subroutine eos_all()
    implicit none

    integer :: i,j
    real :: eint, ekin,pres,velx,vely
    

    ! dens-eint mode:
    ! This assumes that conservative vars are known
    ! (but not prim vars yet) to evaluate pres
    ! == inputs : dens & eint
    ! == outputs: pres
    ! This routine fills out pressure also the GC regions

    !I think this is deprecated and not needed anymore
    
!!$    do i =  gr_i0(XDIM),gr_imax(XDIM)
!!$       do j = gr_i0(YDIM),gr_imax(YDIM)
!!$          ! updated velx, ekin, eint
!!$          velx = gr_U(MOMX_VAR,i,j)/gr_U(DENS_VAR,i,j)
!!$          vely = gr_U(MOMY_VAR,i,j)/gr_U(DENS_VAR,i,j)
!!$          ekin = 0.5*(velx*velx+vely*vely)*gr_U(DENS_VAR,i,j)
!!$          eint = max(gr_U(ENER_VAR,i,j) - ekin,sim_smallPres)
!!$
!!$          ! now let's get a new pressure
!!$          pres = (gr_V(GAME_VAR,i,j)-1.)*eint
!!$          gr_V(PRES_VAR,i,j) = max(sim_smallPres,pres)
!!$       end do
!!$    end do
  end subroutine eos_all


  
  subroutine eos_cell(dens,eint,game,pres)
    implicit none
    real, intent(IN) :: dens,eint,game
    real, intent(OUT):: pres
    
    ! ideal gas law
    ! eint = pres/dens/(game-1)
    pres = max((game-1.)*dens*eint,sim_smallPres)

    
  end subroutine eos_cell

end module eos
