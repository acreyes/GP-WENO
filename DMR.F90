module DMR
#include "definition.h"
  
  use sim_data, only: sim_bcT, sim_time
  use grid_data, only: gr_xCoord, gr_yCoord, gr_ny, gr_nx, gr_imax, gr_ngc, &
                       gr_ibeg, gr_iend


contains

  subroutine DMR_IN(bufferL)
    implicit none
    real, dimension(NUMB_VAR, gr_ngc, gr_ny), intent(INOUT) :: bufferL
    
    !just leave as the IC
    return
  end subroutine DMR_IN

  subroutine DMR_OUT(bufferR)
    implicit none
    real, dimension(NUMB_VAR, gr_ngc, gr_ny), intent(INOUT) :: bufferR

    !just leave as IC
    return
  end subroutine DMR_OUT

  subroutine DMR_bot(bufferB)
    implicit none
    real, dimension(NUMB_VAR, gr_nx, gr_ngc), intent(INOUT) :: bufferB
    
    real, dimension(NUMB_VAR, gr_nx, gr_ngc) :: V_dat

    integer :: i, j

    V_dat(:,:,:) = bufferB(:,:,:)
    do i = gr_ibeg(XDIM), gr_iend(XDIM)
       if (gr_xCoord(i) > 1./6.) then
          do j = 0, gr_ngc-1
             bufferB(:,i-gr_ngc,j+1) = V_dat(:,i-gr_ngc,gr_ngc-j)
             bufferB(VELY_VAR,i-gr_ngc,j+1) = -1.*V_dat(VELY_VAR,i-gr_ngc,gr_ngc-j)
          end do
       else
          bufferB(DENS_VAR,i-gr_ngc,:) = 8.
          bufferB(VELX_VAR,i-gr_ngc,:) = 7.1447096
          bufferB(VELY_VAR,i-gr_ngc,:) = -4.125
          bufferB(PRES_VAR,i-gr_ngc,:) = 116.5
       end if
    end do

    return
  end subroutine DMR_bot

  subroutine DMR_top(bufferT)
    implicit none
    real, dimension(NUMB_VAR, gr_nx, gr_ngc), intent(INOUT) :: bufferT

    integer :: i,j
    real    ::  xmin, uLt, x, y, rt3

    rt3 = sqrt(3.)
    xmin = 1./6. + 10.*sim_bcT/(.5*rt3) + 1./rt3
    !xmin = 1./6. + sqrt(3.) + 7.1447096*sim_time
    !uLt = 7.1447096*sim_time
    !print *, sim_time, sim_bcT

    do i = gr_ibeg(XDIM), gr_iend(XDIM)
       do j = 1, gr_ngc
          x = gr_xCoord(i)
          y = gr_yCoord(j+gr_iend(YDIM))
          if (y-1. > (x-xmin)*rt3) then
          !if ( x < xmin) then
             !print *, i , x, this_image(), gr_xCoord(gr_ibeg(XDIM))
             !still in shock, use left state
             bufferT(DENS_VAR,i-gr_ngc,j) = 8.
             bufferT(VELX_VAR,i-gr_ngc,j) = 7.1447096
             bufferT(VELY_VAR,i-gr_ngc,j) = -4.125
             bufferT(PRES_VAR,i-gr_ngc,j) = 116.5
          else
             !outside of shock
             bufferT(DENS_VAR,i-gr_ngc,j) = 1.4
             bufferT(VELX_VAR,i-gr_ngc,j) = 0.
             bufferT(VELY_VAR,i-gr_ngc,j) = 0.
             bufferT(PRES_VAR,i-gr_ngc,j) = 1.
          end if
       end do
    end do

    return
  end subroutine DMR_top

end module DMR
  
