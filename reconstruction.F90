module reconstruction
  
#include "definition.h"

  use grid_data, only: gr_vL, gr_vR, gr_imax
  use sim_interfaces

contains

  subroutine recon1D_ij(ibeg, iend, sm, sp, dir, V, recon)
    !reconstruct in direction "dir" using 1D stencil ranging from i+sm to i+sp
    !do so over indices ibeg:iend and jbeg:jend
    integer, dimension(NDIM), intent(IN) :: ibeg, iend
    integer, intent(IN) :: sm, sp, dir
    
    real, dimension(NUMB_VAR, gr_imax(XDIM), gr_imax(YDIM), gr_imax(ZDIM)), intent(IN) :: V
    procedure (recon1D), pointer, intent(IN) :: recon

    integer :: i, j, k, si, sj, sk, var, Nx, i0, j0, k0, si_max, sj_max, sk_max, scntr

    real, dimension(sp-sm+1,NUMB_VAR) :: stencil

    Nx = sp-sm+1


    do i = ibeg(XDIM), iend(XDIM)
       do j = ibeg(YDIM), iend(YDIM)
          StencilDo: do k = ibeg(ZDIM), iend(ZDIM)
             !make 1D stencil centered on i,j,k
             !with extent +/- sp/sm in dir
             select case(dir)
             case(XDIM)
                i0 = i + sm
                si_max = i + sp
                j0 = j
                sj_max = j
                k0 = k
                sk_max = k
             case(YDIM)
                i0 = i
                si_max = i
                j0 = j + sm
                sj_max = j + sp
                k0 = k
                sk_max = k
             case(ZDIM)
                i0 = i
                si_max = i
                j0 = j
                sj_max = j
                k0 = k + sm
                sk_max = k + sp
             end select
#ifdef BDRY_VAR
             do si = i0, si_max
                do sj = j0, sj_max
                   do sk = k0, sk_max
                      !if (.true.) then
                      if (V(BDRY_VAR,si,sj,sk) == 1.0) then
                         !we have a physical boundary on the stencil
                         !revert to first order
                         gr_vL(:,i,j,k,dir) = V(:,i,j,k)
                         gr_vR(:,i,j,k,dir) = V(:,i,j,k)
                         cycle StencilDo
                      end if
                   end do
                end do
             end do
#endif             

             do var = 1, NUMB_VAR
                scntr = 0
                do si = i0, si_max
                   do sj = j0, sj_max
                      do sk = k0, sk_max
                         scntr = scntr + 1
                         stencil(scntr,var) = V(var, si, sj,sk)
                      end do !sk
                   end do !sj
                end do    !si
             end do       !var

             call recon(dt, Nx, stencil, gr_vL(:,i,j,k,dir), gr_vR(:,i,j,k,dir), dir)
          end do StencilDo
       end do
    end do
    
  end subroutine recon1D_ij

end module reconstruction
