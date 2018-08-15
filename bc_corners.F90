subroutine bc_corners(V)
#include "definition.h"

  use block_data
  use grid_data, only: gr_imax, gr_ngc, gr_iend, gr_i0, gr_nx, gr_ny, gr_ibeg
  use sim_data , only: sim_bcTypex, sim_bcTypey     
  
  implicit none

  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM)), intent(INOUT) :: V

  integer :: dest, gc, ibeg, iend, jbeg, jend, i0, j0, imax, jmax

  real, dimension(NUMB_VAR,gr_ngc,gr_ngc) :: loc_buffTL, loc_buffTR, loc_buffBR, loc_buffBL

  !needs to be updated for 3D

!!$  i0   = gr_i0(  XDIM)
!!$  imax = gr_imax(XDIM)
!!$  ibeg = gr_ibeg(XDIM)
!!$  iend = gr_iend(XDIM)
!!$
!!$  j0   = gr_i0(  YDIM)
!!$  jmax = gr_imax(YDIM)
!!$  jbeg = gr_ibeg(YDIM)
!!$  jend = gr_iend(YDIM)
!!$
!!$  !load halo data into image buffers
!!$  ! 1D buffers. I'm loading guard cells for use in corners
!!$  ! should work since bc_corners is called after filling 1D direction GCs
!!$  bl_bufCL(:, :, :, 1) = V(:, i0    :ibeg-1, jbeg         :jbeg+gr_ngc-1)
!!$  bl_bufCR(:, :, :, 1) = V(:, iend+1:imax  , jbeg         :jbeg+gr_ngc-1)
!!$  bl_bufCL(:, :, :, 2) = V(:, i0    :ibeg-1, jend-gr_ngc+1:jend         )
!!$  bl_bufCR(:, :, :, 2) = V(:, iend+1:imax  , jend-gr_ngc+1:jend         )
!!$  
!!$  bl_bufCB(:, :, :, 1) = V(:, ibeg         :ibeg+gr_ngc-1, j0    :jbeg-1)
!!$  bl_bufCT(:, :, :, 1) = V(:, ibeg         :ibeg+gr_ngc-1, jend+1:jmax  )
!!$  bl_bufCB(:, :, :, 2) = V(:, iend-gr_ngc+1:iend         , j0    :jbeg-1)
!!$  bl_bufCT(:, :, :, 2) = V(:, iend-gr_ngc+1:iend         , jend+1:jmax  )
!!$
!!$  !corner buffers
!!$  bl_buffBL(:,:,:) = V(:, i0           :i0+gr_ngc-1, j0           :j0+gr_ngc-1)
!!$  bl_buffTL(:,:,:) = V(:, i0           :i0+gr_ngc-1, jend-gr_ngc+1:jend       )
!!$  bl_buffBR(:,:,:) = V(:, iend-gr_ngc+1:iend       , j0           :j0+gr_ngc-1)
!!$  bl_buffTR(:,:,:) = V(:, iend-gr_ngc+1:iend       , jend-gr_ngc+1:jend       )
!!$
!!$  Sync All
!!$
!!$  !do block-block corners
!!$  dest = bl_cornerBC(1)
!!$  if (dest > 0) then
!!$     loc_buffTL(:,:,:) = bl_buffBR(:,:,:)[dest]
!!$  end if
!!$
!!$  dest = bl_cornerBC(2)
!!$  if (dest > 0) then
!!$     loc_buffTR(:,:,:) = bl_buffBL(:,:,:)[dest]
!!$  end if
!!$
!!$  dest = bl_cornerBC(3)
!!$  if (dest > 0) then
!!$     loc_buffBR(:,:,:) = bl_buffTL(:,:,:)[dest]
!!$  end if
!!$
!!$  dest = bl_cornerBC(4)
!!$  if (dest > 0) then
!!$     loc_buffBL(:,:,:) = bl_buffTR(:,:,:)[dest]
!!$  end if
!!$
!!$  !now we do the domain edge corners
!!$  if (bl_i == 1 .or. bl_i == bl_iProcs) then
!!$     !left or right domain boundaries
!!$     dest = bl_cornerBC(1)
!!$     if (dest < 0) then
!!$        loc_buffTL(:,:,:) = bl_bufCL(:,:,:,1)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(2)
!!$     if (dest < 0) then
!!$        loc_buffTR(:,:,:) = bl_bufCR(:,:,:,1)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(3)
!!$     if (dest < 0) then
!!$        loc_buffBR(:,:,:) = bl_bufCR(:,:,:,2)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(4)
!!$     if (dest < 0) then
!!$        loc_buffBL(:,:,:) = bl_bufCL(:,:,:,2)[-dest]
!!$     end if
!!$  end if
!!$
!!$  if (bl_j == 1 .or. bl_j == bl_jProcs)then
!!$     !top or bottom domain boundaries
!!$     dest = bl_cornerBC(1)
!!$     if (dest < 0) then
!!$        loc_buffTL(:,:,:) = bl_bufCT(:,:,:,2)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(2)
!!$     if (dest < 0) then
!!$        loc_buffTR(:,:,:) = bl_bufCT(:,:,:,1)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(3)
!!$     if (dest < 0) then
!!$        loc_buffBR(:,:,:) = bl_bufCB(:,:,:,1)[-dest]
!!$     end if
!!$
!!$     dest = bl_cornerBC(4)
!!$     if (dest < 0) then
!!$        loc_buffBL(:,:,:) = bl_bufCB(:,:,:,2)[-dest]
!!$     end if
!!$  end if
!!$
!!$
!!$  !now we take care to hanlde the domain corners
!!$  if (sim_bcTypex == 'periodic' .and. sim_bcTypey == 'periodic') then
!!$     if (bl_i == 1) then
!!$        if (bl_j == 1) then
!!$           !BL corner
!!$           loc_buffBL(:,:,:) = bl_buffTR(:,:,:)[bl_grid(bl_iProcs,bl_jProcs)]
!!$        end if
!!$        if (bl_j == bl_jProcs) then
!!$           !TL corner
!!$           loc_buffTL(:,:,:) = bl_buffBR(:,:,:)[bl_grid(bl_iProcs,1)]
!!$        end if
!!$     end if
!!$     if (bl_i == bl_iProcs) then
!!$        if (bl_j == 1) then
!!$           !BR corner
!!$           loc_buffBR(:,:,:) = bl_buffTL(:,:,:)[bl_grid(1,bl_jProcs)]
!!$        end if
!!$        if (bl_j == bl_jProcs) then
!!$           !TR corner
!!$           loc_buffTR(:,:,:) = bl_buffBL(:,:,:)[bl_grid(1,1)]
!!$        end if
!!$     end if
!!$     
!!$  end if
!!$
!!$  !now we fill corner guard cells
!!$  V(:,i0    :ibeg-1,j0    :jbeg-1) = loc_buffBL(:,:,:)
!!$  V(:,i0    :ibeg-1,jend+1:jmax  ) = loc_buffTL(:,:,:)
!!$  V(:,iend+1:imax  ,j0    :jbeg-1) = loc_buffBR(:,:,:)
!!$  V(:,iend+1:imax  ,jend+1:jmax  ) = loc_buffTR(:,:,:)

  return
end subroutine bc_corners
