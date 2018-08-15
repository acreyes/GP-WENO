!!!!!!!!!!!!
! here we are to initialize the block information used to facilitate
! the communications necessary for parallelization

subroutine block_init()

#include "definition.h"

  use block_data
  use sim_data
  use grid_data
  use bc, only: bc_init

  implicit none

  integer :: i, j,k, nBlockx, nBlocky, nBlockz, leftover

  !image information
  bl_nProcs = num_images()
  bl_ID     = this_image()
  
  if (bl_iProcs*bl_jProcs*bl_kProcs .ne. bl_nprocs) then
     print *, 'iProcs*jProcs*bl_kProcs must be number of procs'
     stop
  end if
 
  ! make grid
  allocate(bl_grid(bl_iProcs, bl_jProcs, bl_kProcs)); bl_grid = 0
  do i = 1, bl_iProcs
     do j = 1, bl_jProcs
        do k = 1, bl_kProcs
           bl_grid(i,j,k) = j + (i-1)*bl_jProcs + (k-1)*bl_jProcs*bl_iProcs
           if (bl_ID == 1) print *, bl_grid(i,j,k), i, j, k
           if (bl_grid(i,j,k) == bl_ID) then
              bl_i = i
              bl_j = j
              bl_k = k
           end if
        end do
     end do
  end do

!!$  if (bl_ID == 1) then
!!$     print *, "grid"
!!$     print *, bl_grid
!!$  end if
  
  !now we assign grid information to block
  nBlockx = gr_nx / bl_iProcs
  nBlocky = gr_ny / bl_jProcs
  nBlockz = gr_nz / bl_kProcs
  
  gr_glb_nx = gr_nx
  gr_glb_ny = gr_ny
  gr_glb_nz = gr_nz
  
  ! the first and the last interior cell index
  allocate(gr_ibeg(NDIM)); allocate(gr_iend(NDIM))
  allocate(gr_i0(  NDIM)); allocate(gr_imax(NDIM))
  gr_ibeg = gr_ngc + 1
!!$  gr_ibeg(YDIM) = gr_ngc + 1
!!$  gr_ibeg(ZDIM) = gr_ngc + 1
  gr_i0 = 1
!!$  gr_i0(XDIM)   = 1 
!!$  gr_i0(YDIM)   = 1
   
  if (bl_i == bl_iProcs) then
     !last proc in i-dir
     leftover = gr_nx - bl_iProcs*nBlockx
     gr_iend(XDIM) = gr_ngc + nBlockx + leftover
     gr_nx = nBlockx + leftover
  else
     gr_iend(XDIM) = gr_ngc + nBlockx
     gr_nx = nBlockx
  end if
  
  if (bl_j == bl_jProcs) then
     leftover = gr_ny - bl_jProcs*nBlocky
     gr_iend(YDIM) = gr_ngc + nBlocky + leftover
     gr_ny = nBlocky + leftover
  else
     gr_iend(YDIM) = gr_ngc + nBlocky
     gr_ny = nBlocky
  end if

  if (bl_k == bl_kProcs) then
     !last proc in k-dir
     leftover = gr_nz - bl_kProcs*nBlockz
     gr_iend(ZDIM) = gr_ngc + nBlockz + leftover
     gr_nz = nBlockz + leftover
  else
     gr_iend(ZDIM) = gr_ngc + nBlockz
     gr_nz = nBlockz
  end if

  gr_imax(XDIM) = gr_iend(XDIM) + gr_ngc
  gr_imax(YDIM) = gr_iend(YDIM) + gr_ngc
  gr_imax(ZDIM) = gr_iend(ZDIM) + gr_ngc

  allocate(bl_nBLock(NDIM))
  bl_nBlock(XDIM) = nBlockx
  bl_nBlock(YDIM) = nBlocky
  bl_nBlock(ZDIM) = nBlockz


  
  !now we figure out the boundary conditions for the block
  !meaning global domain boundary or block-block boundary
!!$                2      
!!$             ----------
!!$            |          |
!!$            |          |         
!!$         1  |          | 3
!!$            |          |
!!$            |          |
!!$             ----------
!!$                 4      

  !BCS: 1: left x
  !     2: right y
  !     3: right x
  !     4: left y
  !     5: left z
  !     6: right z

  
  !convert bc string to integer code
  !add_back
  call bc_init(sim_bcTypex, sim_xBC)
  call bc_init(sim_bcTypey, sim_yBC)
  call bc_init(sim_bcTypez, sim_zBC)

  
  allocate(bl_BC(2*NDIM)); bl_BC = 0

  !this takes care of domain boundaries
  if (bl_i == 1        ) bl_BC(1) = sim_xBC !left face
  if (bl_i == bl_iProcs) bl_BC(3) = sim_xBC !right face
  
  if (bl_j == 1        ) bl_BC(4) = sim_yBC !bottom face
  if (bl_j == bl_jProcs) bl_BC(2) = sim_yBC !top face

  if (bl_k == 1        ) bl_BC(5) = sim_zBC
  if (bl_k == bl_kProcs) bl_BC(6) = sim_zBC

  !now we check for block-block boundaries
  if (bl_BC(1) == 0) bl_BC(1) = bl_grid(bl_i-1,bl_j,bl_k) !left face
  if (bl_BC(3) == 0) bl_BC(3) = bl_grid(bl_i+1,bl_j,bl_k) !right face

  if (bl_BC(4) == 0) bl_BC(4) = bl_grid(bl_i,bl_j-1,bl_k) !bottom face
  if (bl_BC(2) == 0) bl_BC(2) = bl_grid(bl_i,bl_j+1,bl_k) ! top face

  if (bl_BC(5) == 0) bl_BC(5) = bl_grid(bl_i,bl_j,bl_k-1) !down face
  if (bl_BC(6) == 0) bl_BC(6) = bl_grid(bl_i,bl_j,bl_k+1) !up face



  !now we allocate the coarray buffers
  allocate(bl_buffL(NUMB_VAR, gr_ngc, gr_ny, gr_nz)[*])
  allocate(bl_buffR(NUMB_VAR, gr_ngc, gr_ny, gr_nz)[*])
  allocate(bl_buffB(NUMB_VAR, gr_nx, gr_ngc, gr_nz)[*])
  allocate(bl_buffT(NUMB_VAR, gr_nx, gr_ngc, gr_nz)[*])
  allocate(bl_buffD(NUMB_VAR, gr_nx, gr_ny, gr_ngc)[*])
  allocate(bl_buffU(NUMB_VAR, gr_nx, gr_ny, gr_ngc)[*])
!!!here I am 2/2/18
  !here we deal with the corner BCs
!!$           1             2
!!$             ----------
!!$            |          |
!!$            |          |         
!!$            |          | 
!!$            |          |
!!$            |          |
!!$             ----------  3
!!$           4            

!!$  if (sim_reconMultiD) then
!!$     allocate(bl_cornerBC(2*NDIM)); bl_cornerBC = 0
!!$     !domain edges
!!$     if (bl_i == 1) then
!!$        bl_cornerBC(1) = -bl_grid(bl_i,bl_j+1)
!!$        bl_cornerBC(4) = -bl_grid(bl_i,bl_j-1)
!!$     end if
!!$     if (bl_i == bl_iProcs) then
!!$        bl_cornerBC(2) = -bl_grid(bl_i,bl_j+1)
!!$        bl_cornerBC(3) = -bl_grid(bl_i,bl_j-1)
!!$     end if
!!$     if (bl_j == 1) then
!!$        bl_cornerBC(3) = -bl_grid(bl_i+1,bl_j)
!!$        bl_cornerBC(4) = -bl_grid(bl_i-1,bl_j)
!!$     end if
!!$     if (bl_j == bl_jProcs) then
!!$        bl_cornerBC(1) = -bl_grid(bl_i-1,bl_j)
!!$        bl_cornerBC(2) = -bl_grid(bl_i+1,bl_j)
!!$     end if
!!$     !block-block corners
!!$     if (bl_cornerBC(1) == 0) bl_cornerBC(1) = bl_grid(bl_i-1,bl_j+1)
!!$     if (bl_cornerBC(2) == 0) bl_cornerBC(2) = bl_grid(bl_i+1,bl_j+1)
!!$     if (bl_cornerBC(3) == 0) bl_cornerBC(3) = bl_grid(bl_i+1,bl_j-1)
!!$     if (bl_cornerBC(4) == 0) bl_cornerBC(4) = bl_grid(bl_i-1,bl_j-1)
!!$     !domain corners
!!$     if (bl_i == 1         .and. bl_j == 1        ) bl_cornerBC(4) = 0
!!$     if (bl_i == 1         .and. bl_j == bl_jProcs) bl_cornerBC(1) = 0
!!$     if (bl_i == bl_iProcs .and. bl_j == 1        ) bl_cornerBC(3) = 0
!!$     if (bl_i == bl_iProcs .and. bl_j == bl_jProcs) bl_cornerBC(2) = 0
!!$
!!$     
!!$     allocate(bl_buffBR(NUMB_VAR, gr_ngc, gr_ngc)[*])
!!$     allocate(bl_buffTR(NUMB_VAR, gr_ngc, gr_ngc)[*])
!!$     allocate(bl_buffBL(NUMB_VAR, gr_ngc, gr_ngc)[*])
!!$     allocate(bl_buffTL(NUMB_VAR, gr_ngc, gr_ngc)[*])
!!$
!!$     allocate(bl_bufCB(NUMB_VAR,gr_ngc,gr_ngc,2)[*])
!!$     allocate(bl_bufCT(NUMB_VAR,gr_ngc,gr_ngc,2)[*])
!!$     allocate(bl_bufCL(NUMB_VAR,gr_ngc,gr_ngc,2)[*])
!!$     allocate(bl_bufCR(NUMB_VAR,gr_ngc,gr_ngc,2)[*])
!!$  end if

  allocate(bl_delT[*])
  allocate(bl_maxSpeed[*])
  return



end subroutine block_init
