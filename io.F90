module io
 
#include "definition.h"
  
  use grid_data, only : gr_xCoord, gr_ibeg, gr_iend &
       ,gr_yCoord, gr_zCoord, gr_V, gr_nx, gr_ny, gr_nz &
       , gr_glb_nx, gr_glb_ny, gr_glb_nz, gr_xbeg &
       , gr_dx, gr_dy, gr_dz
  
  use sim_data, only : sim_name, sim_hdf5
  use block_data
  implicit none

  
  integer, save :: nCounter
  
contains

  subroutine io_writeOutput(t,nstep,ioCounter)

    implicit none
    real, intent(IN) :: t
    integer, intent(IN) :: nstep, ioCounter
    if (.true.) then
       call io_writeHDF5(t,nstep,ioCounter)
    else
       call io_writeASCII(t,nstep,ioCounter)
    end if
  end subroutine io_writeOutput

  subroutine io_writeHDF5(t,nstep,ioCounter)

    use HDF5
    implicit none

    real, intent(IN) :: t
    integer, intent(IN) :: nstep,ioCounter

    

    integer :: i,j,k,nVar,nCell,dest
    character(len=50) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend, kbeg, kend
    real, dimension(gr_glb_nx) :: xCoord, xread
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(gr_glb_nz) :: zCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny, gr_glb_nz) :: V, V_read
    integer, dimension(NDIM)   :: img_size
    
    integer, allocatable :: nBlock(:)[:]
    real,    allocatable :: img_V(:,:,:,:)[:], loc_V(:,:,:,:)

    character(len=50) :: dset_prim, dset_x, dset_y, dset_z
    integer(HID_T)    :: file_id, dspace_id, dset_id
    integer           :: error, rank_V, rank_XYZ
    integer(HSIZE_T),allocatable, dimension(:) :: dims_V, dims_XYZ
    !make image buffers for grid
    allocate(img_V(NUMB_VAR, gr_nx, gr_ny, gr_nz)[*])
    allocate(nBlock(NDIM)[*])
    nBlock(XDIM) = gr_nx
    nBlock(YDIM) = gr_ny
    nBlock(ZDIM) = gr_nz
    img_V(:, :, :, :) = gr_V(:, gr_ibeg(XDIM):gr_iend(XDIM), &
                                gr_ibeg(YDIM):gr_iend(YDIM), gr_ibeg(ZDIM):gr_iend(ZDIM))

    Sync All

    if (bl_ID == 1) then
       !I am root
       dset_prim = "prim_vars"
       dset_x    = "xCoord"
       dset_y    = "yCoord"
       dset_z    = "zCoord"
       do i = 1, gr_glb_nx
          xCoord(i) = (real(i)-0.5)*gr_dx + gr_xbeg(XDIM)
       end do
       do j = 1, gr_glb_ny
          yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
       end do
       do k = 1, gr_glb_nz
          zCoord(k) = (real(k)-0.5)*gr_dz + gr_xbeg(ZDIM)
       end do
       do i = 1, bl_iProcs
          do j = 1, bl_jProcs
             do k = 1, bl_kProcs
                dest = bl_grid(i,j,k)
                img_size(:) = nBlock(:)[dest]

                ibeg = (i-1)*bl_nBlock(XDIM) + 1
                jbeg = (j-1)*bl_nBlock(YDIM) + 1
                kbeg = (k-1)*bl_nBlock(ZDIM) + 1
                iend = ibeg -1 + img_size(XDIM)
                jend = jbeg -1 + img_size(YDIM)
                kend = kbeg -1 + img_size(ZDIM)

                allocate(loc_V(NUMB_VAR,img_size(XDIM),img_size(YDIM),img_size(ZDIM)))
                loc_V(:,:,:,:) = img_V(:,:,:,:)[dest]
                V(:, ibeg:iend, jbeg:jend, kbeg:kend) = loc_V(:,:,:,:)
                deallocate(loc_V)
             end do !k
          end do !j
       end do !i
       
       ! convert conter number to character
       write(cCounter,910) ioCounter + 10000

       ! file name for ascii output
       !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
       ofile = trim(sim_name)//'_'//cCounter//'.slug'
       allocate(dims_V(NDIM+1))
       allocate(dims_XYZ(1))
       rank_XYZ = 1
       dims_V = (/ NSYS_VAR,  gr_glb_nx, gr_glb_ny, gr_glb_nz/)
       rank_V = NDIM + 1
       
       !hdf5 fortran interface
       call h5open_f(error)
       !create data file
       call h5fcreate_f(ofile, H5F_ACC_TRUNC_F, file_id, error)

       !datatspace for primitive vars
       call h5screate_simple_f(rank_V, dims_V, dspace_id, error)
       !dataset for primitive vars
       call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write prim vars
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, V(DENS_VAR:NSYS_VAR,:,:,:), dims_V, error)
       !close dataset
       call h5dclose_f(dset_id, error)
       !close dataspace
       call h5sclose_f(dspace_id, error)

#ifdef BDRY_VAR       
       !save physical boundary info
       dims_V = (/ 1, gr_glb_nx, gr_glb_ny, gr_glb_nz /)
       dset_prim = "bdry_var"
       call h5screate_simple_f(rank_V, dims_V, dspace_id, error)
       !dataset for primitive vars
       call h5dcreate_f(file_id, dset_prim, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write prim vars
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, V(BDRY_VAR,:,:,:), dims_V, error)
       !close dataset
       call h5dclose_f(dset_id, error)
       !close dataspace
       call h5sclose_f(dspace_id, error)
#endif       

       !dataspace for x-coordinates
       dims_XYZ = (/gr_glb_nx/)
       
       call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
       call h5dcreate_f(file_id, dset_x, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xCoord, dims_XYZ, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)

       !dataspace for y-coordinates
       dims_XYZ = gr_glb_ny
       call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
       call h5dcreate_f(file_id, dset_y, h5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, yCoord, dims_XYZ, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)

       !dataspace for z-coordinates
       dims_XYZ = gr_glb_nz
       call h5screate_simple_f(rank_XYZ, dims_XYZ, dspace_id, error)
       call h5dcreate_f(file_id, dset_z, h5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
       !write data
       call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, zCoord, dims_XYZ, error)

       call h5dclose_f(dset_id,error)
       call h5sclose_f(dspace_id,error)
       



       !close hdf5 file and interface
       call h5fclose_f(file_id,error)
       call h5close_f(error)

       deallocate(dims_V)
       deallocate(dims_XYZ)

      

    end if
    Sync All
    deallocate(img_V)
    deallocate(nBlock)

100 format()
910 format(i5)
920 format(1x,f16.8,1x,f16.8,1x,NUMB_VAR f32.16)
    
    close(20)

  end subroutine io_writeHDF5

  subroutine io_writeASCII(t,nstep,ioCounter)
    implicit none

    real, intent(IN) :: t
    integer, intent(IN) :: nstep,ioCounter
    
    integer :: i,j,nVar,nCell,dest
    character(len=50) :: ofile
    character(len=5)  :: cCounter

    integer :: ibeg, jbeg, iend, jend
    real, dimension(gr_glb_nx) :: xCoord
    real, dimension(gr_glb_ny) :: yCoord
    real, dimension(NUMB_VAR, gr_glb_nx, gr_glb_ny) :: V
    integer, dimension(NDIM)   :: img_size
    
    integer, allocatable :: nBlock(:)[:]
    real,    allocatable :: img_V(:,:,:)[:], loc_V(:,:,:)

    !make image buffers for grid

    !needs to be updated for 3D
!!$    allocate(img_V(NUMB_VAR, gr_nx, gr_ny)[*])
!!$    allocate(nBlock(NDIM)[*])
!!$    
!!$    nBlock(XDIM) = gr_nx
!!$    nBlock(YDIM) = gr_ny
!!$    img_V(:, :, :) = gr_V(:, gr_ibeg(XDIM):gr_iend(XDIM), gr_ibeg(YDIM):gr_iend(YDIM))
!!$
!!$    Sync All
!!$
!!$    if (bl_ID == 1) then
!!$       !I am root
!!$       do i = 1, gr_glb_nx
!!$          xCoord(i) = ( real(i)-0.5 )*gr_dx + gr_xbeg(XDIM)
!!$       end do
!!$       do j = 1, gr_glb_ny
!!$          yCoord(j) = (real(j)-0.5)*gr_dy + gr_xbeg(YDIM)
!!$       end do
!!$       do i = 1, bl_iProcs
!!$          do j = 1, bl_jProcs
!!$             dest = bl_grid(i,j)
!!$             img_size(:) = nBlock(:)[dest]
!!$
!!$             ibeg = (i-1)*bl_nBlock(XDIM) + 1
!!$             jbeg = (j-1)*bl_nBlock(YDIM) + 1
!!$             iend = ibeg -1 + img_size(XDIM)
!!$             jend = jbeg -1 + img_size(YDIM)
!!$             
!!$             allocate(loc_V(NUMB_VAR,img_size(XDIM),img_size(YDIM)))
!!$             loc_V(:,:,:) = img_V(:,:,:)[dest]
!!$             V(:, ibeg:iend, jbeg:jend) = loc_V(:,:,:)
!!$             deallocate(loc_V)
!!$          end do
!!$       end do
!!$       
!!$       ! convert conter number to character
!!$       write(cCounter,910) ioCounter + 10000
!!$
!!$       ! file name for ascii output
!!$       !ofile = 'slug_'//trim(sim_name)//'_'//cCounter//'.dat'
!!$       ofile = trim(sim_name)//'_'//cCounter//'.dat'
!!$       
!!$       open(unit=20,file=ofile,status='unknown')
!!$       do i=1, gr_glb_nx
!!$          do j = 1, gr_glb_ny
!!$             write(20,920)xCoord(i),ycoord(j),(V(nVar,i,j),nVar=1,EINT_VAR)
!!$          end do
!!$          write(20,100)
!!$       end do
!!$
!!$    end if
!!$    Sync All
!!$    deallocate(img_V)
!!$    deallocate(nBlock)
!!$
!!$100 format()
!!$910 format(i5)
!!$920 format(1x,f16.8,1x,f16.8,1x,NUMB_VAR f32.16)
!!$    
!!$    close(20)
  end subroutine io_writeASCII
  
end module io
