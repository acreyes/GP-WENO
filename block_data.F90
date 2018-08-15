module block_data

#include "definition.h"

  implicit none

  integer,                              save :: bl_iProcs, bl_jProcs,bl_kProcs, bl_ID, bl_nprocs
  integer,                              save :: bl_i, bl_j, bl_k
  integer, allocatable, dimension(:  ), save :: bl_nBlock, bl_BC, bl_cornerBC
  integer, allocatable, dimension(:,:,:), save :: bl_grid

  real, allocatable, dimension(:,:,:,:), codimension[:], save :: bl_buffL , bl_buffR , bl_buffT , bl_buffB, &
                                                                 bl_buffU, bl_buffD
  real, allocatable, dimension(:,:,:,:), codimension[:], save :: bl_buffBL, bl_buffTL, bl_buffTR, bl_buffBR
  real, allocatable, dimension(:,:,:,:,:), codimension[:], save :: bl_bufCL, bl_bufCR, bl_bufCT, bl_bufCB

  real, allocatable, codimension[:] :: bl_delT, bl_maxSpeed
!  real, allocatable, dimension(:,:,:), codimension[:] :: bl_V


end module block_data
  
