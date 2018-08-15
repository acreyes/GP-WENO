module gp_data

#include "definition.h"

  implicit none

  !GP pars
  character(len=MAX_STRING_LENGTH), save :: gp_quad, gp_kernel
  real, save :: gp_el, gp_eldel, gpM_radius , gp_Xdel, gp_Ydel
  integer, save :: gp_radius, gp_Npts, gp_cntrPt

  !GP data
  real, allocatable, dimension(:), save :: gp_evals, gpM_evals
  real, allocatable, dimension(:,:), save :: gp_evecs, gp_Pvecs, gpM_Pvecs, gpM_evecs
  integer, allocatable, dimension(:,:), save :: gp_stencil, gp_1Dstencil


  !these will be truncated from quadruple precision
!  real(KIND=16), save :: gp_Xdel, gp_Ydel
  
  real(KIND=8), allocatable, dimension(:) :: gp_v
  real(KIND=8), allocatable, dimension(:,:) :: gp_w, gp_z , gp_vk
  real(KIND=8), allocatable, dimension(:,:,:) :: gp_zk, gp_Zvecs

  real(KIND=16), allocatable, dimension(:) :: gp4_v
  real(KIND=16), allocatable, dimension(:,:) :: gp4_w, gp4_z , gp4_vk
  real(KIND=16), allocatable, dimension(:,:,:) :: gp4_zk, gp4_Zvecs

  !multi D stencils
  real(KIND=8), allocatable, dimension(:) :: gpM_v
  real(KIND=8), allocatable, dimension(:,:) :: gpM_w, gpM_z

  real(KIND=16), allocatable, dimension(:) :: gpM4_v
  real(KIND=16), allocatable, dimension(:,:) :: gpM4_w, gpM4_z

  !flux interpolations
  real(KIND=8), allocatable, dimension(:) :: gpF_Z2, gpF_Z4



end module gp_data
