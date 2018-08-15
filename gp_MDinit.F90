subroutine gp_MDinit()

#include "definition.h"

  use linalg
  use GP
  use gp_data
  use grid_data, only: gr_dx, gr_dy

  implicit none

  real, dimension(gp_Npts, gp_Npts) :: C, L, W, C2
  real, dimension(2,gp_Npts) :: Tx, Ty, Zx, Zy
  real, dimension(gp_Npts) :: u
  real, dimension(NDIM, gp_Npts) :: stencil

  integer :: N, Nx, Ny, i, j, icntr, i0, j0
  real    :: Xdel, Ydel



  !initialize vars
  C = 0.
  !T = 0.
  u = 1.

  !3x3 stencil = 9 sample points
  Nx = 3
  i0 = 2
  Ny = 3
  j0 = 2
  N = gp_Npts

  

  !take care of \el
  if (gp_el == 0. .and. gp_eldel == 0.) then
     print *, 'both \el and \el / \Delta are zero. exiting'
     stop
  elseif (gp_eldel == 0.) then
     gp_Xdel = gp_el/gr_dx
     gp_Ydel = gp_el/gr_dy
  elseif (gp_el == 0.) then
     gp_el = gp_eldel*gr_dx
     gp_Xdel = gp_eldel
     gp_Ydel = gp_eldel
  end if
  Xdel = gp_Xdel
  Ydel = gp_Ydel
  !double precision vars
  allocate(gpM_z(4,N))
  allocate(gpM_v(N))
  !quadruple precision vars
  allocate(gpM4_z(4,N))
  allocate(gpM4_v(N))
  
  !copy gp_stencil to local stencil var so that it is quad. prec
  do i = 1, N
     stencil(:,i) = REAL(gp_stencil(:,i))
  end do

  !compute covariance matrix and prediction vectors
  do i = 1, N
     do j = 1, N
        C(i,j) = SE(stencil(XDIM,i), stencil(XDIM,j), Xdel)* &
                 SE(stencil(YDIM,i), stencil(YDIM,j), Ydel)
     end do
     Tx(1,i) = SE(-0.5, stencil(XDIM,i), Xdel)*SE(stencil(YDIM,i), 0., Ydel)
     Ty(1,i) = SE(-0.5, stencil(YDIM,i), Ydel)*SE(stencil(XDIM,i), 0., Xdel)
     Tx(2,i) = SE( 0.5, stencil(XDIM,i), Xdel)*SE(stencil(YDIM,i), 0., Ydel)
     Ty(2,i) = SE( 0.5, stencil(YDIM,i), Ydel)*SE(stencil(XDIM,i), 0., Xdel)
  end do


  !solve for weights
  call chol(C, N, L)
  call solve_Axb(C, gpM4_v, u, L, N)
  call solve_CZT(C, Zx, Tx, L, N)
  call solve_CZT(C, Zy, Ty, L, N)

  !truncate to double precision from quad
  gpM_z(1:2,:) = Zx
  gpM_z(3:4,:) = Zy

end subroutine gp_MDinit
