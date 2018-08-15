subroutine gp_MakeStencil()
#include "definition.h"

  use gp_data

  implicit none

  integer, allocatable, dimension(:,:) :: stencil

  real, dimension(5) :: weno5_pts

  integer :: i, j, Ncntr, s, scntr, pt
  real :: R2, x, y

  !round gp_R to ceiling
  if (gpM_radius == 0.0) then
     gpM_radius = REAL(gp_radius)
  else
     gp_radius = CEILING(gpM_radius)
  end if

  allocate(stencil(2,(2*gp_radius+1)**2))
  R2 =  gpM_radius*gpM_radius
  Ncntr = 0
  do i = -gp_radius, gp_radius
     do j = -gp_radius, gp_radius
        x = real(i)
        y = real(j)
        if (x**2 + y**2 .le. R2) then
           Ncntr = Ncntr + 1
           stencil(XDIM,Ncntr) = i
           stencil(YDIM,Ncntr) = j
           if (i == 0 .and. j == 0) then
              gp_cntrPt = Ncntr
           end if
        end if
     end do
  end do
  gp_Npts = Ncntr
  allocate(gp_stencil(NDIM,gp_Npts))
  gp_stencil(:,:) = stencil(:,1:gp_Npts)
  deallocate(stencil)


  !WENO things
  !need to figure out the 1D stencil coordinates on the multiD stencil
  weno5_pts = (/-2., -1., 0., 1., 2.  /)
  allocate(gp_1Dstencil(5,NDIM)); gp_1Dstencil = 0

  do i = 1, 5 !loop over 1D stencil points
     pt = weno5_pts(i)
     !loop over gp stencil
     do j = 1, gp_Npts
        if (gp_stencil(XDIM,j)==pt .and. gp_stencil(YDIM,j) == 0) then
           gp_1Dstencil(i,XDIM) = j
        end if
        if (gp_stencil(YDIM,j)==pt .and. gp_stencil(XDIM,j) == 0) then
           gp_1Dstencil(i,YDIM) = j
        end if
     end do
  end do

end subroutine gp_MakeStencil
