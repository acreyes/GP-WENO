subroutine gr_GPinit()
  !This subroutine is to initialize all of the grid related GP parameters
  !this includes:
  !****gr_gp_stencil
  !****gr_gp_stencilPts
  !
  !We will accomplish this by starting with a box that is 2rx2r (r=gr_radius)
  !and finding the stencil points that fall within the stencil radius for GP

#include "definition.h"

  use grid_data

  implicit none

  real, allocatable, dimension(:,:) :: stencil
  integer :: i,j, R, stencilPts
  real    :: x,y

  !allocate the initial stencil as the max box that circumscribes the actual stencil
  R = FLOOR(gr_radius)
  gr_Tcells = R
  allocate(stencil((2*R+1)**2,2)); stencil = 0. !(stencilPt, x, y)

  stencilPts = 0
  do i = -R,R
     do j = -R,R
        x = REAL(i)
        y = REAL(j)
        if (x**2 + y**2 .le. gr_radius**2) then
           stencilPts = StencilPts + 1
           stencil(StencilPts, :) = (/ x, y /)
        end if
     end do
  end do


  !1D stencil:
!!$  gr_Tcells = 0
!!$  StencilPts = 2*R+1
!!$  do i = -R, R
!!$     y = 0.
!!$     x = REAL(i)
!!$     stencil(i+R+1,:) = (/x, y /)
!!$end do

  
  
  !now we know how many points are in the stencil and can allocate all of our gp grid vars
  gr_GP_stencilPts = stencilPts
  allocate(gr_GP_stencil(StencilPts, 2))
  gr_GP_stencil(:,:) = stencil(1:StencilPts,:)
  
!!$  gr_GP_stencil = 0.
!!$  do i = -2, 2
!!$     gr_GP_stencil(i+3, 1) = REAL(i)
!!$  end do

  allocate(gr_GPv(StencilPts   )); gr_GPv = 0.
  allocate(gr_GPZ(2, StencilPts)); gr_GPZ = 0.
  
end subroutine gr_GPinit
