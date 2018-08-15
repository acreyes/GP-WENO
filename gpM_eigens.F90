subroutine gpM_eigens()

#include "definition.h"

  use gp_data

  implicit none

  real, dimension(gp_Npts,gp_Npts) :: C
  integer :: N, LDA, LWORK, INFO, i, j, k, m

  real, allocatable, dimension(:) :: W, WORK

  N = gp_Npts

  LDA = N
  Lwork = 66*N
  allocate(W(N))
  allocate(WORK(LWORK))

  do i = 1, N
     do j = 1, N
        C(i,j) = SE(REAL(gp_stencil(XDIM,i)),REAL(gp_stencil(XDIM,j)),gp_Xdel)*&
                 SE(REAL(gp_stencil(YDIM,i)),REAL(gp_stencil(YDIM,j)),gp_Ydel)
     end do
  end do

  call DSYEV('V', 'L', N, C, LDA, W, WORK, LWORK, INFO)

  if (INFO < 0) then
     print *, "error in dsyev, info != 1"
     stop
  end if
  !allocate eigensystem vars for GP
  allocate(gpM_evals(N  )); gpM_evals = W
  allocate(gpM_evecs(N,N)); gpM_evecs = C
  allocate(gpM_Pvecs(N,N))

  do i = 1, N
     gpM_Pvecs(:,i) = C(:,i)/sqrt(W(i))
  end do

contains
  function SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/eldel)**2 )
    return
  end function SE



  !deallocate(W)
  !deallocate(WORK)

end subroutine gpM_eigens
