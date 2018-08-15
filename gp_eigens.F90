subroutine gp_eigens()

#include "definition.h"

  use gp_data

  implicit none

  real, allocatable, dimension(:,:) :: C
  integer :: N, LDA, LWORK, INFO, i, j, k, m

  real, allocatable, dimension(:) :: W, WORK

  N = gp_radius+1

  LDA = N
  LWORK = 66*N
  allocate(C(N,N))
  allocate(W(N))
  allocate(WORK(LWORK))


  do i = 1, N
     do j = 1, N
        C(i,j) = SE(REAL(i),REAL(j))
     end do
  end do

  call DSYEV('V', 'L', N, C, LDA, W, WORK, LWORK, INFO)

  if (INFO < 0) then
     print *, "error in dsyev, info != 1"
     stop
  end if
  !allocate eigensystem vars for GP
  allocate(gp_evals(N  )); gp_evals = W
  allocate(gp_evecs(N,N)); gp_evecs = C
  allocate(gp_Pvecs(N,N))

  
!!$  do i = 1, N
!!$     do j = 1, N
!!$        gp_Pvecs(j,i) = dot_product(C(:,i), gp_Zvecs(j,:,1))/sqrt(W(i))
!!$     end do
!!$  end do

  do i = 1, N
     gp_Pvecs(:,i) = C(:,i)/sqrt(W(i))
  end do

contains
  function SE(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/gp_eldel)**2 )
    return
  end function SE

end subroutine gp_eigens
