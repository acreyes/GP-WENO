subroutine sim_GPinit()
  !want to initialize all the GP parameters that do not depend on the data, and only on the grid geometries.
  !I am using sig/del = 2
  !This includes calculating the following:
  !****** the covariance matrix C 
  !****** the cross-correlation matrix T 
  !****** The corresponding vectors and matrices to be used in actual calculations (see eqs 30-32) v & Z
  !******  ****** C.v  = u  (gr_GPv)
  !******  ****** C.Z* = T* (gr_GPZ)
#include "definition.h"

  use grid_data
  use linalg
  use GP

  implicit none

  real, dimension(gr_GP_StencilPts,gr_GP_StencilPts) :: C, L, W, Vt, Sigmai
  real, dimension(2,gr_GP_StencilPts) :: T
  real, dimension(gr_GP_StencilPts,2) :: B
  real, dimension(gr_GP_StencilPts)   ::  u, S
  integer :: i,j,N
  real :: small
  integer :: INFO
  integer, dimension(gr_GP_StencilPts) :: IPIV

  !initialize
  N = gr_GP_StencilPts
  C = 0.
  T = 0.
  u = 1.
  
  

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1, gr_GP_StencilPts
     do j = 1, gr_GP_StencilPts
        C(i,j) = intg_kernel(gr_GP_Stencil(i,:), gr_GP_Stencil(j,:) )
     end do
     T(1:2,i) = intg_predvec(gr_GP_Stencil(i,:))
  end do

!!$  do i = 1,5
!!$     do j = 1,i-1
!!$        !off diaganol terms
!!$        C(i,j) = intg_kernel(stencil(i),stencil(j))
!!$        C(j,i) = C(i,j)
!!$     end do
!!$     C(i,i) = intg_kernel(stencil(i),stencil(i))
!!$     !whiel we're here lets take care of T
!!$     T(1:2, i) = intg_predvec(stencil(i))
!!$  end do
  

  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
!!$  call chol(C, gr_GP_StencilPts, L)
!!$  call solve_Axb(C, gr_GPv, u, L, gr_GP_StencilPts)
!!$  call solve_CZT(C, gr_GPZ, T, L, gr_GP_StencilPts)

  !LU factorization method
  L = C
  call dgetrf(N, N, C, N, IPIV, INFO)
  gr_GPv = 1.
  call dgetrs('N', N, 1, C, N, IPIV, gr_GPv, N, INFO)
  gr_GPZ = T
  call dgetrs('N', N, 1, C, N, IPIV, T(1,:), N, INFO)
  call dgetrs('N', N, 1, C, N, IPIV, T(2,:), N, INFO)
  gr_GPZ = T
end subroutine sim_GPinit
