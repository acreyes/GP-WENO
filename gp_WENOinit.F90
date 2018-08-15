subroutine gp_WENOinit()
  !want to initialize all the GP parameters that do not depend on the data, and only on the grid geometries.
  !I am using sig/del = 2
  !This includes calculating the following:
  !****** the covariance matrix C 
  !****** the cross-correlation matrix T 
  !****** The corresponding vectors and matrices to be used in actual calculations (see eqs 30-32) v & Z
  !******  ****** C.v  = u  (gr_GPv)
  !******  ****** C.Z* = T* (gr_GPZ)
#include "definition.h"

  use grid_data, only: gr_dx
  use gp_data
  use linalg
  use GP

  implicit none


  real, dimension(2*gp_radius+1,2*gp_radius+1) :: C, L
  real, dimension(2,2*gp_radius+1) :: T, Z
  real, dimension(2*gp_radius+1)   :: stencil, u

  real, dimension(gp_radius+1,gp_radius+1) :: Ck, Lk
  real, dimension(2,gp_radius+1) :: Tk, w
  real, dimension(2,gp_radius+1,gp_radius+1) :: Zk

  real, dimension(2*gp_radius+2, gp_radius+1) :: Zmat
  real, dimension(2*gp_radius+2) :: Zvec
  real, dimension(gp_radius+1) :: ul
  real, dimension(2*gp_radius+1) :: un

  integer :: i,j,N,R,m,LR,COL,ROW
  real :: Xdel



  !initialize
  R = gp_radius
  N = 2*R+1
  C = 0.
  T = 0.
  u = 1.

  if (gp_el == 0.) then
     gp_el = gr_dx*gp_eldel
  elseif(gp_eldel == 0.) then
     gp_eldel = gp_el/gr_dx
  end if

  Xdel = gp_el/gr_dx
  gp_Xdel = Xdel
  
  do i = 1, N
     stencil(i) = REAL(i-gp_radius-1)
  end do

  !first thing is to calculate the covariance matrix according to eq. 15
  !since C is symmetric only bother with one side of the diaganol
  do i = 1,N
     do j = 1,N
        C(i,j) = SE(stencil(i),stencil(j),Xdel)
     end do
     T(1, i) = SE(-0.5, stencil(i), Xdel)
     T(2, i) = SE( 0.5, stencil(i), Xdel)
  end do

  !now we need to solve the linear eqns for v & Z (see eqs 30-32)
  call chol(C, N, L)
  call solve_CZT(C, Z, T, L, N)
  
  !now we take care of the ENO stencils
  N = R+1
  do m = 1, N
     do i = 1, N
        stencil(i) = REAL(i-1-R+m-1)
     end do
     do i = 1, N
        do j = 1, N
           Ck(i,j) = SE(stencil(i),stencil(j),Xdel)
        end do
        Tk(1,i) = SE(-0.5, stencil(i), Xdel)
        Tk(2,i) = SE( 0.5, stencil(i), Xdel)
     end do
     call chol(Ck,N,Lk)
     call solve_CZT(Ck, Zk(:,:,m), Tk, Lk, N)
     !no need for Z vectors as in FVM, since we are interpolating ptwise values
  end do

  

  !now we solve for GP-WENO linear weights
  ROW = 2*R+1
  COL = R+1
  ul = 1.
  un = 1.

  do LR = 1, 2
     Zmat = 0.
     do m = 1, COL
        Zmat(m:m+R,m) = Zk(LR,:,m)
        Zmat(2*R+2,m) = dot_product(Zk(LR,:,m), ul)
     end do

     Zvec(1:ROW) = Z(LR,:)
     Zvec(2*R+2) = dot_product(Z(LR,:), un)

     call LSTSQ(ROW+1,COL,Zmat, w(LR,:), Zvec)
  end do

  !allocate gp vars and truncate to double precision
  allocate(gp_Z (2,2*R+1  )); gp_Z  = Z
  allocate(gp_Zk(2,R+1,R+1)); gp_Zk = Zk
  allocate(gp_w (2,R+1    )); gp_w  = w


end subroutine gp_WENOinit
