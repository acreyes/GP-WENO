module linalg

contains

  subroutine chol(A, N, L)
    !Calculate the cholesky decomposition (L) of matrix A that is NxN
    implicit none
    integer             , intent(IN)  :: N
    real, dimension(N,N), intent(IN)  :: A

    real, dimension(N,N), intent(OUT) :: L
    real, dimension(N,N) :: R
    real, dimension(N) :: D

    real :: small
    integer :: i, j, k

    small = 0.
    
    !initialize L
    D = 0.
    L = 0.
    R = 0.

    !!!!!LDL decomposition !!!!!!!!!
    !D(1) = A(1,1)
    do i = 1,N
       do j = 1,i
          if (D(j) == 0.) then
             D(j) = A(j,j)+small
             do k = 1, j-1
                D(j) = D(j) - D(k)*R(j,k)**2
             end do
          end if
          R(i,j) = A(i,j)
          do k = 1, j-1
             R(i,j) = R(i,j) - R(i,k)*R(j,k)*D(k)
          end do
          R(i,j) = R(i,j)/D(j)
       end do
    end do

    do i = 1, N
       do j = 1,i
          L(i,j) = R(i,j)*SQRT(D(j))
       end do
    end do

    
    !stop
  end subroutine chol

  subroutine LSTSQ(M, N, A, x, b)
    !compute the least squares solution to the overdetermined system (ONLY!) Ax=b
    !using the QR factorization of A
    !done using Householder reflections
    implicit none
    integer, intent(IN) :: M, N
    real, dimension(M,N), intent(INOUT) :: A
    real, dimension(M)  , intent(INOUT) :: b

    real, dimension(N)  , intent(INOUT) :: x

    real, dimension(M) :: QTb
    real, dimension(M,1) :: v
    real, dimension(N,N) :: R

    integer :: i, j
    real    :: v2, gam

    R = 0.
    QTb = 0.
    V = 0.
    !we're going to compute the relevant items from the QR factorization
    !R is going to be stored in A
    !we are also going to calculate Q^T*b for the least squares solve
    do i = 1, N
       v = 0.
       v(i:M,1) = A(i:M,i)
       !vi = xi + sign(xi)*x^2
       v(i,1  ) = A(i,i) + SIGN( SQRT(dot_product( A(i:M,i),A(i:M,i) )), A(i,i) )

       v2 = SQRT(dot_product(v(i:M,1),v(i:M,1)))
       v(:,1) = v(:,1)/v2

       
       A(i:M,i:N) = A(i:M,i:N) - 2.*MATMUL( v(i:M,:), MATMUL( TRANSPOSE(v(i:M,:)), A(i:M,i:N) ) )

       
       gam = -2.*dot_product(v(i:M,1), b(i:M))
       b(i:M) = b(i:M) + gam*v(i:M,1)
       
    end do

    !now we just to do a back-substitution to solve Rx=Q^T*b
    do i = N, 1, -1
       x(i) = b(i)
       do j = i+1, N
          x(i) = x(i) - A(i,j)*x(j)
       end do
       x(i) = x(i)/A(i,i)
    end do


  end subroutine LSTSQ

  subroutine SVD(A, U, S, Vt, N)
    !gives the singular value decomposition of A (NxN) = U.S.V*
    !uses the DGESVD subroutine in LAPACK
    !on mac this is included in the accelerate framework (LDFLAG -framework accelerate)
    implicit none

    integer, intent(IN) :: N
    real, dimension(N,N), intent(IN) :: A
    real, dimension(N,N), intent(INOUT) :: U, Vt
    real, dimension(N), intent(INOUT) :: S

    integer :: info, lwork
    !integer ( kind = 4 ), parameter :: lwork 
    real ( kind = 8 ) work(3*N**2)
    lwork = 3*N**2
    call DGESVD('A', 'A', N, N, A, N, S, U, N, Vt, N, work, lwork, info)
    print *, info
    !print *, work
   
    return
  end subroutine SVD

  subroutine multUSV(U, S, V, C, N)
    implicit none
    integer, intent(IN) :: N
    real, dimension(N,N), intent(IN) :: U,V
    real, dimension(N,N), intent(INOUT) :: C
    real, dimension(N), intent(IN) :: s
    real, dimension(N,N) :: D, E

    integer :: i, j

    D = 0.
    do i = 1,N
       D(i,i) = S(i)
    end do
    
    do i = 1,N
       do j = 1,N
          E(i,j) = dot_product(U(i,:), D(:,j))
       end do
    end do

    do i = 1,N
       do j = 1,N
          C(i,j) = dot_product(E(i,:), V(i,:))
       end do
    end do
    return
    
  end subroutine multUSV

  subroutine solve_Axb(A, x, b, L, N)
    !solve matrix equation Ax = b for vector x, A is NxN given the cholesky decomposition L of A
    implicit none
    integer, intent(IN) :: N
    real, dimension(N  ), intent(IN) :: b
    real, dimension(N,N), intent(IN) :: A, L

    real, dimension(N) , intent(OUT) :: x

    real, dimension(N)   :: y
    integer :: i, j

    !initialize
    x = 0.
    y = 0.
    
    !first we need to solve Ly=b using forward sub
    y(1) = b(1)/L(1,1)
    do i = 2,N
       y(i) = b(i)
       do j = 1, i-1
          y(i) = y(i) - L(i,j)*y(j)
       end do
       y(i) = y(i)/L(i,i)
    end do

    !now solve L*x = y
    x(N) = y(N)/L(N,N)
    do i = N-1, 1, -1
       x(i) = y(i)
       do j =i+1, N
          x(i) = x(i) - L(j,i)*x(j)
       end do
       x(i) = x(i)/L(i,i)
    end do

  end subroutine solve_Axb

  subroutine solve_CZT(C, Z, T, L, N)
    !subroutine to solve the matrix equation C.Z* = T* for the matrix Z
    !note that the matrix equation is for the transpose of Z and T
    !C shoul be NxN, while Z and T are 2xN
    implicit none
    integer,                 intent(IN ) :: N
    real   , dimension(N,N), intent(IN ) :: C, L
    real   , dimension(2,N), intent(IN ) :: T

    real, dimension(2,N), intent(OUT) :: Z

    !solve for the first column of Z*
    call solve_Axb(C, Z(1, 1:N), T(1, 1:N), L, N)
    !solve fore the second column of Z*
    call solve_Axb(C, Z(2, 1:N), T(2, 1:N), L, N)
  end subroutine solve_CZT





end module linalg
