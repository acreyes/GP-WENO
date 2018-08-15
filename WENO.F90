module WENO

#include "definition.h"

contains

  subroutine gpM_beta(V, Npts, beta)
    use gp_data, only: gpM_Pvecs
    implicit none

    integer              , intent(IN   ) :: Npts
    real, dimension(Npts), intent(IN   ) :: V
    real                 , intent(INOUT) :: beta

    integer :: i

    beta = 0.
    do i = 1, Npts
       beta = beta + dot_product(V(:), gpM_Pvecs(:,i))**2
    end do

    return
  end subroutine gpM_beta

  subroutine gp_betas(V, R, beta)
    use gp_data , only:  gp_Pvecs, gp_Zvecs
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer i, N, s
    !real, dimension(R+1) :: gamms, kxp, f, gamms2
    !real, dimension(R+1) :: weno_beta
    !real :: gamm

    N = R+1
    !loop over stencils
    do s = 1 ,N
       !loop over eigen values
       beta(s) = dot_product(V(s:s+R), gp_Pvecs(:,1))**2
       do i = 2, N
          beta(s) = beta(s) + dot_product(V(s:s+R), gp_Pvecs(:,i))**2
       end do
    end do

    return
  end subroutine gp_betas



  subroutine betas(V, R, beta)
    !subroutine to calculate the smoothness-indicators for a WENO scheme on a 2R+1 point stencil
    implicit none

    integer,                   intent(IN   ) :: R
    real   , dimension(2*R+1), intent(IN   ) :: V
    real   , dimension(  R+1), intent(INOUT) :: beta

    integer :: i

    i = R+1
    !print *, 'betas'

    select case(R)
    case(1)
       beta(1) = (V(i  )-V(i-1))**2
       beta(2) = (V(i+1)-V(i  ))**2
    case(2)
       beta(1) = 13./12.*(V(i-2) - 2.*V(i-1) + V(i  ) )**2 + 0.25*(   V(i-2) - 4.*V(i-1) + 3.*V(i  ) )**2
       beta(2) = 13./12.*(V(i-1) - 2.*V(i  ) + V(i+1) )**2 + 0.25*(   V(i-1)             -    V(i+1) )**2
       beta(3) = 13./12.*(V(i  ) - 2.*V(i+1) + V(i+2) )**2 + 0.25*(3.*V(i  ) - 4.*V(i+1) +    V(i+2) )**2
    case(3)
       beta(1) = V(i-3)*(  547.*V(i-3) -  3882.*V(i-2) + 4642.*V(i-1) - 1854.*V(i  )) + &
                 V(i-2)*( 7043.*V(i-2) - 17246.*V(i-1) + 7042.*V(i  )               ) + &
                 V(i-1)*(11003.*V(i-1) -  9402.*V(i  )                              ) + 2107.*V(i  )**2
       beta(2) = V(i-2)*(  267.*V(i-2) -  1642.*V(i-1) + 1602.*V(i  ) -  494.*V(i+1)) + &
                 V(i-1)*( 2843.*V(i-1) -  5966.*V(i  ) + 1922.*V(i+1)               ) + &
                 V(i  )*( 3443.*V(i  ) -  2522.*V(i+1)                              ) +  547.*V(i+1)**2
       beta(3) = V(i-1)*(  547.*V(i-1) -  2522.*V(i  ) + 1922.*V(i+1) -  494.*V(i+2)) + &
                 V(i  )*( 3443.*V(i  ) -  5966.*V(i+1) + 1602.*V(i+2)               ) + &
                 V(i+1)*( 2843.*V(i+1) -  1642.*V(i+2)                              ) +  267.*V(i+2)**2
       beta(4) = V(i  )*( 2107.*V(i  ) -  9402.*V(i+1) + 7042.*V(i+2) - 1854.*V(i+3)) + &
                 V(i+1)*(11003.*V(i+1) - 17246.*V(i+2) + 4642.*V(i+3)               ) + &
                 V(i+2)*( 7043.*V(i+2) -  3882.*V(i+3)                              ) +  547.*V(i+3)**2
    case DEFAULT
       beta = 1.
    end select

    return
  end subroutine betas

end module WENO
