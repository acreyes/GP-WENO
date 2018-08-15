module eigensystem

#include "definition.h"
  
  use grid_data

contains

  subroutine eigenvalues(V,lambda,dir)
    implicit none

    real, dimension(NUMB_VAR), intent(IN)  :: V
    real, dimension(NUMB_WAVE),intent(OUT) :: lambda
    integer, intent(IN) :: dir
    
    real :: a, u
    integer :: VEL

    if     (dir == XDIM) then
       VEL = VELX_VAR
    elseif (dir == YDIM) then
       VEL = VELY_VAR
    else
       VEL = VELZ_VAR
    end if
    
    ! sound speed
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))!;print*,a,V(GAMC_VAR),V(PRES_VAR),V(DENS_VAR)
    u = V(VEL)
!!$    if (dir == 2) then
!!$       print *, a, u, u-a, u+a
!!$    end if
    
    lambda(SHOCKLEFT) = u - a
    lambda(SLOWWLEFT) = u
    lambda(CTENTROPY) = u
    lambda(SLOWWRGHT) = u
    lambda(SHOCKRGHT) = u + a
    
    return
  end subroutine eigenvalues

  subroutine eigenvectors(V,conservative,reig,leig,dir)
    implicit none
    integer,                                intent(IN ) :: dir
    real   , dimension(NUMB_VAR),           intent(IN ) :: V
    logical,                                intent(IN ) :: conservative
    real   , dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig, leig

    !some convenient parameters
    real :: d, P, a, ahinv, dinv, Na, a2inv, v1, v2, v3,  H, ekin, g, dot
    integer :: VEL1_VAR, VEL2_VAR, VEL3_VAR, varL

    VEL1_VAR = 1; VEL2_VAR = 2
    if (dir == XDIM) then
       VEL1_VAR = VELX_VAR
       VEL2_VAR = VELY_VAR
       VEL3_VAR = VELZ_VAR
    elseif (dir == YDIM) then
       VEL1_VAR = VELY_VAR
       VEL2_VAR = VELZ_VAR
       VEL3_VAR = VELX_VAR
    elseif (dir == ZDIM) then
       VEL1_VAR = VELZ_VAR
       VEL2_VAR = VELX_VAR
       VEL3_VAR = VELY_VAR
    else
       print *, 'invalid dir to eigenvectors, stopping'
       stop
    end if
    
    d = V(DENS_VAR)
    P = V(PRES_VAR)
    a = SQRT(V(GAMC_VAR)*P/d)
    ahinv = 0.5/a
    dinv = 1./d
    Na = 0.5/(a*a)
    a2inv = 1./(a*a)
    v1 = V(VEL1_VAR)
    v2 = V(VEL2_VAR)
    v3 = V(VEL3_VAR)
    ekin = 0.5*(v1*v1+v2*v2+v3*v3)
    g = V(GAME_VAR) - 1.
    H = ekin+(P/g+P)*dinv


    if (conservative) then
       !eigenvectors for conservative variables

       !right eigenvectors first
       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VEL1_VAR,SHOCKLEFT) = v1 - a
       reig(VEL2_VAR,SHOCKLEFT) = v2
       reig(VEL3_var,SHOCKLEFT) = v3
       reig(PRES_VAR,SHOCKLEFT) = H - a*v1

       reig(DENS_VAR,SLOWWLEFT) = 0.
       reig(VEL1_VAR,SLOWWLEFT) = 0.
       reig(VEL2_VAR,SLOWWLEFT) = 1.
       reig(VEL3_VAR,SLOWWLEFT) = 0.
       reig(PRES_VAR,SLOWWLEFT) = v2

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VEL1_VAR,CTENTROPY) = v1
       reig(VEL2_VAR,CTENTROPY) = v2
       reig(VEL3_VAR,CTENTROPY) = v3
       reig(PRES_VAR,CTENTROPY) = ekin

       reig(DENS_VAR,SLOWWRGHT) = 0.
       reig(VEL1_VAR,SLOWWRGHT) = 0.
       reig(VEL2_VAR,SLOWWRGHT) = 0.
       reig(VEL3_VAR,SLOWWRGHT) = 1.
       reig(PRES_VAR,SLOWWRGHT) = v3

       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VEL1_VAR,SHOCKRGHT) = v1 + a
       reig(VEL2_VAR,SHOCKRGHT) = v2
       reig(VEL3_VAR,SHOCKRGHT) = v3
       reig(PRES_VAR,SHOCKRGHT) = H + v1*a



       leig(DENS_VAR,SHOCKLEFT) = g*ekin+v1*a
       leig(VEL1_VAR,SHOCKLEFT) = -g*v1-a
       leig(VEL2_VAR,SHOCKLEFT) = -g*v2
       leig(VEL3_VAR,SHOCKLEFT) = -g*v3
       leig(PRES_VAR,SHOCKLEFT) = g

       leig(DENS_VAR:PRES_VAR,SHOCKLEFT) = Na*leig(DENS_VAR:PRES_VAR,SHOCKLEFT)

       leig(DENS_VAR,SLOWWLEFT) = -v2
       leig(VEL1_VAR,SLOWWLEFT) = 0.
       leig(VEL2_VAR,SLOWWLEFT) = 1.
       leig(VEL3_VAR,SLOWWLEFT) = 0.
       leig(PRES_VAR,SLOWWLEFT) = 0.

       leig(DENS_VAR,CTENTROPY) = 1.-Na*g*2.*ekin
       leig(VEL1_VAR,CTENTROPY) = g*v1*a2inv
       leig(VEL2_VAR,CTENTROPY) = g*v2*a2inv
       leig(VEL3_VAR,CTENTROPY) = g*v3*a2inv
       leig(PRES_VAR,CTENTROPY) = -g*a2inv

       leig(DENS_VAR,SLOWWRGHT) = -v3
       leig(VEL1_VAR,SLOWWRGHT) = 0.
       leig(VEL2_VAR,SLOWWRGHT) = 0.
       leig(VEL3_VAR,SLOWWRGHT) = 1.
       leig(PRES_VAR,SLOWWRGHT) = 0.
       
       leig(DENS_VAR,SHOCKRGHT) = g*ekin-v1*a
       leig(VEL1_VAR,SHOCKRGHT) = -g*v1+a
       leig(VEL2_VAR,SHOCKRGHT) = -g*v2
       leig(VEL3_VAR,SHOCKRGHT) = -g*v3
       leig(PRES_VAR,SHOCKRGHT) = g

       leig(DENS_VAR:PRES_VAR,SHOCKRGHT) = Na*leig(DENS_VAR:PRES_VAR,SHOCKRGHT)

    else
       !primitive eigenvectors

       reig(DENS_VAR,SHOCKLEFT) = d
       reig(VEL1_VAR,SHOCKLEFT) = -a
       reig(VEL2_VAR,SHOCKLEFT) = 0.
       reig(VEL3_VAR,SHOCKLEFT) = 0.
       reig(PRES_VAR,SHOCKLEFT) = d*a*a

       reig(DENS_VAR,SLOWWLEFT) = 0.
       reig(VEL1_VAR,SLOWWLEFT) = 0.
       reig(VEL2_VAR,SLOWWLEFT) = 1.
       reig(VEL3_VAR,SLOWWLEFT) = 0.
       reig(PRES_VAR,SLOWWLEFT) = 0.

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VEL1_VAR,CTENTROPY) = 0.
       reig(VEL2_VAR,CTENTROPY) = 0.
       reig(VEL3_VAR,CTENTROPY) = 0.
       reig(PRES_VAR,CTENTROPY) = 0.

       reig(DENS_VAR,SLOWWRGHT) = 0.
       reig(VEL1_VAR,SLOWWRGHT) = 0.
       reig(VEL2_VAR,SLOWWRGHT) = 0.
       reig(VEL3_VAR,SLOWWRGHT) = 1.
       reig(PRES_VAR,SLOWWRGHT) = 0.

       reig(DENS_VAR,SHOCKRGHT) = d
       reig(VEL1_VAR,SHOCKRGHT) = a
       reig(VEL2_VAR,SHOCKRGHT) = 0.
       reig(VEL3_VAR,SHOCKRGHT) = 0.
       reig(PRES_VAR,SHOCKRGHT) = d*a*a


       leig(DENS_VAR,SHOCKLEFT) = 0.
       leig(VEL1_VAR,SHOCKLEFT) = -ahinv
       leig(VEL2_VAR,SHOCKLEFT) = 0.
       leig(VEL3_VAR,SHOCKLEFT) = 0.
       leig(PRES_VAR,SHOCKLEFT) = dinv*Na

       leig(DENS_VAR,SLOWWLEFT) = 0.
       leig(VEL1_VAR,SLOWWLEFT) = 0.
       leig(VEL2_VAR,SLOWWLEFT) = 1.
       leig(VEL3_VAR,SLOWWLEFT) = 0.
       leig(PRES_VAR,SLOWWLEFT) = 0.

       leig(DENS_VAR,CTENTROPY) = 1.
       leig(VEL1_VAR,CTENTROPY) = 0.
       leig(VEL2_VAR,CTENTROPY) = 0.
       leig(VEL3_VAR,CTENTROPY) = 0.
       leig(PRES_VAR,CTENTROPY) = -a2inv

       leig(DENS_VAR,SLOWWRGHT) = 0.
       leig(VEL1_VAR,SLOWWRGHT) = 0.
       leig(VEL2_VAR,SLOWWRGHT) = 0.
       leig(VEL3_VAR,SLOWWRGHT) = 1.
       leig(PRES_VAR,SLOWWRGHT) = 0.
       
       leig(DENS_VAR,SHOCKRGHT) = 0.
       leig(VEL1_VAR,SHOCKRGHT) = ahinv
       leig(VEL2_VAR,SHOCKRGHT) = 0.
       leig(VEL3_VAR,SHOCKRGHT) = 0.
       leig(PRES_VAR,SHOCKRGHT) = dinv*Na

    end if

    do varL = 1, NUMB_WAVE
       dot = dot_product(leig(:,varL),reig(:,varL))
       if (abs(dot-1.) > 1.e-4) then
          print *, 'not unity', dot, varL
       end if
    end do

    do varL = 1, NUMB_WAVE
       if (varL < NUMB_WAVE) then
          dot = dot_product(leig(:,varL),reig(:,varL+1))
          if (abs(dot) > 1.e-4) then
             print *, 'not 0', dot, varL, varL+1
          end if
       else
          dot = dot_product(leig(:,varL),reig(:,varL-1))
          if (abs(dot) > 1.e-4 ) then
             print *, 'not 0',dot, varL, varL-1
          end if
       end if
    end do

  end subroutine eigenvectors
  
  subroutine right_eigenvectors(V,conservative,reig, dir)
    implicit none
    integer, intent(IN) :: dir
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: reig

    real :: a, d, ekin, hdai, H, e, p, v1, v2,g
    integer :: VEL1_VAR, VEL2_VAR

    VEL1_VAR = 1; VEL2_VAR = 2
    if (dir == XDIM) then
       VEL1_VAR = VELX_VAR
       VEL2_VAR = VELY_VAR
    elseif (dir == YDIM) then
       VEL1_VAR = VELY_VAR
       VEL2_VAR = VELX_VAR
    end if
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    g = V(GAMC_VAR) - 1.
    v1 = V(VEL1_VAR)
    v2 = V(VEL2_VAR)
    d = V(DENS_VAR)
    p = V(PRES_VAR)
    ekin = 0.5*(v1**2+v2**2)
    H = ekin + a**2/g
    hdai = 0.5/(a*d)
    
    if (conservative) then
       !! Conservative eigenvector
       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VEL1_VAR,SHOCKLEFT) = v1 - a
       reig(VEL2_VAR,SHOCKLEFT) = v2
       reig(PRES_VAR,SHOCKLEFT) = H - a*v1

       reig(DENS_VAR,SLOWWLEFT) = 0.
       reig(VEL1_VAR,SLOWWLEFT) = 0.
       reig(VEL2_VAR,SLOWWLEFT) = 1.
       reig(PRES_VAR,SLOWWLEFT) = v2

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VEL1_VAR,CTENTROPY) = v1
       reig(VEL2_VAR,CTENTROPY) = v2
       reig(PRES_VAR,CTENTROPY) = ekin

!!$       I'll want this for 3D
!!$       reig(DENS_VAR,SLOWWRGHT) = 0.
!!$       reig(VEL1_VAR,SLOWWRGHT) = 0.
!!$       reig(VEL2_VAR,SLOWWRGHT) = 1
!!$       reig(PRES_VAR,SLOWWRGHT) = v
       
       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VEL1_VAR,SHOCKRGHT) = v1 + a
       reig(VEL2_VAR,SHOCKRGHT) = v2
       reig(PRES_VAR,SHOCKRGHT) = H + v1*a
       !reig(:,SHOCKRGHT) = hdai*reig(:,SHOCKRGHT)
       
    else
       !! Primitive eigenvector
       reig(DENS_VAR,SHOCKLEFT) = 1.
       reig(VEL1_VAR,SHOCKLEFT) = -hdai
       reig(VEL2_VAR,SHOCKLEFT) = 0.
       reig(PRES_VAR,SHOCKLEFT) = 0.5

       reig(DENS_VAR,SLOWWLEFT) = 0.
       reig(VEL1_VAR,SLOWWLEFT) = 0.
       reig(VEL2_VAR,SLOWWLEFT) = 1
       reig(PRES_VAR,SLOWWLEFT) = 0.

       reig(DENS_VAR,CTENTROPY) = 1.
       reig(VEL1_VAR,CTENTROPY) = 0.
       reig(VEL2_VAR,CTENTROPY) = 0.
       reig(PRES_VAR,CTENTROPY) = 0.

       reig(DENS_VAR,SHOCKRGHT) = 1.
       reig(VEL1_VAR,SHOCKRGHT) = hdai
       reig(VEL2_VAR,SHOCKRGHT) = 0.
       reig(PRES_VAR,SHOCKRGHT) = 0.5
              
    endif
    
    return
  end subroutine right_eigenvectors


  subroutine left_eigenvectors(V,conservative,leig,dir)
    implicit none
    integer, intent(IN) :: dir
    real, dimension(NUMB_VAR), intent(IN)  :: V
    logical :: conservative
    real, dimension(NSYS_VAR,NUMB_WAVE), intent(OUT) :: leig

    real :: a, d, g, ekin, ha2i, agi,  p, v1, v2, ad
    integer :: VEL1_VAR, VEL2_VAR
    VEL1_VAR = 1; VEL2_VAR = 2
    if (dir == XDIM) then
       VEL1_VAR = VELX_VAR
       VEL2_VAR = VELY_VAR
    elseif (dir == YDIM) then
       VEL1_VAR = VELY_VAR
       VEL2_VAR = VELX_VAR
    end if
    
    ! sound speed, and others
    a = sqrt(V(GAMC_VAR)*V(PRES_VAR)/V(DENS_VAR))
    v1 = V(VEL1_VAR)
    v2 = V(VEL2_VAR)
    d = V(DENS_VAR)
    p = V(PRES_VAR)
    g = V(GAMC_VAR) - 1.
    ad = a*d
    ha2i = 0.5/(a*a)
    agi = a/g
    ekin = 0.5*(v1**2+v2**2)
    
    if (conservative) then
!!$       !! Conservative eigenvector
       leig(DENS_VAR,SHOCKLEFT) = ekin + agi*v1
       leig(VEL1_VAR,SHOCKLEFT) = -agi - v1
       leig(VEL2_VAR,SHOCKLEFT) = -v2
       leig(PRES_VAR,SHOCKLEFT) = 1.

       leig(DENS_VAR,SLOWWLEFT) = -2.*v2*a*agi
       leig(VEL1_VAR,SLOWWLEFT) = 0.
       leig(VEL2_VAR,SLOWWLEFT) = 2*a*agi
       leig(PRES_VAR,SLOWWLEFT) = 0.

       leig(DENS_VAR,CTENTROPY) = 2*a*agi - 2.*ekin
       leig(VEL1_VAR,CTENTROPY) = 2.*v1
       leig(VEL2_VAR,CTENTROPY) = 2.*v2
       leig(PRES_VAR,CTENTROPY) = -2.
       
       leig(DENS_VAR,SHOCKRGHT) = ekin - v1*agi
       leig(VEL1_VAR,SHOCKRGHT) = -v1 + agi
       leig(VEL2_VAR,SHOCKRGHT) = -v2
       leig(PRES_VAR,SHOCKRGHT) = 1.

       leig(:,:) = leig(:,:)*g/(2.*a*a)

!!$       leig(DENS_VAR,SHOCKLEFT) = ha2i*(g*ekin + v1*a)
!!$       leig(VEL1_VAR,SHOCKLEFT) = -ha2i*(g*v1+a)
!!$       leig(VEL2_VAR,SHOCKLEFT) = -ha2i*g*v2
!!$       leig(PRES_VAR,SHOCKLEFT) = g*ha2i
!!$
!!$       leig(DENS_VAR,SLOWWLEFT) = -v2
!!$       leig(VEL1_VAR,SLOWWLEFT) = 0.
!!$       leig(VEL2_VAR,SLOWWLEFT) = 1.
!!$       leig(PRES_VAR,SLOWWLEFT) = 0.
!!$
!!$       leig(DENS_VAR,CTENTROPY) = 1.-ha2i*g*ekin*2.
!!$       leig(VEL1_VAR,CTENTROPY) = g*v1/(a**2)
!!$       leig(VEL2_VAR,CTENTROPY) = g*v2/(a**2)
!!$       leig(PRES_VAR,CTENTROPY) = -g/(a**2)
!!$       
!!$       leig(DENS_VAR,SHOCKRGHT) = ha2i*(g*ekin-v1*a)
!!$       leig(VEL1_VAR,SHOCKRGHT) = -ha2i*(g*v1-a)
!!$       leig(VEL2_VAR,SHOCKRGHT) = -ha2i*g*v2
!!$       leig(PRES_VAR,SHOCKRGHT) = ha2i*g
    else
       !! Primitive eigenvector
       leig(DENS_VAR,SHOCKLEFT) = 0.
       leig(VEL1_VAR,SHOCKLEFT) = -ad
       leig(VEL2_VAR,SHOCKLEFT) = 0.
       leig(PRES_VAR,SHOCKLEFT) = 1.

       leig(DENS_VAR,SLOWWLEFT) = 0.
       leig(VEL1_VAR,SLOWWLEFT) = 0.
       leig(VEL2_VAR,SLOWWLEFT) = 1.
       leig(PRES_VAR,SLOWWLEFT) = 0.

       leig(DENS_VAR,CTENTROPY) = 1.
       leig(VEL1_VAR,CTENTROPY) = 0.
       leig(VEL2_VAR,CTENTROPY) = 0.
       leig(PRES_VAR,CTENTROPY) = -ha2i
       
       leig(DENS_VAR,SHOCKRGHT) = 0.
       leig(VEL1_VAR,SHOCKRGHT) = ad
       leig(VEL2_VAR,SHOCKRGHT) = 0.
       leig(PRES_VAR,SHOCKRGHT) = 1.
       

    endif
    
    return
  end subroutine left_eigenvectors


  
end module eigensystem
