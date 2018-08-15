module GP
#include "definition.h"

  use gp_data

  real :: rt2 = SQRT(2.)



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! Quadrature Rules for Integrated Kernel!!!!!!

  function intg_kernel(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f
    f = 0.
    if (gp_quad == 'exact') then
       f = quad_exact(x, y, eldel)
!!$    elseif (gp_quad == 'midpt') then
!!$       f = K(x, y)
!!$    elseif (gp_quad == 'trap') then
!!$       f = 1./4. * ( &
!!$            K(x + 0.5, y + 0.5) + K(x - 0.5, y + 0.5) + &
!!$            K(x + 0.5, y - 0.5) + K(x - 0.5, y - 0.5) )
!!$    elseif (gp_quad == 'simpson') then
!!$       f = 1./36. * ( &
!!$            K(x+0.5,y+0.5) + 4.*K(x,y+0.5) + K(x-0.5,y+0.5) + 4.*( &
!!$            K(x+0.5,y    ) + 4.*K(x,y    ) + K(x-0.5,y    ) ) + &
!!$            K(x+0.5,y-0.5) + 4.*K(x,y-0.5) + K(x-0.5,y-0.5) )
    
    end if
    return
  end function intg_kernel

  function quad_cross(x, t, eldel) result(f)
    implicit none
    real, intent(IN) :: x, t, eldel
    real :: f
    f = 0.
    if (gp_quad == 'midpt') then
       !f = K(x, t)
!!$    elseif (gp_quad == 'trap') then
!!$       f = 0.5*( K(x-0.5, t) + K(x+0.5, t) )
!!$    elseif (gp_quad == 'simpson') then
!!$       f = 1./6. * ( K(x-0.5, t) + 4.*K(x,t) + K(x+0.5,t) )
    elseif (gp_quad == 'exact') then
       f = eldel*SQRT(.5*PI)*int_egrand(x,t,eldel)
    end if
    return
  end function quad_cross

  function intg_predvec(x, eldel) result(T)
    implicit none
    real, intent(IN) :: x, eldel
    real, dimension(2) :: T

    T(1) = quad_cross(x, -0.5, eldel)
    T(2) = quad_cross(x,  0.5, eldel)
    return
  end function intg_predvec

  

  function K(x, y) result(f)
    implicit none
    real, intent(IN) :: x, y
    real :: f

    f = 0.
    if (gp_kernel == 'matern') then
       !f = matern(x,y)
    elseif (gp_kernel == 'SE') then
       !f = SE(x,y)
    elseif (gp_kernel == 'RQ') then
       !f = RQ(x,y)
    elseif (gp_kernel == 'NN') then
       !f = NN(x,y)
    elseif (gp_kernel == 'GB') then
       !f = gibbs(x,y)
    end if
    return
  end function K
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Kernel Functions !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    

!!!!!!!!!!!! Squared Exponential!!!!!!!!!!!!!!!!!!
  
 
  function SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, r
    r = abs(x-y)
    f = EXP( -0.5*(r/eldel)**2 )
    return
  end function SE

!!$  function SE_der_cov_dxy(x, y) result(f)
!!$    !this is the covariance of the derivative at x and the data at y use SE kernel
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f
!!$    f = SE(x,y)*(y - x)*gr_dx/(gp_el**2)
!!$    return
!!$  end function SE_der_cov_dxy
  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!! SE Kernel !!!!!!!!!!!!!!!!!!!

  function d2_SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, ell

    f = SE(x,y,eldel)/(eldel**2)*(((x-y)/eldel)**2-1.)
    return
    
  end function d2_SE

  function d4_SE(x, y, eldel) result(f)
    implicit none
    real, intent(IN) :: x, y, eldel
    real :: f, ell, xmy

    xmy = (x - y)/eldel

    f = SE(x,y,eldel)/(eldel**4)*( (xmy)**4 - 6.*(xmy)**2 + 3. )
    return
  end function d4_SE

  function int_egrand(x, t, eldel) result(f)
    implicit none
    real, intent(IN) :: x, t, eldel
    real :: f
    !eldel = gp_eldel*SQRT(2.)
    !rt2 = SQRT(2.)

    f = ERF( (x + .5 - t)/(rt2*eldel)) - ERF( (x - .5 - t)/(rt2*eldel) )
    !f = ERF( (x + .5*gr_dx - t)/eldel) - ERF( (x - .5*gr_dx - t)/eldel )

  end function int_egrand

  function quad_exact(x1,x2,eldel) result(Integ)
    
    !exact quadrature, only good for SE kernel
    real, intent(IN) :: x1, x2, eldel

    real :: Integ, yxp, yxn, yxm
    !eldel = gp_eldel*SQRT(2.)

    yxp = (x1 - x2 + 1.)/(rt2*eldel)
    yxn = (x1      -x2)/(rt2*eldel)
    yxm = (x1 - x2 -1.)/(rt2*eldel)
    
    
    Integ = SQRT(PI)*(eldel)**2 *( yxp*ERF(yxp) + yxm*ERF(yxm) &
         - 2.*( yxn*ERF(yxn) + 1./SQRT(PI) *EXP(-yxn**2) ) &
         + 1./SQRT(PI) * ( EXP(-yxp**2) + exp(-yxm**2) ) )
    return
  end function quad_exact

  function quad_mid(x1, x2) result(Integ)
    !midpoint quadrature rule
    implicit none
    real, intent(IN) :: x1, x2
    real :: Integ, eldel
    eldel = gp_eldel*SQRT(2.)

    !Integ = 0.5*SQRT(PI)*eldel*int_egrand(x2, x1, eldel)
    return
  end function quad_mid

  function quad_simps(x1,x2) result(Integ)
    !simpson's quadrature rule
    implicit none
    real, intent(IN) :: x1, x2
    real :: Integ, eldel
    eldel = gp_eldel*SQRT(2.)

    !Integ = 0.5*SQRT(PI)*eldel*( int_egrand(x2,x1-0.5) + 4.*int_egrand(x2,x1) + int_egrand(x2,x1+0.5) )/6.
    return
  end function quad_simps

  function quad_trap(x1,x2) result(Integ)
    !trapezoidal quadrature rule
    implicit none
    real, intent(IN) :: x1, x2
    real :: Integ, eldel
    eldel = gp_eldel*SQRT(2.)

    !Integ = 0.25*SQRT(PI)*eldel*( int_egrand(x2,x1-0.5) + int_egrand(x2,x1+0.5) )
    return
  end function quad_trap
    
  function int_SEcov(x1, x2, eldel) result(Integ)
    !integrates the covariance between cells centered at x1 & x2 in units of x/delta (see eq 15)
    !cell ranges from x-1/2 to x+1/2
    implicit none
    real, intent(IN) :: x1, x2, eldel
    real :: Integ, dN, t, fa, fb, yxp, yxn, yxm
    integer :: i,N
    Integ = 0.
    if (gp_quad == 'exact') then
       Integ = quad_exact(x1,x2,eldel)
    elseif (gp_quad == 'midpt') then
       Integ = quad_mid(x1,x2)
    elseif (gp_quad == 'trap') then
       Integ = quad_trap(x1,x2)
    elseif (gp_quad == 'simpson') then
       Integ = quad_simps(x1,x2)
    end if
    return
  end function int_SEcov

  function cross_cor(x,eldel) result(T)
    !returns the cross-correlation between the left and right states and the cell centered at x
    !see eq 24
    implicit none
    real, intent(IN) :: x, eldel
    real, dimension(2) :: T
    !eldel = gp_eldel*SQRT(2.)
    !gp_eldel = gp_el/gr_dx
    T(1) = int_egrand(x, -.5, eldel)
    T(2) = int_egrand(x, 0.5, eldel)
    T = T*.5*eldel*SQRT(PI)

  end function cross_cor


!!$  !!!!!!!!!!!!!!!! Matern Kernel!!!!!!!!!!!!!!!!!!!!
!!$
!!$  function matern(x,y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f
!!$
!!$    f = 0.
!!$
!!$    if (gp_matern_nu == 0.5) then
!!$       f = mat_1h(x,y)
!!$    elseif (gp_matern_nu == 1.5) then
!!$       f = mat_3h(x,y)
!!$    elseif (gp_matern_nu == 2.5) then
!!$       f = mat_5h(x,y)
!!$    elseif (gp_matern_nu == 3.5) then
!!$       f = mat_7h(x,y)
!!$    elseif (gp_matern_nu == 4.5) then
!!$       f = mat_9h(x,y)
!!$    elseif (gp_matern_nu == 5.5) then
!!$       f = mat_11h(x,y)
!!$    end if
!!$    return
!!$  end function matern
!!$  
!!$!!!!!! half integer matern kernels!!!!!!!!!!!!!!!!
!!$  
!!$
!!$  function mat_1h(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r
!!$    r = abs(x-y)
!!$    f = EXP(-r/gp_eldel)
!!$    return
!!$  end function mat_1h
!!$
!!$  function mat_3h(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, rt3li
!!$    r = abs(x-y)
!!$    rt3li = SQRT(3.)/gp_eldel
!!$    f = (1 + rt3li*r)*EXP(-rt3li*r)
!!$    return
!!$  end function mat_3h
!!$
!!$  function mat_5h(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, rt5li
!!$    r = abs(x-y)
!!$    rt5li = SQRT(5.)/gp_eldel
!!$    f = (1 + rt5li*r + 5./3. * (r/gp_eldel)**2)*EXP(-rt5li*r)
!!$    return
!!$  end function mat_5h
!!$
!!$  function mat_7h(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, rt7li
!!$    r = abs(x-y)
!!$    rt7li = SQRT(7.)/gp_eldel
!!$
!!$    f = ( 1 + rt7li*r + 14./5.*(r/gp_eldel)**2 + 7.*SQRT(7.)/15.*(r/gp_eldel)**3 )*EXP(-rt7li*r)
!!$
!!$    return
!!$  end function mat_7h
!!$
!!$  function mat_9h(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, rli
!!$    r = abs(x-y)
!!$    rli = r/gp_eldel
!!$    
!!$    f = ( 1. + 3.*rli + 27./7.*rli**2 + 18./7.*rli**3 + 27./35.*rli**4 )*EXP(-3.*rli)
!!$    
!!$    return
!!$  end function mat_9h
!!$
!!$  function mat_11h(x,y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, rli
!!$    r = abs(x-y)
!!$    rli = r/gp_eldel
!!$
!!$    f = ( 1. + SQRT(11.)*rli + 44./9.*rli**2 + 11./9.*SQRT(11.)*rli**3 + 121./63.*rli**4 + &
!!$         121./945.*SQRT(11.)*rli**5 )*EXP(-SQRT(11.)*rli)
!!$
!!$    return
!!$  end function mat_11h
!!$  
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!!!!!! Exponential Kernel !!!!!!!!!!
!!$  
!!$  function intg_exp(x, t) result(f)
!!$    !prediction vector for exponential kernel
!!$    implicit none
!!$    real, intent(IN) :: x,t
!!$    real :: f
!!$    f = 0.
!!$    if (x < t) then
!!$       f = gp_eldel*(EXP((x+0.5)/gp_eldel) - EXP((x-0.5)/gp_eldel) )*EXP(-t/gp_eldel)
!!$    elseif (x > t) then
!!$       f = gp_eldel*(EXP(-(x-0.5)/gp_eldel) - EXP(-(x+0.5)/gp_eldel) )*EXP(t/gp_eldel)
!!$    end if
!!$    return
!!$  end function intg_exp
!!$
!!$  function EXP_exact(x1,x2) result(Integ)
!!$    !exact quadrature for exponential kernel
!!$    implicit none
!!$    real, intent(IN) :: x1, x2
!!$    real :: Integ, x, y
!!$    Integ = 0.
!!$    if (x1 .ne. x2) then
!!$       x = MAX(x1, x2)
!!$       y = MIN(x1, x2)
!!$
!!$       Integ = gp_eldel**2 *( EXP(-(x-0.5)/gp_eldel) - EXP(-(x+0.5)/gp_eldel) )*&
!!$            ( EXP((y+0.5)/gp_eldel) - EXP((y-0.5)/gp_eldel) )
!!$    elseif (x1 == x2) then
!!$       x = x1
!!$       y = x2
!!$       Integ = gp_eldel*( &
!!$            2. + gp_eldel*(&
!!$            EXP( (y-0.5)/gp_eldel)*( EXP(-(x+0.5)/gp_eldel) - EXP(-(x-0.5)/gp_eldel) ) - &
!!$            EXP(-(y+0.5)/gp_eldel)*( EXP( (x+0.5)/gp_eldel) - EXP( (x-0.5)/gp_eldel) )   &
!!$            )&
!!$            )
!!$    end if
!!$    return
!!$  end function EXP_exact
!!$
!!$  function int_EXPcov(x1, x2) result(Integ)
!!$
!!$    implicit none
!!$    real, intent(IN) :: x1, x2
!!$    real :: Integ
!!$    Integ = 0.
!!$    if (gp_quad == 'exact') then
!!$       Integ = EXP_exact(x1,x2)
!!$    end if
!!$    return
!!$  end function int_EXPcov
!!$
!!$  function cross_EXP(x) result(T)
!!$    !cross correlation for the prediction using the exponential kernel
!!$    implicit none
!!$    real, intent(IN) :: x
!!$    real, dimension(2) :: T
!!$
!!$    T(1) = intg_EXP(x, -0.5)
!!$    T(2) = intg_EXP(x, 0.5)
!!$    return
!!$  end function cross_EXP

!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!!!!!!!!!!!!!! Rational Quad. !!!!!!!!!!!!!!!!!!
!!$
!!$  function RQ(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r
!!$    r = abs(x-y)
!!$    f = (1. + r**2/(2.*gp_RQ_alpha*gp_eldel**2))**(-gp_RQ_alpha)
!!$    return
!!$  end function RQ

!!$  !!!!!!!!!!! Neural Network (NN) !!!!!!!!!!!!!!!!!!
!!$
!!$  function NN(x,y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, sig0, xx, yy, xy
!!$
!!$    sig0 = 1.
!!$    !print *, sig0, x, y
!!$    xx = sig0**2 + x*x/(gp_eldel**2)
!!$    yy = sig0**2 + y*y/(gp_eldel**2)
!!$    xy = sig0**2 + x*y/(gp_eldel**2)
!!$    f = 2./PI * ASIN(2.*xy/SQRT((1.+2.*xx)*(1.+2.*yy)))
!!$    return
!!$  end function NN
!!$
!!$!!!!!!!!!!!! Gibbs Covariance !!!!!!!!!!!!!!!!!!!!
!!$
!!$  function lscale(x) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x
!!$    real :: f, a, r
!!$
!!$    a = 0.8
!!$    r = x - 1.
!!$    f = 1. - 0.6*EXP(-0.5*(r/a)**2)
!!$    f = f*gp_eldel
!!$    return
!!$  end function lscale
!!$
!!$  function gibbs(x, y) result(f)
!!$    implicit none
!!$    real, intent(IN) :: x, y
!!$    real :: f, r, lx, ly
!!$
!!$    lx = lscale(x)
!!$    ly = lscale(y)
!!$    r = abs(x-y)
!!$
!!$    f = exp( -(r**2)/(lx**2+ly**2) )
!!$    f = f*SQRT( (2.*lx*ly)/(lx**2+ly**2) )
!!$
!!$    return
!!$  end function gibbs

  
end module GP
