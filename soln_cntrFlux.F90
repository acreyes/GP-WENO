subroutine soln_cntrFlux(F, Npts, Fiph)
  ! here we do the flux correction like Del Zanna
  ! Npts should be equal to 5 only

#include "definition.h"

  implicit none
  
  integer, intent(IN) :: Npts
  real, dimension(Npts), intent(IN) :: F
  real, intent(OUT) :: Fiph

  integer :: var
  real :: fac2, fac4, D2, D4, Fiph_fac
  real, dimension(2) :: coeff2
  real, dimension(4) :: coeff4


  fac2 = -1./24.
  fac4 = 1./480.
  coeff2 = 4.*(/1., 1. /)
  coeff4 = (/8., -72., -72., 8. /)/3.
  Fiph_fac = 1. + fac2*4.*(-2.) + fac4*128./3.
  
  D2 = dot_product(coeff2, F(2:3))
  D4 = dot_product(coeff4, F(1:4))

  Fiph = Fiph_fac*F(5) + fac2*D2 + fac4*D4

  return
end subroutine soln_cntrFlux
  
