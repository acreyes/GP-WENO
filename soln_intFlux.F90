subroutine soln_intFlux(F, Npts, Fiph)
  ! here we do the flux correction like Del Zanna
  ! Npts should be equal to 5 only

#include "definition.h"

  implicit none
  
  integer, intent(IN) :: Npts
  real, dimension(Npts), intent(IN) :: F
  real, intent(OUT) :: Fiph

  integer :: var
  real :: fac2, fac4, D2, D4
  real, dimension(3) :: coeff3
  real, dimension(5) :: coeff5


  fac2 = -1./24.
  fac4 = 3./640.
  coeff3 = (/1., -2., 1. /)
  coeff5 = (/1., -4., 6., -4., 1. /)
  
  D2 = dot_product(coeff3, F(2:4))
  D4 = dot_product(coeff5, F(:))

  Fiph = F(3) + fac2*D2 + fac4*D4

  return
end subroutine soln_intFlux
  
