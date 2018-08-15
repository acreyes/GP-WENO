subroutine averageState(vL,vR,vAvg)

#include "definition.h"  

  implicit none
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NUMB_VAR), intent(OUT) :: vAvg  !average state

  vAvg(:) = .5*(vL(:) + vR(:))
  
  return
end subroutine averageState
