subroutine soln_FOG(dt, Nx, V, vL, vR, dir)

#include "definition.h"
  integer, intent(IN) :: dir, Nx
  real   , intent(IN) :: dt

  real, dimension(NUMB_VAR), intent(INOUT) :: vL, vR
  real, dimension(Nx, NUMB_VAR), intent(IN) ::  V

  vL(:) = V(1,:)
  vR(:) = V(1,:)

  return
end subroutine soln_FOG
