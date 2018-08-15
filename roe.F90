subroutine roe(vL,vR,Flux,dir)

#include "definition.h"  

  use grid_data
  use primconsflux, only : prim2flux,prim2cons
  use eigensystem

  implicit none
  integer, intent(IN) :: dir
  real, dimension(NUMB_VAR), intent(IN) :: vL,vR !prim vars
  real, dimension(NSYS_VAR), intent(OUT):: Flux 

  real, dimension(NSYS_VAR)  :: FL,FR,uL,uR
  real, dimension(NUMB_VAR)  :: vAvg
  real, dimension(NUMB_WAVE) :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  logical :: conservative
  real, dimension(NSYS_VAR) :: vec, sigma
  integer :: kWaveNum
  
  ! set the initial sum to be zero
  sigma(DENS_VAR:ENER_VAR) = 0.
  vec(DENS_VAR:ENER_VAR)   = 0.
  
  ! we need conservative eigenvectors
  conservative = .true.
  
  call averageState(vL,vR,vAvg)
  call eigenvalues(vAvg,lambda,dir)
  call eigenvectors(vAvg,conservative,reig,leig,dir)
  !call left_eigenvectors (vAvg,conservative,leig,dir)
  !call right_eigenvectors(vAvg,conservative,reig,dir)
  
  call prim2flux(vL,FL,dir)
  call prim2flux(vR,FR,dir)
  call prim2cons(vL,uL)
  call prim2cons(vR,uR)


  do kWaveNum = 1, NUMB_WAVE
     ! STUDENTS: PLEASE FINISH THIS ROE SOLVER
     vec(DENS_VAR:ENER_VAR) = uR(DENS_VAR:ENER_VAR) - uL(DENS_VAR:ENER_VAR)
     sigma(DENS_VAR:ENER_VAR) = sigma(DENS_VAR:ENER_VAR) + dot_product(leig(DENS_VAR:ENER_VAR,kWaveNum), vec(DENS_VAR:ENER_VAR))*ABS(lambda(kWaveNum))*reig(DENS_VAR:ENER_VAR, kWaveNum)
          
  end do
  
  ! numerical flux
  Flux(DENS_VAR:ENER_VAR) = 0.5*(FL(DENS_VAR:ENER_VAR) + FR(DENS_VAR:ENER_VAR)) - 0.5*sigma


  return
end subroutine roe
