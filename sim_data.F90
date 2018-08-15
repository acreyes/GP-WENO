module sim_data
 
#include "definition.h"

  implicit none

  !! numerics
  real, save :: sim_cfl, sim_tmax, sim_outputIntervalTime, sim_WENeps, sim_dt, sim_time, sim_bcT
  integer, save :: sim_order, sim_nStep, sim_Torder, sim_mval
  character(len=MAX_STRING_LENGTH), save :: sim_name, sim_limiter,sim_riemann, sim_WENO
  logical, save :: sim_charLimiting, sim_RK, sim_fixDt, sim_nlim, sim_reconMultiD, sim_intFlux, sim_hdf5

  !! ICs
  real, save :: sim_gamma
  real, save :: sim_densL,sim_velxL,sim_presL
  real, save :: sim_densR,sim_velxR,sim_presR
  real, save :: sim_shockLoc, sim_x0, sim_y0, sim_boost
  real, save :: sim_smallPres
  real, dimension(NSYS_VAR), save :: sim_Q1, sim_Q2, sim_Q3, sim_Q4
  character(len=MAX_STRING_LENGTH), save :: sim_icType

  !! BCs
  character(len=MAX_STRING_LENGTH), save :: sim_bcTypeX, sim_bcTypeY, sim_bcTypeZ
  integer                         , save :: sim_xBC, sim_yBC, sim_zBC

  !! IO
  integer, save :: sim_ioNfreq
  real,    save :: sim_ioTfreq

  !! GP
  character(len=MAX_STRING_LENGTH), save :: sim_quad, sim_gp_kernel
  real, save :: sim_sigdel, sim_sigma, sim_matern_nu, sim_RQ_alpha
  logical, save :: sim_GPFlux, sim_DongwookFlux

end module sim_data
