subroutine grid_init()

#include "definition.h"

  use grid_data
  use sim_data
  use read_initFile

  implicit none

  

  ! allocate cell coordinates
  allocate(gr_xCoord(gr_nx+2*gr_ngc)); gr_xCoord = 0.0
  allocate(gr_yCoord(gr_ny+2*gr_ngc)); gr_yCoord = 0.0
  allocate(gr_zCoord(gr_nz+2*gr_ngc)); gr_zCoord = 0.0

  

  ! grid delta
  gr_dx = (gr_xend(XDIM) - gr_xbeg(XDIM))/gr_glb_nx
  gr_dy = (gr_xend(YDIM) - gr_xbeg(YDIM))/gr_glb_ny
  gr_dz = (gr_xend(ZDIM) - gr_xbeg(ZDIM))/gr_glb_nz

  ! allocate grid variables
  allocate(gr_U(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM))); gr_U = 0.
  allocate(gr_V(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM))); gr_V = 0.
  allocate(gr_W(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM))); gr_W = 0. !is this even used???

  !These still need to be retooled for multiD
  ! allocate grid Riemann states
  allocate(gr_vL(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM)); gr_vL = 0.
  allocate(gr_vR(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM)); gr_vR = 0.

  ! allocate grid fluxes
  allocate(gr_flux(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM)); gr_flux = 0. !(var,i,j,dir)
  allocate(gr_ptFluxes(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM),NDIM)); gr_ptFluxes = 0.
!!$  allocate(gr_vP(NUMB_VAR,gr_nx,gr_ny,NDIM)); gr_vP = 0.
!!$  allocate(gr_vM(NUMB_VAR,gr_nx,gr_ny,NDIM)); gr_vM = 0.

  ! allocate grid eigensystem
  allocate(gr_maxalphas(NUMB_WAVE,NDIM)); gr_maxalphas = 0.
!!$  allocate(gr_eigval(NUMB_WAVE,gr_imax)); gr_eigval = 0.
!!$  allocate(gr_leigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_leigvc = 0.
!!$  allocate(gr_reigvc(NSYS_VAR,NUMB_WAVE,gr_imax)); gr_reigvc = 0.
!!$
!!$
!!$  ! allocate GP variables
  !lets instead do this in sim_GPinit so that we can allocate according to stencil size
!!$  allocate(gr_GPv(2*gr_radius+1   )); gr_GPv = 0.
!!$  allocate(gr_GPZ(2, 2*gr_radius+1)); gr_GPZ = 0.
  gr_Tcells = 0 !this is the number of transverse cells used in the reconstructions stencil
  return
end subroutine grid_init
