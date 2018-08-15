subroutine soln_reconstruct(dt, V)

#include "definition.h"  

  use gp_data, only: gp_radius, gp_Npts, gp_stencil
  use grid_data
  use sim_data
  use eigensystem
  use primconsflux
  use sim_interfaces
  use reconstruction


  implicit none
  real, intent(IN) :: dt
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)), intent(IN) :: V

  real, allocatable, dimension(:,:) :: recon_stencil
  real, dimension(NUMB_WAVE)          :: lambda
  real, dimension(NSYS_VAR,NUMB_WAVE) :: reig, leig
  integer, dimension(NDIM) :: ibeg, iend
  integer :: i, j, si, sj, scntr, var, Nx, im, ip, dir, i0, j0, is, js, ig, jg, R

  !reconstruction pointers
  procedure (recon1D)          :: soln_FOG , soln_WENO, soln_gpWENO
  procedure (recon1D), pointer :: recon

  if (sim_reconMultiD) then
     !always false
  else
     R = 2
     !figure out reconstruction to use
     select case(sim_order)
     case(1)
        !FOG
        Nx = 1
        im = 0
        ip = 0
        recon => soln_FOG
     case(5)
        !WENO
        Nx = 5
        im = -2
        ip =  2
        recon => soln_WENO
     case(10)
        Nx = 2*gp_radius +1
        im = -gp_radius
        ip =  gp_radius
        !R = gp_radius
        recon => soln_gpWENO
     case default
        !default to WENO reconstruction
        Nx = 5
        im = -2
        ip =  2
        recon => soln_WENO
     end select
     do dir = XDIM, NDIM
        ibeg = gr_ibeg
        iend = gr_iend
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
        !take care of guard cell regions
        !(i) left GC's
        ibeg(dir) = gr_ibeg(dir) - (R+1)
        iend(dir) = gr_ibeg(dir) - 1
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
        !(ii) right GC's
        ibeg(dir) = gr_iend(dir)+1
        iend(dir) = gr_iend(dir)+(R+1)
        call recon1D_ij(ibeg, iend, im, ip, dir, V, recon)
     end do
  end if

  
  return
end subroutine soln_reconstruct
