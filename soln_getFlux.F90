subroutine soln_getFlux(V)

#include "definition.h"  

  use grid_data
  use sim_data
  use primconsflux
  use sim_interfaces
  use gp_data, only: gp_radius

  implicit none
  real, dimension(NUMB_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)), intent(IN) :: V

  
  procedure (rmn_slvr) :: hll, hllc, roe
  procedure (rmn_slvr), pointer :: RP
  procedure (num_flux) :: soln_intFlux, soln_cntrFlux, soln_gpFlux
  procedure (num_flux), pointer :: NumFlux
  integer :: i, j, k, var, dir, l, m, n, R, BDRY_VEL
  integer :: cntr, Npts
  integer, dimension(NDIM) :: ibeg, iend, im
  real    :: D2, D4, Fiph_fac, fac2, fac4
  real, dimension(2) :: coeff2
  real, dimension(3) :: coeff3
  real, dimension(4) :: coeff4
  real, dimension(5) :: coeff5, flux_vec
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) ::  intFlux
  real, dimension(NSYS_VAR,gr_imax(XDIM),gr_imax(YDIM),gr_imax(ZDIM)) ::  cntrFlux

  real, dimension(NUMB_VAR) :: rpL, rpR

  if (sim_order == 10) then
     R = gp_radius
  else
     R = 2
  end if
  
  if (sim_intFlux) then
     if (sim_GPFlux .and. sim_order == 10) then
        NumFlux => soln_gpFlux
        Npts    =  5!2*gp_radius + 1 
     else
        Npts = 5
        NumFlux => soln_intFlux
     end if
  else
     Npts = 5
     NumFlux => soln_cntrFlux
  end if
  
  if (sim_riemann == 'hll') then
     RP => hll
  elseif (sim_riemann == 'hllc') then
     RP => hllc
  elseif (sim_riemann == 'roe') then
     RP => roe
  else
     RP => hll
  end if

  if (sim_intFlux) then
!  if (.false.) then
     !here we correct the fluxes using the interface fluxes
     do dir = XDIM, NDIM
        ibeg = gr_ibeg
        iend = gr_iend
        im = 0
        im(dir) = 1
        ibeg(dir) = gr_ibeg(dir)-2
        iend(dir) = gr_iend(dir)+3
        BDRY_VEL = VELX_VAR + (dir-1)

        do i = ibeg(XDIM), iend(XDIM)
           do j = ibeg(YDIM), iend(YDIM)
              do k = ibeg(ZDIM), iend(ZDIM)

                 !store left & right riemann states
                 rpL(:) = gr_vR(:,i-im(XDIM),j-im(YDIM),k-im(ZDIM),dir)
                 rpR(:) = gr_vL(:,i   ,j   ,k   ,dir)

#ifdef BDRY_VAR
                 if (gr_V(BDRY_VAR,i,j,k) > 0.0 .and. gr_V(BDRY_VAR,i-im(XDIM),j-im(YDIM),k-im(ZDIM)) < 0.0) then
                    rpR = rpL
                    rpR(BDRY_VEL) = -rpL(BDRY_VEL)
                 end if
                 if (gr_V(BDRY_VAR,i,j,k) < 0.0 .and. gr_V(BDRY_VAR,i-im(XDIM),j-im(YDIM),k-im(ZDIM)) > 0.0) then
                    rpL = rpR
                    rpL(BDRY_VEL) = - rpR(BDRY_VEL)
                 end if
#endif                 
                 
                 call RP(rpL, rpR, intFlux(:,i,j,k), dir)
           end do
           end do
        end do

        do i = gr_ibeg(XDIM), gr_iend(XDIM) + im(XDIM)
           do j = gr_ibeg(YDIM), gr_iend(YDIM) + im(YDIM)
              StencilDo: do k = gr_ibeg(ZDIM), gr_iend(ZDIM) + im(ZDIM)
#ifdef BDRY_VAR
                 do l = i -2*im(XDIM), i + 2*im(XDIM)
                    do m = j-2*im(YDIM), j + 2*im(YDIM)
                       do n = k-2*im(ZDIM), k+2*im(ZDIM)
                          if (gr_V(BDRY_VAR,l,m,n) == 1.0) then
                             !physical boundary on flux stencil
                             gr_flux(:,i,j,k,dir) = intFlux(:,i,j,k)
                             cycle StencilDo
                          end if
                       end do
                    end do
                 end do
#endif                 
                 do var = DENS_VAR, ENER_VAR
                    cntr = 0
                    do l = i-2*im(XDIM), i + 2*im(XDIM)
                       do m = j-2*im(YDIM), j + 2*im(YDIM)
                          do n = k-2*im(ZDIM), k+2*im(ZDIM)
                             cntr = cntr + 1
                             flux_vec(cntr) = intFlux(var,l,m,n)
                          end do
                       end do
                    end do
                    call NumFlux(flux_vec, Npts, gr_flux(var,i,j,k,dir))
                 end do
              end do StencilDo
            end do
        end do
     end do

     !now we get y-fluxes
    
  else
     !correct fluxes using cell center fluxes
     !calculate interface fluxes
     do dir = XDIM, NDIM
        
        ibeg = gr_ibeg
        iend = gr_iend
        im = 0
        iend(dir) = gr_iend(dir) + 1
        im(dir) = 1
        
        do i = ibeg(XDIM), iend(XDIM)
           do j = ibeg(YDIM), iend(YDIM)
              do k = ibeg(ZDIM), iend(ZDIM)
                 call RP(gr_vR(  DENS_VAR:GAME_VAR,i-im(XDIM),j-im(YDIM),k-im(ZDIM),dir), &
                         gr_vL(  DENS_VAR:GAME_VAR,i         ,j         ,k         ,dir), &
                         intFlux(DENS_VAR:ENER_VAR,i         ,j         ,k             ), dir)
              end do
           end do
        end do

        do i = gr_i0(XDIM), gr_imax(XDIM)
           do j = gr_i0(YDIM), gr_imax(YDIM)
              do k = gr_i0(ZDIM), gr_imax(ZDIM)
                 call prim2flux(V(:,i,j,k),cntrFlux(:,i,j,k),dir)
              end do
           end do
        end do

        

        do i = gr_ibeg(XDIM), gr_iend(XDIM)+im(XDIM)
           do j = gr_ibeg(YDIM), gr_iend(YDIM)+im(YDIM)
              do k = gr_ibeg(ZDIM), gr_iend(ZDIM) + im(ZDIM)
                 do var = DENS_VAR, ENER_VAR
                    cntr = 0
                    do l = i-2*im(XDIM), i + im(XDIM)
                       do m = j-2*im(YDIM), j + im(YDIM)
                          do n = k-2*im(ZDIM), k + im(ZDIM)
                             cntr = cntr + 1
                             flux_vec(cntr) = cntrFlux(var,l,m,n)
                          end do
                       end do
                    end do
                    flux_vec(5) = intFlux(var,i,j,k)
                    call NumFlux(flux_vec, Npts, gr_flux(var,i,j,k,dir))
                 end do
              end do
           end do
        end do

        
        
     end do

  end if

  


  return
end subroutine soln_getFlux
