subroutine sim_initBlock()
 
#include "definition.h"
  
  use sim_data
  use grid_data, only : gr_V,gr_U,gr_ngc,gr_xbeg &
                       ,gr_i0,gr_imax,gr_xCoord,gr_dx, gr_yCoord, gr_dy, gr_zCoord, gr_dz
  use block_data, only   : bl_i, bl_j, bl_k, bl_nBlock, bl_ID, bl_grid  
  use primconsflux, only : prim2cons
  use bc
  
  implicit none

  integer :: i,j,k,in,jn,kn
  real :: ekin, eint, x, y,z, ranx, rany, small, r2, beta, T, dr2, E, gamm
  real :: rt3, xmin, x0, ymin, sig, xx, yy, zz, dr3, loc, thkns, lngth
  logical :: stepL, stepR


  
  small = 0.01
  sig = 0.05/SQRT(2.)
  call RANDOM_SEED()
  
  ! generate x-coordinate
  do i = gr_i0(XDIM),gr_imax(XDIM)
     gr_xCoord(i) = (real(i-gr_ngc + (bl_i-1)*bl_nBlock(XDIM))-0.5)*gr_dx + gr_xbeg(XDIM)
  end do

  ! generate y-coordinates
  do i = gr_i0(YDIM),gr_imax(YDIM)
     gr_yCoord(i) = (real(i-gr_ngc+(bl_j-1)*bl_nBlock(YDIM))-0.5)*gr_dy + gr_xbeg(YDIM)
  end do

  !generate z-coordinates
  do i = gr_i0(ZDIM),gr_imax(ZDIM)
     gr_zCoord(i) = (real(i-gr_ngc+(bl_k-1)*bl_nBlock(ZDIM)) - 0.5)*gr_dz + gr_xbeg(ZDIM)
  end do
  
  !some handy constants
  gamm = sim_gamma
  !these are used for sedov
  dr2 = (3.5*MIN(gr_dx, gr_dy, gr_dz))
  dr3 = (3.5*MIN(gr_dx, gr_dy, gr_dz))**3
  E = 1.

  rt3 = sqrt(3.)
  xmin = 1./6.
  x0 = xmin + 1./rt3

#ifdef BDRY_VAR
  !initialize boundary vars everywhere 
  gr_V(BDRY_VAR,:,:,:) = -1.0
#endif  

  do i = gr_i0(XDIM),gr_imax(XDIM)
     xx = gr_xCoord(i)
     do j = gr_i0(YDIM),gr_imax(YDIM)
        yy = gr_yCoord(j)
        do k = gr_i0(ZDIM),gr_imax(ZDIM)
           zz = gr_zCoord(k)
           if (sim_icType == 'vortex') then
              beta = 5.
              x = xx - 5.
              y = yy - 5.
              z = zz - 5.
              r2 = x**2 + y**2
              T = 1. - (sim_gamma-1.)*beta*beta*EXP(1.-r2)/(8.*sim_gamma*PI*PI)

              gr_V(DENS_VAR,i,j,k) = (T)**(1./(sim_gamma-1.))
              gr_V(VELX_VAR,i,j,k) = 1. - y*beta/(2.*PI)*EXP(0.5*(1.-r2))
              gr_V(VELY_VAR,i,j,k) = 1. + x*beta/(2.*PI)*EXP(0.5*(1.-r2))
              gr_V(VELZ_VAR,i,j,k) = 1.
              gr_V(PRES_VAR,i,j,k) = gr_V(DENS_VAR,i,j,k)*T

           elseif (sim_icType == 'sedov') then
              r2 = xx**2 + yy**2 + zz**2
              if (sqrt(r2) < dr2) then
                 gr_V(PRES_VAR,i,j,k) = (gamm-1.)*3.*E/(4.*PI*dr3)
              else
                 gr_V(PRES_VAR,i,j,k) = 1.e-5
              end if

              gr_V(DENS_VAR,i,j,k) = 1.
              gr_V(VELX_VAR,i,j,k) = 0.
              gr_V(VELY_VAR,i,j,k) = 0.
              gr_V(VELZ_VAR,i,j,k) = 0.

           elseif (sim_icType == 'sedovChamber') then
              loc = .25
              thkns = 0.025
              lngth = .15
              x = xx
              y = yy
              z = zz

              r2 = xx**2 + yy**2 + zz**2
              if (sqrt(r2) < dr2) then
                 gr_V(PRES_VAR,i,j,k) = (gamm-1.)*3.*E/(4.*PI*dr3)
              else
                 gr_V(PRES_VAR,i,j,k) = 1.e-5
              end if

              gr_V(DENS_VAR,i,j,k) = 1.
              gr_V(VELX_VAR,i,j,k) = 0.
              gr_V(VELY_VAR,i,j,k) = 0.
              gr_V(VELZ_VAR,i,j,k) = 0.

#ifdef BDRY_VAR
              if (abs(z) < loc) then
                 if ( (abs(x-loc) .le. thkns) .and. (abs(y) .le. lngth) ) then
                    gr_V(BDRY_VAR,i,j,k) = 1.
                    gr_V(DENS_VAR,i,j,k) = 0.
                 elseif ( (abs(x+loc) .le. thkns) .and. (abs(y) .le. lngth) ) then
                    gr_V(BDRY_VAR,i,j,k) = 1.
                    gr_V(DENS_VAR,i,j,k) = 0.
                 elseif ( (abs(y+loc) .le. thkns) .and. (abs(x) .le. lngth) ) then
                    gr_V(BDRY_VAR,i,j,k) = 1.
                    gr_V(DENS_VAR,i,j,k) = 0.
                 elseif ( (abs(y-loc) .le. thkns) .and. (abs(x) .le. lngth) ) then
                    gr_V(BDRY_VAR,i,j,k) = 1.
                    gr_V(DENS_VAR,i,j,k) = 0.
                 else
                    gr_V(BDRY_VAR,i,j,k) = -1.0
                 endif
              end if
#endif

           elseif (sim_icType == 'windtunnel2D') then
              if ((xx .ge. 0.6) .and. (zz .le. 0.2)) then
                 !this is the step
                 gr_V(DENS_VAR,i,j,k) = 10.4
                 gr_V(PRES_VAR,i,j,k) = 1.0
                 gr_V(VELX_VAR,i,j,k) = 0.0
                 gr_V(VELY_VAR,i,j,k) = 0.0
                 gr_V(VELZ_VAR,i,j,k) = 0.0
#ifdef BDRY_VAR                 
                 gr_V(BDRY_VAR,i,j,k) = 1.
                 gr_V(DENS_VAR,i,j,k) = 0.
#endif
              else
                 gr_V(DENS_VAR,i,j,k) = 1.4
                 gr_V(PRES_VAR,i,j,k) = 1.0
                 gr_V(VELX_VAR,i,j,k) = 3.0
                 gr_V(VELY_VAR,i,j,k) = 0.0
                 gr_V(VELZ_VAR,i,j,k) = 0.0
              end if

           elseif (sim_icType == 'windtunnel3D') then
              stepL = ( xx .ge. 0.6) .and. (abs(yy) > 0.5)
              stepR = ( xx .ge. 0.5) .and. (abs(yy) .le. 0.5)
              !stepL = ( xx .ge. 0.6) .and. (yy > 0.)
              !stepR = ( xx .ge. 1.2) .and. (yy .le. 0.)
              !if ( (zz .le. 0.2) .and. (stepL  .or. stepR)) then
              if ((zz .le. 0.2) .and. stepR ) then
                 !this is the step
                 gr_V(DENS_VAR,i,j,k) = 10.4
                 gr_V(PRES_VAR,i,j,k) = 1.0
                 gr_V(VELX_VAR,i,j,k) = 0.0
                 gr_V(VELY_VAR,i,j,k) = 0.0
                 gr_V(VELZ_VAR,i,j,k) = 0.0
#ifdef BDRY_VAR                 
                 gr_V(BDRY_VAR,i,j,k) = 1.
                 gr_V(DENS_VAR,i,j,k) = 0.
#endif
              else
                 gr_V(DENS_VAR,i,j,k) = 1.4
                 gr_V(PRES_VAR,i,j,k) = 1.0
                 gr_V(VELX_VAR,i,j,k) = 3.0
                 gr_V(VELY_VAR,i,j,k) = 0.0
                 gr_V(VELZ_VAR,i,j,k) = 0.0
              end if
              
           end if

           gr_V(GAMC_VAR,i,j,k) = sim_gamma
           gr_V(GAME_VAR,i,j,k) = sim_gamma
           gr_V(EINT_VAR,i,j,k) = gr_V(PRES_VAR,i,j,k)/(gr_V(GAME_VAR,i,j,k)-1.)/gr_V(DENS_VAR,i,j,k)

           !initialize conservative vars
           call prim2cons(gr_V(:,i,j,k), gr_U(DENS_VAR:ENER_VAR,i,j,k))
        end do
     end do
  end do

  call bc_apply(gr_V)

  
!!$  do i = gr_i0(XDIM),gr_imax(XDIM)
!!$     do j = gr_i0(YDIM),gr_imax(YDIM)
!!$        if (sim_icType == 'gaussx') then
!!$           !advect gaussian in the x-direction
!!$           x = gr_xCoord(i)-0.5
!!$           y = gr_yCoord(j)-0.5
!!$           gr_V(DENS_VAR,i,j) = 1. + EXP(-100.*(y**2))
!!$           gr_V(VELX_VAR,i,j) = 0.
!!$           gr_V(VELY_VAR,i,j) = 1.
!!$           gr_V(PRES_VAR,i,j) = 1./sim_gamma
!!$
!!$        elseif (sim_icType == 'gaussxy') then
!!$           x = gr_xCoord(i)-0.5
!!$           y = gr_yCoord(j)-0.5
!!$           
!!$           gr_V(DENS_VAR,i,j) = 1. + EXP(-(x**2+y**2)/(0.1)**2)
!!$           gr_V(VELX_VAR,i,j) = 1.
!!$           gr_V(VELY_VAR,i,j) = 1.
!!$           gr_V(PRES_VAR,i,j) = 1./sim_gamma
!!$
!!$        elseif (sim_icType == 'DMR') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           ymin = (x-xmin)*rt3
!!$           if (y > ymin) then
!!$              !in the shock region
!!$              gr_V(DENS_VAR,i,j) = 8.
!!$              gr_V(VELX_VAR,i,j) = 7.1447096
!!$              gr_V(VELY_VAR,i,j) = -4.125
!!$              gr_V(PRES_VAR,i,j) = 116.5
!!$           else
!!$              !outside shock
!!$              gr_V(DENS_VAR,i,j) = 1.4
!!$              gr_V(VELX_VAR,i,j) = 0.
!!$              gr_V(VELY_VAR,i,j) = 0.
!!$              gr_V(PRES_VAR,i,j) = 1.
!!$           end if
!!$           
!!$        elseif (sim_icType == '2DRP') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           if (x > sim_x0 .and. y .ge. sim_y0) then
!!$              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q1(:)
!!$           elseif (x .le. sim_x0 .and. y .ge. sim_y0) then
!!$              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q2(:)
!!$           elseif ( x .le. sim_x0 .and. y < sim_y0) then
!!$              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q3(:)
!!$           elseif (x > sim_x0 .and. y < sim_y0) then
!!$              gr_V(DENS_VAR:PRES_VAR,i,j) = sim_Q4(:)
!!$           end if
!!$           
!!$        elseif (sim_icType == 'shocky') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           if (y < 0.5) then
!!$              gr_V(DENS_VAR,i,j) = 1.
!!$              gr_V(PRES_VAR,i,j) = 1.
!!$           else
!!$              gr_V(DENS_VAR,i,j) = 0.125
!!$              gr_V(PRES_VAR,i,j) = 0.1
!!$           end if
!!$           gr_V(VELX_VAR,i,j) = 0.
!!$           gr_V(VELY_VAR,i,j) = 0.
!!$        elseif (sim_icType == 'shockx') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           if (x < 0.5) then
!!$              gr_V(DENS_VAR,i,j) = 1.
!!$              gr_V(PRES_VAR,i,j) = 1.
!!$           else
!!$              gr_V(DENS_VAR,i,j) = 0.125
!!$              gr_V(PRES_VAR,i,j) = 0.1
!!$           end if
!!$           gr_V(VELX_VAR,i,j) = 0.
!!$           gr_V(VELY_VAR,i,j) = 0.
!!$        elseif (sim_icType == 'shockxy') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           r2 = x**2 + y**2
!!$           if (r2 < 0.25) then
!!$              gr_V(DENS_VAR,i,j) = 1.
!!$              gr_V(PRES_VAR,i,j) = 1.
!!$           else
!!$              gr_V(DENS_VAR,i,j) = 0.125
!!$              gr_V(PRES_VAR,i,j) = 0.1
!!$           end if
!!$           gr_V(VELX_VAR,i,j) = 0.
!!$           gr_V(VELY_VAR,i,j) = 0.
!!$
!!$        elseif (sim_icType == 'vortex') then
!!$           beta = 5.
!!$           x = gr_xCoord(i) - 5.
!!$           y = gr_yCoord(j) - 5.
!!$           r2 = x**2+y**2
!!$           T = 1. - (sim_gamma-1.)*beta*beta*EXP(1.-r2)/(8.*sim_gamma*PI*PI)
!!$           
!!$           gr_V(DENS_VAR,i,j) = (T)**(1./(sim_gamma-1.))
!!$           gr_V(VELX_VAR,i,j) = 1. - y*beta/(2.*PI)*EXP(0.5*(1.-r2))
!!$           gr_V(VELY_VAR,i,j) = 1. + x*beta/(2.*PI)*EXP(0.5*(1.-r2))
!!$           gr_V(PRES_VAR,i,j) = gr_V(DENS_VAR,i,j)*T
!!$
!!$        elseif (sim_icType == 'sedov') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$           r2 = x**2 + y**2
!!$           if (sqrt(r2) < sqrt(dr2)) then
!!$              gr_V(PRES_VAR,i,j) = (gamm-1.)*E/(PI*dr2)
!!$           else
!!$              gr_V(PRES_VAR,i,j) = 1.e-5
!!$           end if
!!$
!!$           gr_V(DENS_VAR,i,j) = 1.
!!$           gr_V(VELX_VAR,i,j) = 0.
!!$           gr_V(VELY_VAR,i,j) = 0.
!!$
!!$        elseif (sim_icType == 'KH') then
!!$           x = gr_xCoord(i)
!!$           y = gr_yCoord(j)
!!$
!!$           if (abs(y) > 0.25) then
!!$              gr_V(VELX_VAR,i,j) = -0.5
!!$              gr_V(DENS_VAR,i,j) = 1.
!!$           else
!!$              gr_V(VELX_VAR,i,j) = 0.5
!!$              gr_V(DENS_VAR,i,j) = 2.
!!$           end if
!!$           gr_V(PRES_VAR,i,j) = 2.5
!!$           gr_V(VELX_VAR,i,j) = gr_V(VELX_VAR,i,j) + sim_boost !+ small*ranx
!!$           gr_V(VELY_VAR,i,j) = sim_boost +  0.1*SIN(4.*PI*x)*(EXP(-(y+0.25)**2/(0.05**2)) + EXP(-(y-0.25)**2/(0.05**2)))!small*rany
!!$        
!!$        end if
!!$
!!$        gr_V(GAMC_VAR,i,j) = sim_gamma
!!$        gr_V(GAME_VAR,i,j) = sim_gamma
!!$        gr_V(EINT_VAR,i,j) = gr_V(PRES_VAR,i,j)/(gr_V(GAME_VAR,i,j)-1.)/gr_V(DENS_VAR,i,j)
!!$        
!!$     end do
!!$  end do
!!$  !initialize cons vars
!!$  do i = gr_i0(XDIM),gr_imax(XDIM)
!!$     do j = gr_i0(YDIM),gr_imax(YDIM)
!!$        call prim2cons(gr_V(:,i,j), gr_U(DENS_VAR:ENER_VAR,i,j))
!!$     end do
!!$  end do
!!$
!!$  call bc_apply(gr_V)
  
  


  
end subroutine sim_initBlock
