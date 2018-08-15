subroutine read_pars()
  !this subroutine is to read all of the 
#include "definition.h"

  use grid_data
  use gp_data
  use sim_data
  use block_data
  use read_initFile
  
  implicit none

  !read sim_data
  call read_initFileInt ('slug.init','sim_order',  sim_order)
  call read_initFileInt ('slug.init','sim_Torder',  sim_Torder)
  call read_initFileInt ('slug.init','sim_nstep',  sim_nStep)
  call read_initFileReal('slug.init','sim_dt',    sim_dt)
  call read_initFileReal('slug.init','sim_cfl',    sim_cfl)
  call read_initFileReal('slug.init','sim_tmax',   sim_tmax)
  call read_initFileReal('slug.init','sim_WENeps',   sim_WENeps)
  call read_initFileReal('slug.init','sim_boost',   sim_boost)
  call read_initFileReal('slug.init','sim_outputIntervalTime',sim_outputIntervalTime)
  call read_initFileChar('slug.init','sim_riemann',sim_riemann)
  call read_initFileChar('slug.init','sim_limiter',sim_limiter)
  call read_initFileChar('slug.init','sim_name',sim_name)
  call read_initFileChar('slug.init','sim_WENO',sim_WENO)
  call read_initFileBool('slug.init','sim_charLimiting',sim_charLimiting)
  call read_initFileBool('slug.init','sim_RK',sim_RK)
  call read_initFileBool('slug.init','sim_fixDt',sim_fixDt)
  call read_initFileBool('slug.init','sim_nlim',sim_nlim)
  sim_GPFlux = .false.
  sim_DongwookFlux = .false.
  
  call read_initFileBool('slug.init','sim_hdf5',sim_hdf5)

  sim_reconMultiD = .false.
  sim_intFlux = .true.

  call read_initFileChar('slug.init','sim_icType',   sim_icType)
  call read_initFileReal('slug.init','sim_densL',    sim_densL)
  call read_initFileReal('slug.init','sim_velxL',    sim_velxL)
  call read_initFileReal('slug.init','sim_presL',    sim_presL)
  call read_initFileReal('slug.init','sim_densR',    sim_densR)
  call read_initFileReal('slug.init','sim_velxR',    sim_velxR)
  call read_initFileReal('slug.init','sim_presR',    sim_presR)
  call read_initFileReal('slug.init','sim_gamma',    sim_gamma)
  call read_initFileReal('slug.init','sim_shockLoc', sim_shockLoc)
  call read_initFileReal('slug.init','sim_smallPres', sim_smallPres)
  call read_initFileReal('slug.init','sim_x0',    sim_x0)
  call read_initFileReal('slug.init','sim_y0',    sim_y0)

  call read_initFileReal('slug.init','sim_dens1',    sim_Q1(DENS_VAR))
  call read_initFileReal('slug.init','sim_velx1',    sim_Q1(VELX_VAR))
  call read_initFileReal('slug.init','sim_vely1',    sim_Q1(VELY_VAR))
  call read_initFileReal('slug.init','sim_pres1',    sim_Q1(PRES_VAR))

  call read_initFileReal('slug.init','sim_dens2',    sim_Q2(DENS_VAR))
  call read_initFileReal('slug.init','sim_velx2',    sim_Q2(VELX_VAR))
  call read_initFileReal('slug.init','sim_vely2',    sim_Q2(VELY_VAR))
  call read_initFileReal('slug.init','sim_pres2',    sim_Q2(PRES_VAR))

  call read_initFileReal('slug.init','sim_dens3',    sim_Q3(DENS_VAR))
  call read_initFileReal('slug.init','sim_velx3',    sim_Q3(VELX_VAR))
  call read_initFileReal('slug.init','sim_vely3',    sim_Q3(VELY_VAR))
  call read_initFileReal('slug.init','sim_pres3',    sim_Q3(PRES_VAR))

  call read_initFileReal('slug.init','sim_dens4',    sim_Q4(DENS_VAR))
  call read_initFileReal('slug.init','sim_velx4',    sim_Q4(VELX_VAR))
  call read_initFileReal('slug.init','sim_vely4',    sim_Q4(VELY_VAR))
  call read_initFileReal('slug.init','sim_pres4',    sim_Q4(PRES_VAR))

  call read_initFileChar('slug.init','sim_bcTypeX',sim_bcTypeX)
  call read_initFileChar('slug.init','sim_bcTypeY',sim_bcTypeY)
  call read_initFileChar('slug.init','sim_bcTypeZ',sim_bcTypeZ)

  call read_initFileReal('slug.init','sim_ioTfreq',  sim_ioTfreq)
  call read_initFileInt ('slug.init','sim_ioNfreq',  sim_ioNfreq)
  call read_initFileInt ('slug.init','sim_mval',  sim_mval)

!!$  call read_initFileChar('slug.init','sim_quad'     ,  sim_quad)
!!$  call read_initFileChar('slug.init','sim_gp_kernel',  sim_gp_kernel)
!!$  call read_initFileReal('slug.init','sim_sigma'    ,  sim_sigma)
!!$  call read_initFileReal('slug.init','sim_sigdel'   ,  sim_sigdel)
!!$  call read_initFileReal('slug.init','sim_matern_nu',  sim_matern_nu)
!!$  call read_initFileReal('slug.init','sim_RQ_alpha',  sim_RQ_alpha)


  !read grid data
  allocate(gr_xend(NDIM)); allocate(gr_xbeg(NDIM))
  call read_initFileInt ('slug.init','gr_nx',   gr_nx)
  call read_initFileInt ('slug.init','gr_ny',   gr_ny)
  call read_initFileInt ('slug.init','gr_nz',   gr_nz)
  call read_initFileInt ('slug.init','gr_ngc',  gr_ngc)
  call read_initFileReal('slug.init','gr_xbeg', gr_xbeg(XDIM))
  call read_initFileReal('slug.init','gr_xend', gr_xend(XDIM))
  call read_initFileReal('slug.init','gr_ybeg', gr_xbeg(YDIM))
  call read_initFileReal('slug.init','gr_yend', gr_xend(YDIM))
  call read_initFileReal('slug.init','gr_zbeg', gr_xbeg(ZDIM))
  call read_initFileReal('slug.init','gr_zend', gr_xend(ZDIM))

  call read_initFileReal('slug.init','gr_radius',   gr_radius)

  !read block data
  call read_initFileInt ('slug.init','bl_iProcs',   bl_iProcs)
  call read_initFileInt ('slug.init','bl_jProcs',   bl_jProcs)
  call read_initFileInt ('slug.init','bl_kProcs',   bl_kProcs)

  !read in GP pars
  gpM_radius = 0.0
  call read_initFileChar('slug.init','gp_quad'  ,  gp_quad  )
  call read_initFileChar('slug.init','gp_kernel',  gp_kernel)
  call read_initFileReal('slug.init','gp_ell'    ,  gp_el    )
  call read_initFileReal('slug.init','gp_eldel' ,  gp_eldel )
  call read_initFileInt ('slug.init','gp_radius',  gp_radius)
  call read_initFileReal('slug.init','gpM_radius'    ,  gpM_radius)
end subroutine read_pars
