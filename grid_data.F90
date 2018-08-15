module grid_data
  implicit none
  real, allocatable, dimension(:), save   :: gr_xCoord, gr_yCoord, gr_zCoord
  
  real, save                            :: gr_dx, gr_dy, gr_dz
  real, allocatable, dimension(:), save ::gr_xbeg, gr_xend
  
  integer, save                            :: gr_nx , gr_ny, gr_nz, gr_ngc, gr_glb_nx, gr_glb_ny, gr_glb_nz
  integer, allocatable, dimension(:), save ::gr_i0, gr_ibeg, gr_iend, gr_imax
  
  real, allocatable, dimension(:,:,:,:)   :: gr_U ! conservative vars
  real, allocatable, dimension(:,:,:,:)   :: gr_V ! primitive vars
  real, allocatable, dimension(:,:,:,:)   :: gr_W ! characteristic vars

  real, allocatable, dimension(:,:,:,:,:)     :: gr_vL   ! left Riemann states
  real, allocatable, dimension(:,:,:,:,:)     :: gr_vR   ! right Riemann states
  real, allocatable, dimension(:,:,:,:,:)     :: gr_flux ! fluxes
  real, allocatable, dimension(:,:,:,:,:) :: gr_ptFluxes !ptwise fluxes at grid points

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !I think these are all only used in the FDM flux formulation and shouldn't be used in this implementation
!!$  real, allocatable, dimension(:,:,:,:,:) :: gr_vP   ! Plus char fluxes 
!!$  real, allocatable, dimension(:,:,:,:,:) :: gr_vM   ! Minus char fluxes
!!$
!!$  real, allocatable, dimension(:,:)     :: gr_eigval ! eigenvalues
!!$  real, allocatable, dimension(:,:,:)   :: gr_reigvc ! right eigenvectors
!!$  real, allocatable, dimension(:,:,:)   :: gr_leigvc ! left  eigenvectors
  real, allocatable, dimension(:,:)     :: gr_maxalphas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  !GP vars
  real, save                               :: gr_radius
  integer, save                            :: gr_gp_stencilPts, gr_Tcells
  real, allocatable, dimension(:  ), save  :: gr_GPv
  real, allocatable, dimension(:,:), save  :: gr_GPZ
  real, allocatable, dimension(:, :), save :: gr_GP_stencil

end module grid_data
