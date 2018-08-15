subroutine grid_finalize()

  use grid_data!, only : gr_xCoord, gr_U, gr_V, gr_W, gr_eigval, gr_leigvc, gr_reigvc

  implicit none

  if (allocated(gr_xCoord) .eqv. .true.) deallocate(gr_xCoord)
  if (allocated(gr_yCoord) .eqv. .true.) deallocate(gr_yCoord)
  if (allocated(gr_i0) .eqv. .true.) deallocate(gr_i0)
  if (allocated(gr_imax) .eqv. .true.) deallocate(gr_imax)
  if (allocated(gr_ibeg) .eqv. .true.) deallocate(gr_ibeg)
  if (allocated(gr_iend) .eqv. .true.) deallocate(gr_iend)
  if (allocated(gr_U) .eqv. .true.) deallocate(gr_U)
  if (allocated(gr_V) .eqv. .true.) deallocate(gr_V)
  if (allocated(gr_W) .eqv. .true.) deallocate(gr_W)

  if (allocated(gr_vR) .eqv. .true.) deallocate(gr_vR) 
  if (allocated(gr_vL) .eqv. .true.) deallocate(gr_vL)
  !if (allocated(gr_vP) .eqv. .true.) deallocate(gr_vP)
  !if (allocated(gr_vM) .eqv. .true.) deallocate(gr_vM)
  if (allocated(gr_flux) .eqv. .true.) deallocate(gr_flux)

!!$  if (allocated(gr_eigval) .eqv. .true.) deallocate(gr_eigval)
!!$  if (allocated(gr_leigvc) .eqv. .true.) deallocate(gr_leigvc)
!!$  if (allocated(gr_reigvc) .eqv. .true.) deallocate(gr_reigvc)

  if (allocated(gr_GPv) .eqv. .true.) deallocate(gr_GPv)
  if (allocated(gr_GPZ) .eqv. .true.) deallocate(gr_GPZ)
  if (allocated(gr_GP_stencil) .eqv. .true.) deallocate(gr_GP_stencil)


  
  return
end subroutine grid_finalize
