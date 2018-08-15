subroutine block_finalize()
  use block_data

  implicit none

  if (allocated(bl_buffL) .eqv. .true.) deallocate(bl_buffL)
  if (allocated(bl_buffR) .eqv. .true.) deallocate(bl_buffR)
  if (allocated(bl_buffT) .eqv. .true.) deallocate(bl_buffT)
  if (allocated(bl_buffB) .eqv. .true.) deallocate(bl_buffB)

  if (allocated(bl_buffTR) .eqv. .true.) deallocate(bl_buffTR)
  if (allocated(bl_buffBR) .eqv. .true.) deallocate(bl_buffBR)
  if (allocated(bl_buffTL) .eqv. .true.) deallocate(bl_buffTL)
  if (allocated(bl_buffBL) .eqv. .true.) deallocate(bl_buffBL)

  if (allocated(bl_bufCL) .eqv. .true.) deallocate(bl_bufCL)
  if (allocated(bl_bufCR) .eqv. .true.) deallocate(bl_bufCR)
  if (allocated(bl_bufCT) .eqv. .true.) deallocate(bl_bufCT)
  if (allocated(bl_bufCB) .eqv. .true.) deallocate(bl_bufCB)

  if (allocated(bl_delT) .eqv. .true.) deallocate(bl_delT)

  return
end subroutine block_finalize
