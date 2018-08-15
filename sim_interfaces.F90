module sim_interfaces
#include "definition.h"
  abstract interface
     subroutine recon1D(dt, Nx, stencil, vL, vR, dir)
       
       integer, intent(IN) :: Nx, dir
       real   , intent(IN) :: dt

       real, dimension(NUMB_VAR), intent(INOUT ) :: vL, vR
       real, dimension(Nx, NUMB_VAR), intent(IN) :: stencil
     end subroutine recon1D
  end interface

  abstract interface
     subroutine reconMD(dt, Npts, stencil, vL, vR)

       integer, intent(IN) :: Npts
       real   , intent(IN) :: dt

       real, dimension(Npts    , NUMB_VAR) :: stencil
       real, dimension(NUMB_VAR, NDIM    ) :: vL, vR
     end subroutine reconMD
  end interface

  abstract interface
     subroutine rmn_slvr(vR, vL, flux, dir)

       integer, intent(IN) :: dir
       real, dimension(NUMB_VAR), intent(IN ) :: vL, vR
       real, dimension(NSYS_VAR), intent(OUT) :: flux

     end subroutine rmn_slvr
  end interface

  abstract interface
     subroutine num_flux(F, Npts, Fiph)

     integer, intent(IN) :: Npts
     real, dimension(NSYS_VAR, Npts), intent(IN) :: F
     real, dimension(NSYS_VAR) :: Fiph

   end subroutine num_flux

  end interface
  
end module sim_interfaces
