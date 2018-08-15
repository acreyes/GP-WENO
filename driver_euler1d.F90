program driver_euler1d

#include "definition.h"
  
  use sim_data
  use grid_data
  use block_data, only : bl_ID
  use gp_data, only : gp_radius
  use io
  use bc
  use eos, only : eos_all

  implicit none

  real :: t,dt, dt_small, dt_init
  integer :: nStep,ioCounter,ioTimeFreqCounter
  real :: ioCheckTime
  
  t = 0.
  sim_time = t
  sim_bcT = t
  nStep = 0
  ioCounter = 0
  ioTimeFreqCounter = 0
  dt_init = 1.e2
  
  
  ! grid_init should be called first before sim_init
  call read_pars()
  call block_init()
  call grid_init()
  call sim_init()

  if (sim_order == 10 .and. .not. sim_reconMultiD) then
     call gp_WENOinit()
     call gp_eigens()
     if (sim_gpFlux) then
        call gp_Fluxinit
     end if
  elseif (sim_reconMultiD .and. sim_order == 10) then
     call gp_MakeStencil
     call gp_MDinit
     call gpM_eigens
     gp_radius = 2
     call gp_WENOinit
     call gp_eigens
  end if
  
  if (bl_id == 1) then
     write(*,*)''
     write(*,*)'================================================='
     write(*,*)'     AMS 260 - CFD, 1D Euler FVM Code            '
     write(*,*)'      Written by Prof. Dongwook Lee              '
     write(*,*)'           Winter Quarter, 2015                  '
     write(*,*)'================================================='
     write(*,*)''


     ! write the initial condition
     write(*,*)''
     write(*,*)'       Initial condition was written!            '
     write(*,*)'================================================='
     write(*,*)'   Steps      Time              dt               '
     write(*,*)'================================================='
     write(*,*)''
  end if

  call io_writeOutput(t, nStep,ioCounter)
  
  do while ( (t < sim_tmax))
     !print *, "1"
     if (sim_fixDt) then
        dt = sim_dt
     else
        call cfl(dt)
     end if

     dt_small = dt_init*2.**nstep
     if (dt_small < dt) dt = dt_small
     
     !check to see if there is a reason to stop
     if (  sim_nlim .and. (nStep .ge. sim_nstep)) then
        exit
     elseif ( abs(t - sim_tmax) .le. dt ) then
        dt = abs(t - sim_tmax)
        !exit
     end if
     !print *, "2"
     sim_time = t
     sim_bcT  = t
     call soln_ReconEvolveAvg(dt)
     !print *, "3"
     call soln_update(dt)
     !print *, "4"

     ! call BC on primitive vars
     ! add_back
     call bc_apply(gr_V)
     !print *, "5"
     
     ! call eos to make sure all through GC regions
     !call eos_all
     !lets remove this because dongwook say its better
     
     ! write outputs every ioNfreq cycle or ioTfreq cycle
     
     ioCheckTime = sim_ioTfreq*real(ioTimeFreqCounter+1)
     if (t-dt < ioCheckTime .and. t>ioCheckTime) then
        if (bl_ID == 1) then
           write(*,*)''
           write(*,*)' Output no.',ioCounter+1, 'has been written      '
           write(*,*)'================================================='
           write(*,*)'   Steps      Time              dt               '
           write(*,*)'================================================='
           write(*,*)''
        end if
        ioCounter = ioCounter + 1
        ioTimeFreqCounter = ioTimeFreqCounter + 1
        call io_writeOutput(t, nStep,ioCounter)
     endif
     !print *, "6"
     if (sim_ioNfreq > 0) then
        if (mod(nStep, sim_ioNfreq) == 0) then
           if (bl_ID == 1) then
              write(*,*)''
              write(*,*)' Output no.',ioCounter+1, 'has been written      '
              write(*,*)'================================================='
              write(*,*)'   Steps      Time              dt               '
              write(*,*)'================================================='
              write(*,*)''
           end if
           ioCounter = ioCounter + 1
           call io_writeOutput(t, nStep,ioCounter)
        endif
     endif
  
     ! update your time and step count
     t = t + dt
     nStep = nStep + 1

     if (bl_ID==1) write(*,900)nstep,t,dt
     if (dt .le. 0.) then
        exit
     end if
  end do


  !! Let's write the final result before exiting
  if (bl_ID == 1) then
     write(*,*)''
     write(*,*)' Final output no.',ioCounter+1, 'has been written'
     write(*,*)'================================================='
     write(*,*)'        The final tmax has reached, bye!         '
     write(*,*)'================================================='
     write(*,*)''
  end if
  call io_writeOutput(t, nStep,ioCounter+1)

  !! finalize and deallocate memories
  call grid_finalize()
  call block_finalize()

900 format(1x,i5,f16.8,1x,f16.8)
  
end program driver_euler1d
