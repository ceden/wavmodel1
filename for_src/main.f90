
program main
 use main_module
 use timing_module   
 implicit none
 include "mpif.h"
 integer :: ierr,iargc,n
 character (len=80) :: arg

 !----------------------------------------------------------------------
 ! initialize parallelisation with MPI
 !----------------------------------------------------------------------
 call mpi_init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, my_pe, ierr)
 call MPI_Comm_size(MPI_COMM_WORLD, n_pes, ierr)
 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,' '
 if (my_pe==0)  print *,' '
 
 if (n_pes>1) then
   if (iargc() < 2)  then
      if (my_pe==0) print*,'ERROR: not enough command line input'
      call MPI_ABORT(mpi_comm_world, 99, IERR)
   endif
   call getarg(1,arg); read(arg,*) n_pes_i
   call getarg(2,arg); read(arg,*) n_pes_k
 else 
  n_pes_i = 1; n_pes_k = 1
 endif
 
 if (n_pes_i*n_pes_k /= n_pes) then
   if (my_pe==0) print*,' ERROR: n_pes_i times n_pes_k not equal number of PEs'
   call MPI_ABORT(MPI_COMM_WORLD, 99, ierr)
 endif
 if (my_pe==0)  print '(a,i4,a,i4)','Using processor grid ',n_pes_i,' x ',n_pes_k

 !----------------------------------------------------------------------
 ! setup the model
 !----------------------------------------------------------------------
 call tic('setup')
 if (my_pe==0)  print *,' setting parameters '
 call set_parameter()
 call pe_decomposition
 
 !----------------------------------------------------------------------
 ! allocate memory
 !----------------------------------------------------------------------
 if (my_pe==0) print *,' allocate '
 call allocate_main_module()
 
 if (my_pe==0) print *,' initialization of diagnostic output  '
 call init_snap_cdf
     
 if (my_pe==0) then
     print*,' '
     print'(a,i4,a,i4)',' nx x nz = ',nx,' x ',nz
     print'(a,f8.2,a,f8.2,a)',' Lx x Lz = ',Lx,' x ',Lz,' m'
     print'(a,f8.5,a)',' Delta t = ',dt,' s'
     print'(a,f8.5,a)',' Delta x = ',dx,' m'
     print'(a,f8.5,a)',' Delta z = ',dz,' m'
     print'(a,f8.5,a)',' runlen  = ',runlen,' s'
     print*,' '
 endif
 
 time = 0.; itt = 0
 taum2 = 1; taum1 = 2; tau = 3

 if (my_pe==0)  print *,' setting initial conditions '
 call set_initial_conditions
 
 if (my_pe==0)  print *,' done setting initial conditions '
 call toc('setup')

 if (my_pe==0)  print *,' starting main loop '

  
 !----------------------------------------------------------------------
 ! reading restart
 !----------------------------------------------------------------------
 call read_restart 
 if (enable_particles) call particles_read_restart
 !----------------------------------------------------------------------
 ! start main integration loop
 !----------------------------------------------------------------------

 call tic('loop')   
 do while (time < runlen) 

   call integrate
     
   call tic('diagnostics')     
   if ( mod(itt,max(1,int(tsmonint/dt)))  == 0 .or. itt == 0)  call diagnose  
   if ( mod(itt,max(1,int( snapint/dt)))  == 0 .or. itt == 0)  call diag_snap 
   if (enable_particles) then
    call integrate_particles 
    if ( mod(itt,max(1,int( snapint/dt)))  == 0 .or. itt == 0)  call write_particles  
    if ( mod(itt,max(1,int( meanint/dt)))  == 0 .and. itt > 0)  call write_mean_particles  
   endif 
   call toc('diagnostics')
   
   tau   = mod(tau,3)+1
   taum1 = mod(taum1,3)+1
   taum2 = mod(taum2,3)+1
   time = time + dt       
   itt = itt + 1 
   
 enddo
 call toc('loop') 

 call write_restart
 if (enable_particles) call particles_write_restart
 
!--------------------------------------------------------------
!     show timing results here
!--------------------------------------------------------------
 !do n = 0,n_pes
  n=0
     call mpi_barrier(MPI_COMM_WORLD, ierr)
     if (my_pe == n) then
        print'(/,a,i4)','Timing summary for PE #',my_pe 
        print'(a,f15.2,a)',' costs for measuring      = ',timing_secs('tictoc'),' s'
        print'(a,f15.2,a)',' setup time summary       = ',timing_secs('setup'),' s'
        print'(a,f15.2,a)',' loop                     = ',timing_secs('loop'),' s'
        print'(a,f15.2,a)',' pressure                 = ',timing_secs('pressure'),' s'
        print'(a,f15.2,a)',' boundary exchange        = ',timing_secs('boundary'),' s'
        print'(a,f15.2,a)',' diagnostics              = ',timing_secs('diagnostics'),' s'
       
       endif
  !enddo
  
 if (my_pe==0) then
  open(10,file='ritt',form='formatted',status='unknown')
  write(10,*) itt
  close(10)
 endif
 
 !----------------------------------------------------------------------
 ! cancel parallelisation and quit
 !----------------------------------------------------------------------
 if (my_pe==0) print'(/a/)','cancelling MPI service'
 call mpi_finalize(ierr)
 
end program main




