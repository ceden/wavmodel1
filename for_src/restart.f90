
subroutine write_restart
 use main_module
 implicit none
 integer :: io=25
 character :: fname*80 

 write(fname,'(a,i9,a)')  'restart_pe=',my_pe,'.dta'
 call replace_space_zero(fname)
 
 if (my_pe==0) print'(a,a)',' writing to restart file ',fname(1:len_trim(fname))
 if (my_pe==0) print'(a,i12)',' at itt = ',itt
  
 open(io,file=fname,form='unformatted',status='unknown')
 write(io,err=10) nx,nz,itt,taum2,taum1,tau,time
 write(io,err=10) is_pe,ie_pe,ks_pe,ke_pe
 write(io,err=10) u,v,w,phi,p,du,dv,dw,dphi
 close(io)
 return
 10 continue
 print'(a)',' Warning: error writing restart file'
end subroutine write_restart


subroutine read_restart
 use main_module
 implicit none
 include 'mpif.h'
 character :: fname*80 
 logical :: file_exists
 integer :: io=25,nx_,nz_,ierr
 integer :: is_,ie_,ks_,ke_
 
 write(fname,'(a,i9,a)')  'restart_pe=',my_pe,'.dta'
 call replace_space_zero(fname)
 inquire ( FILE=fname, EXIST=file_exists )
 if (.not. file_exists) then
  if (my_pe==0) then
         print'(a,a)',' found no restart file ',fname(1:len_trim(fname))
         print'(a)',' proceeding with initial conditions'
  endif
  return
 endif
 if (my_pe==0) print'(a,a)',' reading restart file ',fname(1:len_trim(fname))
 
 open(io,file=fname,form='unformatted',status='old',err=10)
 read(io,err=10) nx_,nz_,itt,taum2,taum1,tau,time

 if (my_pe==0) print*,' time = ',time,' itt = ',itt
 runlen = runlen + time
 
 if (my_pe==0) print*,' adjusting run length to = ',runlen

 if (nx/=nx_ .or. nz/= nz_) then 
       if (my_pe==0) then
        print*,' read from restart dimensions: ',nx_,nz_
        print*,' does not match dimensions   : ',nx,nz
       endif
       call MPI_ABORT(mpi_comm_world, 99, IERR)
 endif
  
 read(io,err=10) is_,ie_,ks_,ke_

 if (is_/=is_pe.or.ie_/=ie_pe.or.ks_/=ks_pe.or.ke_/=ke_pe) then
       if (my_pe==0) then
        print*,' read from restart PE boundaries: ',is_,ie_,ks_,ke_
        print*,' which does not match           : ',is_pe,ie_pe,ks_pe,ke_pe
       endif
       call MPI_ABORT(mpi_comm_world, 99, IERR)
 endif

 read(io,err=10) u,v,w,phi,p,du,dv,dw,dphi  
 close(io)
 return
 10 continue
 print'(a)',' Warning: error reading restart file'
 call MPI_ABORT(mpi_comm_world, 99, IERR)
end subroutine read_restart



subroutine replace_space_zero(name)
      implicit none
      character (len=*) :: name
      integer  :: i
      do i=1,len_trim(name)
          if (name(i:i)==' ')name(i:i)='0'
      enddo
end subroutine replace_space_zero


