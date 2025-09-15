


subroutine pe_decomposition
  use main_module   
  implicit none
  integer :: n,tag=0,iloc(20),ierr
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status

! ----------------------------------
!      domain decomposition for each PE
! ----------------------------------
   if (n_pes>1) then
      if (n_pes_i*n_pes_k /= n_pes ) &
        call halt_stop(' n_pes_i x n_pes_k not equal number of PEs')

      i_blk = (nx-1)/n_pes_i + 1    ! i-extent of each block
      k_blk = (nz-1)/n_pes_k + 1

      my_blk_i = mod(my_pe,n_pes_i)+1! number of PE in i-dir.
      my_blk_k = (my_pe)/(n_pes_i) + 1 ! number of PE in k-dir.

      is_pe = (my_blk_i-1)*i_blk + 1 ! start index in i-dir of this PE
      ie_pe = min(my_blk_i*i_blk,nx)
      ks_pe = (my_blk_k-1)*k_blk + 1
      ke_pe = min(my_blk_k*k_blk,nz)

! ----------------------------------
!     last block might have been truncated
! ----------------------------------   
      i_blk = ie_pe-is_pe+1      
      k_blk = ke_pe-ks_pe+1      
! ----------------------------------
! ----------------------------------
!     check for incorrect domain decomposition
! ----------------------------------

      if (my_blk_i==n_pes_i .and. is_pe>ie_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in i-direction'
       print*,' choose other number of PEs in i-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_k==n_pes_k .and. ks_pe>ke_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in k-direction'
       print*,' choose other number of PEs in k-direction'
       call halt_stop(' in pe_decomposition')
      endif
   else
       n_pes_i = 1; n_pes_k = 1
       i_blk = nx;  k_blk = nz
       my_blk_i = 1 ; my_blk_k = 1
       ks_pe = 1; ke_pe = nz
       is_pe = 1; ie_pe = nx
   endif
! ----------------------------------
!      print out the PE decomposition, let PE 0 talk
! ----------------------------------
   do n=0,n_pes-1
     if (n==0) then
       iloc(1:4) = (/is_pe,ie_pe,ks_pe,ke_pe/)
     else
       if (my_pe==n) call mpi_send((/is_pe,ie_pe,ks_pe,ke_pe/),4,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)    
       if (my_pe==0) call mpi_recv(iloc,4,mpi_integer,n,tag,MPI_COMM_WORLD,status,ierr)
     endif
     if (my_pe==0) print'(a,i4,a,i4,a,i4,a,i4,a,i4)', &
           'domain of PE #',n,' i=',iloc(1),':',iloc(2),' k=',iloc(3),':',iloc(4)
   enddo

   if (my_pe==0) print*,' '

   do n=0,n_pes-1
     if (n==0) then
       iloc(1:2) = (/my_blk_i,my_blk_k/)
     else
       if (my_pe==n) call mpi_send((/my_blk_i,my_blk_k/),2,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
       if (my_pe==0) call mpi_recv(iloc,2,mpi_integer,n,tag,MPI_COMM_WORLD,status,ierr)
     endif
     if (my_pe==0) print'(a,i4,a,i4,a,i4)', &
          ' pe#',n,' : my_blk_i=',iloc(1) ,' my_blk_k=',iloc(2) 
   enddo

   if (my_pe==0) print*,' '
   
end subroutine pe_decomposition


subroutine barrier
!--------------------------------------------------------------
!--------------------------------------------------------------
      implicit none     
      integer :: ierr
      include "mpif.h"
      call mpi_barrier(MPI_COMM_WORLD,ierr)
end subroutine barrier


subroutine halt_stop(string)
!--------------------------------------------------------------
!     controlled stop, should not be called from python
!--------------------------------------------------------------
      implicit none
      character*(*) :: string
      integer :: ierr,code,my_pe
      include "mpif.h"
      call mpi_comm_rank(MPI_COMM_WORLD,my_pe,ierr)
      print*,' global pe #',my_pe,' : ',string
      print*,' global pe #',my_pe,' aborting '
      code=99
      call MPI_ABORT(mpi_comm_world, code, IERR)
end subroutine halt_stop


subroutine pe0_bcast(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      real(real_type), intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_real8,0,MPI_COMM_WORLD,ierr)
end subroutine pe0_bcast




subroutine global_max(x)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real(real_type),intent(inout)    :: x
      real(real_type)    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MAX,mpi_comm_world,ierr)
      x = x_sym2
 end subroutine global_max



subroutine global_sum(x)
!--------------------------------------------------------------
!     Do a sum of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real(real_type),intent(inout)    :: x
      real(real_type)    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_SUM,mpi_comm_world,ierr)
      x = x_sym2
end subroutine global_sum




subroutine global_max_int(x)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MAX,mpi_comm_world,ierr)
      x = x_sym2
 end subroutine global_max_int
 

subroutine global_min_int(x)
!--------------------------------------------------------------
!     Get the min of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MIN,mpi_comm_world,ierr)
      x = x_sym2
end subroutine global_min_int




 subroutine global_max_int2(x,len)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x(len)
      integer,intent(in) :: len
      integer    :: x_sym(len),x_sym2(len),ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,len,mpi_integer,MPI_MAX,mpi_comm_world,ierr)
      x = x_sym2
 end subroutine global_max_int2




 subroutine global_max2(x,len)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x(len)
      integer,intent(in) :: len
      real*8    :: x_sym(len),x_sym2(len)
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,len,mpi_real8,MPI_MAX,mpi_comm_world,ierr)
      x = x_sym2
 end subroutine global_max2




subroutine border_exchg_3D(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  call border_exchg_in_z(a)
  call border_exchg_in_x(a)
end subroutine border_exchg_3D


subroutine border_exchg_in_x(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! only in x direction, for periodic boundaries
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  integer  ::  tag=0, ierr,i,len,east,west
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm
  
  my_comm = MPI_COMM_WORLD     
  west = my_pe-1
  if (my_blk_i==1      ) west = my_pe+ n_pes_i-1
  east = my_pe+1
  if (my_blk_i==n_pes_i) east = my_pe- (n_pes_i-1)
  
  if ( n_pes_i > 1) then
     len=(ke_pe-ks_pe+1+2*onx)
     do i=1,onx
       call mpi_send(a(is_pe+i-1,:),len,mpi_real8,west,tag,my_comm,ierr)
       call mpi_recv(a(ie_pe+i  ,:),len,mpi_real8,east,tag,my_comm,status,ierr)      
     enddo
     do i=1,onx
       call mpi_send(a(ie_pe-i+1,:),len,mpi_real8,east,tag,my_comm,ierr)
       call mpi_recv(a(is_pe-i  ,:),len,mpi_real8,west,tag,my_comm,status,ierr)
     enddo
  else
    do i=1,onx
     a(nx+i,:) = a(i     ,:)
     a(1-i ,:) = a(nx-i+1,:) 
    enddo
  endif
end subroutine border_exchg_in_x


subroutine border_exchg_in_z(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! only in vertical direction, no periodic boundaries
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(inout)  :: a(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  integer  ::  tag=0, ierr,k,len,top,bottom
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm
  
  my_comm = MPI_COMM_WORLD
  bottom = my_pe - n_pes_i
  top    = my_pe + n_pes_i   
  if (n_pes_k > 1) then
     len=(ie_pe-is_pe+1+2*onx)
     do k=1,onx
       if (my_blk_k /=1 )       call mpi_send(a(:,ks_pe+k-1),len,mpi_real8,bottom,tag,my_comm,ierr)
       if (my_blk_k /= n_pes_k) call mpi_recv(a(:,ke_pe+k  ),len,mpi_real8,top,tag,my_comm,status,ierr)
      enddo
     do k=1,onx
       if (my_blk_k /= n_pes_k) call mpi_send(a(:,ke_pe-k+1),len,mpi_real8,top,tag,my_comm,ierr)
       if (my_blk_k /=1 )       call mpi_recv(a(:,ks_pe-k),len,mpi_real8,bottom,tag,my_comm,status,ierr)
     enddo
  endif
end subroutine border_exchg_in_z



subroutine retrieve_in_z(a,b)
!--------------------------------------------------------------
! PE with my_blk_k == 1 receives full array b with length nz in z
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(in)  :: a(ks_pe:ke_pe)
  real(real_type), intent(out) :: b(nz)
  integer  ::  tag=0, ierr, n , ks,ke,kk(2)
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm, sender, receiver

   my_comm = MPI_COMM_WORLD    
   if (my_blk_k == 1) then
    b(ks_pe:ke_pe) = a(ks_pe:ke_pe) ! local copy
   endif
 
   do n = 2, n_pes_k
     sender   = my_pe + n_pes_i*(n-1)
     receiver = my_pe - n_pes_i*(n-1)
     if (my_blk_k == 1) then
      call mpi_recv(kk,2            ,mpi_real8,sender,tag,my_comm,status,ierr)
      ks=kk(1);ke=kk(2)
      call mpi_recv(b(ks:ke),ke-ks+1,mpi_real8,sender,tag,my_comm,status,ierr)
     else if (my_blk_k == n) then
       kk = (/ks_pe,ke_pe/)
       call mpi_send(kk,2           ,mpi_real8,receiver,tag,my_comm,ierr)
       call mpi_send(a,ke_pe-ks_pe+1,mpi_real8,receiver,tag,my_comm,ierr)
     endif
   enddo 
end subroutine retrieve_in_z   
   
   
subroutine distribute_in_z(a,b) 
!--------------------------------------------------------------
! PE with my_blk_k == 1 distributes full array of length nz in z
!--------------------------------------------------------------
  use main_module   
  implicit none
  real(real_type), intent(in)  :: a(nz)
  real(real_type), intent(out) :: b(ks_pe:ke_pe)
  integer  ::  tag=0, ierr, n , ks,ke,kk(2)
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
  integer :: my_comm, sender, receiver

   my_comm = MPI_COMM_WORLD      
 
   if (my_blk_k == 1) then
    b(ks_pe:ke_pe) = a(ks_pe:ke_pe) ! local copy
   endif
   
   do n = 2, n_pes_k
     sender   = my_pe + n_pes_i*(n-1)
     receiver = my_pe - n_pes_i*(n-1)
     if (my_blk_k == 1) then
      call mpi_recv(kk,2            ,mpi_real8,sender,tag,my_comm,status,ierr)
      ks=kk(1);ke=kk(2)
      call mpi_send(a(ks:ke),ke-ks+1,mpi_real8,sender,tag,my_comm,ierr)      
     else if (my_blk_k == n) then
       kk = (/ks_pe,ke_pe/)
       call mpi_send(kk,2           ,mpi_real8,receiver,tag,my_comm,ierr)
       call mpi_recv(b,ke_pe-ks_pe+1,mpi_real8,receiver,tag,my_comm,status,ierr)
     endif
   enddo 
end subroutine distribute_in_z


