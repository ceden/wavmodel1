
subroutine diagnose()
 use main_module
 implicit none 
 integer :: i,k
 real(real_type) :: umax,wmax,cfl_h,cfl_v

 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe
     umax = max(umax,abs(u(i,k)))
     wmax = max(wmax,abs(w(i,k)))
  enddo
 enddo
 call global_max(umax)
 call global_max(wmax)

 if (umax /= umax) then
   if (my_pe==0) print*,' ERROR: in diagnose detected NaN, aborting at t=',time,' itt=',itt
   call halt_stop('in diagnose')
 endif
 
 cfl_h = umax*dt/dx
 cfl_v = wmax*dt/dz
 if (my_pe==0) print'(a,f8.2,a,i8,a,i5,a,f5.2,a,f5.2)',' diagnosing at t=',time,' itt=',itt, &
            ' solver iterations = ',itts,' max hor. CFL = ', cfl_h,' max ver. CFL = ', cfl_v 
 
end subroutine diagnose



#ifdef with_netcdf_parallel 

 
 subroutine init_snap_cdf

 use main_module
 implicit none
 
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,id,xdim,zdim,tdim,nc_mode,i,k
 character :: fname*80,name*32
 real(real_type), parameter :: spval = -1.0d33
 real(real_type) :: x(nx),z(nz)
  
 fname = 'snap.cdf' 
 if (my_pe==0)  print*,'preparing file ',fname(1:len_trim(fname))

 nc_mode = IOR(nf_clobber,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_CREATE_par(fname,nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)

 iret = nf_def_dim(ncid,'x', nx,xdim)
 iret = nf_def_dim(ncid,'z', nz,zdim)  
 iret = nf_def_dim(ncid,'Time', nf_unlimited,tdim)
 
 iret = nf_def_var(ncid, 'x',NF_DOUBLE,1,xdim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('x'),'x')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m')
  
 iret = nf_def_var(ncid, 'z',NF_DOUBLE,1,zdim,id) 
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('z'),'z')
 iret = nf_put_att_text(ncid,id,'units',len_trim('m'),'m') 

 iret = nf_def_var(ncid, 'Time',NF_DOUBLE,1,tdim,id)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('time'),'time')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('seconds'),'seconds')

 id  = ncvdef (ncid,'u',NF_DOUBLE,3,(/xdim, zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
 id  = ncvdef (ncid,'v',NF_DOUBLE,3,(/xdim, zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
 id  = ncvdef (ncid,'w',NF_DOUBLE,3,(/xdim, zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  
 id  = ncvdef (ncid,'phi',NF_DOUBLE,3,(/xdim, zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
  
 id  = ncvdef (ncid,'p',NF_DOUBLE,3,(/xdim, zdim,tdim/),iret)
 iret = nf_put_att_text(ncid,id,'long_name',len_trim('pressure'),'pressure')
 iret = nf_put_att_text(ncid,id,'unit',len_trim('m^2/s^2'),'m^2/s^2')
 
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'rho_o',NF_DOUBLE,1,rho_o) 
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'rho_a',NF_DOUBLE,1,rho_a)
 iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'g',NF_DOUBLE,1,g)

 iret= nf_enddef(ncid)
 
 do i=1,nx
   x(i)=(i-1)*dx + dx/2.
 enddo     
 iret=nf_inq_varid(ncid,'x',id)
 iret= nf_put_vara_double(ncid,id,is_pe,ie_pe-is_pe+1,x(is_pe:ie_pe)) 

 do k=1,nz
   z(k)=(k-1)*dz + dz/2. - Lz/2.
 enddo   
 iret=nf_inq_varid(ncid,'z',id)
 iret= nf_put_vara_double(ncid,id,ks_pe,ke_pe-ks_pe+1,z(ks_pe:ke_pe)) 

 iret= nf_close(ncid) 
 if (my_pe==0) print*,' done' 
end subroutine init_snap_cdf
 
subroutine diag_snap
 use main_module
 implicit none 
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,id,nc_mode,start(3),count(3),ilen,k
 real(real_type) :: aloc(is_pe:ie_pe,ks_pe:ke_pe)
 character :: fname*80

 fname = 'snap.cdf' 
 if (my_pe==0)  print*,'writing to file ',fname(1:len_trim(fname)),' at itt=',itt

 nc_mode = IOR(nf_write,nf_mpiio)
 nc_mode = IOR(nc_mode,nf_netcdf4 )
 nc_mode = IOR(nc_mode,nf_classic_model )
 
 iret = NF_open_par(fname,nc_mode,MPI_COMM_WORLD,MPI_INFO_NULL, ncid) 
 if (my_pe==0.and.iret/=0) print*,'at opening:',nf_strerror(iret)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)
 
 iret=nf_inq_dimid(ncid,'Time',id)
 iret=nf_inq_dimlen(ncid,id,ilen)
 ilen=ilen+1
 
 iret=nf_inq_varid(ncid,'Time',id) 
 iret = nf_var_par_access(ncid, id, nf_collective)
 iret= nf_put_vara_double(ncid,id,ilen,1,time) 
 
 start = (/is_pe,        ks_pe        ,ilen/)
 count = (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/)
 
 iret=nf_inq_varid(ncid,'u',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = u(is_pe:ie_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'v',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = v(is_pe:ie_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'w',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = w(is_pe:ie_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'phi',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = phi(is_pe:ie_pe,ks_pe:ke_pe)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)

 iret=nf_inq_varid(ncid,'p',id)
 iret = nf_var_par_access(ncid, id, nf_collective)
 aloc = p(is_pe:ie_pe,ks_pe:ke_pe,tau)
 iret= nf_put_vara_double(ncid,id,start,count,aloc)
      
 iret= nf_close(ncid) 

 end subroutine diag_snap
  
  

#else 


subroutine init_snap_cdf
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,i,xdim,zdim,timeid,timedim,id 
 real(real_type) :: x(nx),z(nz)

 call barrier()
 if (my_pe==0) then
 
  print*,'creating file snap.cdf'
  !iret = NF_CREATE ('snap.cdf',IOR(nf_write,nf_classic_model ),ncid) 
  !iret = NF_CREATE ('snap.cdf',IOR(nf_write,NF_NETCDF4  ),ncid) 
  iret = NF_CREATE ('snap.cdf',or(NF_CLASSIC_MODEL,nf_64bit_offset ),ncid)
  if (iret/=0) print*,'at opening:',nf_strerror(iret)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xdim  = ncddef(ncid, 'x', nx, iret)
  zdim  = ncddef(ncid, 'z', nz, iret)  
  Timedim = ncddef(ncid,'Time', nf_unlimited, iret)
  timeid  = ncvdef (ncid,'Time', NF_DOUBLE,1,(/timedim/),iret)
  iret = nf_put_att_text(ncid,timeid,'long_name',len_trim('time'),'time')
  iret = nf_put_att_text(ncid,timeid,'unit',len_trim('seconds'),'seconds')
  
  id  = ncvdef (ncid,'x',NF_DOUBLE,1,(/xdim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('x'),'x')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m'),'m')
  id  = ncvdef (ncid,'z',NF_DOUBLE,1,(/zdim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('z'),'z')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m'),'m')
  id  = ncvdef (ncid,'u',NF_DOUBLE,3,(/xdim, zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  id  = ncvdef (ncid,'v',NF_DOUBLE,3,(/xdim, zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  id  = ncvdef (ncid,'w',NF_DOUBLE,3,(/xdim, zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('velocity'),'velocity')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s'),'m/s')
  id  = ncvdef (ncid,'phi',NF_DOUBLE,3,(/xdim,zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('buoyancy'),'buoyancy')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m/s^2'),'m/s^2')
  id  = ncvdef (ncid,'p',NF_DOUBLE,3,(/xdim, zdim,Timedim/),iret)
  iret = nf_put_att_text(ncid,id,'long_name',len_trim('pressure'),'pressure')
  iret = nf_put_att_text(ncid,id,'unit',len_trim('m^2/s^2'),'m^2/s^2')
 
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'rho_o',NF_DOUBLE,1,rho_o) 
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'rho_a',NF_DOUBLE,1,rho_a)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'g',NF_DOUBLE,1,g)
  iret= nf_enddef(ncid) 
  
  do i=1,nx
   x(i)=(i-1)*dx + dx/2.
  enddo 
  iret=nf_inq_varid(ncid,'x',id)
  iret= nf_put_vara_double(ncid,id,(/1/),(/nx/),x)  
  
  do i=1,nz
   z(i)=(i-1)*dz + dz/2. - Lz/2.
  enddo 
  iret=nf_inq_varid(ncid,'z',id)
  iret= nf_put_vara_double(ncid,id,(/1/),(/nz/),z)   
  iret= nf_close(ncid) 

  print*,'done creating file snap.cdf'
 endif
 call barrier
end subroutine init_snap_cdf


subroutine diag_snap
 use main_module
 
 implicit none
 include "netcdf.inc"
 include "mpif.h"
 integer :: ncid,iret,n
 integer :: tdimid,ilen,timeid,id
 integer :: tag=1,ist(3),isz(3),ien(3)
 integer, dimension(MPI_STATUS_SIZE) :: Status
 real(real_type), allocatable :: a(:,:)

 if (my_pe==0) then
   print*,'writing to file snap.cdf at t=',time
   iret=nf_open('snap.cdf',NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',tdimid)
   iret=nf_inq_dimlen(ncid, tdimid,ilen)
   iret=nf_inq_varid(ncid,'Time',timeid)
   ilen=ilen+1
   iret= nf_put_vara_double(ncid,timeid,(/ilen/),(/1/),(/time/))
   
   iret=nf_inq_varid(ncid,'u',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,ks_pe,ilen/), &
                            (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/),u(is_pe:ie_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'v',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,ks_pe,ilen/), &
                            (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/),u(is_pe:ie_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'w',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,ks_pe,ilen/), &
                            (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/),w(is_pe:ie_pe,ks_pe:ke_pe))
   iret=nf_inq_varid(ncid,'phi',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,ks_pe,ilen/), &
                            (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/),phi(is_pe:ie_pe,ks_pe:ke_pe))   
   iret=nf_inq_varid(ncid,'p',id)
   iret= nf_put_vara_double(ncid,id,(/is_pe,ks_pe,ilen/), &
                            (/ie_pe-is_pe+1,ke_pe-ks_pe+1,1/),p(is_pe:ie_pe,ks_pe:ke_pe,tau)) 
                                                                                    
 endif

 do n=1,n_pes-1
  call mpi_barrier(MPI_COMM_WORLD, iret)
  if (my_pe==n) then
        call mpi_send((/is_pe,ks_pe/),2,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send((/ie_pe,ke_pe/),2,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send((/ie_pe-is_pe+1,ke_pe-ks_pe+1/),2,mpi_integer,0,tag,mpi_comm_world,iret)
        call mpi_send(u(is_pe:ie_pe,ks_pe:ke_pe), &               
          (ie_pe-is_pe+1)*(ke_pe-ks_pe+1),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(v(is_pe:ie_pe,ks_pe:ke_pe), &               
          (ie_pe-is_pe+1)*(ke_pe-ks_pe+1),mpi_real8,0,tag,mpi_comm_world,iret)          
        call mpi_send(w(is_pe:ie_pe,ks_pe:ke_pe), &
          (ie_pe-is_pe+1)*(ke_pe-ks_pe+1),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(phi(is_pe:ie_pe,ks_pe:ke_pe), &
          (ie_pe-is_pe+1)*(ke_pe-ks_pe+1),mpi_real8,0,tag,mpi_comm_world,iret)
        call mpi_send(p(is_pe:ie_pe,ks_pe:ke_pe,tau), &          
          (ie_pe-is_pe+1)*(ke_pe-ks_pe+1),mpi_real8,0,tag,mpi_comm_world,iret)
     
  else if (my_pe==0) then
        call mpi_recv(ist,2,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(ien,2,mpi_integer,n,tag,mpi_comm_world,status,iret)
        call mpi_recv(isz,2,mpi_integer,n,tag,mpi_comm_world,status,iret)
        allocate(a(ist(1):ien(1),ist(2):ien(2)) )
        
        call mpi_recv(a,isz(1)*isz(2),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'u',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ilen/),(/isz(1),isz(2),1/),a)
               
        call mpi_recv(a,isz(1)*isz(2),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'v',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ilen/),(/isz(1),isz(2),1/),a)
               
        call mpi_recv(a,isz(1)*isz(2),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'w',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ilen/),(/isz(1),isz(2),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'phi',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ilen/),(/isz(1),isz(2),1/),a)
        
        call mpi_recv(a,isz(1)*isz(2),mpi_real8,n,tag,mpi_comm_world,status,iret)
        iret=nf_inq_varid(ncid,'p',id)
        iret= nf_put_vara_double(ncid,id,(/ist(1),ist(2),ilen/),(/isz(1),isz(2),1/),a)
              
        deallocate(a)
  endif
 enddo  
 if (my_pe == 0) iret= nf_close(ncid) 
end subroutine diag_snap

#endif