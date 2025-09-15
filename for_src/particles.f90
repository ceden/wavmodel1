

module particles_module
!=======================================================================
! module for particles
!=======================================================================
 implicit none 
 integer :: nptraj  = -1 ! number of particles
 integer :: nr_prop =  5 ! number of properties to store
 real*8, allocatable :: pxz(:,:) ! position of particles
 real*8, allocatable :: puwp(:,:) ! velocity and thickness of particles
 real*8, allocatable :: pxz_m(:,:), puwp_m(:,:) ! mean position, etc
 integer, allocatable :: pik(:,:) ! index of tracer box
 logical, allocatable :: particle_active(:)
 real*8, allocatable  :: xt(:),xu(:),zt(:),zu(:)
 real*8 :: n_sum = 0.
end module particles_module


subroutine allocate_particles(n_in)
!=======================================================================
! allocate variables for particles
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 integer, intent(in) ::  n_in
 nptraj = n_in
 allocate( pik(2,nptraj) ); pik = 0
 allocate( pxz(2,nptraj) ); pxz = 0.0
 allocate( puwp(nr_prop,nptraj) ); puwp = 0.0
 allocate( pxz_m(2,nptraj) ); pxz_m = 0.0
 allocate( puwp_m(nr_prop,nptraj) ); puwp_m = 0.0
 allocate( particle_active( nptraj) ); particle_active = .false.
 allocate( xt(is_pe-2:ie_pe+2), xu(is_pe-2:ie_pe+2) ); xt=0;xu=0
 allocate( zt(ks_pe-2:ke_pe+2), zu(ks_pe-2:ke_pe+2) ); zt=0;zu=0
end subroutine allocate_particles



subroutine init_particles
!=======================================================================
! initialize everything, user defined function to seed particles is
! already called using allocate_particles
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 integer :: i,k,n

 do i=is_pe-2,ie_pe+2
  xt(i) = (i-1)*dx+dx/2
  xu(i) = i*dx
 enddo
 
 do k=ks_pe-2,ke_pe+2
  zt(k) = (k-1)*dz+dz/2
  zu(k) = k*dz
 enddo
 
 ! find positions of particles
 do n=1,nptraj
  if ( pxz(1,n)>xu(is_pe-1) .and. pxz(1,n)<=xu(ie_pe) .and. &
       pxz(2,n)>zu(ks_pe-1) .and. pxz(2,n)<=zu(ke_pe) ) then
   pik(1,n) = minloc(  (pxz(1,n) - xt(is_pe:ie_pe) )**2  ,1)+is_pe-1
   pik(2,n) = minloc(  (pxz(2,n) - zt(ks_pe:ke_pe) )**2  ,1)+ks_pe-1
  else
   pik(:,n) = -99
  endif 
 enddo
 pxz_m=0; puwp_m=0; n_sum = 0
 
 call particle_pe_domain()
 call particle_distribute()
 call init_write_particles()
 call init_write_mean_particles()
end subroutine init_particles




subroutine particle_pe_domain()
!=======================================================================
!  is particle inside domain of this pe?
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 integer :: n
 do n=1,nptraj
  if (pik(1,n)>=is_pe .and. pik(1,n)<=ie_pe .and. pik(2,n)>=ks_pe .and. pik(2,n)<=ke_pe ) then
        particle_active(n) = .true. 
  else
        particle_active(n) = .false.
  endif
 enddo
end subroutine particle_pe_domain


subroutine particle_distribute
!=======================================================================
!  distribute particles to all pes
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 integer :: n
 do n=1,nptraj
  if (.not.particle_active(n)) then
    pxz(:,n)  = -1d12
    pik(:,n)  = -1
    puwp(:,n) = -1d12
    pxz_m(:,n) = -1e12
    puwp_m(:,n) = -1e12
  endif
 enddo
 call global_max_int2(pik, 2*nptraj)
 call global_max2(    pxz, 2*nptraj)
 call global_max2(    puwp,nr_prop*nptraj)
 call global_max2(    pxz_m, 2*nptraj)
 call global_max2(    puwp_m, nr_prop*nptraj)
end subroutine particle_distribute



subroutine integrate_particles
!=======================================================================
!       integrate particles
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 integer :: i,k,n
 real*8 :: xe,xw,zn,zs,xezs,xezn,xwzs,xwzn
 !real*8 :: uu,ww

 call particle_pe_domain()
 do n=1,nptraj
     if (particle_active(n)) then
          
!-----------------------------------------------------------------------
!      pij gives tracer box of the particle, 
!      find u-box and distances to borders, account for free slip
!      and interpolate u on particle position
!-----------------------------------------------------------------------
       i  = pik(1,n);
       xe = (xu(i) - pxz(1,n)); xw = (pxz(1,n)-(xu(i)-dx )   )

       k  = pik(2,n);
       if (pxz(2,n) > zt(k) )  k=k+1
       zn = (zt(k) - pxz(2,n)); zs = (pxz(2,n) - zt(k-1))

       xezs = xe*zs/(dx*dz); xwzs = xw*zs/(dx*dz) 
       xezn = xe*zn/(dx*dz); xwzn = xw*zn/(dx*dz)
       puwp(1,n) = u(i-1,k)*xezs + u(i,k)*xwzs + u(i-1,k-1)*xezn + u(i,k-1)*xwzn
       !uu = um1(i-1,k)*xezs + um1(i,k)*xwzs + um1(i-1,k-1)*xezn + um1(i,k-1)*xwzn 

!-----------------------------------------------------------------------
!      find v-box and distances to borders, account for free slip
!      and interpolate w on particle position
!-----------------------------------------------------------------------
       i  = pik(1,n); 
       if (pxz(1,n) > xt(i) ) then
        i=i+1
        xe = (xt(i-1)+dx - pxz(1,n)); xw = (pxz(1,n)-xt(i-1))
       else
        xe = (xt(i) - pxz(1,n)); xw = (pxz(1,n)-(xt(i)-dx) )
       endif

       k  = pik(2,n);
       zn = (zu(k) - pxz(2,n)); zs = (pxz(2,n) - zu(k-1))

       xezs = xe*zs/(dx*dz); xwzs = xw*zs/(dx*dz)
       xezn = xe*zn/(dx*dz); xwzn = xw*zn/(dx*dz)      
       puwp(2,n) = w(i-1,k)*xezs + w(i,k)*xwzs + w(i-1,k-1)*xezn + w(i,k-1)*xwzn 
       !ww = wm1(i-1,k)*xezs + wm1(i,k)*xwzs + wm1(i-1,k-1)*xezn + wm1(i,k-1)*xwzn   
       
!-----------------------------------------------------------------------
!      interpolate also phi and v on particle position
!-----------------------------------------------------------------------
       i  = pik(1,n); 
       if (pxz(1,n) > xt(i) ) then
        i=i+1
        xe = (xt(i-1)+dx - pxz(1,n)); xw = (pxz(1,n)-xt(i-1))
       else
        xe = (xt(i) - pxz(1,n)); xw = (pxz(1,n)-(xt(i)-dx) )
       endif

       k  = pik(2,n);
       if (pxz(2,n) > zt(k) )  k=k+1
       zn = (zt(k) - pxz(2,n)); zs = (pxz(2,n) - zt(k-1))
     
       xezs = xe*zs/(dx*dz); xwzs = xw*zs/(dx*dz)
       xezn = xe*zn/(dx*dz); xwzn = xw*zn/(dx*dz)  
       puwp(3,n) =   v(i-1,k)*xezs +   v(i,k)*xwzs +   v(i-1,k-1)*xezn +   v(i,k-1)*xwzn 
       puwp(4,n) = phi(i-1,k)*xezs + phi(i,k)*xwzs + phi(i-1,k-1)*xezn + phi(i,k-1)*xwzn 
       puwp(5,n) = p(i-1,k,tau)*xezs + p(i,k,tau)*xwzs + p(i-1,k-1,tau)*xezn + p(i,k-1,tau)*xwzn         
!-----------------------------------------------------------------------
!      integrate the particle trajectory forward for one time step
!-----------------------------------------------------------------------
       pxz(1,n) = pxz(1,n) + dt*puwp(1,n)  ! Euler forward
       pxz(2,n) = pxz(2,n) + dt*puwp(2,n)
       !pxz(1,n) = pxz(1,n) + dt*(puwp(1,n) + uu)/2. ! Runge-Kutta 2. order
       !pxz(2,n) = pxz(2,n) + dt*(puwp(2,n) + ww)/2. 
       
       pxz_m(:,n) = pxz_m(:,n) + pxz(:,n)
       puwp_m(:,n) = puwp_m(:,n) + puwp(:,n)
       
!-----------------------------------------------------------------------
!      update index of bounding tracer volume
!-----------------------------------------------------------------------
       i  = pik(1,n); k  = pik(2,n)
       if (pxz(1,n) >= xu(i)) then
            pik(1,n) = i + 1
       else if (pxz(1,n) < (xu(i)-dx)  ) then
            pik(1,n) = i - 1
       endif

       if (pxz(2,n) >= zu(k)) then
            pik(2,n) = k + 1
       else if (pxz(2,n) < zu(k-1)) then
            pik(2,n) = k - 1
       endif   
!-----------------------------------------------------------------------
!      periodic boundary conditions
!-----------------------------------------------------------------------
       if (pik(1,n)>nx) then
          pik(1,n)=pik(1,n)-nx; pxz(1,n)=pxz(1,n)-nx*dx
       endif
       if (pik(1,n)<1) then
          pik(1,n)=pik(1,n)+nx; pxz(1,n)=pxz(1,n)+nx*dx
       endif

     endif ! particle_active
 enddo ! nptraj
 call particle_distribute
 n_sum = n_sum + 1
end subroutine integrate_particles





subroutine init_write_mean_particles
!=======================================================================
! initialize netcdf output
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 include "netcdf.inc"
 integer :: ncid,tdim,pdim,id,iret

 if (my_pe==0) then
   ncid = nccre ('mean_particles.cdf', NCCLOB, iret)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   tdim = ncddef(ncid, 'Time', nf_unlimited, iret)
   pdim = ncddef(ncid, 'number', max(1,nptraj) , iret)   
   id  = ncvdef (ncid,'Time', NCFLOAT,1,tdim,iret)
   iret = nf_put_att_text (ncid,id,'long_name',4,'time')
   iret = nf_put_att_int  (ncid,id,'time_origin',nf_int,1,0)  
   
   id  = ncvdef (ncid,'mpos_x', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',10,'mean x position') 
   id  = ncvdef (ncid,'mpos_z', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',10,'mean y position')         
   id  = ncvdef (ncid,'um', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'mean velocity') 
   id  = ncvdef (ncid,'vm', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'mean velocity')  
   id  = ncvdef (ncid,'wm', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'mean velocity') 
   id  = ncvdef (ncid,'phim', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',9,'mean thickness')
   id  = ncvdef (ncid,'pm', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',9,'mean thickness')
  
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dt',NF_DOUBLE,1,dt)   
   iret = nf_close (ncid)
 endif
end subroutine init_write_mean_particles



subroutine init_write_particles
!=======================================================================
! initialize netcdf output
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 include "netcdf.inc"
 integer :: ncid,tdim,pdim,id,iret

 if (my_pe==0) then
   ncid = nccre ('particles.cdf', NCCLOB, iret)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   tdim = ncddef(ncid, 'Time', nf_unlimited, iret)
   pdim = ncddef(ncid, 'number', max(1,nptraj) , iret)   
   id  = ncvdef (ncid,'Time', NCFLOAT,1,tdim,iret)
   iret = nf_put_att_text (ncid,id,'long_name',4,'time')
   iret = nf_put_att_int  (ncid,id,'time_origin',nf_int,1,0)  
   
   id  = ncvdef (ncid,'pos_x', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',10,'x position') 
   id  = ncvdef (ncid,'pos_z', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',10,'y position')   
   id  = ncvdef (ncid,'u', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'velocity') 
   id  = ncvdef (ncid,'v', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'velocity')  
   id  = ncvdef (ncid,'w', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',8,'velocity') 
   id  = ncvdef (ncid,'phi', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',9,'thickness')
   id  = ncvdef (ncid,'p', NCFLOAT,2,(/pdim,tdim/),iret)
   iret = nf_put_att_text (ncid,id,'long_name',9,'thickness')
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',NF_DOUBLE,1,dx)
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dz',NF_DOUBLE,1,dz)
   iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dt',NF_DOUBLE,1,dt)  
   iret = nf_close (ncid)
 endif
end subroutine init_write_particles


subroutine write_mean_particles
!=======================================================================
! write to netcdf file
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 include "netcdf.inc"
 integer :: ncid,tdim,tid,iret,ilen,id
 real*8 :: fxa
 
 if (n_sum>0) then
  pxz_m(:,:)  = pxz_m(:,:)/n_sum
  puwp_m(:,:) = puwp_m(:,:)/n_sum
 endif
 if (my_pe==0) then
  print*,'writing mean particles to file mean_particles.cdf'
  iret=nf_open('mean_particles.cdf',NF_WRITE,ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  iret=nf_inq_dimid(ncid,'Time',tdim)
  iret=nf_inq_dimlen(ncid, tdim,ilen)
  iret=nf_inq_varid(ncid,'Time',tid)
  ilen=ilen+1
  fxa = itt*dt
  iret= nf_put_vara_double(ncid,tid,(/ilen/),1,fxa)  
  iret=nf_inq_varid(ncid,'mpos_x',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),pxz_m(1,1:nptraj))
  iret=nf_inq_varid(ncid,'mpos_z',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),pxz_m(2,1:nptraj))
  iret=nf_inq_varid(ncid,'um',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp_m(1,1:nptraj))
  iret=nf_inq_varid(ncid,'vm',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp_m(3,1:nptraj))
  iret=nf_inq_varid(ncid,'wm',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp_m(2,1:nptraj))  
  iret=nf_inq_varid(ncid,'phim',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp_m(4,1:nptraj)) 
  iret=nf_inq_varid(ncid,'pm',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp_m(5,1:nptraj)) 
  call ncclos (ncid, iret)
 endif 
 puwp_m = 0.0; pxz_m = 0.0; n_sum = 0.0
end subroutine write_mean_particles


subroutine write_particles
!=======================================================================
! write to netcdf file
!=======================================================================
 use main_module
 use particles_module
 implicit none 
 include "netcdf.inc"
 integer :: ncid,tdim,tid,iret,ilen,id
 real*8 :: fxa

 if (my_pe==0) then
  print*,'writing particles to file particles.cdf'
  iret=nf_open('particles.cdf',NF_WRITE,ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  iret=nf_inq_dimid(ncid,'Time',tdim)
  iret=nf_inq_dimlen(ncid, tdim,ilen)
  iret=nf_inq_varid(ncid,'Time',tid)
  ilen=ilen+1
  fxa = itt*dt
  iret= nf_put_vara_double(ncid,tid,(/ilen/),1,fxa)
  
  iret=nf_inq_varid(ncid,'pos_x',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),pxz(1,1:nptraj))
  iret=nf_inq_varid(ncid,'pos_z',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),pxz(2,1:nptraj))
  iret=nf_inq_varid(ncid,'u',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp(1,1:nptraj))
  iret=nf_inq_varid(ncid,'v',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp(3,1:nptraj))
  iret=nf_inq_varid(ncid,'w',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp(2,1:nptraj))  
  iret=nf_inq_varid(ncid,'phi',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp(4,1:nptraj))
  iret=nf_inq_varid(ncid,'p',id)
  iret= nf_put_vara_double(ncid,id,(/1,ilen/),(/nptraj,1/),puwp(5,1:nptraj))
  call ncclos (ncid, iret)
 endif
end subroutine write_particles


subroutine particles_read_restart
!=======================================================================
! read unfinished averages from file
!=======================================================================
 use main_module
 use particles_module
 implicit none
 character (len=80) :: filename
 logical :: file_exists
 integer :: io=33,nptraj_

 filename = 'particles_restart.dta'
 inquire ( FILE=filename, EXIST=file_exists )
 if (.not. file_exists) then
      if (my_pe==0) then
         print'(a,a,a)',' file ',filename(1:len_trim(filename)),' not present'
         print'(a)',' reading no restart for particles'
      endif
      return
 endif
 if (my_pe==0) print'(2a)',' reading particles from ',filename(1:len_trim(filename))
 
 open(io,file=filename,form='unformatted',status='old',err=10)
 read(io,err=10) nptraj_
 if (nptraj_ /= nptraj) then
       if (my_pe==0) then
        print*,' read number of particles ',nptraj_
        print*,' which does not match ',nptraj
       endif
       goto 10
 endif
 read(io,err=10) pik,pxz,puwp,pxz_m,puwp_m,n_sum
 close(io)

 call particle_pe_domain()
 call particle_distribute()
 return
 10 continue
 print'(a)',' Warning: error reading file'
end subroutine particles_read_restart



subroutine particles_write_restart
!=======================================================================
! write unfinished averages to restart file
!=======================================================================
 use main_module
 use particles_module
 implicit none
 character (len=80) :: filename
 integer :: io=33

 if (my_pe==0) then
  filename = 'particles_restart.dta'
  print'(a,a)',' writing particles to ',filename(1:len_trim(filename))  
  open(io,file=filename,form='unformatted',status='unknown')   
  write(io,err=10) nptraj
  write(io,err=10) pik,pxz,puwp,pxz_m,puwp_m,n_sum
  close(io)
 endif 
 
 return
 10 continue
 print'(a)',' Warning: error writing file'
end subroutine particles_write_restart



 
