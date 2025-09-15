
module config_module
 implicit none
 integer :: fac = 16
 real*8 :: om_glob,k_glob, f_glob
end module config_module

subroutine set_parameter
 use main_module 
 use config_module
 implicit none
 
 nx=50*fac; nz = 50*fac
 Lx=1.; Lz = 1.
 dx=Lx/nx; dz=Lz/nz; dt = 0.0025/fac
 
 k_glob = 4*pi
 f_glob = sqrt(9.81*k_glob)/50
 om_glob = sqrt(f_glob**2 + 9.81*k_glob*tanh(k_glob*Lz/2.)) 
 
 congr_eps=1e-3
 max_iterations = 25000
 runlen =  2*pi/om_glob
 snapint = runlen/100.
 meanint = 2*pi/om_glob
 tsmonint = dt
 enable_dst3_advection = .true.
 enable_multidim_advection = .true.
 enable_AB3_time_stepping  = .true.
 
 phi_diff = 1e-4
 gamma = 0.6
 enable_particles = .true.
 
 Coriolis = f_glob
 !nu_o = 1.3e-6
 !nu_a = 1.6e-5 
end subroutine set_parameter 



subroutine set_initial_conditions
 use main_module  
 use config_module
 implicit none
 integer :: i,k,n
 real(real_type) :: H,xt(nx),zt(nz)
 real(real_type) :: a_loc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: zeta,zeta0

 
  do i=1,nx
   xt(i)=(i-1)*dx + dx/2.
  enddo  
  do i=1,nz
   zt(i)=(i-1)*dz + dz/2. - Lz/2.
  enddo 
  
  H=Lz/2.
  zeta0 = 0.01  

  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe  
   
    zeta = zeta0*cos(k_glob*xt(i))
    zeta = zeta + k_glob*zeta0**2/2*cos(2*k_glob*xt(i))
    
    !u(i,k) = real(k_glob*zeta0*om_glob*exp(cmplx(0,1)*(k_glob*xt(i))) )
    !u(i,k) = u(i,k)*cosh(k_glob*(zt(k)+H))/(k_glob*sinh(k_glob*H))
    
    u(i,k) = zeta0*om_glob*exp(k_glob*zt(k))*cos(k_glob*xt(i))
    u(i,k) = u(i,k) + k_glob*om_glob*zeta0**2*exp(2*k_glob*zt(k))*cos(2*k_glob*xt(i))
    if (zeta <= zt(k) ) u(i,k)=0d0

    
    !w(i,k) = real(-k_glob*cmplx(0,1)*zeta0*om_glob*exp(cmplx(0,1)*(k_glob*xt(i))) )
    !w(i,k) = w(i,k)*sinh(k_glob*(zt(k)+H))/(k_glob*sinh(k_glob*H))
    
    w(i,k) = zeta0*om_glob*exp(k_glob*zt(k))*sin(k_glob*xt(i))
    w(i,k) = w(i,k) + k_glob*om_glob*zeta0**2*exp(2*k_glob*zt(k))*sin(2*k_glob*xt(i))
    if ( zeta <= zt(k) ) w(i,k)=0d0

    phi(i,k) = -tanh( (zt(k)-zeta)/(4*dz) )   
   
   enddo 
  enddo
  if (my_blk_k == n_pes_k) w(is_pe:ie_pe,nz) = 0d0  


if (.false. ) then
  a_loc=0d0
  do n=1,2
   call border_exchg_3D(u)
   do k=ks_pe,ke_pe
    do i=is_pe,ie_pe
     a_loc(i,k) =              u(i  ,k)+u(i  ,min(nz,k+1))+u(i  ,max(1,k-1))
     a_loc(i,k) = a_loc(i,k) + u(i+1,k)+u(i+1,min(nz,k+1))+u(i+1,max(1,k-1))
     a_loc(i,k) = a_loc(i,k) + u(i-1,k)+u(i-1,min(nz,k+1))+u(i-1,max(1,k-1))
     a_loc(i,k) = a_loc(i,k)/9.0
    enddo
   enddo
   u = a_loc
  enddo
  
  a_loc=0d0
  do n=1,2 
   call border_exchg_3D(w)   
   do k=ks_pe,ke_pe
    do i=is_pe,ie_pe
     a_loc(i,k) =              w(i  ,k)+w(i  ,min(nz-1,k+1))+w(i  ,max(1,k-1))
     a_loc(i,k) = a_loc(i,k) + w(i+1,k)+w(i+1,min(nz-1,k+1))+w(i+1,max(1,k-1))
     a_loc(i,k) = a_loc(i,k) + w(i-1,k)+w(i-1,min(nz-1,k+1))+w(i-1,max(1,k-1))
     a_loc(i,k) = a_loc(i,k)/9.0
    enddo   
   
    !a_loc(:,k) = (w(:,k)+w(:,min(nz-1,k+1))+w(:,max(1,k-1)))/3.0
   enddo
   w = a_loc
   if (my_blk_k == n_pes_k) w(is_pe:ie_pe,nz) = 0d0 
  enddo
  
endif

  call border_exchg_3D(u)
  call border_exchg_3D(w)
  if (my_blk_k == n_pes_k) w(is_pe:ie_pe,nz) = 0d0 
  call border_exchg_3D(phi)
  
  
  if (enable_particles) call set_particles
  
end subroutine set_initial_conditions


subroutine set_particles
 use main_module  
 use config_module
 use particles_module 
 implicit none
 integer :: n,pmax=4000
 real(real_type) :: fxa
 
  call allocate_particles(pmax)
  do n=1,pmax
   call random_number(fxa)
   pxz(1,n) = fxa*dx*nx
   call random_number(fxa)
   pxz(2,n) = fxa*dz*nz
   call pe0_bcast(pxz(:,n),2)
  enddo 
  call init_particles
end subroutine set_particles  



subroutine set_forcing 
end subroutine set_forcing
