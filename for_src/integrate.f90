
subroutine integrate
 use main_module
 use timing_module   
 implicit none
 integer :: i,k
 real(real_type) :: diff(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: conc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: mu(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: nu(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: rho(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: a_tri(ks_pe:ke_pe),b_tri(ks_pe:ke_pe),c_tri(ks_pe:ke_pe)
 real(real_type) :: a_tri_all(nz),b_tri_all(nz),c_tri_all(nz),d_tri_all(nz),x_tri_all(nz)
 
 ! allow for user defined forcing
 call set_forcing()

 ! concentration
 conc = phi/2.0+0.5
 if (my_blk_k == 1)       conc(:,1)  = conc(:,2)
 if (my_blk_k == n_pes_k) conc(:,nz) = conc(:,nz-1) 

 ! friction
 if (nu_o/=0d0 .or. nu_a/=0d0)  nu = (1+phi)/2. * nu_o + (1-phi)/2. * nu_a

 ! chem. potential
 mu = min( max(12*conc**2 - 12*conc + 2,-1.),2.)
 diff = phi_diff*(1-gamma*phi**2)*mu

 ! coefficient for pressure solver, harmonic interpolation for 1/rho 
 rho = rho_o*(1+phi)/2 + rho_a*(1-phi)/2   
 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe
   irho_ip(i,k) = r1 + 0.5*(phi(i+1,k) + phi(i,k) )*r2
   irho_im(i,k) = r1 + 0.5*(phi(i-1,k) + phi(i,k) )*r2
   irho_kp(i,k) = r1 + 0.5*(phi(i,k+1) + phi(i,k) )*r2
   irho_km(i,k) = r1 + 0.5*(phi(i,k-1) + phi(i,k) )*r2
  enddo
 enddo 
 call tic('boundary')
 call border_exchg_3D(irho_ip) 
 call border_exchg_3D(irho_im) 
 call border_exchg_3D(irho_kp) 
 call border_exchg_3D(irho_km)  
 call toc('boundary')


 ! fluxes, gravity, Coriolis
 call u_tendency 
 if (coriolis /=0d0 .or. enable_v_velocity ) call v_tendency
 if (coriolis /=0d0) then   
   do k=ks_pe,ke_pe
    do i=is_pe,ie_pe
     du(i,k,tau) = du(i,k,tau) + (v(i  ,k) + v(i+1,k))*0.5*coriolis
     dv(i,k,tau) = dv(i,k,tau) - (u(i-1,k) + u(i  ,k))*0.5*coriolis
    enddo
   enddo
 endif 
 call w_tendency
 call phi_tendency 
 
 
 if (phi_diff>0d0) then
  ! Cahn-Hilliard magic for phi 
  !for k=1
  ! phi^n_ik =  phi^n-1_ik +  dt/dz**2*( (a)^k+1 *phi^n_i,k+1 - (a)^k+1  phi^n_ik  )
  ! for  k=2:nz-1
  ! phi^n_ik =  phi^n-1_ik +  dt/dz**2*
  ! ( (a)^k+1 *phi^n_i,k+1 - ((a)^k+1 + (a)^k-1) phi^n_ik + (a)^k-1 phi^n_i,k-1 )
  ! for k=nz
  ! phi^n_ik =  phi^n-1_ik +  dt/dz**2* ( -  (a)^k-1 phi^n_ik + (a)^k-1 phi^n_i,k-1 )
  do i=is_pe,ie_pe
    a_tri(ks_pe:ke_pe) =   - dt/dz**2*0.5*(diff(i,ks_pe:ke_pe)+diff(i,ks_pe-1:ke_pe-1) )        
    b_tri(ks_pe:ke_pe) = 1+dt/dz**2*( 0.5*(diff(i,ks_pe:ke_pe)+diff(i,ks_pe-1:ke_pe-1)) &  
                                     +0.5*(diff(i,ks_pe:ke_pe)+diff(i,ks_pe+1:ke_pe+1)) )        
    c_tri(ks_pe:ke_pe) =   - dt/dz**2*0.5*(diff(i,ks_pe:ke_pe)+diff(i,ks_pe+1:ke_pe+1))
        
    if (my_blk_k == 1) then
          a_tri(1) = 0.
          b_tri(1) = 1+dt/dz**2*0.5*(diff(i,2)+diff(i,1))
    else if   (my_blk_k == n_pes_k) then
          c_tri(nz) = 0.
          b_tri(nz) = 1+dt/dz**2*0.5*(diff(i,nz)+diff(i,nz-1) )
    endif
                        
    call retrieve_in_z(a_tri,a_tri_all)
    call retrieve_in_z(b_tri,b_tri_all)
    call retrieve_in_z(c_tri,c_tri_all)
    call retrieve_in_z(phi(i,ks_pe:ke_pe),d_tri_all)
    call solve_tridiag(a_tri_all,b_tri_all,c_tri_all,d_tri_all,x_tri_all,nz)
    call distribute_in_z(x_tri_all,phi(i,ks_pe:ke_pe))
  enddo
 endif
 
 if (nu_o/=0d0 .or. nu_a/=0d0) then
  do k=ks_pe-1,ke_pe
   do i=is_pe-1,ie_pe  
    fe(i,k) =                                       nu(i+1,k)*(u(i+1,k)-u(i,k))/dx
    ft(i,k) =  0.25*(nu(i+1,k)+nu(i,k)+nu(i+1,k+1)+nu(i,k+1))*(u(i,k+1)-u(i,k))/dz
   enddo
  enddo 
  if (my_blk_k == n_pes_k) ft(:,nz-1:nz+onx)=0d0
  if (my_blk_k == 1      ) ft(:,1-onx:0)    =0d0
  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe  
     u(i,k) = u(i,k) + dt*( (fe(i,k)-fe(i-1,k))/dx + (ft(i,k)-ft(i,k-1))/dz  )
   enddo
  enddo  
  
  if (coriolis /=0d0 .or. enable_v_velocity ) then
   do k=ks_pe-1,ke_pe
    do i=is_pe-1,ie_pe  
     fe(i,k) =  0.5*(nu(i+1,k)+nu(i,k))*(v(i+1,k)-v(i,k))/dx
     ft(i,k) =  0.5*(nu(i,k+1)+nu(i,k))*(v(i,k+1)-v(i,k))/dz
    enddo
   enddo 
   if (my_blk_k == n_pes_k) ft(:,nz-1:nz+onx)=0d0
   if (my_blk_k == 1      ) ft(:,1-onx:0)    =0d0
   do k=ks_pe,ke_pe
    do i=is_pe,ie_pe  
     v(i,k) = v(i,k) + dt*( (fe(i,k)-fe(i-1,k))/dx + (ft(i,k)-ft(i,k-1))/dz  )
    enddo
   enddo      
  endif
  
  do k=ks_pe-1,ke_pe
   do i=is_pe-1,ie_pe  
    fe(i,k) =  0.25*(nu(i+1,k)+nu(i,k)+nu(i+1,k+1)+nu(i,k+1))*(w(i+1,k)-w(i,k))/dx
    ft(i,k) =                                       nu(i,k+1)*(w(i,k+1)-w(i,k))/dz
   enddo
  enddo 
  if (my_blk_k == n_pes_k) ft(:,nz-2:nz+onx)=0d0
  if (my_blk_k == 1      ) ft(:,1-onx:0)    =0d0

  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe  
    w(i,k) = w(i,k) + dt*( (fe(i,k)-fe(i-1,k))/dx + (ft(i,k)-ft(i,k-1))/dz  )
   enddo
  enddo  
  if (my_blk_k == n_pes_k) w(:,nz-1:nz+onx)=0d0
 endif
 


 if (enable_AB3_time_stepping) then  
  if (itt==0) then
   u = u + dt*du(:,:,tau) 
   if (coriolis /=0d0 .or. enable_v_velocity) v = v + dt*dv(:,:,tau) 
   w = w + dt*dw(:,:,tau) 
   phi = phi + dt*dphi(:,:,tau) 
  elseif (itt==1) then     
   u = u + dt*( 1.5*du(:,:,tau) - 0.5*du(:,:,taum1)) 
   if (coriolis /=0d0 .or. enable_v_velocity) v = v + dt*( 1.5*dv(:,:,tau) - 0.5*dv(:,:,taum1)) 
   w = w + dt*( 1.5*dw(:,:,tau) - 0.5*dw(:,:,taum1)) 
   phi = phi + dt*( 1.5*dphi(:,:,tau) - 0.5*dphi(:,:,taum1)) 
  else
   u = u + dt*( AB3_a*du(:,:,tau) + AB3_b*du(:,:,taum1) + AB3_c*du(:,:,taum2)) 
   if (coriolis /=0d0 .or. enable_v_velocity ) &
    v = v + dt*( AB3_a*dv(:,:,tau) + AB3_b*dv(:,:,taum1) + AB3_c*dv(:,:,taum2)) 
   w = w + dt*( AB3_a*dw(:,:,tau) + AB3_b*dw(:,:,taum1) + AB3_c*dw(:,:,taum2)) 
   phi = phi + dt*( AB3_a*dphi(:,:,tau) + AB3_b*dphi(:,:,taum1) + AB3_c*dphi(:,:,taum2)) 
  endif 
 else
  u = u + dt*du(:,:,tau) 
  if (coriolis /=0d0 .or. enable_v_velocity) v = v + dt*dv(:,:,tau) 
  w = w + dt*dw(:,:,tau) 
  phi = phi + dt*dphi(:,:,tau) 
 endif

 call tic('boundary')
 call border_exchg_3D(u) 
 if (coriolis /=0d0 .or. enable_v_velocity) call border_exchg_3D(v) 
 call border_exchg_3D(w)
 call border_exchg_3D(phi)
 call toc('boundary')

 call tic('pressure')
 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe
   fe(i,k) =  ((u(i,k)-u(i-1,k))/dx + (w(i,k)-w(i,k-1))/dz ) /dt   
   precond(i,k) =  -dx**2/irho_ip(i,k) -dx**2/irho_im(i,k) &
                   -dz**2/irho_kp(i,k) -dz**2/irho_km(i,k)
  enddo
 enddo  
 
 p(:,:,tau) = 2*p(:,:,taum1) - p(:,:,taum2)
 call congrad(fe,p(:,:,tau))
 call toc('pressure')
 
 call tic('boundary')
 call border_exchg_3D(p(:,:,tau)) 
 call toc('boundary')
 
 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe
   u(i,k) = u(i,k) - dt* (p(i+1,k,tau) - p(i,k,tau))/dx*irho_ip(i,k)
   w(i,k) = w(i,k) - dt* (p(i,k+1,tau) - p(i,k,tau))/dz*irho_kp(i,k)
  enddo
 enddo
 if (my_blk_k == n_pes_k) w(is_pe:ie_pe,nz) = 0d0  
 
 call tic('boundary')
 call border_exchg_3D(u) 
 call border_exchg_3D(w) 
 call toc('boundary')
   
end subroutine integrate






subroutine solve_tridiag(a,b,c,d,x,n)
      implicit none
      integer, parameter :: real_type = KIND(1.d0)
!---------------------------------------------------------------------------------
!        a - sub-diagonal (means it is the diagonal below the main diagonal)
!        b - the main diagonal
!        c - sup-diagonal (means it is the diagonal above the main diagonal)
!        d - right part
!        x - the answer
!        n - number of equations
!---------------------------------------------------------------------------------
        integer,intent(in) :: n
        real(real_type) ,dimension(n),intent(in) :: a,b,c
        real(real_type) ,dimension(n),intent(in ) :: d
        real(real_type) ,dimension(n),intent(out) :: x
        real(real_type) ,dimension(n) :: dp
        real(real_type) ,dimension(n) :: cp
        real(real_type)  :: m,fxa
        integer i
 
! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           fxa = 1D0/m
           cp(i) = c(i)*fxa
           dp(i) = (d(i)-dp(i-1)*a(i))*fxa
         enddo
! initialize x
         x(n) = dp(n) 
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)  
        end do
end subroutine solve_tridiag







