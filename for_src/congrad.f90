

subroutine congrad(b,x)
 use main_module
 use timing_module   
 implicit none
 real(real_type), intent(inout)  :: b(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), intent(inout)  :: x(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), allocatable :: Ap(:,:),r(:,:),z(:,:),p1(:,:)
 real(real_type) :: alpha,rsold,rsnew,dot,mymax,fxa
 real(real_type) :: step,step1,estimated_error,convergence_rate,smax
 
 allocate( r(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ); r=0.
 allocate(Ap(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ); Ap=0.
 allocate( z(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ); z=0.
 allocate( p1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ); p1=0.
 
 call apply_op(x,r)
 r(is_pe:ie_pe,ks_pe:ke_pe) = b(is_pe:ie_pe,ks_pe:ke_pe) - r(is_pe:ie_pe,ks_pe:ke_pe) 
 call apply_pre(r,z)  
 p1 = z
 call tic('boundary')
 call border_exchg_3D(p1)
 call toc('boundary')
 rsold = dot(z,r)
 do itts=1,max_iterations
    call apply_op(p1,Ap)    
    fxa = dot(p1,Ap)
    if (fxa == 0.) then
      estimated_error = 0.
      return 
    endif  
    alpha = rsold / fxa
    x = x + alpha * p1
    r = r - alpha * Ap
    call apply_pre(r,z) 
    rsnew = dot(z,r)
    p1 = z + (rsnew / rsold) * p1
    call tic('boundary')
    call border_exchg_3D(p1)
    call toc('boundary')
    rsold = rsnew
    smax = mymax(p1)
    step = abs(alpha) * smax  
    if (itts == 1) then
          step1 = step
          estimated_error = step
          if (step < congr_eps) return 
    else if (step < congr_eps) then          
          convergence_rate = exp(log(step/step1)/(itts-1))
          estimated_error = step*convergence_rate/(1.0-convergence_rate)
          if (estimated_error < congr_eps) return       
    endif      
 enddo     
 deallocate(r,Ap,z,p1)
end subroutine congrad
  
  
  
subroutine apply_op(p1,Ap1)
 use main_module
 implicit none
 real(real_type), intent(inout)  ::  p1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), intent(inout)  :: Ap1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 integer :: i,k
 real(real_type) :: px_p,px_m,pz_p,pz_m
 
 if (my_blk_k == 1      )  p1(:,0)    = p1(:,1)   
 if (my_blk_k == n_pes_k)  p1(:,nz+1) = p1(:,nz)  
 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe   
    px_p = ( p1(i+1,k) - p1(  i,k) )/dx*irho_ip(i,k) 
    px_m = ( p1(  i,k) - p1(i-1,k) )/dx*irho_im(i,k)
    pz_p = ( p1(i,k+1) - p1(i,k  ) )/dz*irho_kp(i,k) 
    pz_m = ( p1(i,k  ) - p1(i,k-1) )/dz*irho_km(i,k)
    Ap1(i,k) =  (px_p - px_m)/dx + (pz_p - pz_m)/dz
  enddo
 enddo 
end subroutine apply_op


subroutine apply_pre(p1,Ap1)
 use main_module
 implicit none
 real(real_type), intent(inout)  ::  p1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), intent(inout)  :: Ap1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 integer :: i,k
  
 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe 
   !Ap1(i,k) =  p1(i,k)*(-dx**2/irho_ip(i,k) -dx**2/irho_im(i,k) &
   !                     -dz**2/irho_kp(i,k) -dz**2/irho_km(i,k) )
   Ap1(i,k) =  p1(i,k)*precond(i,k)
  enddo
 enddo 
end subroutine apply_pre


function dot(p1,p2)
  use main_module   
  implicit none
  real(real_type), intent(in)  ::  p1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type), intent(in)  ::  p2(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type) :: dot
  integer :: i,k   
  dot=0d0
  do k=ks_pe,ke_pe
    do i=is_pe,ie_pe
      dot = dot + p1(i,k)*p2(i,k)
    enddo
  enddo
  call global_sum(dot)
end function dot


function mymax(p1)
  use main_module   
  implicit none
  real(real_type), intent(in)  ::  p1(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type) :: mymax
  integer :: i,k
  mymax=0d0
  do k=ks_pe,ke_pe
    do i=is_pe,ie_pe
      mymax = max( abs(p1(i,k)), mymax )  
    enddo
  enddo
  call global_max(mymax) 
end function mymax
