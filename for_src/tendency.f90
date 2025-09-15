


subroutine u_tendency
 ! u- tendency
 use main_module
 use timing_module   
 implicit none
 real(real_type) :: uu(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: ww(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: a_loc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
   
 if (enable_multidim_advection) then
  uu=0d0;ww=0d0
  uu(is_pe:ie_pe,ks_pe:ke_pe) = 0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe+1:ie_pe+1,ks_pe:ke_pe))
  ww(is_pe:ie_pe,ks_pe:ke_pe) = 0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe+1:ie_pe+1,ks_pe:ke_pe))
  call tic('boundary')
  call border_exchg_3D(uu) 
  call border_exchg_3D(ww)  
  call toc('boundary')

  if (enable_dst3_advection) then
    call    dst3_x(u,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_upwind3_advection) then
    call upwind3_x(u,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_superbee_advection) then  
    call superbee_x(u,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = u(is_pe:ie_pe,ks_pe:ke_pe) &
             - dt*( (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
                    -u(is_pe:ie_pe,ks_pe:ke_pe)* &
                    ( uu(is_pe:ie_pe,ks_pe:ke_pe) -  uu(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx )
  call tic('boundary')
  call border_exchg_3D(a_loc)
  call toc('boundary')
  
  if (enable_dst3_advection) then
    call    dst3_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe)) 
  else if (enable_upwind3_advection) then
    call upwind3_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_superbee_advection) then  
    call superbee_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
       
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = a_loc(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz &
                   -u(is_pe:ie_pe,ks_pe:ke_pe)* &
                    ( ww(is_pe:ie_pe,ks_pe:ke_pe) - ww(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz )
                                        
  du(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
            (a_loc(is_pe:ie_pe,ks_pe:ke_pe) - u(is_pe:ie_pe,ks_pe:ke_pe))/dt  

 else
  if (enable_dst3_advection) then
   call    dst3_x(u,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe+1:ie_pe+1,ks_pe:ke_pe)) )
   call    dst3_z(u,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe+1:ie_pe+1,ks_pe:ke_pe)) )   
  else if (enable_upwind3_advection) then
   call upwind3_x(u,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe+1:ie_pe+1,ks_pe:ke_pe)) )
   call upwind3_z(u,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe+1:ie_pe+1,ks_pe:ke_pe)) ) 
  else if (enable_superbee_advection) then  
   call superbee_x(u,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe+1:ie_pe+1,ks_pe:ke_pe)) )
   call superbee_z(u,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe+1:ie_pe+1,ks_pe:ke_pe)) ) 
  else
    call halt_stop('no advection scheme chosen')
  endif 
  du(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
           - (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe) )/dx &
           - (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1) )/dz  
 endif
end subroutine u_tendency




subroutine w_tendency
 ! w tendency
 use main_module
 use timing_module   
 implicit none
 real(real_type) :: uu(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: ww(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type) :: a_loc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 
 if (enable_multidim_advection) then 
  uu=0d0;ww=0d0
  uu(is_pe:ie_pe,ks_pe:ke_pe) = 0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe:ie_pe,ks_pe+1:ke_pe+1))
  ww(is_pe:ie_pe,ks_pe:ke_pe) = 0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe:ie_pe,ks_pe+1:ke_pe+1))
  call tic('boundary')
  call border_exchg_3D(uu) 
  call border_exchg_3D(ww)  
  call toc('boundary')
  
  if (enable_dst3_advection) then
   call    dst3_x(w,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_upwind3_advection) then
   call upwind3_x(w,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_superbee_advection) then  
   call superbee_x(w,uu(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif  
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = w(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
                   -w(is_pe:ie_pe,ks_pe:ke_pe)* &
                   ( uu(is_pe:ie_pe,ks_pe:ke_pe) -  uu(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx )
  call tic('boundary')
  call border_exchg_3D(a_loc)
  call toc('boundary')
  
  if (enable_dst3_advection) then                                        
   call    dst3_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe)) 
  else if (enable_upwind3_advection) then 
   call upwind3_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe)) 
  else if (enable_superbee_advection) then   
   call superbee_z(a_loc,ww(is_pe:ie_pe,ks_pe:ke_pe)) 
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = a_loc(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz &
                   -w(is_pe:ie_pe,ks_pe:ke_pe)* &
                    ( ww(is_pe:ie_pe,ks_pe:ke_pe) - ww(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz )
                                        
  dw(is_pe:ie_pe,ks_pe:ke_pe,tau) = -g  &
            + (a_loc(is_pe:ie_pe,ks_pe:ke_pe) - w(is_pe:ie_pe,ks_pe:ke_pe))/dt  
 else
  if (enable_dst3_advection) then
    call    dst3_x(w,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe:ie_pe,ks_pe+1:ke_pe+1)) ) 
    call    dst3_z(w,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe:ie_pe,ks_pe+1:ke_pe+1)) )    
  else if (enable_upwind3_advection) then 
    call upwind3_x(w,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe:ie_pe,ks_pe+1:ke_pe+1)) ) 
    call upwind3_z(w,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe:ie_pe,ks_pe+1:ke_pe+1)) ) 
  else if (enable_superbee_advection) then   
    call superbee_x(w,0.5*(u(is_pe:ie_pe,ks_pe:ke_pe) + u(is_pe:ie_pe,ks_pe+1:ke_pe+1)) ) 
    call superbee_z(w,0.5*(w(is_pe:ie_pe,ks_pe:ke_pe) + w(is_pe:ie_pe,ks_pe+1:ke_pe+1)) )   
  else
    call halt_stop('no advection scheme chosen')
  endif
  dw(is_pe:ie_pe,ks_pe:ke_pe,tau) =  -g  &
          - (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
          - (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz 
 endif
  
 if (my_blk_k == n_pes_k) dw(is_pe:ie_pe,nz,tau) = 0d0            
end subroutine w_tendency




subroutine v_tendency
 ! v tendency
 use main_module
 use timing_module   
 implicit none
 real(real_type) :: a_loc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)

 if (enable_multidim_advection) then 
  if (enable_dst3_advection) then
   call    dst3_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_upwind3_advection) then  
   call upwind3_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_superbee_advection) then 
   call superbee_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = v(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
                   -v(is_pe:ie_pe,ks_pe:ke_pe)* &
                   ( u(is_pe:ie_pe,ks_pe:ke_pe) -  u(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx )
  call tic('boundary')
  call border_exchg_3D(a_loc)
  call toc('boundary')
  
  if (enable_dst3_advection) then                                        
   call    dst3_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe))  
  else if (enable_upwind3_advection) then  
   call upwind3_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe)) 
  else if (enable_superbee_advection) then  
   call superbee_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = a_loc(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz &
                   -v(is_pe:ie_pe,ks_pe:ke_pe)* &
                    ( w(is_pe:ie_pe,ks_pe:ke_pe) - w(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz )
                                        
  dv(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
            (a_loc(is_pe:ie_pe,ks_pe:ke_pe) - v(is_pe:ie_pe,ks_pe:ke_pe))/dt  
 else
  if (enable_dst3_advection) then
   call    dst3_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
   call    dst3_z(v,w(is_pe:ie_pe,ks_pe:ke_pe))     
  else if (enable_upwind3_advection) then  
   call upwind3_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
   call upwind3_z(v,w(is_pe:ie_pe,ks_pe:ke_pe))  
  else if (enable_superbee_advection) then   
   call superbee_x(v,u(is_pe:ie_pe,ks_pe:ke_pe))
   call superbee_z(v,w(is_pe:ie_pe,ks_pe:ke_pe))    
  else
    call halt_stop('no advection scheme chosen')
  endif
  dv(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
          - (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
          - (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz  
 endif
end subroutine v_tendency




subroutine phi_tendency
 ! phi tendency
 use main_module
 use timing_module   
 implicit none
 real(real_type) :: a_loc(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)

 if (enable_multidim_advection) then 
  if (enable_dst3_advection) then
   call    dst3_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_upwind3_advection) then  
   call upwind3_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
  else if (enable_superbee_advection) then 
   call superbee_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = phi(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
                   -phi(is_pe:ie_pe,ks_pe:ke_pe)* &
                   ( u(is_pe:ie_pe,ks_pe:ke_pe) -  u(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx )
  call tic('boundary')
  call border_exchg_3D(a_loc)
  call toc('boundary')
  
  if (enable_dst3_advection) then                                        
   call    dst3_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe))  
  else if (enable_upwind3_advection) then  
   call upwind3_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe)) 
  else if (enable_superbee_advection) then  
   call superbee_z(a_loc,w(is_pe:ie_pe,ks_pe:ke_pe))
  else
    call halt_stop('no advection scheme chosen')
  endif
  a_loc(is_pe:ie_pe,ks_pe:ke_pe) = a_loc(is_pe:ie_pe,ks_pe:ke_pe) &
            - dt*( (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz &
                   -phi(is_pe:ie_pe,ks_pe:ke_pe)* &
                    ( w(is_pe:ie_pe,ks_pe:ke_pe) - w(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz )
                                        
  dphi(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
            (a_loc(is_pe:ie_pe,ks_pe:ke_pe) - phi(is_pe:ie_pe,ks_pe:ke_pe))/dt  
 else
  if (enable_dst3_advection) then
   call    dst3_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
   call    dst3_z(phi,w(is_pe:ie_pe,ks_pe:ke_pe))     
  else if (enable_upwind3_advection) then  
   call upwind3_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
   call upwind3_z(phi,w(is_pe:ie_pe,ks_pe:ke_pe))  
  else if (enable_superbee_advection) then   
   call superbee_x(phi,u(is_pe:ie_pe,ks_pe:ke_pe))
   call superbee_z(phi,w(is_pe:ie_pe,ks_pe:ke_pe))    
  else
    call halt_stop('no advection scheme chosen')
  endif
  dphi(is_pe:ie_pe,ks_pe:ke_pe,tau) =  &
          - (fe(is_pe:ie_pe,ks_pe:ke_pe) - fe(is_pe-1:ie_pe-1,ks_pe:ke_pe))/dx &
          - (ft(is_pe:ie_pe,ks_pe:ke_pe) - ft(is_pe:ie_pe,ks_pe-1:ke_pe-1))/dz  
 endif
end subroutine phi_tendency


