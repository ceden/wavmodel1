


subroutine superbee_x(var,uu)
!---------------------------------------------------------------------------------
! from MITgcm
!---------------------------------------------------------------------------------
 use main_module  
 use timing_module  
 implicit none
 real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), intent(in) :: uu(is_pe:ie_pe,ks_pe:ke_pe)
 integer :: i,k
 real*8 :: Rjp,Rj,Rjm,uCFL=0.5,Cr
 ! Suberbee      Limiter(Cr)=max(0.,max(min(1.,2*Cr),min(2.,Cr)))
 ! Sweby         Limiter(Cr)=max(0.,max(min(1.,1.5*Cr),min(1.5.,Cr)))
 real*8 :: Limiter
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 
 !Limiter(Cr)=max(0.D0,max(min(1.D0,1.5D0*Cr), min(1.5D0,Cr))) 

 do k=ks_pe,ke_pe
  do i=is_pe,ie_pe
    uCFL = ABS( uu(i,k)*dt/dx )
    Rjp=(var(i+2,k)-var(i+1,k))
    Rj =(var(i+1,k)-var(i  ,k))
    Rjm=(var(i  ,k)-var(i-1,k))
    IF (Rj.NE.0.) THEN
          IF (uu(i,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
    ELSE
          IF (uu(i,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
    ENDIF
    Cr=Limiter(Cr)
    fe(i,k) = uu(i,k)*(var(i+1,k)+var(i,k))*0.5d0 -ABS(uu(i,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
  enddo
 enddo
 call tic('boundary')
 call border_exchg_3D(fe)
 call toc('boundary')
end subroutine superbee_x


subroutine superbee_z(var,ww)
!---------------------------------------------------------------------------------
! from MITgcm
!---------------------------------------------------------------------------------
 use main_module  
 use timing_module  
 implicit none
 real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
 real(real_type), intent(in) :: ww(is_pe:ie_pe,ks_pe:ke_pe)
 integer :: i,k
 real*8 :: Rjp,Rj,Rjm,uCFL=0.5,Cr
 ! Suberbee      Limiter(Cr)=max(0.,max(min(1.,2*Cr),min(2.,Cr)))
 ! Sweby         Limiter(Cr)=max(0.,max(min(1.,1.5*Cr),min(1.5.,Cr)))
 real*8 :: Limiter
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 
 !Limiter(Cr)=max(0.D0,max(min(1.D0,1.5D0*Cr), min(1.5D0,Cr))) 

 do k=ks_pe,ke_pe     
  do i=is_pe,ie_pe
    Rjp=(var(i,k+2)-var(i,k+1))
    Rj =(var(i,k+1)-var(i,k  ))
    Rjm=(var(i,k  )-var(i,k-1))
    uCFL = ABS( ww(i,k)*dt/dz )
    IF (Rj.NE.0.) THEN
          IF (ww(i,k).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
    ELSE
          IF (w(i,k).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
    ENDIF
    Cr=Limiter(Cr)
    ft(i,k) = ww(i,k)*(var(i,k+1)+var(i,k))*0.5d0-ABS(w(i,k))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
  enddo
 enddo
 
 if (my_blk_k == 1) then
    k=1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
 endif
  
 if (my_blk_k == n_pes_k) then
    ft(is_pe:ie_pe,nz)=0.0
    k=nz-1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
 endif
    
 call tic('boundary')
 call border_exchg_3D(ft)
 call toc('boundary')

end subroutine superbee_z




subroutine dst3_x(var,uu) 
!---------------------------------------------------------------------------------
! from MITgcm
!   Compute advective Flux of Tracer using 3rd Order DST Sceheme with flux limiting               
!---------------------------------------------------------------------------------
  use main_module  
  use timing_module  
  implicit none
  real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type), intent(in) :: uu(is_pe:ie_pe,ks_pe:ke_pe)
  integer :: i,k
  real(real_type) :: Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM
  real(real_type), parameter :: oneSixth=1.D0/6.D0
  real(real_type), parameter :: thetaMax = 1.D+20 
  
  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe
     uCFL = ABS( uu(i,k)*dt/dx )
     Rjp=(var(i+2,k)-var(i+1,k))
     Rj =(var(i+1,k)-var(i  ,k))
     Rjm=(var(i  ,k)-var(i-1,k))
     d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
     d1=(1d0 -uCFL*uCFL)*oneSixth
     IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
     ELSE
          thetaP=Rjm/Rj
     ENDIF
     IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
     ELSE
          thetaM=Rjp/Rj
     ENDIF
     psiP=d0+d1*thetaP
     psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
     psiM=d0+d1*thetaM
     psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))

     fe(i,k)= 0.5*(uu(i,k)+ABS(uu(i,k))) *( var(i  ,k) + psiP*Rj )  &
             +0.5*(uu(i,k)-ABS(uu(i,k))) *( var(i+1,k) - psiM*Rj )
   enddo
  enddo

  call tic('boundary')
  call border_exchg_3D(fe)
  call toc('boundary')
end subroutine dst3_x



subroutine dst3_z(var,ww) 
!---------------------------------------------------------------------------------
! from MITgcm
!   Compute advective Flux of Tracer using 3rd Order DST Sceheme with flux limiting               
!---------------------------------------------------------------------------------
  use main_module  
  use timing_module  
  implicit none
  real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type), intent(in) :: ww(is_pe:ie_pe,ks_pe:ke_pe)
  integer :: i,k
  real(real_type) :: Rjp,Rj,Rjm,uCFL,d0,d1,thetaP,thetaM,psiP,psiM
  real(real_type), parameter :: oneSixth=1.D0/6.D0
  real(real_type), parameter :: thetaMax = 1.D+20 

  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe
     Rjp=(var(i,k+2)-var(i,k+1))
     Rj =(var(i,k+1)-var(i,k  ))
     Rjm=(var(i,k  )-var(i,k-1))
     uCFL = ABS( ww(i,k)*dt/dz )
     d0=(2d0 -uCFL)*(1d0 -uCFL)*oneSixth
     d1=(1d0 -uCFL*uCFL)*oneSixth
     IF ( ABS(Rj)*thetaMax .LE. ABS(Rjm) ) THEN
             thetaP=SIGN(thetaMax,Rjm*Rj)
     ELSE
          thetaP=Rjm/Rj
     ENDIF
     IF ( ABS(Rj)*thetaMax .LE. ABS(Rjp) ) THEN
          thetaM=SIGN(thetaMax,Rjp*Rj)
     ELSE
          thetaM=Rjp/Rj
     ENDIF
     psiP=d0+d1*thetaP
     psiP=MAX(0d0,MIN(MIN(1d0,psiP), thetaP*(1d0 -uCFL)/(uCFL+1d-20) ))
     psiM=d0+d1*thetaM
     psiM=MAX(0d0,MIN(MIN(1d0,psiM), thetaM*(1d0 -uCFL)/(uCFL+1d-20) ))
         
     ft(i,k)= 0.5*(ww(i,k)+ABS(ww(i,k))) *( var(i,k  ) + psiP*Rj )  &
             +0.5*(ww(i,k)-ABS(ww(i,k))) *( var(i,k+1) - psiM*Rj )
   enddo
  enddo

  if (my_blk_k == 1) then
    k=1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
  endif
  
  if (my_blk_k == n_pes_k) then
    ft(is_pe:ie_pe,nz)=0.0
    k=nz-1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
  endif
    
  call tic('boundary')
  call border_exchg_3D(ft)
  call toc('boundary')
end subroutine dst3_z








subroutine upwind3_x(var,uu)
!---------------------------------------------------------------------------------
! upwind third order from MITgcm
!---------------------------------------------------------------------------------
  use main_module  
  use timing_module  
  implicit none
  real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type), intent(in) :: uu(is_pe:ie_pe,ks_pe:ke_pe)
      
  integer :: i,k
  real*8 :: Rjp,Rj,Rjm,Rjjp,Rjjm
  real*8, parameter :: oneSixth=1.D0/6.D0

  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe
     Rjp=(var(i+2,k)-var(i+1,k))
     Rj =(var(i+1,k)-var(i  ,k))
     Rjm=(var(i  ,k)-var(i-1,k))
     Rjjp=Rjp-Rj; Rjjm=Rj-Rjm
     fe(i,k) = uu(i,k)*(var(i+1,k)+var(i,k) -oneSixth*( Rjjp+Rjjm )  )*0.5d0   &
                         +ABS(uu(i,k))*0.5d0*oneSixth*( Rjjp-Rjjm )
   enddo
  enddo
  
  call tic('boundary')
  call border_exchg_3D(fe)
  call toc('boundary')
end subroutine upwind3_x


subroutine upwind3_z(var,ww)
!---------------------------------------------------------------------------------
! upwind third order from MITgcm
!---------------------------------------------------------------------------------
  use main_module  
  use timing_module  
  implicit none
  real(real_type), intent(in) :: var(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx)
  real(real_type), intent(in) :: ww(is_pe:ie_pe,ks_pe:ke_pe)      
  integer :: i,k
  real*8 :: Rjp,Rj,Rjm,Rjjp,Rjjm
  real*8, parameter :: oneSixth=1.D0/6.D0

  do k=ks_pe,ke_pe
   do i=is_pe,ie_pe
     Rjp=(var(i,k+2)-var(i,k+1))
     Rj =(var(i,k+1)-var(i,k  ))
     Rjm=(var(i,k  )-var(i,k-1))
     Rjjp=Rjp-Rj; Rjjm=Rj-Rjm
     ft(i,k) = ww(i,k)*(var(i,k+1)+var(i,k)   -oneSixth*( Rjjp+Rjjm )   )*0.5d0   &
                           +ABS(ww(i,k))*0.5d0*oneSixth*( Rjjp-Rjjm )
   enddo
  enddo
  
  if (my_blk_k == 1) then
    k=1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
  endif
  
  if (my_blk_k == n_pes_k) then
    ft(is_pe:ie_pe,nz)=0.0
    k=nz-1
    ft(is_pe:ie_pe,k) =  ww(is_pe:ie_pe,k)*0.5*(var(is_pe:ie_pe,k)+var(is_pe:ie_pe,k+1))
  endif
    
  call tic('boundary')
  call border_exchg_3D(ft)
  call toc('boundary')
end subroutine upwind3_z




