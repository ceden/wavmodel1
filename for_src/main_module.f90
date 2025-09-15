

module main_module
 implicit none

 integer, parameter :: real_type = KIND(1.d0)

!---------------------------------------------------------------------------------
! some fixed parameter
!---------------------------------------------------------------------------------
  real(real_type), parameter :: version = 0.1
  real(real_type) :: pi = 3.14159265358979323846264338327950588 ! will be set below
!---------------------------------------------------------------------------------
! switches for configuration
!---------------------------------------------------------------------------------
  logical :: enable_upwind3_advection  = .false.
  logical :: enable_dst3_advection     = .false.
  logical :: enable_superbee_advection = .false.
  logical :: enable_multidim_advection = .false.
  logical :: enable_AB3_time_stepping  = .false.
  logical :: enable_particles          = .false.
  logical :: enable_v_velocity         = .false.
!---------------------------------------------------------------------------------
! model parameter
!---------------------------------------------------------------------------------
  integer:: nx,nz    ! number of grid points   
  real(real_type) :: Lx,Lz    ! extent of domain in m
  real(real_type) :: dx,dz,dt ! extent of grid cell in m, time step
  real(real_type),parameter :: g=9.81, rho_a=1.,rho_o = 1000.
  real(real_type),parameter :: r1 = (rho_a + rho_o)/(2*rho_a*rho_o)
  real(real_type),parameter :: r2 = (rho_a - rho_o)/(2*rho_a*rho_o)
  real(real_type) :: phi_diff = 0d0, gamma = 0.5, Coriolis = 0d0
  real(real_type) :: nu_o=0d0, nu_a=0d0
!---------------------------------------------------------------------------------
! time stepping parameter
!---------------------------------------------------------------------------------
  integer ::  taum2 = 1, taum1 = 2, tau = 3 !  time step index
 
  real(real_type), parameter :: AB3_a=  23d0/12d0 , AB3_b = -16d0/12d0, AB3_c = 5d0/12d0
  integer :: itt = 0                         ! current time step number
  real(real_type) :: time,runlen = 0.        ! current time in s, length of integration in s

!---------------------------------------------------------------------------------
! model fields
!---------------------------------------------------------------------------------
  real(real_type), allocatable :: u(:,:),du(:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: v(:,:),dv(:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: w(:,:),dw(:,:,:)    ! velocity and tendency
  real(real_type), allocatable :: phi(:,:),dphi(:,:,:)    ! order variables
  real(real_type), allocatable :: p(:,:,:)                ! pressure
  
  real(real_type), allocatable :: fe(:,:), ft(:,:)
  real(real_type), allocatable :: irho_ip(:,:),irho_im(:,:),irho_kp(:,:),irho_km(:,:)
  real(real_type), allocatable :: precond(:,:)
!---------------------------------------------------------------------------------
!     Parallel domain setup
!---------------------------------------------------------------------------------
 
  integer :: n_pes     ! total number of processors
  integer :: my_pe     ! index of this processor from 0 to n_pes-1
  integer :: n_pes_i   ! total number of processors in x direction
  integer :: n_pes_k
  integer :: i_blk,k_blk
  integer :: my_blk_i  ! index of this processor in x direction from 1 to n_pes_i
  integer :: my_blk_k
  integer :: is_pe     ! start index of grid points in x direction of this processor
  integer :: ie_pe     ! end index of grid points in x direction of this processor
  integer :: ks_pe
  integer :: ke_pe
  integer :: onx=2     ! number of overlapping points in x and y direction
 
 !---------------------------------------------------------------------------------
! conjugate gradient solver related parameter
!---------------------------------------------------------------------------------
  integer :: max_iterations=5000,itts ! maximal and actual iterations of Poisson solver
  real(real_type)  :: congr_eps=1e-4 ! epsilon cut off criterion for solver
 
!---------------------------------------------------------------------------------
! diagnostic setup
!---------------------------------------------------------------------------------
 real(real_type) :: tsmonint = 0. , snapint = 0., meanint = 0.
end module main_module



subroutine  allocate_main_module
 use main_module
 implicit none
 
 pi = acos(0d0)*2d0
 
 allocate(u(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(v(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(w(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(phi(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate(p(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 u=0;v=0;w=0;phi=0;p=0
 allocate(du(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(dv(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(dw(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 allocate(dphi(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx,3) )
 dphi=0;du=0;dv=0;dw=0
 allocate( fe(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 allocate( ft(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) )
 fe=0;ft=0
 allocate( irho_ip(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ) 
 allocate( irho_im(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ) 
 allocate( irho_kp(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ) 
 allocate( irho_km(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ) 
 irho_ip = 0;irho_im=0; irho_kp=0;irho_km=0
 allocate( precond(is_pe-onx:ie_pe+onx,ks_pe-onx:ke_pe+onx) ) 
 precond= 0
end subroutine  allocate_main_module
