module Godunov
  use GeneralUtilities
  use Riemann
  implicit none
  real(8), parameter :: max_dt = 1.d0
!  real(8), dimension(7), parameter :: dxi_a = [1.d0, .5d0, .25d0, .2d0, &
!       2.d0, 4.d0, 5.d0]
!  real(8), dimension(7), parameter :: deta_a = [1.d0, .5d0, .25d0, .2d0, &
!       2.d0, 4.d0, 5.d0]
!  real(8), dimension(7), parameter :: dzeta_a = [1.d0, .5d0, .25d0, .2d0, &
!       2.d0, 4.d0, 5.d0]
  real(8) :: dxi, deta, dzeta, dxi_inv, deta_inv, dzeta_inv, dV_inv
!  real(8), parameter :: dxi   = 1.d0
!  real(8), parameter :: deta  = 1.d0
!  real(8), parameter :: dzeta = 1.d0
!  real(8), parameter :: dxi_inv   = 1.d0/dxi
!  real(8), parameter :: deta_inv  = 1.d0/deta
!  real(8), parameter :: dzeta_inv = 1.d0/dzeta
!  real(8), parameter :: dV_inv = dxi_inv*deta_inv*dzeta_inv
!  integer :: update_type = 1 ! 1 = FV, 2 = HUI3D

  interface primtocons
     module procedure primtoconsarray
     module procedure primtoconspoint
  end interface primtocons

  interface constoprim
     module procedure constoprimarray
     module procedure constoprimpoint
  end interface constoprim

contains
!  integer elemental function gt0(x)
!    implicit none
!    real(8), intent(in) :: x
!    gt0 = ishft( int(sign(1.d0,x) + 1) , -1 )
!  end function gt0

  ! Computes specific energy, given the array [ p, rho, u, v, w ]
  real(8) function energy_func(in)
    implicit none
    real(8), dimension(:), intent(in) :: in
    energy_func = 0.5d0*(in(3)**2+in(4)**2+in(5)**2) + in(1)/(in(2)*gamma2)
  end function energy_func

  subroutine primtoconsarray(main)
    ! Assumes the structure of prim(:) is :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ]-pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Returns the structure of cons(:) :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(4) :: en_shape
    real(8), allocatable, dimension(:,:,:) :: energy
    en_shape = shape(main)
    allocate(energy(en_shape(2), en_shape(3), en_shape(4)))
    energy = main(21,:,:,:)*( main(1,:,:,:)*gamma1 + 0.5d0*main(2,:,:,:)&
         *( main(3,:,:,:)**2 + main(4,:,:,:)**2 + main(5,:,:,:)**2 ) )
    main(1,:,:,:) = main(2,:,:,:)*main(21,:,:,:)
    main(2,:,:,:) = main(1,:,:,:)*main(3,:,:,:)
    main(3,:,:,:) = main(1,:,:,:)*main(4,:,:,:)
    main(4,:,:,:) = main(1,:,:,:)*main(5,:,:,:)
    main(5,:,:,:) = energy
    deallocate(energy)
  end subroutine primtoconsarray

  subroutine primtoconspoint(main)
    ! Assumes the structure of prim(:) is :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ] - pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Returns the structure of cons(:) :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:), intent(inout) :: main
    real(8) :: energy, J
    J = Jacobian(main(6:14))
    energy = J*( main(1)*gamma1 + 0.5d0*main(2)&
         *( main(3)**2 + main(4)**2 + main(5)**2 ) )
    main(1) = main(2)*J
    main(2) = main(1)*main(3)
    main(3) = main(1)*main(4)
    main(4) = main(1)*main(5)
    main(5) = energy
  end subroutine primtoconspoint
  
  subroutine constoprimarray(main)
    ! Returns the structure of prim(:) :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ]-pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Assume the structure of cons(:) is :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(4) :: p_shape
    real(8), allocatable, dimension(:,:,:) :: temp1, temp2, p
    allocate( temp1(p_shape(2), p_shape(3), p_shape(4)) )
    allocate( temp2(p_shape(2), p_shape(3), p_shape(4)) )
    allocate(     p(p_shape(2), p_shape(3), p_shape(4)) )
    temp1 = 1.d0/main(21,:,:,:)
    temp2 = 1.d0/main(1,:,:,:)
    p = gamma2*temp1*( main(5,:,:,:) - .5d0*temp2&
         *( main(2,:,:,:)**2 + main(3,:,:,:)**2 + main(4,:,:,:)**2 )&
         )
    main(5,:,:,:) = main(4,:,:,:)*temp2
    main(4,:,:,:) = main(3,:,:,:)*temp2
    main(3,:,:,:) = main(2,:,:,:)*temp2
    main(2,:,:,:) = main(1,:,:,:)*temp1
    main(1,:,:,:) = p
    deallocate( temp1, temp2, p )
  end subroutine constoprimarray

  subroutine constoprimpoint(main)
    ! Returns the structure of prim(:) :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ]-pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Assume the structure of cons(:) is :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:), intent(inout) :: main
    real(8) :: temp1, temp2, p, J
    J = Jacobian(main(6:14))
    temp1 = 1.d0/J
    temp2 = 1.d0/main(1)
    p = gamma2*temp1*( main(5) - .5d0*temp2&
         *( main(2)**2 + main(3)**2 + main(4)**2 )&
         )
    main(5) = main(4)*temp2
    main(4) = main(3)*temp2
    main(3) = main(2)*temp2
    main(2) = main(1)*temp1
    main(1) = p
  end subroutine constoprimpoint

  function invnorm3(in)
    real(8), dimension(3), intent(in) :: in
    real(8) :: invnorm3
    invnorm3 = 1.d0/sqrt( in(1)**2 + in(2)**2 + in(3)**2 )
  end function invnorm3

  subroutine grid_coords(grad, normal, tangential1, tangential2)
! Compute an orthonormal coordinate system given an initial, 
! unnormalized vector.
    real(8), dimension(:), intent(in) :: grad
    real(8), dimension(3), intent(out) :: normal, tangential1, tangential2
    real(8) :: temp1, temp3
    real(8), dimension(3) :: temp2, temp4
    temp1 = invnorm3(grad)
    normal = grad*temp1
    if( grad(2)**2 + grad(3)**2 < EPS )then
       tangential1 = (/ 0.d0, 1.d0, 0.d0 /)
       tangential2 = (/ 0.d0, 0.d0, 1.d0 /)
    else
       temp2 = (/ 0.d0, -grad(3), grad(2) /)
       temp3 = invnorm3(temp2)
       temp4 = (/ grad(2)**2 + grad(3)**2, -grad(1)*grad(2), -grad(1)*grad(3) /)
       tangential1 = temp2*temp3
       tangential2 = temp1*temp3*temp4
    end if
  end subroutine grid_coords

!!$  subroutine prim_update_HUI3D(main,dt_out,dt_in,CFL,nx,ny,nz,options)
!!$    implicit none
!!$    !f2py intent(in,out) :: main
!!$    !f2py integer :: nx, ny, nz
!!$    !f2py intent(in) :: nx, ny, nz
!!$    real(8), dimension(21,-1*options(201):nx+options(201)-1,&
!!$         -1*options(201):ny+options(201)-1,-1*options(201):nz+options(201)-1),&
!!$         intent(inout) :: main
!!$    real(8), dimension(21,0:nx-1,0:ny-1,0:nz-1) :: main_temp
!!$    real(8), intent(out) :: dt_out
!!$    real(8), intent(in), optional :: dt_in
!!$    real(8), intent(in), optional :: CFL
!!$    integer, intent(in) :: nx,ny,nz
!!$    ! Options values are used to activate specific routine options.
!!$    ! - Options(202) controls the spatial order of accuracy. 
!!$    ! - Options(203) controls grid motion.
!!$    ! - Options(204) controls the boundary conditions 
!!$    integer, dimension(:), intent(in) :: options
!!$    integer :: i, j, k, m, n, im, jm, km, ip, jp, kp, case_no
!!$    real(8) :: area, dv_inv, max_wave_speed_temp
!!$    real(8), dimension(5) :: cons, prim
!!$    real(8), dimension(5) :: left_interface, right_interface
!!$    real(8), dimension(5) :: left_flux, right_flux
!!$    real(8), dimension(3) :: grid_vel, grid_pos
!!$    real(8), dimension(3,3) :: row_ops_mat, vels_transform
!!$    real(8), dimension(9) :: metric, metric_inverse
!!$    real(8), dimension(21) :: temp, center, left, right, geom_avg
!!$    integer, parameter :: splitting_type = 1
!!$    real(8) :: dt, junk
!!$    real(8), dimension(4) :: riemann_middle_states
!!$    real(8), dimension(5) :: riemann_wave_speeds
!!$    integer :: spatial_order
!!$    integer :: grid_motion
!!$    spatial_order = options(202)
!!$    grid_motion = options(203)
!!$    select case(options(204))
!!$    case(1)
!!$       bc_func => 
!!$    dt = dt_in
!!$    dt_out = dt_in
!!$    dv_inv = 1.d0
!!$    do n = 1, 3
!!$       if( n .eq. 1 )then
!!$          case_no = 1
!!$          area = deta*dzeta
!!$          im =-1; ip = 1; jm = 0; jp = 0; km = 0; kp = 0
!!$       elseif( n .eq. 2 )then
!!$          case_no = 2
!!$          area = dxi*dzeta
!!$          jm =-1; jp = 1; im = 0; ip = 0; km = 0; kp = 0
!!$       elseif( n .eq. 3 )then
!!$          case_no = 3
!!$          area = dxi*deta
!!$          km =-1; kp = 1; im = 0; ip = 0; jm = 0; jp = 0
!!$       else
!!$          write(*,*) "Error in prim_update, invalid value for n"
!!$          stop
!!$       end if
!!$
!!$       call bc_func(main,size(main,2),size(main,3),size(main,4))
!!$       do k = 0, nz-1
!!$          do j = 0, ny-1
!!$             do i = 0, nx-1
!!$                row_ops_mat = row_ops_mat_func(case_no)
!!$                center = main(:,i,j,k)
!!$                left = main(:,i+im,j+jm,k+km)
!!$
!!$                if(spatial_order .eq. 2.and.(&
!!$                     (n==1.and.i>0.and.i<nx-1).or.&
!!$                     (n==2.and.j>0.and.j<ny-1).or.&
!!$                     (n==3.and.k>0.and.k<nz-1)))&
!!$                     call MUSCL_HUI(main(1:5,i+2*im,j+2*jm,k+2*km),&
!!$                     main(1:5,i+im,j+jm,k+km),main(1:5,i,j,k),&
!!$                     main(1:5,i+ip,j+jp,k+kp),left(1:5),center(1:5))
!!$
!!$                metric = .5d0*(left(6:14)+center(6:14))
!!$                metric_inverse = MetricInverse(metric)
!!$                grid_vel = .5d0*(left(15:17)+center(15:17))
!!$
!!$                vels_transform = matmul(&
!!$                     MetrictoMatrix(metric_inverse),row_ops_mat)
!!$                do m = 1, 3
!!$                   vels_transform(m,:) = vels_transform(m,:)&
!!$                        /sqrt(sum(vels_transform(m,:)**2))
!!$                end do
!!$
!!$                left(3:5) = matmul(vels_transform,left(3:5))
!!$                center(3:5) = matmul(vels_transform,center(3:5))
!!$                grid_vel = matmul(vels_transform,grid_vel)
!!$                geom_avg = 0d0
!!$                geom_avg(6:14) = center(6:14)
!!$                geom_avg(15:17) = grid_vel
!!$                call riemann_solve(left,center,n,1,[0d0],left_interface,&
!!$                     max_wave_speed_temp)
!!$
!!$                vels_transform = matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric))
!!$                do m = 1, 3
!!$                   vels_transform(m,:) = vels_transform(m,:)&
!!$                        /sqrt(sum(vels_transform(m,:)**2))
!!$                end do
!!$
!!$                left_interface(3:5) = matmul(vels_transform,left_interface(3:5))
!!$
!!$                center = main(:,i,j,k)
!!$                right = main(:,i+ip,j+jp,k+kp)
!!$                if(spatial_order .eq. 2.and.(&
!!$                     (n==1.and.i>0.and.i<nx-1).or.&
!!$                     (n==2.and.j>0.and.j<ny-1).or.&
!!$                     (n==3.and.k>0.and.k<nz-1)))&
!!$                     call MUSCL_HUI(main(1:5,i+im,j+jm,k+km),&
!!$                     main(1:5,i,j,k),main(1:5,i+ip,j+jp,k+kp),&
!!$                     main(1:5,i+2*ip,j+2*jp,k+2*kp),center(1:5),right(1:5))
!!$
!!$                metric = .5d0*(center(6:14)+right(6:14))
!!$                metric_inverse = MetricInverse(metric)
!!$                grid_vel = .5d0*(center(15:17)+right(15:17))
!!$
!!$                vels_transform = matmul(&
!!$                     MetrictoMatrix(metric_inverse),row_ops_mat)
!!$                do m = 1, 3
!!$                   vels_transform(m,:) = vels_transform(m,:)&
!!$                        /sqrt(sum(vels_transform(m,:)**2))
!!$                end do
!!$
!!$                right(3:5) = matmul(vels_transform,right(3:5))
!!$                center(3:5) = matmul(vels_transform,center(3:5))
!!$                grid_vel = matmul(vels_transform,grid_vel)
!!$                geom_avg = 0d0
!!$                geom_avg(6:14) = center(6:14)
!!$                geom_avg(15:17) = grid_vel
!!$                call riemann_solve(center,right,n,1,[0d0],right_interface,&
!!$                     max_wave_speed_temp)
!!$                
!!$                vels_transform = matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric))
!!$                do m = 1, 3
!!$                   vels_transform(m,:) = vels_transform(m,:)&
!!$                        /sqrt(sum(vels_transform(m,:)**2))
!!$                end do
!!$
!!$                right_interface(3:5)=matmul(vels_transform,right_interface(3:5))
!!$
!!$                junk = sum(matmul(matmul(MetrictoMatrix(metric_inverse),&
!!$                     row_ops_mat),matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric)))&
!!$                     -reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[3,3])&
!!$                     **2)
!!$                if(junk > 1d-14)then
!!$                   write(*,*) "Inverse not working!!"
!!$                   write(*,*) junk
!!$                   stop
!!$                end if
!!$                     
!!$                center = main(:,i,j,k)
!!$                
!!$                if(grid_motion .eq. 1)then
!!$                   if(n==1)then
!!$                      center(6:8) = center(6:8) + .25d0*dt*area*dv_inv*&
!!$                           (right_interface(3:5)-left_interface(3:5))
!!$                   else if(n==2)then
!!$                      center(9:11) = center(9:11) + .25d0*dt*area*dv_inv*&
!!$                           (right_interface(3:5)-left_interface(3:5))
!!$                   else if(n==3)then
!!$                      center(12:14) = center(12:14) + .25d0*dt*area*dv_inv*&
!!$                           (right_interface(3:5)-left_interface(3:5))
!!$                   end if
!!$                end if
!!$                left_flux  = flux( left_interface,center,n)
!!$                right_flux = flux(right_interface,center,n)
!!$                call primtocons(center)
!!$                center(1:5) = center(1:5) - dt*area*dv_inv*&
!!$                     (right_flux-left_flux)
!!$                call constoprim(center)
!!$                if(grid_motion .eq. 0)then
!!$                   center(15:17) = 0d0
!!$                elseif(grid_motion .eq. 1)then
!!$                   center(15:17) = center(3:5)*.25d0
!!$                end if
!!$                center(18:20) = center(18:20) + dt*center(15:17)
!!$                center(21) = Jacobian(center(6:14))
!!$
!!$                main_temp(:,i,j,k) = center
!!$             end do
!!$          end do
!!$       end do
!!$       main(:,0:nx-1,0:ny-1,0:nz-1) = main_temp
!!$       
!!$!       main(15:17,0:nx-1,0:ny-1,0:nz-1) = main(3:5,0:nx-1,0:ny-1,0:nz-1)*.0d0
!!$!       main(18:20,0:nx-1,0:ny-1,0:nz-1) = main(18:20,0:nx-1,0:ny-1,0:nz-1)&
!!$!            *dt*area*dv_inv*main(15:17,0:nx-1,0:ny-1,0:nz-1)
!!$    end do
!!$  end subroutine prim_update_HUI3D

  subroutine prim_update_FV(main,dt_out,dt_in,CFL,nx,ny,nz,opts)
! Advance the solution using the integral form of the equations
! This subroutine assumes that main is the full array of primitive variables. 
! main must also include the boundary cell values. That is, main contains a 
! 0-index and an n + 1 index containing the contents prescribed by the 
! boundary conditions.
    implicit none
    integer, dimension(:), intent(in) :: opts
!!$    type(prim_update_FV_options), intent(in) :: opts
!!$    real(8), dimension(21,&
!!$         -1*opts%ghost_pts:nx+opts%ghost_pts-1,&
!!$         -1*opts%ghost_pts:ny+opts%ghost_pts-1,&
!!$         -1*opts%ghost_pts:nz+opts%ghost_pts-1),intent(inout) :: main
    real(8), dimension(21,-1*opts(101):nx+opts(101)-1,&
         -1*opts(101):ny+opts(101)-1,-1*opts(101):nz+opts(101)-1),&
         intent(inout) :: main
!f2py intent(in,out) :: main
    real(8), intent(out) :: dt_out
    real(8), intent(in), optional :: dt_in
    real(8), intent(in), optional :: CFL
    integer, intent(in) :: nx,ny,nz
    integer :: spatial_order
    integer :: grid_motion
    integer :: time_step_scheme
    integer, dimension(2) :: splitting
    integer :: i, j, k
!    real(8), dimension(14,3,nx+1,ny+1,nz+1) :: fluxes
    real(8), dimension(14,0:nx,0:ny-1,0:nz-1) :: fluxx
    real(8), dimension(14,0:nx-1,0:ny,0:nz-1) :: fluxy
    real(8), dimension(14,0:nx-1,0:ny-1,0:nz) :: fluxz
    real(8), dimension(14,nx,ny,nz) :: temp
    real(8) :: max_wave_speed, dt, max_dt_grid
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8), dimension(3) :: GradX, GradY, GradZ
    real(8), dimension(9) :: geom_avg
    real(8), dimension(3,3) :: dXidX, dXdXi
    real(8), dimension(5) :: interface_vars
    real(8), dimension(21) :: StateL, StateR
    real(8), dimension(14) :: junk
    splitting = opts(1:2)
    spatial_order = opts(102)
    grid_motion = opts(103)
    time_step_scheme = opts(104)
    ! Set module values for dxi, deta, dzeta, based on opts.
    dxi = dxi_a(opts(3))
    deta = deta_a(opts(4))
    dzeta = dzeta_a(opts(5))
    dxi_inv = 1.d0/dxi; deta_inv = 1.d0/deta; dzeta_inv = 1.d0/dzeta
    dV_inv = dxi_inv*deta_inv*dzeta_inv
    
  ! Riemann_solve expects the left and right states to express velocities in
  ! grid-oriented components: normal, tangential, tangential.
    fluxx = 0d0; fluxy = 0d0; fluxz = 0d0
!! Grid motion has been moved to grid_motion_driver.f90
!    select case(grid_motion)
!    case(0)
!       main(15:17,:,:,:) = 0d0
!    case(1)
!       main(15:17,:,:,:) = .999d0*main(3:5,:,:,:)
!    case(2)
!       main(15:17,:,:,:) = .5
!    case default
!       write(*,*) "Invalid grid motion specification!"
!       stop
!    end select
    if(splitting(1).eq. 1 .or.(splitting(1).eq. 2 .and.splitting(2).eq. 2 ))&
         then
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx
                StateL = main(:,i-1,j,k)
                StateR = main(:,i,j,k)
                if(spatial_order .eq. 2 .and. i > 0 .and. i < nx)&
                     call MUSCL_HUI(main(1:5,i-2,j,k),main(1:5,i-1,j,k),&
                     main(1:5,i,j,k),main(1:5,i+1,j,k),StateL(1:5),StateR(1:5))
                call compute_fluxes_FV(StateL,StateR,fluxx(:,i,j,k),1,&
                     max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
             end do
          end do
       end do
    end if
    if(splitting(1).eq. 1 .or.(splitting(1).eq. 2 .and.splitting(2).eq. 2 ))&
         then
       do k = 0, nz-1
          do j = 0, ny
             do i = 0, nx-1
                StateL = main(:,i,j-1,k)
                StateR = main(:,i,j,k)
                if(spatial_order .eq. 2 .and. j > 0 .and. j < ny)&
                     call MUSCL_HUI(main(1:5,i,j-2,k),main(1:5,i,j-1,k),&
                     main(1:5,i,j,k),main(1:5,i,j+1,k),StateL(1:5),StateR(1:5))
                call compute_fluxes_FV(StateL,StateR,fluxy(:,i,j,k),2,&
                     max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
             end do
          end do
       end do
    end if
    if(splitting(1).eq. 1 .or.(splitting(1).eq. 2 .and.splitting(2).eq. 3))then
       do k = 0, nz
          do j = 0, ny-1
             do i = 0, nx-1
                StateL = main(:,i,j,k-1)
                StateR = main(:,i,j,k)
                if(spatial_order .eq. 2 .and. k > 0 .and. k < nz)&
                     call MUSCL_HUI(main(:,i,j,k-2),main(:,i,j,k-1),&
                     main(:,i,j,k),main(:,i,j,k+1),StateL,StateR)
                call compute_fluxes_FV(StateL,StateR,fluxz(:,i,j,k),3,&
                     max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
             end do
          end do
       end do
    end if
    max_wave_speed = max(max_wave_speed,EPS)
    select case(time_step_scheme)
    case(0)
       write(*,*) "Using given timestep", time_step_scheme
       dt = dt_in
    case(1)
       dt = min(CFL/max_wave_speed,dt_in)
    case default
       write(*,*) "Error in Godunov prim_update, invalid time step flag!"
       stop
    end select
    dt_out = dt
!    write(*,*) "Got here -- prim_update"
    call primtocons(main(:,0:nx-1,0:ny-1,0:nz-1))
    main(1:14,0:nx-1,0:ny-1,0:nz-1) = main(1:14,0:nx-1,0:ny-1,0:nz-1) - (&
         (fluxx(:,1:nx,:,:)-fluxx(:,0:nx-1,:,:))*deta*dzeta + &
         (fluxy(:,:,1:ny,:)-fluxy(:,:,0:ny-1,:))*dxi*dzeta + &
         (fluxz(:,:,:,1:nz)-fluxz(:,:,:,0:nz-1))*dxi*deta &
         )*dt*dV_inv
! Update extra variables
  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1
           main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
        end do
     end do
  end do
    call constoprim(main(:,0:nx-1,0:ny-1,0:nz-1))
! Update grid position
    main(18:20,0:nx-1,0:ny-1,0:nz-1) = main(18:20,0:nx-1,0:ny-1,0:nz-1) + &
         main(15:17,0:nx-1,0:ny-1,0:nz-1)*dt
! Update extra variables
  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1
           main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
        end do
     end do
  end do
  end subroutine prim_update_FV

  subroutine compute_fluxes_FV(inL, inR, flux_vec, case_no,&
      max_wave_speed,dt,dV_in,debug_flag)
    implicit none
    real(8), dimension(21), intent(in) :: inL, inR
    real(8), dimension(:), intent(out) :: flux_vec
    integer, intent(in) :: case_no
    real(8), intent(inout) :: max_wave_speed
    real(8), intent(in) :: dt
    real(8), dimension(3),intent(in) :: dV_in
    logical, intent(in), optional :: debug_flag

    real(8), dimension(21) :: StateL, StateR
    real(8), dimension(21) :: geom_avg
    real(8) :: dA, dV_inv
    real(8), dimension(5) :: interface_vars
    real(8), dimension(3) :: GradXi, GradEta, GradZeta
    real(8), dimension(3) :: GradX, GradY, GradZ
    real(8), dimension(3,3) :: dX_dXi_u, dXi_dX_u
    real(8) :: temp_wave_speed
    integer, dimension(3,3) :: row_ops_mat
    real(8), dimension(9) :: metric, metric_inverse
    real(8), dimension(3,3) :: vels_transform
    real(8), dimension(3) :: grid_vel
    integer :: m
    StateL = inL
    StateR = inR
    geom_avg = .5d0*(inL+inR)
    flux_vec = 0.d0
    dV_inv = 1.d0/(product(dV_in))
    row_ops_mat = row_ops_mat_func(case_no)
    metric = geom_avg(6:14)
    metric_inverse = MetricInverse(metric)

    vels_transform = matmul(MetrictoMatrix(metric_inverse),row_ops_mat)
    do m = 1, 3
       vels_transform(m,:) = vels_transform(m,:)&
            /sqrt(sum(vels_transform(m,:)**2))
    end do
    
    StateL(3:5) = matmul(vels_transform,StateL(3:5))
    StateR(3:5) = matmul(vels_transform,StateR(3:5))
    StateL(15:17) = matmul(vels_transform,StateL(15:17))
    StateR(15:17) = matmul(vels_transform,StateR(15:17))
    
    call riemann_solve(StateL,StateR,case_no,1,[0d0],interface_vars,&
         temp_wave_speed)

    vels_transform = matmul(transpose(row_ops_mat),MetrictoMatrix(metric))
    do m = 1, 3
       vels_transform(m,:) = vels_transform(m,:)&
            /sqrt(sum(vels_transform(m,:)**2))
    end do
    interface_vars(3:5) = matmul(vels_transform, interface_vars(3:5))
!    geom_avg(15:17) = .25d0*interface_vars(3:5)
    flux_vec(1:5) = flux(interface_vars,geom_avg,case_no)
    select case(case_no)
    case(1)
       flux_vec(6:8) = -geom_avg(15:17)
       dA = dV_in(2)*dV_in(3)
    case(2)
       flux_vec(9:11) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(3)
    case(3)
       flux_vec(12:14) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(2)
    case default
       write(*,*) "Invalid case_no in compute_fluxes -- case_no = ",case_no
       stop
    end select
!    flux_vec(6:14)=0d0
    max_wave_speed = max(max_wave_speed,temp_wave_speed)
  end subroutine compute_fluxes_FV

  function row_ops_mat_func(case_no)
    implicit none
    integer, dimension(3,3) :: row_ops_mat_func, row_ops_mat
    integer, intent(in) :: case_no
        select case(case_no)
    case(1)
       row_ops_mat = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    case(2)
       row_ops_mat = reshape([0,1,0,1,0,0,0,0,1],[3,3])
    case(3)
       row_ops_mat = reshape([0,0,1,0,1,0,1,0,0],[3,3])
    case default
       write(*,*) "Invalid case_no in compute_fluxes -- case_no = ",case_no
       stop
    end select
    row_ops_mat_func = row_ops_mat
  end function row_ops_mat_func

  function flux(in,geom_avg,case_no)
    implicit none
    real(8), dimension(:), intent(in) :: in, geom_avg
    integer, intent(in) :: case_no
    real(8), dimension(5) :: flux
    real(8) :: D, grads(3,3), grad(3), J, Jinv
!!$    write(*,*) "Got this far"
    J = Jacobian(geom_avg(6:14))
    Jinv = 1.d0/J
    call ComputationalGrads(geom_avg(6:14),J,grads(:,1),grads(:,2),grads(:,3))
    grad = grads(:,case_no)
    D = sum(([ in(3), in(4), in(5) ] - &
         [ geom_avg(15), geom_avg(16), geom_avg(17) ])*grad)
    flux = [&
         in(2)*J*D, &
         in(2)*J*D*in(3) + J*grad(1)*in(1), &
         in(2)*J*D*in(4) + J*grad(2)*in(1), &
         in(2)*J*D*in(5) + J*grad(3)*in(1), &
         in(2)*J*D*energy_func(in(1:5)) + J*sum(grad*in(3:5))*in(1) ]
!!$    write(*,*) "J = ", J
!!$    write(*,*) "D = ", D
!!$    write(*,*) "e = ", energy_func(in(1:5))
!!$    write(*,*) "v = ", in(3:5)
!!$    write(*,*) "V = ", geom_avg(15:17)
!!$    write(*,*) "grad = ", grad
!!$    write(*,*) "Energy term = ", J*in(2)*D*energy_func(in(1:5))
!!$    write(*,*) "Case No. = ", case_no
  end function flux

  subroutine MUSCL_HUI(leftleft, left, right, rightright, outleft, outright)
    implicit none
    real(8), dimension(5), intent(in) :: leftleft, left, right, rightright
    real(8), dimension(5), intent(out) :: outleft, outright
    integer :: n
    
    do n = 1, 5
       outleft(n)  =  left(n) + 5d-1*( left(n) - leftleft(n))*minmod(&
            (right(n)-left(n))/( left(n) - leftleft(n)))
       outright(n) = right(n) - 5d-1*(rightright(n)-right(n))*minmod(&
            (right(n)-left(n))/(rightright(n)-right(n)))
    end do
    outleft  =  left + 5d-1*( left - leftleft)*&
         minmod((right-left)/( left - leftleft))
    outright = right - 5d-1*(rightright-right)*&
         minmod((right-left)/(rightright-right))
    
  contains
    elemental function minmod(in)
      implicit none
      real(8), intent(in) :: in
      real(8) :: minmod
      minmod = max(0d0,min(1d0,in))
    end function minmod
  end subroutine MUSCL_HUI

  subroutine FreeExitConditions(main,nx,ny,nz)
    implicit none
    real(8), dimension(21,nx,ny,nz) :: main
    integer, intent(in) :: nx, ny, nz

!!$    nx = size(main,2)
!!$    ny = size(main,3)
!!$    nz = size(main,4)
    main(:, 1,:,:) = main(:,   2,:,:)
    main(:,nx,:,:) = main(:,nx-1,:,:)
    main(:,:, 1,:) = main(:,:,   2,:)
    main(:,:,ny,:) = main(:,:,ny-1,:)
    main(:,:,:, 1) = main(:,:,:,   2)
    main(:,:,:,nz) = main(:,:,:,nz-1)
  end subroutine FreeExitConditions

end module Godunov
