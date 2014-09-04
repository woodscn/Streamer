module UCS_tester
  use grid_motion_driver_mod
  use Godunovdriver
  use TimeAdvancementStuff
  ! Check that the combined code correctly solves the 2-D riemann problem
  ! with various schemes for grid motion.
  real(8), dimension(21) :: top_state
  real(8), dimension(21) :: bottom_state
  real(8), dimension(21) :: base_state
  real(8), dimension(11) :: SRS_Exact_init
contains

  integer function UCS_test_main()
    implicit none
!    UCS_test_main = Steady2DRiemann(100,5000,.00005d0)
    UCS_test_main = ObliqueShock(100,5000,.0001d0)
  end function UCS_test_main

  function Steady2DRiemann(ny,nt,dt)
    implicit none
    integer :: Steady2DRiemann
    integer, intent(in) :: ny, nt
    real(8), intent(in) :: dt
    real(8), dimension(:,:,:,:), allocatable :: main, main2
    real(8), dimension(:,:,:), allocatable :: new_column
    integer :: nx, nz
    integer :: i, j, k, l, m, n
    real(8) :: dx, dy, dz, dxi, deta, dzeta, dt_out
    integer, dimension(1000) :: opts
    logical :: check_remove_column

    
    ! Set simulation options
    opts = 0
    opts(1:5) = [2,0,1,1,1]
    opts(6:7) = [1,4]
    opts(101) = 1
    opts(102) = 2
    opts(104) = 1
    ! Initialize the main data array
    nx = 1; nz = 1
    dx = .6d0/((ny*6)/10); dy = 1d0/(ny); dz = dx
    dxi = 1d0; deta = 1d0; dzeta = 1d0
    top_state = 0d0; bottom_state = 0d0; base_state = 0d0
    top_state(1:5) = [.25d0,.5d0,7d0*sqrt(.25d0/.5d0*1.4d0),0d0,0d0]
    bottom_state(1:5) = [1d0,1d0,2.4d0*sqrt(1.4d0),0d0,0d0]
    base_state(6:21) = [dx/dxi,0d0,0d0,0d0,dy/deta,0d0,0d0,0d0,dz/dzeta,&
         0d0,0d0,0d0,0d0,0d0,0d0,dx*dy*dz/(dxi*deta*dzeta)]
    allocate(main(21,nx+2,ny+2,nz+2))
    do k = 2, nz+1
       do j = 2, ny+1
          do i = 2, nx+1
             if(j>ny/2)then
                main(:,i,j,k) = top_state+base_state
             else
                main(:,i,j,k) = bottom_state+base_state
             end if
             main(18,i,j,k) = (i-2)*dx
             main(19,i,j,k) = (j-2)*dy
             main(20,i,j,k) = (k-2)*dz
          end do
       end do
    end do
    call write_files_matlab(main(:,2:nx+1,2:ny+1,2:nz+1),0d0,nx,ny,&
         first_flag=.true.)
    do n = 1, nt
       do m = 1, 3
          opts(2) = m
          ! Apply boundary conditions
          call SteadyRiemannBCs(main)
          ! Set grid motion
          call grid_motion_driver(main,opts)
!          main(15:17,:,:,:) = main(3:5,:,:,:)*.999
          ! Advance flow variables
          call prim_update(main,dt_out,dt,.25d0,nx,ny,nz,opts)
!       if(CheckCreateColumn(main(:,2,2:ny+1,2:nz+1),&
!            main(:,1,2:ny+1,2:nz+1),nx,ny))then
       end do
       if(maxval(main(18,2,2:ny+1,2:nz+1))>.01d0)then
          write(*,*) "Creating a column"
          allocate(main2(21,nx+2,ny+2,nz+2))
          main2 = main
          deallocate(main)
          nx = nx + 1
          allocate(main(21,nx+2,ny+2,nz+2))
          main(:,2:nx+2,:,:) = main2
          main(:,1,:,:) = main2(:,1,:,:)
          deallocate(main2)
          call SteadyRiemannBCs(main)
!          allocate(new_column(21,ny,nz))
!          call CreateColumn(main(:,3,2:ny+1,2:nz+1),&
!               main(:,1,2:ny+1,2:nz+1),ny,nz,new_column)
          call SteadyRiemannBCs(main)
!          main(:,2,2:ny+1,2:nz+1) = main(:,1,2:ny+1,2:nz+1)
!          main(6,2,2:ny+1,2:nz+1) = main(18,2,2:ny+1,2:nz+1)/dxi

          call SteadyRiemannBCs(main)
!          deallocate(new_column)
       end if
       check_remove_column = .false.
       do k = 2, nz+1
          do j = 2, ny+1
             check_remove_column = (check_remove_column .or. &
                  main(18,nx+1,j,k)>.8d0)
          end do
       end do
       if(check_remove_column)then
          write(*,*) "Removing a column"
          allocate(main2(21,nx+2,ny+2,nz+2))
          main2 = main
          deallocate(main)
          nx = nx - 1
          allocate(main(21,nx+2,ny+2,nz+2))
          main(:,:,:,:) = main2(:,1:nx+2,:,:)
          deallocate(main2)
       end if
       if(mod(n,10)==0)call write_files_matlab(main(:,2:nx+1,2:ny+1,2),&
            0d0,nx,ny,first_flag=.false.)
    end do
    Steady2DRiemann = 1
  end function Steady2DRiemann

  subroutine SteadyRiemannBCs(main)
    implicit none
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer :: j, k

    main(:,1,:,:) = main(:,2,:,:)
    do j = 2,size(main,3)-1
       do k = 2, size(main,4)-1
          if(j > (size(main,3)-2)/2+1)then
             main(:,1,j,k) = top_state+base_state
          else
             main(:,1,j,k) = bottom_state+base_state
          end if
          main(18,1,j,k) = 0d0
          main(19:20,1,j,k) = main(19:20,2,j,k)
       end do
    end do
    main(:,size(main,2),:,:) = main(:,size(main,2)-1,:,:)
    main(:,:,1,:) = main(:,:,2,:)
    main(:,:,size(main,3),:) = main(:,:,size(main,3)-1,:)
    main(:,:,:,1) = main(:,:,:,2)
    main(:,:,:,size(main,4)) = main(:,:,:,size(main,4)-1)
  end subroutine SteadyRiemannBCs

!!$  function SteadyRiemannSolutionInit(pL,dL,ML,alphaL,pR,dR,MR,alphaR,gamma)
!!$    implicit none
!!$    real(8), intent(in) :: pL, dL, ML, alphaL, pR, dR, MR, alphaR, gamma
!!$    integer SteadyRiemannSolutionInit
!!$    real(8) :: pstar, alphastar, dstarL, dstarR, MstarL, MstarR
!!$    real(8), dimension(5) :: angles
!!$
!!$
!!$    call SRS_star_state(pL,dL,ML,alphaL,pR,dR,MR,alphaR,gamma,&
!!$         pstar,alphastar,dstarL,dstarR,MstarL,MstarR)
!!$
!!$    if(pstar/pL<=1d0)then
!!$       angles(2) = SRS_left_fan_angle(MstarL,alphastar)
!!$       angles(1) = SRS_left_fan_angle(ML,alphaL)
!!$    else
!!$       angles(2) = 0d0
!!$       angles(1) = SRS_left_shock_angle(pstar,pL,alphaL,ML)
!!$    end if
!!$    if(pstar/pR<=1d0)then
!!$       angles(5) = SRS_right_fan_angle(MstarR,alphastar)
!!$       angles(4) = SRS_right_fan_angle(MR,alphaR)
!!$    else
!!$       angles(5) = SRS_right_shock_angle(pstar,pR,alphaR,MR)
!!$       angles(4) = 0d0
!!$    end if
!!$
!!$    SRS_Exact_Init = [pstar,alphastar,dstarL,dstarR,MstarL,MstarR,&
!!$         angles(1),angles(2),angles(3),angles(4),angles(5)]
!!$
!!$    SteadyRiemannSolutionInit = 0
!!$
!!$  end function SteadyRiemannSolution
!!$
!!$  function SRS_star_state(pL,dL,ML,alphaL,pR,dR,MR,alphaR,gamma,&
!!$       pstar,alphastar,dstarL,dstarR,MstarL,MstarR)
!!$
!!$  end function SRS_star_state
!!$
!!$  function SRS_left_fan_angle(M,alpha)
!!$
!!$  end function SRS_left_fan_angle
!!$
!!$  function SRS_left_shock_angle(pstar,p_0,alpha,M)
!!$
!!$  end function SRS_left_shock_angle
!!$  
!!$  

  function ObliqueShock(ny,nt,dt)
    implicit none
    integer :: ObliqueShock
    integer, intent(in) :: ny, nt
    real(8), intent(in) :: dt
    real(8), dimension(:,:,:,:), allocatable :: main, main2
    real(8), dimension(:,:,:), allocatable :: new_column
    integer :: nx, nz
    integer :: i, j, k, l, m, n
    real(8) :: dx, dy, dz, dxi, deta, dzeta, dt_out
    integer, dimension(1000) :: opts
    logical :: check_remove_column

    
    ! Set simulation options
    opts = 0
    opts(1:5) = [2,0,1,1,1]
    opts(6:7) = [1,4]
    opts(101) = 1
    opts(102) = 2
    opts(104) = 1
    ! Initialize the main data array
    nx = 1; nz = 1
    dx = 1d0/ny; dy = dx; dz = dx
    dxi = 1d0; deta = 1d0; dzeta = 1d0
    top_state = 0d0; bottom_state = 0d0; base_state = 0d0
    base_state = [1d0,1d0,1.8d0*sqrt(1.4d0),0d0,0d0,&
         dx/dxi,0d0,0d0,0d0,dy/deta,0d0,0d0,0d0,dz/dzeta,&
         0d0,0d0,0d0,0d0,0d0,0d0,dx*dy*dz/(dxi*deta*dzeta)]
    allocate(main(21,nx+2,ny+2,nz+2))
    do k = 2, nz+1
       do j = 2, ny+1
          do i = 2, nx+1
             if(j>ny/2)then
                main(:,i,j,k) = top_state+base_state
             else
                main(:,i,j,k) = bottom_state+base_state
             end if
             main(15:17,:,:,:) = main(3:5,:,:,:)*.999d0
             main(18,i,j,k) = (i-2)*dx
             main(19,i,j,k) = (j-1.5)*dy
             main(20,i,j,k) = (k-2)*dz
          end do
       end do
    end do
    call write_files_matlab(main(:,2:nx+1,2:ny+1,2:nz+1),0d0,nx,ny,&
         first_flag=.true.)
    do n = 1, nt
       do m = 1, 3
          opts(2) = m
          ! Apply boundary conditions
          call ObliqueShockBCs(main,1.8d0,1.4d0)
          ! Set grid motion
          call grid_motion_driver(main,opts)
!          main(15:17,:,:,:) = main(3:5,:,:,:)*.999
          ! Advance flow variables
          call prim_update(main,dt_out,dt,.25d0,nx,ny,nz,opts)
!       if(CheckCreateColumn(main(:,2,2:ny+1,2:nz+1),&
!            main(:,1,2:ny+1,2:nz+1),nx,ny))then
       end do
       if(maxval(main(18,2,2:ny+1,2:nz+1))>.01d0)then
          write(*,*) "Creating a column"
          allocate(main2(21,nx+2,ny+2,nz+2))
          main2 = main
          deallocate(main)
          nx = nx + 1
          allocate(main(21,nx+2,ny+2,nz+2))
          main(:,2:nx+2,:,:) = main2
          main(:,1,:,:) = main2(:,1,:,:)
          deallocate(main2)
          call ObliqueShockBCs(main,1.8d0,1.4d0)
!          allocate(new_column(21,ny,nz))
!          call CreateColumn(main(:,3,2:ny+1,2:nz+1),&
!               main(:,1,2:ny+1,2:nz+1),ny,nz,new_column)
!          call SteadyRiemannBCs(main)
!          main(:,2,2:ny+1,2:nz+1) = main(:,1,2:ny+1,2:nz+1)
!          main(6,2,2:ny+1,2:nz+1) = main(18,2,2:ny+1,2:nz+1)/dxi
!          call SteadyRiemannBCs(main)
!          deallocate(new_column)
       end if
       check_remove_column = .false.
       do k = 2, nz+1
          do j = 2, ny+1
             check_remove_column = (check_remove_column .or. &
                  main(18,nx+1,j,k)>1d0)
          end do
       end do
       if(check_remove_column)then
          write(*,*) "Removing a column"
          allocate(main2(21,nx+2,ny+2,nz+2))
          main2 = main
          deallocate(main)
          nx = nx - 1
          allocate(main(21,nx+2,ny+2,nz+2))
          main(:,:,:,:) = main2(:,1:nx+2,:,:)
          deallocate(main2)
       end if
       if(mod(n,10)==0)call write_files_matlab(main(:,2:nx+1,2:ny+1,2),&
            0d0,nx,ny,first_flag=.false.)
    end do
    ObliqueShock = 1
  end function ObliqueShock

  subroutine ObliqueShockBCs(main,mach,gamma)
    implicit none
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer i,j,k
    real(8) :: theta, mach, gamma
    intent(in) :: mach, gamma
    real(8), parameter :: beta = PI*.25d0
    real(8) :: dx, dy, dz, dxi, deta, dzeta

    dx = .01; dy = .01; dz = .01
    dxi = 1; deta = 1; dzeta = 1
    ! Outflow condition
    main(:,size(main,2),:,:) = main(:,size(main,2)-1,:,:)
    ! Front & back conditions
    main(:,:,:,1) = main(:,:,:,2)
    main(:,:,:,size(main,4)) = main(:,:,:,size(main,4)-1)
    ! Top condition
    main(:,:,size(main,3),:) = main(:,:,size(main,3)-1,:)
    ! Bottom condition
    theta = atan(2d0/tan(beta)*(mach**2*sin(beta)**2-1)&
         /(mach**2*(gamma+cos(2d0*beta))+2d0))
    do i = 2, size(main,2)
       j = i
       if(main(18,i,2,2)>.1d0)then
          exit
       end if
    end do
    call wall_reflect(main(:,2:j-1,1,2),main(:,2:j-1,2,2),0d0)
    call wall_reflect(main(:,j:size(main,2),1,2),main(:,j:size(main,2),2,2)&
         ,theta)
    ! Inflow condition
    do i = 1,size(main,3)
       main(:,1,i,2) = [1d0,1d0,mach*sqrt(gamma),0d0,0d0,&
         dx/dxi,0d0,0d0,   0d0,dy/deta,0d0,   0d0,0d0,dz/dzeta,&
         .999d0*sqrt(gamma)*1.8d0,0d0,0d0,&
         0d0,0d0,0d0,dx*dy*dz/(dxi*deta*dzeta)]
       main(19,1,i,:) = (i-1.5)*dy
    end do
  end subroutine ObliqueShockBCs

  subroutine wall_reflect(bound,interior,angle)
    implicit none
    real(8), dimension(:,:), intent(in) :: interior
    real(8), dimension(:,:), intent(out) :: bound
    real(8), intent(in) :: angle

    bound(2,:) = interior(2,:)
    bound(1,:)   = interior(1,:)
    bound(3,:)   = cos(2.*angle)*interior(3,:) + sin(2.*angle)*interior(4,:)
    bound(4,:)   = sin(2.*angle)*interior(3,:) - cos(2.*angle)*interior(4,:)
    bound(15,:)  = cos(2.*angle)*interior(15,:)+ sin(2.*angle)*interior(16,:)
    bound(16,:)  = sin(2.*angle)*interior(15,:)- cos(2.*angle)*interior(16,:)
    if(abs(angle)<=.5*PI)then
       bound(6,:)   = interior(6,:)
       bound(10,:)   = (interior(10,:)*(interior(6,:)+bound(6,:)&
            *sin(angle)**2) - sin(angle)&
            *cos(angle)*(interior(6,:)*interior(9,:)+interior(7,:)&
            *interior(10,:)+bound(6,:)*interior(9,:)))&
            /(bound(6,:)*(cos(angle)**2 - &
            sin(angle)**2)-interior(6,:)&
            *sin(angle)**2+interior(7,:)*sin(angle)&
            *cos(angle))
       bound(9,:) = tan(angle)*(interior(10,:)+bound(10,:))&
            - interior(9,:)
       bound(7,:) = tan(angle)*(interior(6,:)+bound(6,:))&
            - interior(7,:)
    elseif(abs(angle)>.5*PI)then
       bound(7,:)   = interior(7,:)
       bound(9,:)   = (interior(9,:)*(interior(7,:)+bound(7,:)&
            *cos(angle)**2) - sin(angle)&
            *cos(angle)*(interior(6,:)*interior(9,:)+interior(7,:)&
            *interior(10,:)+bound(7,:)*interior(10,:)))&
            /(bound(7,:)*(sin(angle)**2 - &
            cos(angle)**2)-interior(7,:)&
            *cos(angle)**2+interior(6,:)*sin(angle)&
            *cos(angle))
       bound(6,:) = 1.d0/tan(angle)*(interior(7,:)+bound(7,:))&
            - interior(6,:)
       bound(10,:) = 1.d0/tan(angle)*(interior(9,:)+bound(9,:))&
            - interior(10,:)
    end if
    bound(14,:) = interior(14,:)
    bound(18:21,:) = interior(18:21,:)
    
  end subroutine wall_reflect

end module UCS_tester
