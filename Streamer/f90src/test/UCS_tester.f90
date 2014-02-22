module UCS_tester
  use grid_motion_driver_mod
  use Godunovdriver
  use TimeAdvancementStuff
  ! Check that the combined code correctly solves the 2-D riemann problem
  ! with various schemes for grid motion.
  real(8), dimension(21) :: top_state
  real(8), dimension(21) :: bottom_state
  real(8), dimension(21) :: base_state
contains

  integer function UCS_test_main()
    implicit none
    UCS_test_main = Steady2DRiemann(100,5000,.0001d0)
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
    opts(1:5) = [1,0,1,1,1]
    opts(6:7) = [0,4]
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
       ! Apply boundary conditions
       call SteadyRiemannBCs(main)
       ! Set grid motion
!       call grid_motion_driver(main,opts)
       ! Advance flow variables
       call prim_update(main,dt_out,dt,.25d0,nx,ny,nz,opts)
!       if(CheckCreateColumn(main(:,2,2:ny+1,2:nz+1),&
!            main(:,1,2:ny+1,2:nz+1),nx,ny))then
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
!          call SteadyRiemannBCs(main)
!          allocate(new_column(21,ny,nz))
!          call CreateColumn(main(:,3,2:ny+1,2:nz+1),&
!               main(:,1,2:ny+1,2:nz+1),ny,nz,new_column)
          call SteadyRiemannBCs(main)
          main(:,2,2:ny+1,2:nz+1) = main(:,1,2:ny+1,2:nz+1)
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
       if(mod(n,10)==0)call write_files_matlab(main(:,2:nx+1,2:ny+1,2),0d0,nx,ny,&
            first_flag=.false.)
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
end module UCS_tester
