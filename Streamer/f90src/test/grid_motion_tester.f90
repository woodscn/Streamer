!program test
!  use grid_motion_tester
!  i= grid_motion_test()
!end program test
module grid_motion_tester
  use GeneralUtilities
  use grid_motion_driver
  integer, parameter :: nx = 103, ny = 103, nz = 3, nt = 100
  real(8), parameter :: dx = 0.01, dy = 0.01, dz = 0.01, dt = 0.01
  ! Test the means for grid control, including grid-angle-preservation,
  ! Jacobian preservation, and whatever else.
contains
  integer function grid_motion_test()
    implicit none
    real(8), dimension(21,nx,ny,nz) :: main_data
    integer, dimension(1000) :: options
    real(8), dimension(nx,ny) :: fluxA, fluxB, fluxL, fluxM
    integer :: n, i, j
    real(8) :: err
    call grid_motion_tester_init_2D(main_data,options)
    do n = 1, nt
       call h_update_driver(main_data,options) ! Update grid velocity
       do i = 1, nx
          do j = 1, ny
             main_data(3:5,i,j,1) = &
                  flow_velocity(n*dt,(i-1)*dx,(j-1)*dy,(0)*dz)
          end do
       end do
       call TwoDGradient(main_data(15,:,:,1),dx,dy,nx,ny,fluxA,fluxL)
       call TwoDGradient(main_data(16,:,:,1),dx,dy,nx,ny,fluxB,fluxM)
       main_data( 6,:,:,1) = main_data( 6,:,:,1) + dt*fluxA
       main_data( 7,:,:,1) = main_data( 7,:,:,1) + dt*fluxB
       main_data( 9,:,:,1) = main_data( 9,:,:,1) + dt*fluxL
       main_data(10,:,:,1) = main_data(10,:,:,1) + dt*fluxM
       main_data(18:19,:,:,:) = main_data(18:19,:,:,:) &
            + dt*main_data(15:16,:,:,:)
    end do
    err = sqrt(sum(sum(sum((&
         main_data(6,:,:,:)*main_data(9,:,:,:)+&
         main_data(7,:,:,:)*main_data(10,:,:,:)&
         )**2,1),1),1))
    if (err < 1d-14)then
       grid_motion_test = 0
    else
       grid_motion_test = 1
    end if
  end function grid_motion_test
  
  function flow_velocity(t,x,y,z)
    real(8), intent(in) :: t,x,y,z
    real(8), dimension(3) :: flow_velocity
    real(8) :: u,v,w
    real(8), parameter :: PI = 3.1415926535897932
    real(8), parameter :: PIi = 1.d0/PI
    real(8), parameter :: atu = 1., axu = 1., ayu = 1., azu = 1.
    real(8), parameter :: atv = 1., axv = 1., ayv = 1., azv = 1.
    real(8), parameter :: atw = 1., axw = 1., ayw = 1., azw = 1.
    real(8), parameter :: btu = 0., bxu = 0., byu = 0., bzu = 0.
    real(8), parameter :: btv = 0., bxv = 0., byv = 0., bzv = 0.
    real(8), parameter :: btw = 0., bxw = 0., byw = 0., bzw = 0.
    u = cos(atu*t+btu)*cos(axu*x+bxu)*cos(ayu*y+byu)*cos(azu*z+bzu)
    v = cos(atv*t+btv)*cos(axv*x+bxv)*cos(ayv*y+byv)*cos(azv*z+bzv)
    w = cos(atw*t+btw)*cos(axw*x+bxw)*cos(ayw*y+byw)*cos(azw*z+bzw)
    flow_velocity = [u,v,w]
  end function flow_velocity

  subroutine grid_motion_tester_init_2D(main,options)
    implicit none
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(:), intent(inout) :: options
    integer i, j
    
    main = 0.d0
    do i = 1, nx
       do j = 1, ny
          main(3:5,i,j,2) = flow_velocity(0.d0,(i-1)*dx,(j-1)*dy,0.d0)
          main(6:14,i,j,2) = [1.,0.,0.,0.,1.,0.,0.,0.,1.]
          main(15:17,i,j,2) = .999*main(3:5,i,j,1)
          main(18,i,j,:) = (i-1)*dx
          main(19,i,j,:) = (j-1)*dy
          main(20,i,j,:) = 0.*dz
       end do
    end do
    options = 0
    options(3:5) = 1
    options(6:7) = [3,3]

  end subroutine grid_motion_tester_init_2D

  function grid_motion_reader(in)
    implicit none
    integer, intent(in) :: in
    integer grid_motion_reader
    select case(in)
    case(0)
       write(*,*) "Grid_motion_test passed "
    case default
       write(*,*) "Error in grid_motion_test"
    end select
  end function grid_motion_reader
end module grid_motion_tester
