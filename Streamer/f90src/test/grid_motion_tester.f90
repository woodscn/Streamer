!program test
!  use grid_motion_tester
!  i= grid_motion_test()
!end program test
module grid_motion_tester
  use GeneralUtilities
  use grid_motion_driver_mod
  real(8), parameter :: xmin = 0d0, xmax = 1d0, ymin = 0d0, ymax = 1d0
  ! Test the means for grid control, including grid-angle-preservation,
  ! Jacobian preservation, and whatever else.
contains
  integer function grid_motion_test()
    implicit none
    integer :: l, m, n, i, j, k
    integer, parameter :: len_nx_a = 2
    integer, parameter, dimension(len_nx_a) :: nx_a = [103,23]
    real(8), dimension(len_nx_a) :: err
    real(8), dimension(len_nx_a) :: dx_a, dy_a
    real(8) :: order

    integer, parameter ::  len_nt_a = 2
    integer, parameter, dimension(len_nt_a) :: nt_a = [11,1001]
    real(8), dimension(len_nt_a) :: dt_a = [.1, .1]
    real(8), dimension(len_nt_a) :: orthogonality
    real(8), dimension(21,53,53,3) :: main_data
    integer, dimension(1000) :: options
    real(8) :: t    
    real(8), dimension(51,51) :: FluxA, FluxB, FluxL, FluxM
    grid_motion_test = 0

    ! Test the spatial convergence of the grid_update_driver routine
    do n = 1, len_nx_a
       call compute_rms_error_2d(nx_a(n),nx_a(n),3,1d0,0d0,1d0,0d0,err(n),&
            dx_a(n), dy_a(n))
    end do
    order = (log(err(2))-log(err(1)))/(log(dx_a(2))-log(dx_a(1)))
    if (order < .5d0)then
       grid_motion_test = 1
    end if

    ! Test the efficacy of grid-angle preservation.
    call grid_motion_tester_init_2D(main_data,options,.02d0,.02d0,53,53)
    do n = 1, len_nt_a
       t = 0.d0
       do m = 1, nt_a(n)
          do i = 1, 53
             do j = 1, 53
                do k = 1, 3
                   main_data(3:5,i,j,k) = flow_velocity(&
                        t,.02d0*(i-2),.02d0*(j-2),.02d0*(k-2))
                end do
             end do
          end do
          call grid_motion_driver(main_data,options)
          main_data(6:17,:,:, 1) = main_data(6:17,:,:, 2)
          main_data(6:17,:,:, 3) = main_data(6:17,:,:, 2)
          main_data(6:17, 1,:,:) = main_data(6:17, 2,:,:)
          main_data(6:17,53,:,:) = main_data(6:17,52,:,:)
          main_data(6:17,:, 1,:) = main_data(6:17,:, 2,:)
          main_data(6:17,:,53,:) = main_data(6:17,:,52,:)
          call TwoDGradient(&
               main_data(15,2:52,2:52,2),.02d0,.02d0,51,51,fluxA,fluxB)
          call TwoDGradient(&
               main_data(16,2:52,2:52,2),.02d0,.02d0,51,51,fluxL,fluxM)
          main_data( 6,2:52,2:52,2) = main_data( 6,2:52,2:52,2)+dt_a(n)*fluxA
          main_data( 7,2:52,2:52,2) = main_data( 7,2:52,2:52,2)+dt_a(n)*fluxB
          main_data( 9,2:52,2:52,2) = main_data( 9,2:52,2:52,2)+dt_a(n)*fluxL
          main_data(10,2:52,2:52,2) = main_data(10,2:52,2:52,2)+dt_a(n)*fluxM
       end do
       orthogonality(n) = sqrt(sum((&
            main_data(6,2:52,2:52,2)*main_data( 9,2:52,2:52,2)+&
            main_data(7,2:52,2:52,2)*main_data(10,2:52,2:52,2))**2)/(51*51))
    end do
    if (maxval(abs(orthogonality))>1d-12)then
       grid_motion_test = 2
    end if
  end function grid_motion_test

  subroutine compute_rms_error_2d(nx,ny,nz,xmax,xmin,ymax,ymin,err,dx,dy)
    integer, intent(in) :: nx, ny, nz
    real(8), intent(in) :: xmax, xmin, ymax, ymin
    real(8), intent(out) :: dx, dy, err
    real(8), dimension(21,nx,ny,nz) :: main_data
    real(8), dimension(nx,ny) :: fluxA, fluxB, fluxL, fluxM
    real(8), dimension(nx,ny) :: g, t, dgdx, dgdy, dtdx, dtdy
    real(8), dimension(nx,ny) :: S2, T2, alpha, beta, gamma
    integer, dimension(1000) :: options

    dx = (xmax-xmin)/(nx-3); dy = (ymax-ymin)/(ny-3);
    call grid_motion_tester_init_2D(main_data,options,dx,dy,nx,ny)
    call grid_motion_driver(main_data,options)
    g = log(sqrt(main_data(15,:,:,2)**2+main_data(16,:,:,2)**2))
    t = atan2(main_data(16,:,:,2),main_data(15,:,:,2))
    call TwoDGradient(g,dx,dy,nx,ny,dgdx,dgdy)
    call TwoDGradient(t,dx,dy,nx,ny,dtdx,dtdy)
    ! Big, long formula to compute the residual.
    S2 = main_data(9,:,:,2)**2+main_data(10,:,:,2)**2
    T2 = main_data(6,:,:,2)**2+main_data(7,:,:,2)**2
    alpha = S2&
         *(main_data(6,:,:,2)*sin(t)-main_data(7,:,:,2)*cos(t))
    beta = T2&
         *(main_data(10,:,:,2)*cos(t)-main_data(9,:,:,2)*sin(t))
    gamma = -S2*dtdx&
         *(main_data(6,:,:,2)*cos(t)+main_data(7,:,:,2)*sin(t))&
         +T2*dtdy&
         *(main_data(9,:,:,2)*cos(t)+main_data(10,:,:,2)*sin(t))
    err = sqrt(sum(&
         (alpha*dgdx+beta*dgdy-gamma)**2)/size(alpha))
  end subroutine compute_rms_error_2d

  function flow_velocity(t,x,y,z)
    real(8), intent(in) :: t,x,y,z
    real(8), dimension(3) :: flow_velocity
    real(8) :: u,v,w
    real(8), parameter :: PI = 3.1415926535897932
    real(8), parameter :: PIi = 1.d0/PI
    real(8), parameter :: atu = 1., axu = 1., ayu = 1., azu = 0.
    real(8), parameter :: atv = 1., axv = 1., ayv = 1., azv = 0.
    real(8), parameter :: atw = 0., axw = 0., ayw = 0., azw = 0.
    real(8), parameter :: btu = 0., bxu = 0., byu = 0., bzu = 0.
    real(8), parameter :: btv = 0., bxv = 0., byv = 0., bzv = 0.
    real(8), parameter :: btw = 0., bxw = 0., byw = 0., bzw = 0.
    u = cos(atu*t+btu)*cos(axu*x+bxu)*cos(ayu*y+byu)*cos(azu*z+bzu)
    v = sin(atv*t+btv)*sin(axv*x+bxv)*sin(ayv*y+byv)*sin(azv*z+bzv)
    w = sin(atw*t+btw)*sin(axw*x+bxw)*sin(ayw*y+byw)*sin(azw*z+bzw)
    flow_velocity = [u,v,w]
  end function flow_velocity

  subroutine grid_motion_tester_init_2D(main,options,dx,dy,nx,ny)
    implicit none
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(:), intent(inout) :: options
    real(8), intent(in) :: dx, dy
    integer, intent(in) :: nx, ny
    integer i, j, k
    
    main = 0.d0
    options = 0
    options(3:5) = 1
    options(6:7) = [4,3]
    do i = 1, nx
       do j = 1, ny
          do k = 1, 3
             main(1:2,i,j,k) = 1.d0
             main(3:5,i,j,k) = flow_velocity(0.d0,(i-2)*dx,(j-2)*dy,(k-2)*1d0)
             main(6:14,i,j,k) = [1.,0.,0.,0.,1.,0.,0.,0.,1.]
             main(15:17,i,j,k) = .5*main(3:5,i,j,2)
             main(18,i,j,k) = (i-2)*dx
             main(19,i,j,k) = (j-2)*dy
             main(20,i,j,k) = (k-2)*1
             main(21,i,j,k) = 1.
          end do
       end do
    end do

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
