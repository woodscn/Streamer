module Riemann_tester
  use Riemann
  real(8), dimension(5) :: riemann_sol
  real(8), dimension(21):: test_geom
  real(8) :: mws_test
  ! Though it's possible to do specific tests of the guessp algorithm,
  ! it seems less-effective to do so, since a bad guess will still converge.
contains
  integer function RieErrorReader(in)
    implicit none
    integer, intent(in) :: in
    select case(in)
    case(0)
       write(*,*) "   All tests passed"
    case(1)
       write(*,*) "   Failure in one of the Riemann Star state checks"
    case default
       write(*,*) "   Unexpected error code"
    end select
  end function RieErrorReader

  integer function RieTester(output_dir)
    implicit none
    character(len=*), intent(in) :: output_dir
    real(8) :: U, V, W
    real(8), dimension(21) :: test_base
    real(8), dimension(21) :: test_1_left
    real(8), dimension(21) :: test_1_right
    real(8), dimension(21) :: test_2_left
    real(8), dimension(21) :: test_2_right
    real(8), dimension(21) :: test_3_left
    real(8), dimension(21) :: test_3_right
    real(8), dimension(21) :: test_4_left
    real(8), dimension(21) :: test_4_right
    real(8), dimension(21) :: test_5_left
    real(8), dimension(21) :: test_5_right
    real(8), dimension(4) :: test_1_sol
    real(8), dimension(4) :: test_2_sol
    real(8), dimension(4) :: test_3_sol
    real(8), dimension(4) :: test_4_sol
    real(8), dimension(4) :: test_5_sol
    real(8), dimension(5) :: test_1_speeds
    real(8), dimension(5) :: test_2_speeds
    real(8), dimension(5) :: test_3_speeds
    real(8), dimension(5) :: test_4_speeds
    real(8), dimension(5) :: test_5_speeds
    real(8), dimension(21) :: moving_test_1_left
    real(8), dimension(21) :: moving_test_1_right
    real(8), dimension(21) :: complicated_grid_left
    real(8), dimension(21) :: complicated_grid_right
    real(8), dimension(4) :: moving_test_1_sol
    real(8), dimension(5) :: moving_test_1_speeds
    integer, parameter :: nx = 1
    real(8), dimension(nx) :: x
    real(8), dimension(5,nx) :: out
    real(8) :: max_wave_speed
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    integer, parameter :: nxfull = 101
    real(8), dimension(nxfull) :: xfull
    real(8), dimension(5,nxfull) :: outfull
    integer :: n, dir
    RieTester = 0
    
    test_base = [0d0,0d0,0d0,0d0,0d0,&
         1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0,&
         U,V,W,0d0,0d0,0d0,0d0]
    test_1_left = test_base; test_1_right = test_base
    test_2_left = test_base; test_2_right = test_base
    test_3_left = test_base; test_3_right = test_base
    test_4_left = test_base; test_4_right = test_base
    test_5_left = test_base; test_5_right = test_base
    
    test_1_left(1:5)  = [1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_1_right(1:5) = [.1d0, .125d0, 0.d0, 0.d0, 0.d0]
    test_2_left(1:5)  = [.4d0, 1.d0, -2.d0, 0.d0, 0.d0]
    test_2_right(1:5) = [.4d0, 1.d0, 2.d0, 0.d0, 0.d0]
    test_3_left(1:5)  = [1.d3, 1.d0, 0.d0, 0.d0, 0.d0]
    test_3_right(1:5) = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_4_left(1:5)  = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_4_right(1:5) = [1.d2, 1.d0, 0.d0, 0.d0, 0.d0]
    test_5_left(1:5)  = [460.894d0, 5.99924d0, 19.5975d0, 0d0, 0d0]
    test_5_right(1:5) = [46.095d0, 5.99242d0, -6.19633d0, 0d0, 0d0]
    test_1_sol = [.30313d0, .92745d0, .42632d0, .26557d0]
    test_1_speeds = [-1.18322d0,-0.0702745d0,0.92745d0,0d0,1.75216d0]
    test_2_sol = [.00189d0, .00000d0, .02185d0, .02185d0]
    test_2_speeds = [-2.74833d0,-0.347992d0,0d0,0.347992d0,2.74833d0]
    test_3_sol = [460.894d0, 19.5975d0, .57506d0, 5.99924d0]
    test_3_speeds = [-37.4166d0,-13.8997d0,19.5975d0,0d0,23.5175d0]
    test_4_sol = [46.095d0, -6.19633d0, 5.99242d0, .57511d0]
    test_4_speeds = [-7.43747d0,0d0,-6.19633d0,4.39658d0,11.8322d0]
    test_5_sol = [1691.64d0, 8.68975d0, 14.2823d0, 31.0426d0]
    test_5_speeds = [0.789631d0,0d0,8.68975d0,0d0,12.2507d0]
    ! Test riemann_solve against the test problems from Toro. This only tests
    ! whether Pstar, Ustar, DstarL and DstarR are correct, and whether the 
    ! computed wave speeds match as computed in Mathematica, and visually
    ! compared (eyeballed) with Toro's plots. This leaves the solutions within
    ! the expansion fan untested.
    test_flag = .true.
!    test_sol = test_1_sol
!      subroutine riemann_solve(left, right, dir, nx, x, out, max_wave_speed,&
!       riemann_middle_states, riemann_wave_speeds)
    do dir = 1, 3
       call riemann_solve(test_1_left,test_1_right,dir,1,[0d0],out,&
            max_wave_speed,riemann_middle_states,riemann_wave_speeds)
       if(.not.(maxval(abs(riemann_middle_states-test_1_sol))<5d-6))&
            RieTester = 1
       if(.not.(maxval(&
            abs((riemann_wave_speeds-test_1_speeds)/test_1_speeds))<5d-5))&
            RieTester = 1
       call riemann_solve(test_2_left,test_2_right,dir,1,[0d0],out,&
            max_wave_speed,riemann_middle_states,riemann_wave_speeds)
       if(.not.(maxval(&
            abs(riemann_middle_states-test_2_sol)/test_2_sol)<5d-3))&
            RieTester = 1
       if(.not.(maxval(abs(riemann_wave_speeds-test_2_speeds))<5d-4))&
            RieTester = 1
       call riemann_solve(test_3_left,test_3_right,dir,1,[0d0],out,&
            max_wave_speed,riemann_middle_states,riemann_wave_speeds)
       if(.not.(maxval(&
            abs(riemann_middle_states-test_3_sol)/test_3_sol)<5d-6))&
            RieTester = 1
       if(.not.(maxval(abs(riemann_wave_speeds-test_3_speeds))<7d-5))&
            RieTester = 1
       call riemann_solve(test_4_left,test_4_right,dir,1,[0d0],out,&
            max_wave_speed,riemann_middle_states,riemann_wave_speeds)
       if(.not.(maxval(&
            abs(riemann_middle_states-test_4_sol)/test_4_sol)<5d-6))&
            RieTester = 1
       if(.not.(maxval(abs(riemann_wave_speeds-test_4_speeds))<5d-5))&
            RieTester = 1
       call riemann_solve(test_5_left,test_5_right,dir,1,[0d0],out,&
            max_wave_speed,riemann_middle_states,riemann_wave_speeds)
       if(.not.(maxval(&
            abs(riemann_middle_states-test_5_sol)/test_5_sol)<5d-6))&
            RieTester = 1
       if(.not.(maxval(abs(riemann_wave_speeds-test_5_speeds))<8d-5))&
            RieTester = 1
    end do
    ! Visual tests are possible to check that the solution within the fan is 
    !correct. This writes the density to a file that can then be plotted.
    do n = 1, nxfull
       xfull(n) = 0d0+1d0/(nxfull-1d0)*(n-1)
    end do
    call riemann_solve(test_1_left,test_1_right,1,nxfull,(xfull-.5d0)/.15d0,&
         outfull,max_wave_speed,riemann_middle_states,riemann_wave_speeds)
    open(unit=9492, file=trim(output_dir)//'riemann_test.dat')
    write(9492,*) outfull(2,:)
    write(9492,*) xfull
    close(9492)

    ! I need to test this with more complex geometric variables
    test_geom = [0d0,0d0,0d0,0d0,0d0,&
         1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0,&
         0d0,0d0,0d0,0d0,0d0,0d0,1d0]

    complicated_grid_left = test_1_left
    complicated_grid_left(6:21) = test_geom(6:21)
    complicated_grid_right = test_1_right
    complicated_grid_right(6:21) = test_geom(6:21)
    do n = 1, nxfull
       xfull(n) = 0d0+1d0/(nxfull-1)*(n-1)
    end do
    call riemann_solve(complicated_grid_left,complicated_grid_right,1,nxfull,&
         (xfull-.5d0)/.15d0,outfull,max_wave_speed,riemann_middle_states,&
         riemann_wave_speeds)
    open(unit=9493,file=trim(output_dir)//'complicated_riemann_test.dat')
    write(9493,*)outfull(2,:)
    write(9493,*)xfull
    close(9493)


!!$    moving_test_1_left = test_1_left
!!$    moving_test_1_left(3) = moving_test_1_left(3) + .3d0
!!$    moving_test_1_left(15) = moving_test_1_left(15) + .3d0
!!$    moving_test_1_right = test_1_right
!!$    moving_test_1_right(3) = moving_test_1_right(3) + .3d0
!!$    moving_test_1_right(15) = moving_test_1_right(15) + .3d0
!!$    call riemann_solve(moving_test_1_left,moving_test_1_right,1,1,[0d0],out,max_wave_speed,&
!!$         moving_test_1_sol,moving_test_1_speeds)
!!$    write(*,*) "Stationary sol = ", riemann_middle_states
!!$    write(*,*) "Moving sol     = ", moving_test_1_sol
!!$    write(*,*) "Stationary speeds = ", riemann_wave_speeds
!!$    write(*,*) "Moving speeds     = ", moving_test_1_speeds
!!$    read(*,*)
!!$    
!!$
  end function RieTester
end module Riemann_tester

