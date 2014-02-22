program test_program
  use generalutilitiestester
  use godunov_tester
  use riemann_tester
  use grid_motion_tester
  use UCS_tester
!  use ode_solvers_tester
  implicit none
  integer :: result, junk
  character(len=32) :: output_dir
  call get_command_argument(1,value=output_dir)
  write(*,*) "General Utilities test:"
  result = GUTest()
  junk = GUErrorReader(result)
  write(*,*) "Riemann tester:"
  result = RieTester(output_dir)
  junk = RieErrorReader(result)
!!$  write(*,*) "Grid Motion Tester:"
!!$  result = grid_motion_test()
!!$  junk = grid_motion_reader(result)
!!$  write(*,*) "Godunov tester:"
!!$  result = GodTester(output_dir)
!!$  junk = GodErrorReader(result)
  write(*,*) "UCS_tester:"
  result = UCS_test_main()

!!$  write(*,*) "ODE_Solvers tester:"
!!$  result = ODEtester()
!!$  junk = ODEErrorReader(result)
  write(*,*) "Done"
end program test_program
