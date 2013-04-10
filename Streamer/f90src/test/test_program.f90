program test_program
  use generalutilitiestester
  use godunov_tester
  use riemann_tester
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
  write(*,*) "Godunov tester:"
  result = GodTester(output_dir)
  junk = GodErrorReader(result)
!!$  write(*,*) "ODE_Solvers tester:"
!!$  result = ODEtester()
!!$  junk = ODEErrorReader(result)
  write(*,*) "Done"
end program test_program
