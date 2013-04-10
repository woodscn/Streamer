module ODE_Solvers
  use Source_functions
contains

  subroutine Forward_Euler(main,dt,options)
    implicit none
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(:), optional :: options
    real(8) :: dt
    real(8), dimension(:,:,:,:), allocatable :: source_array
    allocate(source_array(size(main,1),size(main,2),size(main,3),size(main,4)))
    call Source(main,source_array,[1])
    main = dt*source_array+main
  end subroutine Forward_Euler

end module ODE_Solvers
