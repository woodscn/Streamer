module grid_motion_driver_mod
  use grid_motion_mod
contains
  subroutine grid_motion_driver(main,options)
    implicit none
    real(8), dimension(:,:,:,:) :: main
    integer, dimension(:) :: options
    !f2py intent(in,out) :: main
    call grid_motion(main,options)
  end subroutine grid_motion_driver
end module grid_motion_driver_mod
