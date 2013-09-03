module Godunovdriver
  use Godunov
  interface
     subroutine bc_func(main,nx,ny,nz)
       implicit none
       real(8), dimension(21,nx,ny,nz),intent(inout) :: main
       integer, intent(in) :: nx, ny, nz
     end subroutine bc_func
  end interface
contains
  subroutine prim_update(main,dt_out,dt_in,CFL,nx,ny,nz,options)
    implicit none
    real(8), dimension(:,:,:,:) :: main
    !f2py intent(in,out) :: main
    integer, intent(in) :: nx, ny, nz
    integer, dimension(:), intent(in) :: options 
    real(8) :: dt_in, CFL, dt_out
    !f2py intent(out) :: dt_out
    ! Options is an array of integers that controls the behavior of prim_update
    ! and its various subroutines. It is organized as follows.
    ! Elements 1-100 (Fortran numbering) control prim_update itself
    ! Elements 101-200 control prim_update_FV
    ! Elements 201-300 control prim_update_HUI3D
    ! Elements greater than 300 are reserved for future expansion
    select case(options(1))
    case(1)
       call prim_update_FV(main,dt_out,dt_in,CFL,nx,ny,nz,options)
    case(2)
       write(*,*) "Error in prim_update: "
       write(*,*) "-Dimensionally split algorithms not yet implemented!"
       stop
!       call prim_update_HUI3D(main,dt_out,dt_in,CFL,nx,ny,nz,options)
    case default
       write(*,*) "Bad update_type value!!!"
       stop
    end select
  end subroutine prim_update
end module Godunovdriver
