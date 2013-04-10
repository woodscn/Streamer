module Godunov_driver
  use Godunov
contains
  subroutine prim_update(main,dt_out,bc_func,bcextent,dt_in,CFL,nx,ny,nz,options)
    implicit none
    real(8), dimension(:,:,:,:) :: main
    !f2py intent(in,out) :: main
    integer :: bcextent, nx, ny, nz
    ! options(3) determines which prim_update is called
    integer, dimension(:), intent(in) :: options 
    real(8) :: dt_in, CFL, dt_out
    external :: bc_func
    !f2py intent (callback) bc_func
    !f2py external bc_func
    !f2py dimension(:,:,:,:) real(8) x
    !f2py 
    select case(options(3))
    case(1)
       call prim_update_FV(main,dt_out,bc_func,bcextent,dt_in,CFL,nx,ny,nz,options)
    case(2)
       call prim_update_HUI3D(main,dt_out,bc_func,bcextent,dt_in,CFL,nx,ny,nz,options)
    case default
       write(*,*) "Bad update_type value!!!"
       stop
    end select
  end subroutine prim_update
end module Godunov_driver
