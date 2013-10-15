! -*- f90 -*-

module FortranNormalVectors
contains
subroutine ApplyReflectionConditions(main_data,patch_id,out,nx,ny,dim,t)
use BoundaryConditionsStuff, only: WallReflect
implicit none
integer, intent(in) :: nx, ny, dim
real(8), intent(out), dimension(21,nx,ny) :: out
real(8), intent(in), dimension(21,nx,ny) :: main_data
integer(8), intent(in) :: patch_id
integer :: i, j
real(8) :: x, y, z, t
intent(in) :: t
real(8), dimension(3) :: normal
do i = 1, size(main_data,2)
do j = 1, size(main_data,3)
x = main_data(18,i,j)
y = main_data(19,i,j)
z = main_data(20,i,j)
select case (patch_id)
case(4323088016_8)
normal=[real(1,8),real(0,8),real(0,8)]
case(4358689552_8)
normal=[real(1,8),real(0,8),real(0,8)]
case(4358689744_8)
normal=[real(-1,8),real(0,8),real(0,8)]
case(4358690320_8)
normal=[real(0,8),real(1,8),real(0,8)]
case(4358689808_8)
normal=[real(0.267949192431000d0,8),real(-1,8),real(0,8)]
case(4358691792_8)
normal=[real(0,8),real(1,8),real(0,8)]
case(4358690960_8)
normal=[real(0,8),real(-1,8),real(0,8)]
end select
out(:,i,j) = WallReflect(main_data(:,i,j), normal,&
reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],[3,3]),&
int(dim,4))
end do
end do
end subroutine ApplyReflectionConditions
end module FortranNormalVectors