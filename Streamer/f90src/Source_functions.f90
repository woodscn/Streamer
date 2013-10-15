module Source_functions
contains  
  function Gravity(in,t_in,g,n)
    use GeneralUtilities
    implicit none
    !f2py depend(n) in, Gravity
    integer, intent(in) :: n
    real(8), dimension(n), intent(in) :: in
    real(8), dimension(n) :: Gravity
    real(8), intent(in) :: t_in
    real(8), dimension(3), intent(in) :: g

    Gravity = 0d0
    Gravity(3:5) = matmul(MetrictoMatrix(MetricInverse(in(6:14))),g*in(2))
  end function Gravity
!!$
!!$  subroutine Laplacian(in,out,dx_in)
!!$    implicit none
!!$    real(8), dimension(:,:,:,:), intent(in) :: in
!!$    real(8), dimension(:,:,:,:), intent(out) :: out
!!$    real(8), dimension(3) :: dx_in
!!$    integer :: i, j, k
!!$
!!$    out = 0d0
!!$    forall(i=2:size(in,2)-1,j=2:size(in,3)-1,k=2:size(in,4)-1)
!!$       out(:,i-1,j-1,k-1) = &
!!$            (in(:,i-1,j,k)-2d0*in(:,i,j,k)+in(:,i+1,j,k))/dx_in(1)**2 + &
!!$            (in(:,i,j-1,k)-2d0*in(:,i,j,k)+in(:,i,j+1,k))/dx_in(2)**2 + &
!!$            (in(:,i,j,k-1)-2d0*in(:,i,j,k)+in(:,i,j,k+1))/dx_in(3)**2
!!$    end forall
!!$  end subroutine Laplacian
end module Source_functions