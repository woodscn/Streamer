module GeneralUtilities
contains
  subroutine ComputationalGrads(Metric,Jac,gradXi,gradEta,gradZeta)
    implicit none
    ! Metric has the form:
    !   1  2  3  4  5  6  7  8  9  
    ! [ A, B, C, L, M, N, P, Q, R ]
    ! Gradients are taken with respect to global Cartesian system
    real(8), intent(in), dimension(9) :: Metric
    real(8), intent(in) :: Jac
    real(8), intent(out), dimension(3) :: gradXi
    real(8), intent(out), dimension(3) :: gradEta
    real(8), intent(out), dimension(3) :: gradZeta
    real(8) :: J

    J = Jac
    
    gradXi   = [ Metric(5)*Metric(9) - Metric(6)*Metric(8),&
                 Metric(6)*Metric(7) - Metric(4)*Metric(9),&
                 Metric(4)*Metric(8) - Metric(5)*Metric(7) ]/J

    gradEta  = [ Metric(3)*Metric(8) - Metric(2)*Metric(9),&
                 Metric(1)*Metric(9) - Metric(3)*Metric(7),&
                 Metric(2)*Metric(7) - Metric(1)*Metric(8) ]/J

    gradZeta = [ Metric(2)*Metric(6) - Metric(3)*Metric(5),&
                 Metric(3)*Metric(4) - Metric(1)*Metric(6),&
                 Metric(1)*Metric(5) - Metric(2)*Metric(4) ]/J
    if(sum((matmul(transpose(reshape(metric,[3,3])),&
         reshape([gradXi,gradEta,gradZeta],[3,3]))&
         - reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2) > 1.d-10)then
       write(*,*) "Failed matrix inverse test!"
       write(*,*) (matmul(transpose(reshape(metric,[3,3])),&
            reshape([gradXi,gradEta,gradZeta],[3,3]))&
            - reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2
       read(*,*)
    end if
    
  end subroutine ComputationalGrads

  ! Computes the determinant of a 3 x 3 matrix using a brute-force method:
  !       | A L P |
  !  J =  | B M Q |
  !       | C N R |
  ! Assumes the structure of in(:) is :
  !   1  2  3  4  5  6  7  8  9 
  ! [ A, B, C, L, M, N, P, Q, R ]
  real(8) function Jacobian(in)
    implicit none
    real(8), dimension(9), intent(in) :: in
    Jacobian = &
         in(1)*in(5)*in(9) - in(1)*in(6)*in(8) + & ! A*M*R - A*N*Q
         in(2)*in(6)*in(7) - in(2)*in(4)*in(9) + & ! B*N*P - B*L*R
         in(3)*in(4)*in(8) - in(3)*in(5)*in(7)     ! C*L*Q - C*M*P
  end function Jacobian

  function GradstoMatrix(Grad1,Grad2,Grad3)
    ! matmul(GradstoMatrix,vector) ==
    ! [ sum(GradstoMatrix(1,:)*vector),
    !   sum(GradstoMatrix(2,:)*vector),
    !   sum(GradstoMatrix(3,:)*vector) ]
    implicit none
    real(8), dimension(3), intent(in) :: Grad1, Grad2, Grad3
    real(8), dimension(3,3) :: GradstoMatrix
    GradstoMatrix = transpose(reshape([Grad1,Grad2,Grad3],[3,3]))
  end function GradstoMatrix
  
  subroutine TwoDGradient(in,dx,dy,nx,ny,gradx,grady)
    implicit none
    real(8), intent(in), dimension(nx,ny) :: in
    real(8), intent(in) :: dx, dy
    integer, intent(in) :: nx, ny
    real(8), intent(out), dimension(nx,ny) :: gradx, grady
    integer :: i, j

    if(nx>1)then
       ! Central differencing where possible
       gradx(2:nx-1,:) = .5d0*(in(3:nx,:)-in(1:nx-2,:))/dx
       ! Forward & backward differencing elsewhere
       gradx(1,:)  = (in( 2,:)-in(   1,:))/dx
       gradx(nx,:) = (in(nx,:)-in(nx-1,:))/dx
    else
       gradx(:,:) = 0.d0
    end if
    if(ny>1)then
       ! Central differencing where possible
       grady(:,2:ny-1) = .5d0*(in(:,3:ny)-in(:,1:ny-2))/dy
       ! Forward & backward differencing elsewhere
       grady(:,1)  = (in(:, 2)-in(:,   1))/dy
       grady(:,ny) = (in(:,ny)-in(:,ny-1))/dy
    else
       grady(:,:) = 0.d0
    end if
  end subroutine TwoDGradient


end module GeneralUtilities
