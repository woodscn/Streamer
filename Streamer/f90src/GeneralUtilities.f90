module GeneralUtilities
  real(8), parameter :: PI = 3.141592653589793
  real(8), parameter :: EPS = 5.d-15
  real(8), parameter :: EPSs = 1d-4
  real(8), parameter :: gamma_const = 1.4d0
  real(8), parameter :: gamma1 = 1.d0/(gamma_const-1.d0)
  real(8), parameter :: gamma2 = (gamma_const-1.d0)
  real(8), parameter :: gamma3 = (gamma_const - 1.d0)/(2.d0*gamma_const)
  real(8), parameter :: gamma4 = 1.d0/gamma3
  real(8), parameter :: gamma5 = (gamma_const-1.d0)/(gamma_const+1.d0)
  real(8), parameter :: gamma6 = 1.d0/(gamma_const+1.d0)
  real(8), parameter :: gamma7 = 1.d0/gamma3
  real(8), dimension(7), parameter :: dxi_a = [1., .5, .25, .2, 2., 4., 5.]
  real(8), dimension(7), parameter :: deta_a = [1., .5, .25, .2, 2., 4., 5.]
  real(8), dimension(7), parameter :: dzeta_a = [1., .5, .25, .2, 2., 4., 5.]
contains
  pure function MetricInverse(Metric)
    implicit none
    ! Given a list of metric variables:
    ! dx/dxi, dy/dxi, dz/dxi, dx/deta, dy/deta, dz/deta, 
    ! dx/dzeta, dy/dzeta, dz/dzeta
    ! ( A, B, C, L, M, N, P, Q, R )
    ! Return the inverse of this list:
    ! dxi/dx, deta/dx, dzeta/dx, dxi/dy, deta/dy, dzeta/dy,
    ! dxi/dz, deta/dz, dzeta/dz.
    real(8), dimension(9), intent(in) :: Metric
    real(8), dimension(9) :: MetricInverse
    real(8) :: J
    
    J = Jacobian(Metric)
    MetricInverse = [&
         Metric(5)*Metric(9) - Metric(6)*Metric(8), &
         Metric(3)*Metric(8) - Metric(2)*Metric(9), &
         Metric(2)*Metric(6) - Metric(3)*Metric(5), &
         Metric(6)*Metric(7) - Metric(4)*Metric(9), &
         Metric(1)*Metric(9) - Metric(3)*Metric(7), &
         Metric(3)*Metric(4) - Metric(1)*Metric(6), &
         Metric(4)*Metric(8) - Metric(5)*Metric(7), &
         Metric(2)*Metric(7) - Metric(1)*Metric(8), &
         Metric(1)*Metric(5) - Metric(2)*Metric(4) ]/J
  end function MetricInverse

  function MetrictoMatrix(Metric)
    implicit none
    ! Given a 9-element, rank-1 array, convert it to a 3x3 matrix,
    ! such that the elements are ordered column-wise:
    ! [ A, B, C, L, M, N, P, Q, R ]
    !               
    !              | |
    !             \   /
    !              \ /
    !               
    !            A, L, P
    !            B, M, Q
    !            C, N, R
    !
    ! This is consistent with the statement:
    ! dx = matmul(MetrictoMatrix(Metric),dxi)
    real(8), dimension(9), intent(in) :: Metric
    real(8), dimension(3,3) :: MetrictoMatrix

    MetrictoMatrix = reshape(Metric,[3,3])
  end function MetrictoMatrix
  
!!$  function NormalizedMetric(in)
!!$    implicit none
!!$    real(8), dimension(9), intent(in) :: in
!!$    real(8), dimension(9) :: NormalizedMetric
!!$    integer :: n
!!$    do n = 1, 3
!!$       NormalizedMetric(n:9:3) = in(n:9:3)/sqrt(sum(in(n:9:3)**2))
!!$    end do
!!$  end function NormalizedMetric

  subroutine ComputationalGrads(metric,jac,grad_xi,grad_eta,grad_zeta)
    implicit none
    ! Metric has the form:
    !   1  2  3  4  5  6  7  8  9  
    ! [ A, B, C, L, M, N, P, Q, R ]
    ! dx/dxi, dy/dxi, dz/dxi, dx/deta, dy/deta, 
    ! dz/deta, dx/dzeta, dy/dzeta, dz/dzeta
    ! 
    ! ComputationalGrads returns:
    ! grad_xi   = [   dxi/dx,   dxi/dy,   dxi/dz ]
    ! grad_eta  = [  deta/dx,  deta/dy,  deta/dz ]
    ! grad_zeta = [ dzeta/dx, dzeta/dy, dzeta/dz ]
    real(8), intent(in), dimension(9) :: metric
    real(8), intent(in) :: jac
    real(8), intent(out), dimension(3) :: grad_xi
    real(8), intent(out), dimension(3) :: grad_eta
    real(8), intent(out), dimension(3) :: grad_zeta
    real(8), dimension(9) :: inv_metric
    
    inv_metric = MetricInverse(metric)
    grad_xi   = inv_metric(1:9:3)
    grad_eta  = inv_metric(2:9:3)
    grad_zeta = inv_metric(3:9:3)
  end subroutine ComputationalGrads

  real(8) pure function Jacobian(in)
    implicit none
    ! Computes the determinant of a 3 x 3 matrix using a brute-force method:
    !       | A L P |
    !  J =  | B M Q |
    !       | C N R |
    ! Assumes the structure of in(:) is :
    !   1  2  3  4  5  6  7  8  9 
    ! [ A, B, C, L, M, N, P, Q, R ]
    real(8), dimension(9), intent(in) :: in
    Jacobian = &
         in(1)*in(5)*in(9) - in(1)*in(6)*in(8) + & ! A*M*R - A*N*Q
         in(2)*in(6)*in(7) - in(2)*in(4)*in(9) + & ! B*N*P - B*L*R
         in(3)*in(4)*in(8) - in(3)*in(5)*in(7)     ! C*L*Q - C*M*P
  end function Jacobian

!!$  function GradstoMatrix(Grad1,Grad2,Grad3)
!!$    ! 
!!$    ! matmul(GradstoMatrix,vector) ==
!!$    ! [ sum(GradstoMatrix(1,:)*vector),
!!$    !   sum(GradstoMatrix(2,:)*vector),
!!$    !   sum(GradstoMatrix(3,:)*vector) ]
!!$    implicit none
!!$    real(8), dimension(3), intent(in) :: Grad1, Grad2, Grad3
!!$    real(8), dimension(3,3) :: GradstoMatrix
!!$    GradstoMatrix = transpose(reshape([Grad1,Grad2,Grad3],[3,3]))
!!$  end function GradstoMatrix
  
  subroutine TwoDGradient(in,dx,dy,nx,ny,gradx,grady)
    implicit none
    ! Compute the 2-dimensional gradient of a matrix using central 
    ! differencing wherever possible. Returns gradient matrices 
    ! the same size and shape as the input matrix. If the input matrix 
    ! has length 1 in either dimension, then it is assumed that the
    ! gradient in that dimension is 0.
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

  function MatrixInverse(in)
    ! Compute the inverse of a 3x3 matrix
    implicit none
    real(8), dimension(3,3), intent(in) :: in
    real(8), dimension(3,3) :: MatrixInverse
    real(8) :: J
    
    J = ( &
         in(1,1)*in(2,2)*in(3,3) + &
         in(1,2)*in(2,3)*in(3,1) + &
         in(1,3)*in(2,1)*in(3,2) - &
         in(1,1)*in(2,3)*in(3,2) - &
         in(1,2)*in(2,1)*in(3,3) - &
         in(1,3)*in(2,2)*in(3,1) )

    MatrixInverse = transpose(reshape( [ & 
         in(3,3)*in(2,2)-in(3,2)*in(2,3) , &
         in(3,2)*in(1,3)-in(3,3)*in(1,2) , &
         in(2,3)*in(1,2)-in(2,2)*in(1,3) , &
         in(3,1)*in(2,3)-in(3,3)*in(2,1) , &
         in(3,3)*in(1,1)-in(3,1)*in(1,3) , &
         in(2,1)*in(1,3)-in(2,3)*in(1,1) , &
         in(3,2)*in(2,1)-in(3,1)*in(2,2) , &
         in(3,1)*in(1,2)-in(3,2)*in(1,1) , &
         in(2,2)*in(1,1)-in(2,1)*in(1,2) ] &
         ,[3,3])/J)
!!$    if(.true. .and. maxval(matmul(in,MatrixInverse)&
!!$            -reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2>1.d-15)then
!!$          write(*,*) "MatrixInverse failed!!"
!!$          write(*,*) matmul(in,MatrixInverse)
!!$          stop
!!$    end if
  end function MatrixInverse

  function vectorProjection(in,normal)
    implicit none
    real(8), intent(in), dimension(3) :: in
    real(8), intent(in), dimension(3) :: normal
    real(8), dimension(3) :: vectorProjection
    vectorProjection = normal*&
         (dot_product(in,normal)/dot_product(normal,normal))
  end function vectorProjection

  real(8) function SoundSpeed(point)
    implicit none
    real(8), intent(in), dimension(21) :: point
    SoundSpeed = sqrt(1.4d0*point(1)/point(2))
  end function SoundSpeed

!!$  logical function pnpoly(npol,xp,yp,x,y)
!!$    ! Check to see if a point lies on the interior of a polygon.
!!$    ! The polygon is given by  x,y points in xp & yp. The point
!!$    ! to test is given by x, y.
!!$    ! Returns true if within the polygon, false otherwise.
!!$
!!$    ! Uses the method of counting the number of times a ray from 
!!$    ! the point crosses the polygon. Original C code given on
!!$    ! http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
!!$    ! and written by Randolph Franklin.
!!$
!!$    ! int pnpoly(int npol, float *xp, float *yp, float x, float y)
!!$    !     {
!!$    !       int i, j, c = 0;
!!$    !       for (i = 0, j = npol-1; i < npol; j = i++) {
!!$    !         if ((((yp[i] <= y) && (y < yp[j])) ||
!!$    !              ((yp[j] <= y) && (y < yp[i]))) &&
!!$    !             (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
!!$    !           c = !c;
!!$    !       }
!!$    !       return c;
!!$    !     }
!!$    implicit none
!!$    integer :: npol
!!$    real(8) :: x, y
!!$    real(8), dimension(npol) :: xp, yp
!!$    integer :: i, j, n
!!$
!!$    pnpoly = .false.
!!$    n = size(xp)
!!$    do i = 1, size(xp)
!!$       if(i==1)then
!!$          j = n
!!$       else
!!$          j = i-1
!!$       end if
!!$       if((((yp(i)<=y).and.(y<yp(j))).or.&
!!$            ((yp(j)<=y).and.(y<yp(i)))).and.&
!!$            (x<(xp(j)-xp(i))*(y-yp(i))/(yp(j)-yp(i))+xp(i)))&
!!$            pnpoly = .not. pnpoly
!!$       write(*,*) "i = ",i,"j = ",j       
!!$    end do
!!$
!!$  end function pnpoly 

end module GeneralUtilities


!!$program test_regime
!!$  use generalutilitiestest
!!$  implicit none
!!$  integer :: result, junk
!!$  result = GUTest()
!!$  junk = GUErrorReader(result)
!!$!  write(*,*) GUErrorReader(result)
!!$end program test_regime
