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

    MetrictoMatrix = reshape([Metric(1),Metric(4),Metric(7),&
         Metric(2),Metric(5),Metric(8),Metric(3),Metric(6),&
         Metric(9)],[3,3])
  end function MetrictoMatrix
  
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

  function GradstoMatrix(Grad1,Grad2,Grad3)
    ! 
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
    if(.true. .and. maxval(matmul(in,MatrixInverse)&
            -reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2>1.d-15)then
          write(*,*) "MatrixInverse failed!!"
          write(*,*) matmul(in,MatrixInverse)
          stop
    end if
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

module GeneralUtilitiesTest
  use GeneralUtilities
  real(8), dimension(:), allocatable :: lstsqx, lstsqy
contains
  integer function GUErrorReader(in)
    implicit none
    integer, intent(in) :: in
    write(*,*) " GeneralUtilities.f90 diagnostic:"
    select case(in)
    case(0)
       write(*,*) "   All tests passed"
    case(1) 
       write(*,*) "   Combined test of metrictomatrix and metricinverse failed"
    case(2)
       write(*,*) "   TwoDGradient does not converge at the expected rate"
    case default
       write(*,*) "   Unexpected error code"
    end select
    GUErrorReader = 0
  end function GUErrorReader

  integer function GUTest()
    ! Test the following routines:
    !   function MetricInverse(Metric)
    !   function MetrictoMatrix(Metric)
    !   subroutine ComputationalGrads(metric,jac,grad_xi,grad_eta,grad_zeta)
    !   function Jacobian(in)
    !   function GradstoMatrix(Grad1,Grad2,Grad3)
    !   subroutine TwoDGradient(in,dx,dy,nx,ny,gradx,grady)

    ! Returns an integer error code.
    ! 0: All tests passed
    ! 1: Combined test of metrictomatrix and metricinverse failed
    ! 2: TwoDGradient does not converge at the expected rate
    implicit none
    real(8), dimension(9) :: metric
    real(8), dimension(3,3) :: matrix
    real(8), dimension(0:3) :: p1, p2, p3
    integer :: out, i, j, k, n, nx, ny
    real(8), dimension(400,400) :: test, grad_test_x, grad_test_y
    real(8), dimension( 50, 50) :: num_test_1x, num_test_1y
    real(8), dimension(100,100) :: num_test_2x, num_test_2y
    real(8), dimension(200,200) :: num_test_3x, num_test_3y
    real(8), dimension(400,400) :: num_test_4x, num_test_4y
    real(8), dimension(4) :: dxes, rmserrors
    real(8), dimension(2) :: fitted_poly
    real(8) :: dx, dy
    real(8) :: x, y
    
    out = 0
    ! Try variable ranges on these random numbers to check to handle 
    ! ill-conditioned matrices
    call random_number(metric)
    call random_number(p1)
    call random_number(p2)
    call random_number(p3)
    matrix=matmul(metrictomatrix(metric),metrictomatrix(metricinverse(metric)))
    if(maxval((matrix-reshape([1.,0.,0.,0.,1.,0.,0.,0.,1.],[3,3]))**2)>1.d-14) out = 1
    nx = size(test,1)
    ny = size(test,2)
    dx = 12./nx
    dy = 13./ny
    do j = 1, ny
       do i = 1, nx
          x = i*dx
          y = j*dy
          test(i,j) = 0.d0
          grad_test_x(i,j) = 0.d0
          grad_test_y(i,j) = 0.d0
          do n = 0, size(p1) - 1
             test(i,j) = test(i,j) + p1(n)*x**n + p2(n)*y**n
             grad_test_x(i,j) = grad_test_x(i,j) + n*p1(n)*x**(n-1)
             grad_test_y(i,j) = grad_test_y(i,j) + n*p2(n)*y**(n-1)
          end do
       end do
    end do
    call TwoDGradient(test(1:nx:8,1:ny:8),dx*8,dy*8,nx/8,ny/8,&
         num_test_1x,num_test_1y)
    call TwoDGradient(test(1:nx:4,1:ny:4),dx*4,dy*4,nx/4,ny/4,&
         num_test_2x,num_test_2y)
    call TwoDGradient(test(1:nx:2,1:ny:2),dx*2,dy*2,nx/2,ny/2,&
         num_test_3x,num_test_3y)
    call TwoDGradient(test,dx,dy,nx,ny,num_test_4x,num_test_4y)

    dxes = [dx*8,dx*4,dx*2,dx]
    rmserrors = [sqrt(sum((num_test_1x-grad_test_x(1:nx:8,1:ny:8))**2)&
         /size(num_test_1x)),&
         sqrt(sum((num_test_2x-grad_test_x(1:nx:4,1:ny:4))**2)&
         /size(num_test_2x)),&
         sqrt(sum((num_test_3x-grad_test_x(1:nx:2,1:ny:2))**2)&
         /size(num_test_3x)),&
         sqrt(sum((num_test_4x-grad_test_x)**2)/size(num_test_4x))]

    fitted_poly = polyfit(log(dxes),log(rmserrors),1)
    if(fitted_poly(2) < 1.) out = 2
    GUTest = out
    
  end function GUTest

  integer function lstsq_init(x,y)
    implicit none
    real(8), dimension(:), intent(in) :: x,y
    integer :: nx
    nx = size(x,1)
    allocate(lstsqx(nx), lstsqy(nx))
    lstsqx = x ; lstsqy = y
    lstsq_init = 0
  end function lstsq_init
  
  integer function lstsq_close()
    leastsq_close = 0
    deallocate(lstsqx,lstsqy)
  end function lstsq_close

  subroutine minpack_function_fitting(xdat,ydat,fcn,x,fvec,fjac,tol,info)
    implicit none
    real(8), dimension(:), intent(in) :: xdat,ydat
    external fcn
    real(8), dimension(:), intent(inout) :: x
!    integer, intent(in) :: m, n, ldfjac
    real(8), dimension(:), intent(out) :: fvec
    real(8), dimension(:,:), intent(out) :: fjac
    real(8), intent(in) :: tol
    integer, intent(out) :: info
    integer, dimension(:), allocatable :: wa, ipvt
    integer :: out, m, n, ldfjac, lwa,iflag
    
iflag = 2
    m = size(xdat,1)
    n = size(x,1)
    ldfjac = m
    allocate(wa(2*(5*n+m)),ipvt(n))
    lwa = size(wa)
    out = lstsq_init(xdat,ydat)
    call exponential_with_y_offset(m,n,x,fvec,fjac,ldfjac,iflag)
    call lmder1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,lwa)
    out = lstsq_close()
    deallocate(wa,ipvt)
  end subroutine minpack_function_fitting

  subroutine exponential_with_y_offset(m,n,x,fvec,fjac,ldfjac,iflag)
    integer, intent(in) :: m
    integer, intent(in) :: n
    integer, intent(in) :: ldfjac
    real(8), dimension(n), intent(in) :: x
    real(8), dimension(m), intent(out) :: fvec
    real(8), dimension(ldfjac,n), intent(out) :: fjac
    integer, intent(inout) :: iflag
    integer i
    ! If iflag = 1 calculate the functions at x and return
    !     this vector in fvec. Do not alter fjac.
    ! If iflag = 2 calculate the jacobian at x and return
    !     this matrix in fjac. Do not alter fvec.
    ! The value of iflag should not be changed unless the
    !     user wants to terminate execution of the least-
    !     squares fit. In this case, set iflag < 0.
    if(.not.m==size(lstsqx,1))then
       write(*,*) "Error in least-squares fit!!"
       write(*,*) "  Incompatible dimensions passed to fcn!"
       stop
    end if
    select case(iflag)
    case(1)
       forall (i = 1: m)
          fvec(i) = x(1)*lstsqx(i)**x(2)+x(3)-lstsqy(i)
       end forall
    case(2)
       forall (i = 1: m)
          fjac(i,:) = [lstsqx(i)**x(2),x(1)*lstsqx(i)**x(2)*log(lstsqx(i)),1d0]
       end forall
    end select
  end subroutine exponential_with_y_offset
!!$interface ConvergenceFit
!!$   subroutine 4DConvergenceFit(numerical,exact,nmax,out)
!!$     implicit none
!!$     type(4d_array_pointer), dimension(:), intent(in) :: numerical
!!$     real(8), dimension(:,:,:,:), intent(in) :: exact
!!$     real(8), dimension(:,:,:,:), intent(out) :: out
!!$   end subroutine 4DConvergenceFit
!!$   subroutine 2DConvergenceFit(numerical,exact,nmax,out)
!!$     type(4d_array_pointer), dimension(:), intent(in) :: numerical
!!$     real(8), dimension(:,:,:,:), intent(in) :: exact
!!$     real(8), dimension(:,:), intent(out) :: out
!!$   end subroutine 2DConvergenceFit
!!$   subroutine 1DConvergenceFit(numerical,exact,nmax,out)
!!$     type(4d_array_pointer), dimension(:), intent(in) :: numerical
!!$     real(8), dimension(:,:,:,:), intent(in) :: exact
!!$     real(8), dimension(:), intent(out) :: out
!!$end interface ConvergenceFit
!!$
!!$subroutine 4DConvergenceFit(numerical,exact,nmax,out)
!!$  implicit none
!!$  type(4d_array_pointer), dimension(:), intent(in) :: numerical
!!$  real(8), dimension(:,:,:,:), intent(in) :: numerical, exact
!!$  integer, intent(in) :: nmax
!!$  real(8), dimension(:,:,:,:), intent(out) :: out
!!$  integer :: n
!!$  real(8), dimension(nmax) :: dxes
!!$
!!$
!!$  do n = 1, nmax
!!$     
!!$  
!!$  
!!$end subroutine 4DConvergenceFit

  ! From Rosetta Code: rosettacode.org/wiki/PolynomialRegression#Fortran
  ! Requires LAPACK library
  ! Verified for single-case linear fit against Mathematica, 27 Aug 2012.
  function polyfit(vx, vy, d)
    implicit none
    integer, intent(in)                   :: d
    integer, parameter                    :: dp = selected_real_kind(15, 307)
    real(dp), dimension(d+1)              :: polyfit
    real(dp), dimension(:), intent(in)    :: vx, vy
 
    real(dp), dimension(:,:), allocatable :: X
    real(dp), dimension(:,:), allocatable :: XT
    real(dp), dimension(:,:), allocatable :: XTX
 
    integer :: i, j
 
    integer     :: n, lda, lwork
    integer :: info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
 
    n = d+1
    lda = n
    lwork = n
 
    allocate(ipiv(n))
    allocate(work(lwork))
    allocate(XT(n, size(vx)))
    allocate(X(size(vx), n))
    allocate(XTX(n, n))
 
    ! prepare the matrix
    do i = 0, d
       do j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       end do
    end do
 
    XT  = transpose(X)
    XTX = matmul(XT, X)
 
    ! calls to LAPACK subs DGETRF and DGETRI
    call DGETRF(n, n, XTX, lda, ipiv, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
    call DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    if ( info /= 0 ) then
       print *, "problem"
       return
    end if
 
    polyfit = matmul( matmul(XTX, XT), vy)
 
    deallocate(ipiv)
    deallocate(work)
    deallocate(X)
    deallocate(XT)
    deallocate(XTX)
 
  end function
 
end module GeneralUtilitiesTest

!!$program test_regime
!!$  use generalutilitiestest
!!$  implicit none
!!$  integer :: result, junk
!!$  result = GUTest()
!!$  junk = GUErrorReader(result)
!!$!  write(*,*) GUErrorReader(result)
!!$end program test_regime
