module GeneralUtilitiesTester
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
 
end module GeneralUtilitiesTester