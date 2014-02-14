module h_update_mod
! This module is responsible for solving the PDE for h. It uses an iterative, 
!   first-order, backward difference scheme. Boundary values must be provided 
!   at h(1,:) and h(:,1).
  use GeneralUtilities
contains

  subroutine g_update(A, B, L, M, u_grid, v_grid, u, v, dxi, deta, h0)
    implicit none
    real(8), intent(in) :: A(:,:) , B(:,:) ,  L(:,:) ,   M(:,:)
    real(8), intent(in) :: u(:,:) , v(:,:) , dxi , deta , h0
    real(8), intent(out):: u_grid(:,:), v_grid(:,:)
    real(8)             :: tol  , err !, dxi , deta
    real(8), allocatable:: dthetadxi(:,:) , dthetadeta(:,:) , temp(:,:)
    real(8), allocatable:: alpha(:,:) , beta(:,:) , gamma(:,:) , g(:,:)
    real(8), allocatable:: S2(:,:) , T2(:,:) , q(:,:) , theta(:,:) 
    integer             :: ny , nx , i , j , n , k, iter

    nx = size(A,1) ; ny = size(A,2) ! Define these size parameters
    allocate( g(nx,ny) , q(nx,ny) , theta(nx,ny) , S2(nx,ny) , T2(nx,ny) )
! The PDE we solve to determine h is written in terms of polar velocity 
! components. 
    call vpolar( u , v , q , theta ) ! Transform from cartesian to polar
    u_grid(:,:) = h0*u
    v_grid(:,:) = h0*v
! The PDE has better continuity properties when written in terms of g 
! rather than h.
    g = log(q*h0)
    S2 = L**2 + M**2 ; T2 = A**2 + B**2 ! Define useful parameters.
    tol = 1.0e-12 ! Define convergence criterion.
!!$! We need to compute the gradient of theta (polar velocity angle). We use a 
!!$!   backward-difference method, for consistency with the iterative scheme.
!!$    allocate( dthetadxi(ny,nx) , dthetadeta(ny,nx) )
!!$    do i = 2 , nx
!!$       do j = 2 , ny
!!$          dthetadxi(j,i) = (theta(j,i)-theta(j,i-1))/dxi ;
!!$          dthetadeta(j,i)= (theta(j,i)-theta(j-1,i))/deta;
!!$       end do
!!$    end do
    call TwoDGradient(theta,dxi,deta,nx,ny,dthetadxi,dthetadeta)

! The PDE can be written in the form alpha*dg/dxi + beta*dg/deta + gamma = 0.
!   We define alpha, beta, and gamma.
    allocate( alpha(nx,ny) , beta(nx,ny) , gamma(nx,ny) )
    alpha = S2*(A*cos(theta)-B*sin(theta))
    beta  = T2*(M*cos(theta)-L*sin(theta))
    gamma = S2*(B*sin(theta)+A*cos(theta))*dthetadxi &
          - T2*(L*cos(theta)+M*sin(theta))*dthetadeta

! Iteratively solve for g using a backward difference based on BC's at g(1,:)
!   and g(:,1). 
    err = 1 ; 
    allocate( temp(nx,ny) )
    temp = g
    do iter = 1, 100
       if (err < tol)then 
          exit
       end if
!       temp(:,1) = temp(:,2) + deta/beta(:,2) * ( gamma(:,2) + alpha(:,2) * 0 )
       do i = 2 , nx
          do j = 2 , ny
             temp(i,j) = (alpha(i,j)*g(i,j-1)/dxi + beta(i,j)*g(i-1,j) &
                   /deta - gamma(i,j) )/(alpha(i,j)/dxi + beta(i,j)/deta)
          end do
       end do
!       temp(:,1) = temp(:,2) ; temp(1,:) = temp(2,:)
       err = maxval(abs(temp-g))
       g   = temp
    end do
    if (iter == 100)then
       write(*,*) "Maximum number of iterations exceeded"
       stop
    end if
!!$! Convert the updated g back to h.
!!$    h = exp(g)/q
!!$    h = h/maxval(h)*h0
! Convert the value of g into grid velocity components.
    u_grid = exp(g)*cos(theta)
    v_grid = exp(g)*sin(theta)
    deallocate( temp )
    deallocate( g , q , theta )
    deallocate( dthetadxi , dthetadeta )
    deallocate( alpha , beta , gamma )

  end subroutine g_update

  subroutine vpolar(u,v,r,theta)
! This subroutine computes polar velocity components based on cartesian 
! components.
    real(8), intent(in) :: u(:,:) , v(:,:)
    real(8), intent(out):: r(:,:) , theta(:,:)

    r     = sqrt(u**2+v**2)
    theta = atan2(v,u)
  end subroutine vpolar

end module h_update_mod

