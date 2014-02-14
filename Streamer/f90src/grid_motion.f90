module grid_motion_mod
! This module is responsible for solving the PDE for h. It uses an iterative, 
!   first-order, backward difference scheme. Boundary values must be provided 
!   at h(1,:) and h(:,1).
  use GeneralUtilities
  real(8), dimension(4), parameter :: h0_array = [0., .25, .5, .999]
contains

  subroutine grid_motion(main,options)
    implicit none
    real(8), dimension(:,:,:,:) :: main
    integer, dimension(:) :: options
    integer :: grid_motion_scheme, nxi, neta, nzeta
    real(8) :: h0, dxi, deta, dzeta

    
    dxi = dxi_a(options(3))
    deta = deta_a(options(4))
    dzeta = dzeta_a(options(5))
    nxi = size(main,2); neta = size(main,3); nzeta = size(main,4)
    grid_motion_scheme = options(6)
    h0 = h0_array(options(7))
    select case(grid_motion_scheme)
    case(0)
       main(15:17,:,:,:) = 0d0
    case(1)
       main(15:17,:,:,:) = h0*main(3:5,:,:,:)
    case(2)
       main(15:17,:,:,:) = h0
!!$    case(3)
!!$       if(size(main,4) .ne. 3)then
!!$          write(*,*) "Incompatible rank for input array!"
!!$          stop
!!$       end if
!!$       allocate(h(neta-2,nxi-2))
!!$       call h_update(&
!!$            transpose(main( 6,2:nxi-1,2:neta-1,2)),&
!!$            transpose(main( 7,2:nxi-1,2:neta-1,2)),&
!!$            transpose(main( 9,2:nxi-1,2:neta-1,2)),&
!!$            transpose(main(10,2:nxi-1,2:neta-1,2)),&
!!$            h,& 
!!$            transpose(main(3,2:nxi-1,2:neta-1,2)),&
!!$            transpose(main(4,2:nxi-1,2:neta-1,2)),& 
!!$            dxi,deta,h0)
!!$       main(15,2:nxi-1,2:neta-1,2) = transpose(h)*main(3,2:nxi-1,2:neta-1,2)
!!$       main(16,2:nxi-1,2:neta-1,2) = transpose(h)*main(4,2:nxi-1,2:neta-1,2)
!!$       deallocate(h)
    case(4)
       if(size(main,4) .ne. 3)then
          write(*,*) "Incompatible rank for input array!"
          stop
       end if
       call g_update(main(6,2:nxi-1,2:neta-1,2),main(7,2:nxi-1,2:neta-1,2),&
            main(9,2:nxi-1,2:neta-1,2),main(10,2:nxi-1,2:neta-1,2),&
            main(15,2:nxi-1,2:neta-1,2),main(16,2:nxi-1,2:neta-1,2),&
            main(3,2:nxi-1,2:neta-1,2),main(4,2:nxi-1,2:neta-1,2),&
            dxi,deta,h0)
       main(15:16,:,:,   1) = main(15:16,:,:,     2)
       main(15:16,:,:,   3) = main(15:16,:,:,     2)
       main(15:16,   1,:,:) = main(15:16,     2,:,:)
       main(15:16, nxi,:,:) = main(15:16, nxi-1,:,:)
       main(15:16,:,   1,:) = main(15:16,:,     2,:)
       main(15:16,:,neta,:) = main(15:16,:,neta-1,:)

    case default
       write(*,*) "Invalid grid motion specification!"
       stop
    end select
!    open(unit=1023,file='iamtestingthisstuff.txt',action='write')
!    write(1023,*) main(15,:,neta-2:neta,2)
  end subroutine grid_motion

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
    allocate(dthetadxi(nx,ny),dthetadeta(nx,ny))
    call vpolar( u , v , q , theta ) ! Transform from cartesian to polar
    u_grid(:,:) = h0*u
    v_grid(:,:) = h0*v
! The PDE has better continuity properties when written in terms of g 
! rather than h.
    g = log(q*h0)
    S2 = L**2 + M**2 ; T2 = A**2 + B**2 ! Define useful parameters.
    tol = 1.0e-12 ! Define convergence criterion.
    call TwoDGradient(theta,dxi,deta,nx,ny,dthetadxi,dthetadeta)

! The PDE can be written in the form alpha*dg/dxi + beta*dg/deta + gamma = 0.
!   We define alpha, beta, and gamma.
    allocate( alpha(nx,ny) , beta(nx,ny) , gamma(nx,ny) )
    alpha = S2*(A*sin(theta)-B*cos(theta))
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
       do i = 2 , nx
          do j = 2 , ny
             temp(i,j) = (alpha(i,j)*g(i-1,j)/dxi + beta(i,j)*g(i,j-1) &
                   /deta - gamma(i,j) )/(alpha(i,j)/dxi + beta(i,j)/deta)
          end do
       end do
       err = maxval(abs(temp-g))
       g   = temp
    end do
    if (iter == 100)then
       write(*,*) "Maximum number of iterations exceeded"
       stop
    end if
! Convert the value of g into grid velocity components.
    u_grid = exp(g)*cos(theta)
    v_grid = exp(g)*sin(theta)
    deallocate(temp)
    deallocate(g, q, theta, S2, T2)
    deallocate(dthetadxi, dthetadeta)
    deallocate(alpha, beta, gamma)

  end subroutine g_update

  subroutine vpolar(u,v,r,theta)
! This subroutine computes polar velocity components based on cartesian 
! components.
    real(8), intent(in) :: u(:,:) , v(:,:)
    real(8), intent(out):: r(:,:) , theta(:,:)

    r     = sqrt(u**2+v**2)
    theta = atan2(v,u)
  end subroutine vpolar

end module grid_motion_mod

