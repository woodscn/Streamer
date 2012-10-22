module Riemann
  use GeneralUtilities
  logical :: verbose = .false.
  logical :: riemann_test_flag 
  logical :: test_flag = .false.
  real(8), dimension(4) :: test_sol
  real(8), dimension(:,:), allocatable :: exact_sol
contains
  
!!$  function fan_Jacobian(J0,psi,u0,u_grid,a0,EV,S,n)
!!$    real(8), intent(in) :: J0, psi, u0, u_grid, a0, EV, S
!!$    integer, intent(in) :: n ! 1 for left fan, 2 for right fan
!!$    
!!$    real(8) :: bn, exp
!!$    integer :: pm 
!!$   
!!$    pm = (-1)**n
!!$    bn = -pm*(gamma_const+1d0)/((gamma_const-1d0)*(u0-ugrid)-pm*2d0*a0)
!!$    exp = (1d0-gamma_const)/(2d0*gamma_const)
!!$
!!$    fan_Jacobian = J0*((psi**exp-bn)/(psi**exp*(1d0-bn)))**(&
!!$         (4d0*a0*EV)/(S*(gamma_const+1d0)))
!!$  end function fan_Jacobian

  subroutine riemann_solve(left, right, dir, nx, x, out, max_wave_speed,&
       riemann_middle_states, riemann_wave_speeds,ierr_out)
    ! Riemann_solve accepts two physical flow states of the form:
    !     [ pressure, mass density, velocity_normal, velocity_tangential_1,
    !       velocity_tangential_2 ]
    !   as well as an average geometric state of the form:
    !     [ pressure, density, v_norm, v_tan1, v_tan2, A, B, C, L, M, N, 
    !       P, Q, R, U, V, W, x, y, z, J ]
    ! Riemann_solve solves the 1-dimensional riemann problem given by 
    !   the left and right states, as well as the grid motion U, V, W.
    !   It returns the state value [ p, rho, v_norm, v_tan1, v_tan2 ] 
    !   at x/t = 0. It also updates the maximum wavespeed encountered,
    !   for use in computing an optimal time step. 
    implicit none
    real(8), dimension(21), intent(in) :: left, right
    integer, intent(in) :: dir, nx
    real(8), dimension(nx), intent(in) :: x
    real(8), dimension(5,nx), intent(out) :: out
    real(8), intent(out) :: max_wave_speed
    real(8), dimension(4), intent(out), optional :: riemann_middle_states
    real(8), dimension(5), intent(out), optional :: riemann_wave_speeds
    integer, intent(out), optional :: ierr_out

    real(8) :: DL, PL, UL, VL, WL, AL
    real(8) :: DR, PR, UR, VR, WR, AR
    real(8) :: Uavg
    real(8), dimension(9) :: met_inv
    real(8) :: J, S
    real(8) :: Pstar, Ustar, DstarL, DstarR
    real(8), dimension(4) :: riemann_middle_states_temp
    real(8), dimension(5) :: riemann_wave_speeds_temp
    real(8), parameter :: tol = 1.d-10
    real(8) :: PsiL, PsiR, temp, fL, fR, dfL, dfR
    integer :: n, ierror

    ierror = 0
    temp = 1d0
    DL =  left(2); PL =  left(1); UL =  left(3); VL =  left(4); WL =  left(5)
    DR = right(2); PR = right(1); UR = right(3); VR = right(4); WR = right(5)
    AL = sqrt(gamma_const*PL/DL); AR = sqrt(gamma_const*PR/DR)
    Uavg = .5d0*(left(15)+right(15))

    met_inv = MetricInverse(.5d0*(left(6:14)+right(6:14)));
    J = Jacobian(.5d0*(left(6:14)+right(6:14)))
    S = sqrt(sum((J*met_inv(dir:9:3))**2))
    
    Pstar = guessp(left,right)
!!$    if(verbose)write(*,*) "Initial guess P = " , Pstar
    temp = 1.d0 ; n = 0
    do while(abs(temp) .gt. tol)
!!$       if(verbose)write(*,*) "Step ",n , "Pstar = " , Pstar; 
       n = n + 1
       if(n .gt. 10)then 
          ierror = 1
          exit
       end if
       PsiL = Pstar/PL
       call u_fun( left,Pstar,fL,dfL)
       PsiR = Pstar/PR
       call u_fun(right,Pstar,fR,dfR)
       temp = ( UR - UL + fR + fL )/( dfL + dfR )
       Pstar = max( Pstar - temp , tol )
!!$       if(verbose)write(*,*) "fL , fR , Pstar = " , fL , fR , Pstar
    end do
    Ustar = .5*(UR+fR+UL-fL)
    DstarL = beta(Pstar/PL)*DL
    DstarR = beta(Pstar/PR)*DR
    riemann_middle_states_temp = [Pstar,Ustar,DstarL,DstarR]
    riemann_wave_speeds_temp = wave_speeds(&
         left,right,dir,riemann_middle_states_temp,Uavg,J,S)
    max_wave_speed = maxval(abs(riemann_wave_speeds_temp))
!!$    if(verbose)write(*,*) Ustar , DstarL , DstarR
    do n = 1, nx
       call sample(x(n),left,right,riemann_middle_states_temp,&
            riemann_wave_speeds_temp,Uavg,J,S,out(:,n))
    end do
    if(present(riemann_middle_states))&
         riemann_middle_states = riemann_middle_states_temp
    if(present(riemann_wave_speeds))&
         riemann_wave_speeds = riemann_wave_speeds_temp
    if(present(ierr_out)) ierr_out = ierror
!!$    if(verbose)write(*,*)"RiemannSolve = ", out
  end subroutine riemann_solve

  pure function wave_speeds(left,right,dir,riemann_middle_states,Uavg,J,S)
    implicit none
    real(8), dimension(21), intent(in) :: left, right
    integer, intent(in) :: dir
    real(8), dimension(4), intent(in) :: riemann_middle_states
    real(8), intent(in) :: Uavg
    real(8), intent(in) :: J, S
    real(8), dimension(5) :: wave_speeds

    real(8), dimension(9) :: met_inv
    real(8) :: Pstar, Ustar, DstarL, DstarR

    wave_speeds = 0d0
    Pstar = riemann_middle_states(1)
    Ustar = riemann_middle_states(2)
    DstarL = riemann_middle_states(3)
    DstarR = riemann_middle_states(4)
    if(Pstar/left(1)>1d0)then
       wave_speeds(1) = S/J*(left(3)-Uavg-sqrt(left(1)*gamma_const/left(2))&
            *sqrt((gamma_const+1d0)/(2d0*gamma_const)*(Pstar/left(1)-1d0)+1d0))
    else
       wave_speeds(1) = S/J*(left(3)-Uavg-sqrt(left(1)*gamma_const/left(2)))
       wave_speeds(2) = S/J*(Ustar-Uavg-sqrt(Pstar*gamma_const/DstarL))
    end if
    wave_speeds(3) = S/J*(Ustar-Uavg)
    if(Pstar/right(1)>1d0)then
       wave_speeds(5) = S/J*(right(3)-Uavg+sqrt(right(1)*gamma_const/right(2))&
            *sqrt((gamma_const+1d0)/(2d0*gamma_const)*(Pstar/right(1)-1d0)+1d0))
    else
       wave_speeds(4) = S/J*(Ustar-Uavg+sqrt(Pstar*gamma_const/DstarR))
       wave_speeds(5) = S/J*(right(3)-Uavg+sqrt(right(1)*gamma_const/right(2)))
    end if
    
  end function wave_speeds

  real(8) pure function guessp(left,right)
    implicit none
    real(8), dimension(:), intent(in) :: left, right
    real(8) :: aL, aR, gL, gR, tol
    aL = sqrt(gamma_const* left(1)/ left(2))
    aR = sqrt(gamma_const*right(1)/right(2))
    tol = 1d-8
    ! Linearised guess
    guessp = .5*(left(1)+right(1))&
         -.125*(right(3)-left(3))*(left(2)+right(2))*(aL+aR)
    if(.not.( guessp .gt. min(left(1),right(1)) .and. &
         guessp .lt. max(left(1),right(1)) &
         .and. max(left(1),right(1))/min(left(1),right(1)) .le. 2.0))then
       if(guessp .lt. min(left(1),right(1)))then
          ! Two-rarefaction solution
          guessp = (&
               (aL+aR-.5*(gamma_const-1.)*(right(3)-left(3)))&
               /(aL/left(1)**((gamma_const-1.)/(2.*gamma_const))&
               +aR/right(1)**((gamma_const-1.)/(2.*gamma_const)))&
               )**(2.*gamma_const/(gamma_const-1.))
       else
          ! Two-shock solution
          gL=sqrt( 2./((gamma_const+1.)* left(2))/(guessp+ left(1)&
               *(gamma_const-1.)/(gamma_const+1.)))
          gR=sqrt( 2./((gamma_const+1.)*right(2))/(guessp+right(1)&
               *(gamma_const-1.)/(gamma_const+1.)))
          guessp = max(tol,(gL*left(1)+gR*right(1)-(right(3)-left(3)))/(gL+gR))
       end if
    end if
    !          guessp = .5*( PL + PR )
  end function guessp

  pure subroutine u_fun(in,pstar,f,df)
    real(8), dimension(:), intent(in) :: in
    real(8), intent(in) :: pstar
    real(8) , intent(out) :: f , df
    real(8) :: A, B, psi, a0
    psi = pstar/in(1)
    a0  = sqrt(gamma_const*in(1)/in(2))
    if( psi .gt. 1. )then
       A = 2.d0/((gamma_const+1.d0)*in(2))
       B = (gamma_const-1.d0)/(gamma_const+1.d0)*in(1)
       f = in(1)*(psi-1.d0)*sqrt(A/(in(1)*psi+B))
       df= sqrt(A/(B+in(1)*psi))*(1.d0-in(1)*(psi-1.d0)/(2.*(B+psi*in(1))))
    else
       f = 2.d0*a0/(gamma_const-1.d0)*(psi**((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
       df= 1.d0/(in(2)*a0)*psi**((gamma_const+1.d0)/(-2.d0*gamma_const))
    end if
  end subroutine u_fun

  pure function beta(psi)
    implicit none
    real(8), intent(in) :: psi
    real(8) :: beta
    if( psi .gt. 1 )then
       beta = ((gamma_const+1.)*psi+gamma_const-1.)/(gamma_const+1.+(gamma_const-1.)*psi)
    else
       beta = psi**(1.d0/gamma_const)
    end if
  end function beta

  pure subroutine sample(x,left,right,riemann_middle_states,riemann_wave_speeds,&
       Uavg,J,S,out)    
    implicit none
    real(8), dimension(:), intent(in) :: left, right, riemann_middle_states
    real(8), dimension(:), intent(in) :: riemann_wave_speeds
    real(8), intent(in) :: x, J, S, Uavg
    real(8), dimension(:) , intent(out) :: out
    real(8) :: PsiL , Pstar, Ustar, DstarL, DstarR
    real(8) :: PsiR , aL , aR , betaL , betaR
    real(8) :: cL , cR , cLT , cLH , cRT , cRH , h
    logical :: test_flag
!!$    test_flag = .false.
!!$    if(verbose) test_flag = .true.
    Pstar  = riemann_middle_states(1); Ustar  = riemann_middle_states(2)
    DstarL = riemann_middle_states(3); DstarR = riemann_middle_states(4)
    PsiL = Pstar/left(1)
    PsiR = Pstar/right(1)
    aL   = sqrt(gamma_const* left(1)/ left(2))
    aR   = sqrt(gamma_const*right(1)/right(2))
    betaL= DstarL/ left(2)
    betaR= DstarR/right(2)
    if( riemann_wave_speeds(3) .gt. x )then
!!$       write(*,*) "The boundary lies to the left of the contact wave"
       out(4) = left(4) ; out(5) = left(5)
       if( PsiL .gt. 1.d0 )then
!!$          if(test_flag) write(*,*) " Left shock"
          cL = riemann_wave_speeds(1)
!!$          if(test_flag) write(*,*) "Left shock speed = " , cL
          if( cL .gt. x )then
!!$             write(*,*) " The boundary lies to the left of the shock"
             out(1) = left(1)
             out(3) = left(3)
             out(2) = left(2)
          else
!!$             write(*,*) " The boundary lies in the left central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarL
          end if
       else
!!$          if(test_flag) write(*,*) "Left rarefaction wave"
          cLT = riemann_wave_speeds(1)
!!$          if(test_flag) write(*,*) "Left rarefaction tail speed = " , cLT
          cLH = riemann_wave_speeds(2)
!!$          if(test_flag) write(*,*) "Left rarefaction head speed = " , cLH
          if( cLT .gt. x )then
!!$             write(*,*) " The boundary lies to the left of the wave"
             out(1) = left(1)
             out(3) = left(3)
             out(2) = left(2)
          elseif( cLH .lt. x )then
!!$             write(*,*) "The boundary lies in the left central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarL
          else
!!$             write(*,*) "The boundary lies within the left expansion wave"
             out(1) = left(1)*( 2.d0*gamma6 + gamma5/aL*(left(3)-Uavg-J/S*x) )**gamma7
             out(3) = left(3) - 2.d0*aL/(gamma_const-1.d0)*((out(1)/left(1))&
                  **((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = left(2)*(out(1)/left(1))**(1.d0/gamma_const)
          end if
       end if
    else
!!$       write(*,*) " The boundary lies to the right of the contact wave"
       out(4) = right(4) ; out(5) = right(5)
       if( PsiR .gt. 1.d0 )then
!!$          if(test_flag) write(*,*) " Right shock"
          cR = riemann_wave_speeds(5)
!!$          if(test_flag) write(*,*) " Right shock speed = " , cR
          if( cR .lt. x )then
!!$             write(*,*) " The boundary lies to the right of the shock"
             out(1) = right(1)
             out(3) = right(3)
             out(2) = right(2)
          else
!!$             write(*,*) " The boundary lies in the right central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarR
          end if
       else
!!$          if(test_flag) write(*,*) " Right rarefaction wave"
          cRT = riemann_wave_speeds(5)
!!$          if(test_flag) write(*,*) "Right rarefaction tail speed = " , cRT
          cRH = riemann_wave_speeds(4)
!!$          if(test_flag) write(*,*) "Right rarefaction head speed = " , cRH
          if( cRT .lt. x )then
!!$             write(*,*) " The boundary lies to the right of the wave"
             out(1) = right(1)
             out(3) = right(3)
             out(2) = right(2)
          elseif( cRH .gt. x )then
!!$             write(*,*) "The boundary lies in the right central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarR
          else
!!$             write(*,*) " The boundary lies within the right expansion wave"
             out(1) = right(1)*( 2.d0*gamma6 - gamma5/aR*(right(3)-Uavg-J/S*x) )**gamma7
             out(3) = right(3) + 2.d0*aR/(gamma_const-1.d0)*((out(1)/right(1))**&
                  ((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = right(2)*(out(1)/right(1))**(1.d0/gamma_const)
          end if
       end if
    end if
  end subroutine sample
end module Riemann

module Riemann_tester
  use Riemann
  real(8), dimension(5) :: riemann_sol
  real(8), dimension(21):: test_geom
  real(8) :: mws_test
  ! Though it's possible to do specific tests of the guessp algorithm,
  ! it seems less-effective to do so, since a bad guess will still converge.
contains
  integer function RieErrorReader(in)
    implicit none
    integer, intent(in) :: in
    select case(in)
    case(0)
       write(*,*) "   All tests passed"
    case(1)
       write(*,*) "   Failure in one of the Riemann Star state checks"
    case default
       write(*,*) "   Unexpected error code"
    end select
  end function RieErrorReader

  integer function RieTester()
    implicit none
    real(8) :: U, V, W
    real(8), dimension(21) :: test_base
    real(8), dimension(21) :: test_1_left
    real(8), dimension(21) :: test_1_right
    real(8), dimension(21) :: test_2_left
    real(8), dimension(21) :: test_2_right
    real(8), dimension(21) :: test_3_left
    real(8), dimension(21) :: test_3_right
    real(8), dimension(21) :: test_4_left
    real(8), dimension(21) :: test_4_right
    real(8), dimension(21) :: test_5_left
    real(8), dimension(21) :: test_5_right
    real(8), dimension(4) :: test_1_sol
    real(8), dimension(4) :: test_2_sol
    real(8), dimension(4) :: test_3_sol
    real(8), dimension(4) :: test_4_sol
    real(8), dimension(4) :: test_5_sol
    real(8), dimension(5) :: test_1_speeds
    real(8), dimension(5) :: test_2_speeds
    real(8), dimension(5) :: test_3_speeds
    real(8), dimension(5) :: test_4_speeds
    real(8), dimension(5) :: test_5_speeds
    integer, parameter :: nx = 1
    real(8), dimension(nx) :: x
    real(8), dimension(5,nx) :: out
    real(8) :: max_wave_speed
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    integer, parameter :: nxfull = 101
    real(8), dimension(nxfull) :: xfull
    real(8), dimension(5,nxfull) :: outfull
    integer :: n    
    RieTester = 0
    
    test_base = [0d0,0d0,0d0,0d0,0d0,&
         1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0,&
         U,V,W,0d0,0d0,0d0,0d0]
    test_1_left = test_base; test_1_right = test_base
    test_2_left = test_base; test_2_right = test_base
    test_3_left = test_base; test_3_right = test_base
    test_4_left = test_base; test_4_right = test_base
    test_5_left = test_base; test_5_right = test_base
    
    test_1_left(1:5)  = [1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_1_right(1:5) = [.1d0, .125d0, 0.d0, 0.d0, 0.d0]
    test_2_left(1:5)  = [.4d0, 1.d0, -2.d0, 0.d0, 0.d0]
    test_2_right(1:5) = [.4d0, 1.d0, 2.d0, 0.d0, 0.d0]
    test_3_left(1:5)  = [1.d3, 1.d0, 0.d0, 0.d0, 0.d0]
    test_3_right(1:5) = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_4_left(1:5)  = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    test_4_right(1:5) = [1.d2, 1.d0, 0.d0, 0.d0, 0.d0]
    test_5_left(1:5)  = [460.894d0, 5.99924d0, 19.5975d0, 0d0, 0d0]
    test_5_right(1:5) = [46.095d0, 5.99242d0, -6.19633d0, 0d0, 0d0]
    test_1_sol = [.30313d0, .92745d0, .42632d0, .26557d0]
    test_1_speeds = [-1.18322d0,-0.0702745d0,0.92745d0,0d0,1.75216d0]
    test_2_sol = [.00189d0, .00000d0, .02185d0, .02185d0]
    test_2_speeds = [-2.74833d0,-0.347992d0,0d0,0.347992d0,2.74833d0]
    test_3_sol = [460.894d0, 19.5975d0, .57506d0, 5.99924d0]
    test_3_speeds = [-37.4166d0,-13.8997d0,19.5975d0,0d0,23.5175d0]
    test_4_sol = [46.095d0, -6.19633d0, 5.99242d0, .57511d0]
    test_4_speeds = [-7.43747d0,0d0,-6.19633d0,4.39658d0,11.8322d0]
    test_5_sol = [1691.64d0, 8.68975d0, 14.2823d0, 31.0426d0]
    test_5_speeds = [0.789631d0,0d0,8.68975d0,0d0,12.2507d0]
    ! Test riemann_solve against the test problems from Toro. This only tests
    ! whether Pstar, Ustar, DstarL and DstarR are correct, and whether the 
    ! computed wave speeds match as computed in Mathematica, and visually
    ! compared (eyeballed) with Toro's plots. This leaves the solutions within
    ! the expansion fan untested.
    test_flag = .true.
    test_geom = [0d0,0d0,0d0,0d0,0d0,&
         1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0,&
         0d0,0d0,0d0,0d0,0d0,0d0,1d0]
!    test_sol = test_1_sol
!      subroutine riemann_solve(left, right, dir, nx, x, out, max_wave_speed,&
!       riemann_middle_states, riemann_wave_speeds)
    call riemann_solve(test_1_left,test_1_right,1,1,[0d0],out,max_wave_speed,&
         riemann_middle_states,riemann_wave_speeds)
    if(maxval(abs(riemann_middle_states-test_1_sol))>5d-6)&
         RieTester = 1
    if(maxval(abs((riemann_wave_speeds-test_1_speeds)/test_1_speeds))>5d-5)&
         RieTester = 1
    call riemann_solve(test_2_left,test_2_right,1,1,[0d0],out,max_wave_speed,&
         riemann_middle_states,riemann_wave_speeds)
    if(maxval(abs(riemann_middle_states-test_2_sol)/test_2_sol)>5d-3)&
         RieTester = 1
    if(maxval(abs(riemann_wave_speeds-test_2_speeds))>5d-4)&
         RieTester = 1
    call riemann_solve(test_3_left,test_3_right,1,1,[0d0],out,max_wave_speed,&
         riemann_middle_states,riemann_wave_speeds)
    if(maxval(abs(riemann_middle_states-test_3_sol)/test_3_sol)>5d-6)&
         RieTester = 1
    if(maxval(abs(riemann_wave_speeds-test_3_speeds))>7d-5)&
         RieTester = 1
    call riemann_solve(test_4_left,test_4_right,1,1,[0d0],out,max_wave_speed,&
         riemann_middle_states,riemann_wave_speeds)
    if(maxval(abs(riemann_middle_states-test_4_sol)/test_4_sol)>5d-6)&
         RieTester = 1
    if(maxval(abs(riemann_wave_speeds-test_4_speeds))>5d-5)&
         RieTester = 1
    call riemann_solve(test_5_left,test_5_right,1,1,[0d0],out,max_wave_speed,&
         riemann_middle_states,riemann_wave_speeds)
    if(maxval(abs(riemann_middle_states-test_5_sol)/test_5_sol)>5d-6)&
         RieTester = 1
    if(maxval(abs(riemann_wave_speeds-test_5_speeds))>8d-5)&
         RieTester = 1

    ! Visual tests are possible to check that the solution within the fan is correct.
    ! This writes the density to a file that can then be plotted.
    do n = 1, nxfull
       xfull(n) = 0d0+1d0/(nxfull-1d0)*(n-1)
    end do
    call riemann_solve(test_1_left,test_1_right,1,nxfull,(xfull-.5d0)/.15d0,outfull,&
         max_wave_speed,riemann_middle_states,riemann_wave_speeds)
    open(unit=9492, file='riemann_test.dat')
    write(9492,*) outfull(2,:)
    write(9492,*) xfull
    close(9492)

  end function RieTester
end module Riemann_tester

module Godunov
  use GeneralUtilities
  use Riemann
  implicit none
  real(8), parameter :: max_dt = 1.d0
  real(8), parameter :: dxi   = 1.d0
  real(8), parameter :: deta  = 1.d0
  real(8), parameter :: dzeta = 1.d0
  real(8), parameter :: dxi_inv   = 1.d0/dxi
  real(8), parameter :: deta_inv  = 1.d0/deta
  real(8), parameter :: dzeta_inv = 1.d0/dzeta
  real(8), parameter :: dV_inv = dxi_inv*deta_inv*dzeta_inv
  integer, parameter :: update_type = 2 ! 1 = FV, 2 = HUI3D

  interface primtocons
     module procedure primtoconsarray
     module procedure primtoconspoint
  end interface primtocons

  interface constoprim
     module procedure constoprimarray
     module procedure constoprimpoint
  end interface constoprim

contains
!  integer elemental function gt0(x)
!    implicit none
!    real(8), intent(in) :: x
!    gt0 = ishft( int(sign(1.d0,x) + 1) , -1 )
!  end function gt0

  ! Computes specific energy, given the array [ p, rho, u, v, w ]
  real(8) function energy_func(in)
    implicit none
    real(8), dimension(:), intent(in) :: in
    energy_func = 0.5d0*(in(3)**2+in(4)**2+in(5)**2) + in(1)/(in(2)*gamma2)
  end function energy_func

  subroutine primtoconsarray(main)
    ! Assumes the structure of prim(:) is :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ] - pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Returns the structure of cons(:) :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(4) :: en_shape
    real(8), allocatable, dimension(:,:,:) :: energy
    en_shape = shape(main)
    allocate(energy(en_shape(2), en_shape(3), en_shape(4)))
    energy = main(21,:,:,:)*( main(1,:,:,:)*gamma1 + 0.5d0*main(2,:,:,:)&
         *( main(3,:,:,:)**2 + main(4,:,:,:)**2 + main(5,:,:,:)**2 ) )
    main(1,:,:,:) = main(2,:,:,:)*main(21,:,:,:)
    main(2,:,:,:) = main(1,:,:,:)*main(3,:,:,:)
    main(3,:,:,:) = main(1,:,:,:)*main(4,:,:,:)
    main(4,:,:,:) = main(1,:,:,:)*main(5,:,:,:)
    main(5,:,:,:) = energy
    deallocate(energy)
  end subroutine primtoconsarray

  subroutine primtoconspoint(main)
    ! Assumes the structure of prim(:) is :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ] - pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Returns the structure of cons(:) :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:), intent(inout) :: main
    real(8) :: energy
    energy = main(21)*( main(1)*gamma1 + 0.5d0*main(2)&
         *( main(3)**2 + main(4)**2 + main(5)**2 ) )
    main(1) = main(2)*main(21)
    main(2) = main(1)*main(3)
    main(3) = main(1)*main(4)
    main(4) = main(1)*main(5)
    main(5) = energy
  end subroutine primtoconspoint
  
  subroutine constoprimarray(main)
    ! Returns the structure of prim(:) :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ] - pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Assume the structure of cons(:) is :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:,:,:,:), intent(inout) :: main
    integer, dimension(4) :: p_shape
    real(8), allocatable, dimension(:,:,:) :: temp1, temp2, p
    allocate( temp1(p_shape(2), p_shape(3), p_shape(4)) )
    allocate( temp2(p_shape(2), p_shape(3), p_shape(4)) )
    allocate(     p(p_shape(2), p_shape(3), p_shape(4)) )
    temp1 = 1.d0/main(21,:,:,:)
    temp2 = 1.d0/main(1,:,:,:)
    p = gamma2*temp1*( main(5,:,:,:) - .5d0*temp2&
         *( main(2,:,:,:)**2 + main(3,:,:,:)**2 + main(4,:,:,:)**2 )&
         )
    main(5,:,:,:) = main(4,:,:,:)*temp2
    main(4,:,:,:) = main(3,:,:,:)*temp2
    main(3,:,:,:) = main(2,:,:,:)*temp2
    main(2,:,:,:) = main(1,:,:,:)*temp1
    main(1,:,:,:) = p
    deallocate( temp1, temp2, p )
  end subroutine constoprimarray

  subroutine constoprimpoint(main)
    ! Returns the structure of prim(:) :
    !   1   2   3  4  5  
    ! [ p, rho, u, v, w ] - pressure, mass density, cartesian velocity components
    ! J is the Jacobian
    ! Assume the structure of cons(:) is :
    !     1       2         3         4       5 
    ! [ rho J, rho J u , rho J v , rho J w , J e ]
    !      e is energy, defined for the ideal gas law by :
    !      .5*rho*(u^2+v^2+w^2) + p/(gamma-1)
    real(8), dimension(:), intent(inout) :: main
    real(8) :: temp1, temp2, p
    temp1 = 1.d0/main(21)
    temp2 = 1.d0/main(1)
    p = gamma2*temp1*( main(5) - .5d0*temp2&
         *( main(2)**2 + main(3)**2 + main(4)**2 )&
         )
    main(5) = main(4)*temp2
    main(4) = main(3)*temp2
    main(3) = main(2)*temp2
    main(2) = main(1)*temp1
    main(1) = p
  end subroutine constoprimpoint

  function invnorm3(in)
    real(8), dimension(3), intent(in) :: in
    real(8) :: invnorm3
    invnorm3 = 1.d0/sqrt( in(1)**2 + in(2)**2 + in(3)**2 )
  end function invnorm3

  subroutine grid_coords(grad, normal, tangential1, tangential2)
! Compute an orthonormal coordinate system given an initial, unnormalized vector.
    real(8), dimension(:), intent(in) :: grad
    real(8), dimension(3), intent(out) :: normal, tangential1, tangential2
    real(8) :: temp1, temp3
    real(8), dimension(3) :: temp2, temp4
    temp1 = invnorm3(grad)
    normal = grad*temp1
    if( grad(2)**2 + grad(3)**2 < EPS )then
       tangential1 = (/ 0.d0, 1.d0, 0.d0 /)
       tangential2 = (/ 0.d0, 0.d0, 1.d0 /)
    else
       temp2 = (/ 0.d0, -grad(3), grad(2) /)
       temp3 = invnorm3(temp2)
       temp4 = (/ grad(2)**2 + grad(3)**2, -grad(1)*grad(2), -grad(1)*grad(3) /)
       tangential1 = temp2*temp3
       tangential2 = temp1*temp3*temp4
    end if
  end subroutine grid_coords

  subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)
    implicit none
    real(8), dimension(:,:,:,:) :: main
    !f2py intent(in,out) :: main
    integer :: bcextent, nx, ny, nz
    real(8) :: dt_in, CFL
    select case(update_type)
    case(1)
       call prim_update_FV(main,bcextent,dt_in,CFL,nx,ny,nz)
    case(2)
       call prim_update_HUI3D(main,bcextent,dt_in,CFL,nx,ny,nz)
    case default
       write(*,*) "Bad update_type value!!!"
       stop
    end select
  end subroutine prim_update

!!$  function compute_interfaces_HUI3D(left_in,right_in,row_ops_mat)
!!$    implicit none
!!$    real(8), dimension(21), intent(in) :: left, right
!!$    integer, dimension(3,3), intent(in) :: row_ops_mat
!!$    real(8), dimension(5) :: compute_interfaces_HUI3D, face
!!$    real(8), dimension(9) :: metric, metric_inverse
!!$    real(8), dimension(3) :: grid_vel
!!$    real(8), dimension(3,3) :: vels_transform, vels_transform_inv
!!$    real(8), dimension(21) :: geom_avg
!!$    metric = .5d0*(left(6:14)+right(6:14))
!!$    metric_inverse = MetricInverse(metric)
!!$    grid_vel = .5d0*(left(15:17)+right(15:17))
!!$    vels_transform = matmul(MetrictoMatrix(metric_inverse),row_ops_mat)
!!$    left(3:5) = matmul(left(3:5),vels_transform)
!!$    center(3:5) = matmul(center(3:5),vels_transform)
!!$    grid_vel = matmul(grid_vel,vels_transform)
!!$    geom_avg = 0d0
!!$    geom_avg(6:14) = center(6:14)
!!$    geom_avg(15:17) = grid_vel
!!$    face = riemann_solve(left,right,geom_avg,max_wave_speed_temp,0.d0)
!!$    vels_transform_inv = matmul(transpose(row_ops_mat),MetrictoMatrix(metric))
!!$    face(3:5) = matmul(vels_transform_inv,face(3:5))
!!$    compute_interfaces_HUI3D = face
!!$
!!$  end function compute_interfaces_HUI3D

  subroutine prim_update_HUI3D(main,bcextent,dt_in,CFL,nx,ny,nz)
    implicit none
    real(8), dimension(21,-1*bcextent:nx+bcextent-1,-1*bcextent:ny+bcextent-1,&
         -1*bcextent:nz+bcextent-1), intent(inout) :: main
    real(8), dimension(21,0:nx-1,0:ny-1,0:nz-1) :: main_temp
    real(8), intent(in), optional :: dt_in
    real(8), intent(in), optional :: CFL
    integer, intent(in) :: nx,ny,nz,bcextent
    integer :: i, j, k, m, n, im, jm, km, ip, jp, kp, case_no
    real(8) :: area, dv_inv, max_wave_speed_temp
    real(8), dimension(5) :: cons, prim
    real(8), dimension(5) :: left_interface, right_interface
    real(8), dimension(5) :: left_flux, right_flux
    real(8), dimension(3) :: grid_vel, grid_pos
    real(8), dimension(3,3) :: row_ops_mat, vels_transform
    real(8), dimension(9) :: metric, metric_inverse
    real(8), dimension(21) :: temp, center, left, right, geom_avg
    integer, parameter :: splitting_type = 1
    real(8) :: dt, junk
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    dt = dt_in
    dv_inv = 1.d0
    do n = 1, 1
       if( mod(n,2) .ne. 0 )then
          case_no = 1
          area = deta*dzeta
!!$          time = dt
          im =-1; ip = 1; jm = 0; jp = 0; km = 0; kp = 0
       else
          case_no = 2
          area = dxi*dzeta
!!$          time = dt
          jm =-1; jp = 1; im = 0; ip = 0; km = 0; kp = 0
       end if
       do k = 0, nz-1
          do j = 0, ny-1
             do i = 0, nx-1
                row_ops_mat = row_ops_mat_func(case_no)
                center = main(:,i,j,k)
                left = main(:,i+im,j+jm,k+km)
                call MUSCL_HUI(main(1:5,i+2*im,j+2*jm,k+2*km),&
                     main(1:5,i+im,j+jm,k+km),main(1:5,i,j,k),&
                     main(1:5,i+ip,j+jp,k+kp),left(1:5),center(1:5))
!!$                metric = .5d0*(left(6:14)+center(6:14))
!!$                metric_inverse = MetricInverse(metric)
!!$                grid_vel = .5d0*(left(15:17)+center(15:17))
!!$                vels_transform = matmul(&
!!$                     MetrictoMatrix(metric_inverse),row_ops_mat)
!!$                do m = 1, 3
!!$                   vels_transform(n,:) = vels_transform(n,:)&
!!$                        /sqrt(sum(vels_transform(n,:)**2))
!!$                end do
!!$                left(3:5) = matmul(vels_transform,left(3:5))
!!$                center(3:5) = matmul(vels_transform,center(3:5))
!!$                grid_vel = matmul(vels_transform,grid_vel)
!!$                geom_avg = 0d0
!!$                geom_avg(6:14) = center(6:14)
!!$                geom_avg(15:17) = grid_vel
                call riemann_solve(left,center,n,1,[0d0],left_interface,&
                     max_wave_speed_temp)
!                left_interface = riemann_solve(left, center, &
!                     .5d0*(left+center), max_wave_speed_temp, 0.d0)
!!$                vels_transform = matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric))
!!$                do m = 1, 3
!!$                   vels_transform(n,:) = vels_transform(n,:)&
!!$                        /sqrt(sum(vels_transform(n,:)**2))
!!$                end do
!!$                left_interface(3:5) = matmul(vels_transform,left_interface(3:5))
                
                center = main(:,i,j,k)
                right = main(:,i+ip,j+jp,k+kp)
                call MUSCL_HUI(main(1:5,i+im,j+jm,k+km),&
                     main(1:5,i,j,k),main(1:5,i+ip,j+jp,k+kp),&
                     main(1:5,i+2*ip,j+2*jp,k+2*kp),center(1:5),right(1:5))
!!$                metric = .5d0*(center(6:14)+right(6:14))
!!$                metric_inverse = MetricInverse(metric)
!!$                grid_vel = .5d0*(center(15:17)+right(15:17))
!!$                vels_transform = matmul(&
!!$                     MetrictoMatrix(metric_inverse),row_ops_mat)
!!$                do m = 1, 3
!!$                   vels_transform(n,:) = vels_transform(n,:)&
!!$                        /sqrt(sum(vels_transform(n,:)**2))
!!$                end do
!!$                right(3:5) = matmul(vels_transform,right(3:5))
!!$                center(3:5) = matmul(vels_transform,center(3:5))
!!$                grid_vel = matmul(vels_transform,grid_vel)
!!$                geom_avg = 0d0
!!$                geom_avg(6:14) = center(6:14)
!!$                geom_avg(15:17) = grid_vel
                call riemann_solve(center,right,n,1,[0d0],right_interface,&
                     max_wave_speed_temp)
!!$                right_interface = riemann_solve(center, right, &
!!$                     .5d0*(right+center), max_wave_speed_temp, 0.d0)
!!$                vels_transform = matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric))
!!$                do m = 1, 3
!!$                   vels_transform(n,:) = vels_transform(n,:)&
!!$                        /sqrt(sum(vels_transform(n,:)**2))
!!$                end do
!!$                right_interface(3:5)=matmul(vels_transform,right_interface(3:5))

!!$                junk = sum(matmul(matmul(MetrictoMatrix(metric_inverse),&
!!$                     row_ops_mat),matmul(transpose(row_ops_mat),&
!!$                     MetrictoMatrix(metric)))&
!!$                     -reshape([1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0],[3,3]))&
!!$                     **2
!!$                if(junk > 1d-14)then
!!$                   write(*,*) "Inverse not working!!"
!!$                   write(*,*) junk
!!$                   stop
!!$                end if
                     
                center = main(:,i,j,k)
!!$                if(n==1)then
!!$                   center(6:8) = center(6:8) + 0d0*dt*area*dv_inv*&
!!$                        (right_interface(3:5)-left_interface(3:5))
!!$                   left_flux = flux(left_interface,center,1)
!!$                   right_flux = flux(right_interface,center,1)
!!$                else if(n==2)then
!!$                   center(9:11) = center(9:11) + 0d0*dt*area*dv_inv*&
!!$                        (right_interface(3:5)-left_interface(3:5))
!!$                   left_flux = flux(left_interface,center,2)
!!$                   right_flux = flux(right_interface,center,2)
!!$                else if(n==3)then
!!$                   center(12:14) = center(12:14) + 0d0*dt*area*dv_inv*&
!!$                        (right_interface(3:5)-left_interface(3:5))
!!$                   left_flux = flux(left_interface,center,3)
!!$                   right_flux = flux(right_interface,center,3)
!!$                end if
                left_flux  = flux( left_interface,center,n)
                right_flux = flux(right_interface,center,n)
                call primtocons(center)
                center(1:5) = center(1:5) - dt*area*dv_inv*&
                     (right_flux-left_flux)
                call constoprim(center)
!!$                center(15:17) = center(3:5)*.0d0
!!$                center(18:20) = center(18:20) + 0d0*dt*area*dv_inv*center(15:17)
!!$                center(21) = Jacobian(center(6:14))
!!$
!!$                if(any(left_flux .ne. right_flux))then
!!$                   write(*,*) "i, j, k, = ",i,j,k
!!$                   write(*,*) "ip = ",i+ip,j+jp,k+kp
!!$                   write(*,*) "im = ",i+im,j+jm,k+km
!!$                   write(*,*) "n = ", n
!!$                   write(*,*) "  left = ", main(1:5,i+im,j+jm,k+km)
!!$                   write(*,*) "center = ", main(1:5,i,j,k) 
!!$                   write(*,*) " right = ", main(1:5,i+ip,j+jp,k+kp)
!!$                   write(*,*) " left_interface = ", left_interface
!!$                   write(*,*) "right_interface = ", right_interface
!!$                   call riemann_solve(main(1:5,i,j,k),main(1:5,i+ip,j+jp,k+kp),&
!!$                        n,1,[0d0],right_interface,max_wave_speed_temp,&
!!$                        riemann_middle_states,riemann_wave_speeds)
!!$                   write(*,*) riemann_middle_states
!!$                   write(*,*) riemann_wave_speeds
!!$                   read(*,*)
!!$                   
!!$                end if
                main_temp(:,i,j,k) = center
                if(abs(main_temp(5,i,j,k)) > 1d-12)then
                   write(*,*) "z-component of velocity detected!!!"
                   read(*,*)
                end if
             end do
          end do
       end do
       main(:,0:nx-1,0:ny-1,0:nz-1) = main_temp
!!$       main(15:17,0:nx-1,0:ny-1,0:nz-1) = main(3:5,0:nx-1,0:ny-1,0:nz-1)*.0d0
!!$       main(18:20,0:nx-1,0:ny-1,0:nz-1) = main(18:20,0:nx-1,0:ny-1,0:nz-1)&
!!$            *dt*area*dv_inv*main(15:17,0:nx-1,0:ny-1,0:nz-1)
    end do
  end subroutine prim_update_HUI3D

  subroutine prim_update_FV(main,bcextent,dt_in,CFL,nx,ny,nz)
! Advance the solution using the integral form of the equations
! This subroutine assumes that main is the full array of primitive variables. main 
! must also include the boundary cell values. That is, main contains a 0-index and 
! an n + 1 index containing the contents prescribed by the boundary conditions.
    implicit none
    real(8), dimension(21,-1*bcextent:nx+bcextent-1,-1*bcextent:ny+bcextent-1,&
         -1*bcextent:nz+bcextent-1), intent(inout) :: main
!f2py intent(in,out) :: main
    real(8), intent(in), optional :: dt_in
    real(8), intent(in), optional :: CFL
    integer, intent(in) :: nx,ny,nz,bcextent
    integer :: i, j, k
!    real(8), dimension(14,3,nx+1,ny+1,nz+1) :: fluxes
    real(8), dimension(14,0:nx,0:ny-1,0:nz-1) :: fluxx
    real(8), dimension(14,0:nx-1,0:ny,0:nz-1) :: fluxy
    real(8), dimension(14,0:nx-1,0:ny-1,0:nz) :: fluxz
    real(8), dimension(14,nx,ny,nz) :: temp
    real(8) :: max_wave_speed, dt, max_dt_grid
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8), dimension(3) :: GradX, GradY, GradZ
    real(8), dimension(9) :: geom_avg
    real(8), dimension(3,3) :: dXidX, dXdXi
    real(8), dimension(5) :: interface_vars
    real(8), dimension(21) :: StateL, StateR
    real(8), dimension(14) :: junk
  ! Riemann_solve expects the left and right states to express velocities in
  ! grid-oriented components: normal, tangential, tangential.
    do k = 0, nz-1
       do j = 0, ny-1
          do i = 0, nx
             call compute_fluxes_FV(main(1:5,i-1,j,k), main(1:5,i,j,k),&
                  .5d0*(main(:,i-1,j,k)+main(:,i,j,k)),fluxx(:,i,j,k),&
                  1,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
          end do
       end do
    end do
    do k = 0, nz-1
       do j = 0, ny
          do i = 0, nx-1
             call compute_fluxes_FV(main(1:5,i,j-1,k), main(1:5,i,j,k),&
                  .5d0*(main(:,i,j-1,k)+main(:,i,j,k)),fluxy(:,i,j,k),&
                  2,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta],debug_flag=.true.)
          end do
       end do
    end do
    do k = 0, nz
       do j = 0, ny-1
          do i = 0, nx-1
             call compute_fluxes_FV(main(:,i,j,k-1), main(:,i,j,k),&
                  .5d0*(main(:,i,j,k-1)+main(:,i,j,k)),fluxz(:,i,j,k),&
                  3,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
          end do
       end do
    end do
    max_wave_speed = max(max_wave_speed,EPS)
    if(dt_in > 0.d0)then
       dt = dt_in
    elseif(CFL>0.d0)then
       dt = min(CFL/max_wave_speed,max_dt)
       dt = min(maxval((main(6,0,0:ny-1,0:nz-1)-main(18,0,0:ny-1,0:nz-1))&
            /main(15,0,0:ny-1,0:nz-1)),dt)
    else
       write(*,*) "Error in Godunov prim_update, no time step information given"
       read(*,*)
       stop
    end if
    call primtocons(main(:,0:nx-1,0:ny-1,0:nz-1))
    main(1:14,0:nx-1,0:ny-1,0:nz-1) = main(1:14,0:nx-1,0:ny-1,0:nz-1) - (&
         (fluxx(:,1:nx,:,:)-fluxx(:,0:nx-1,:,:))*deta*dzeta + &
         (fluxy(:,:,1:ny,:)-fluxy(:,:,0:ny-1,:))*dxi*dzeta + &
         (fluxz(:,:,:,1:nz)-fluxz(:,:,:,0:nz-1))*dxi*deta &
         )*dt*dV_inv

    call constoprim(main(:,0:nx-1,0:ny-1,0:nz-1))
! Update grid velocity
!!$    call grid_velocity_update(main)
!    main(15:17,:,:,:) = main(3:5,:,:,:)*.25d0
! Update grid position
    main(18:20,0:nx-1,0:ny-1,0:nz-1) = main(18:20,0:nx-1,0:ny-1,0:nz-1) + &
         main(15:17,0:nx-1,0:ny-1,0:nz-1)*dt
! Update extra variables
  do k = 0, nz-1
     do j = 0, ny-1
        do i = 0, nx-1
           main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
        end do
     end do
  end do
  end subroutine prim_update_FV

  subroutine compute_fluxes_FV(inL, inR, geom_avg, flux_vec, case_no,&
       max_wave_speed,dt,dV_in,debug_flag)
    implicit none
    real(8), dimension(:), intent(in) :: inL, inR
    real(8), dimension(21) :: StateL, StateR
    real(8), dimension(:) :: geom_avg
    real(8), dimension(3),intent(in) :: dV_in
    real(8) :: dA, dV_inv
    real(8),intent(in) :: dt
    integer, intent(in) :: case_no
    logical, intent(in), optional :: debug_flag
    real(8), dimension(:), intent(out) :: flux_vec
    real(8), intent(inout) :: max_wave_speed
    real(8), dimension(5) :: interface_vars
    real(8), dimension(3) :: GradXi, GradEta, GradZeta
    real(8), dimension(3) :: GradX, GradY, GradZ
    real(8), dimension(3,3) :: dX_dXi_u, dXi_dX_u
    real(8) :: temp_wave_speed
    integer, dimension(3,3) :: row_ops_mat
    StateL = inL
    StateR = inR
    flux_vec = 0.d0
    dV_inv = 1.d0/(product(dV_in))
    call ComputationalGrads(geom_avg(6:14),Jacobian(geom_avg(6:14)),&
         GradXi,GradEta,gradZeta)
    GradX = geom_avg(6:12:3)
    GradY = geom_avg(7:13:3)
    GradZ = geom_avg(8:14:3)
    ! Create normalized transformation matrices
    ! These matrices are given as:
    ! dXi_dX = | GradXi(1)  GradEta(1)  GradZeta(1) |
    !          | GradXi(2)  GradEta(2)  GradZeta(2) |
    !          | GradXi(3)  GradEta(3)  GradZeta(3) |
    ! However, when computing fluxes for the eta and zeta directions,
    ! these must be reordered to fit with the riemann solver, which 
    ! requires that vector components be given as | normal, tangential, tangential |.
    ! dX_dXi begins as the matrix inverse of dXi_dX, but this re-ordering,
    ! corresponding to elementary column operations, affects this inverse
    ! property. In order to mirror this effect, dX_dXi must also be 
    ! re-ordered, though with row operations instead. 
    dXi_dX_u = GradstoMatrix(GradXi/sqrt(sum(GradXi**2)),&
         GradEta/sqrt(sum(GradEta**2)),GradZeta/sqrt(sum(GradZeta**2)))
    dX_dXi_u = GradstoMatrix(GradX/sqrt(sum(GradX**2)),&
         GradY/sqrt(sum(GradY**2)),GradZ/sqrt(sum(GradZ**2)))
    dX_dXi_u = reshape([geom_avg(6),geom_avg(9),geom_avg(12),&
         geom_avg(7),geom_avg(10),geom_avg(13),&
         geom_avg(8),geom_avg(11),geom_avg(14)],[3,3])
    dXi_dX_u = MetrictoMatrix(MetricInverse(geom_avg(6:14)))
    row_ops_mat = row_ops_mat_func(case_no)
    select case(case_no)
    case(1)
       flux_vec(6:8) = -geom_avg(15:17)
       dA = dV_in(2)*dV_in(3)
    case(2)
       flux_vec(9:11) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(3)
    case(3)
       flux_vec(12:14) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(2)
    case default
       write(*,*) "Invalid case_no in compute_fluxes -- case_no = ",case_no
       stop
    end select
    dXi_dX_u = matmul(dXi_dX_u,row_ops_mat)
    dX_dXi_u = matmul(transpose(row_ops_mat),dX_dXi_u)
    if(maxval(matmul(dXi_dX_u,dX_dXi_u)-reshape([1,0,0,0,1,0,0,0,1],[3,3]))>1d-13)then
       write(*,*) "Inverses not working!"
       write(*,*) 
       write(*,*) dXi_dX_u
       write(*,*) dX_dXi_u
       write(*,*) matmul(dXi_dX_u,dX_dXi_u)-reshape([1,0,0,0,1,0,0,0,1],[3,3])
       read(*,*)
    end if
    StateL(3:5) = matmul(dXi_dX_u,StateL(3:5))
    StateR(3:5) = matmul(dXi_dX_u,StateR(3:5))
    geom_avg(15:17) = matmul(dXi_dX_u,geom_avg(15:17))
!!$    interface_vars = riemann_solve(StateL,StateR,geom_avg,temp_wave_speed)
    call riemann_solve(StateL,StateR,case_no,1,[0d0],interface_vars,&
         temp_wave_speed)
    interface_vars(3:5) = matmul(dX_dXi_u,interface_vars(3:5))
!!$    write(*,*) "Interface = "
!!$    write(*,*) interface_vars
    geom_avg(15:17) = matmul(dX_dXi_u,geom_avg(15:17))
!!$    geom_avg(6:14) = geom_avg(6:14)-flux_vec(6:14)*dt*dV_inv
!!$    flux_vec(6:14) = 0.d0
    flux_vec(1:5) = flux(interface_vars,geom_avg,case_no)
    max_wave_speed = max(max_wave_speed,temp_wave_speed)
  end subroutine compute_fluxes_FV

  function row_ops_mat_func(case_no)
    implicit none
    integer, dimension(3,3) :: row_ops_mat_func, row_ops_mat
    integer, intent(in) :: case_no
        select case(case_no)
    case(1)
       row_ops_mat = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    case(2)
       row_ops_mat = reshape([0,1,0,1,0,0,0,0,1],[3,3])
    case(3)
       row_ops_mat = reshape([0,0,1,1,0,0,0,1,0],[3,3])
    case default
       write(*,*) "Invalid case_no in compute_fluxes -- case_no = ",case_no
       stop
    end select
    row_ops_mat_func = row_ops_mat
  end function row_ops_mat_func

  function flux(in,geom_avg,case_no)
    implicit none
    real(8), dimension(:), intent(in) :: in, geom_avg
    integer, intent(in) :: case_no
    real(8), dimension(5) :: flux
    real(8) :: D, grads(3,3), grad(3), J, Jinv
!!$    write(*,*) "Got this far"
    J = geom_avg(21)
    J = Jacobian(geom_avg(6:14))
    Jinv = 1.d0/geom_avg(21)
    call ComputationalGrads(geom_avg(6:14),J,grads(:,1),grads(:,2),grads(:,3))
    grad = grads(:,case_no)
    D = sum(([ in(3), in(4), in(5) ] - &
         [ geom_avg(15), geom_avg(16), geom_avg(17) ])*grad)
    flux = [&
         in(2)*J*D, &
         in(2)*J*D*in(3) + J*grad(1)*in(1), &
         in(2)*J*D*in(4) + J*grad(2)*in(1), &
         in(2)*J*D*in(5) + J*grad(3)*in(1), &
         in(2)*J*D*energy_func(in(1:5)) + J*sum(grad*in(3:5))*in(1) ]
!!$    write(*,*) "J = ", J
!!$    write(*,*) "D = ", D
!!$    write(*,*) "e = ", energy_func(in(1:5))
!!$    write(*,*) "v = ", in(3:5)
!!$    write(*,*) "V = ", geom_avg(15:17)
!!$    write(*,*) "grad = ", grad
!!$    write(*,*) "Energy term = ", J*in(2)*D*energy_func(in(1:5))
!!$    write(*,*) "Case No. = ", case_no
  end function flux

  subroutine MUSCL_HUI(leftleft, left, right, rightright, outleft, outright)
    implicit none
    real(8), dimension(5), intent(in) :: leftleft, left, right, rightright
    real(8), dimension(5), intent(out) :: outleft, outright
    integer :: n
    
    do n = 1, 5
       outleft(n)  =  left(n) + 5d-1*( left(n) - leftleft(n))*minmod(&
            (right(n)-left(n))/( left(n) - leftleft(n)))
       outright(n) = right(n) - 5d-1*(rightright(n)-right(n))*minmod(&
            (right(n)-left(n))/(rightright(n)-right(n)))
    end do
    outleft  =  left + 5d-1*( left - leftleft)*&
         minmod((right-left)/( left - leftleft))
    outright = right - 5d-1*(rightright-right)*&
         minmod((right-left)/(rightright-right))
    
  contains
    elemental function minmod(in)
      implicit none
      real(8), intent(in) :: in
      real(8) :: minmod
      minmod = max(0d0,min(1d0,in))
    end function minmod
  end subroutine MUSCL_HUI

end module Godunov

module Godunov_tester
  ! Test the following routines:
  ! function riemann_solve(left, right, geom_avg, max_wave_speed, verbose_flag, t_out)
  ! function energy_func(in)
  ! subroutine primtocons(main,nx,ny,nz)
  ! subroutine constoprim(main,nx,ny,nz)
  ! function invnorm3(in)
  ! subroutine grid_coords(grad, normal, tangential1, tangential2)
  ! subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)
  ! subroutine compute_fluxes(inL, inR, geom_avg, flux_vec, case_no,&
  !      max_wave_speed,dt,dV_in,debug_flag)
  ! function flux(in,geom_avg,case_no)
  use Godunov
  use GeneralUtilitiesTest
contains
  integer function GodErrorReader(in)
    integer, intent(in) :: in
    select case(in)
    case(0) 
       write(*,*) "   All tests passed"
    case(1)
       write(*,*) "   Primitive-Conservative mutual inverse test failed"
    case(2)
       write(*,*) "   Grid_coords failed to return an orthonormal system"
    case(3)
       write(*,*) "   Flux fails provided test problem"
    case default
       write(*,*) "   Unexpected error code"
    end select
    GodErrorReader = 0
  end function GodErrorReader

  integer function GodConvergenceTester1D(left,right,t_out)
    implicit none
    real(8), dimension(5), intent(in) :: left, right
    real(8), intent(in) :: t_out
    integer, parameter :: nx = 400
    real(8), dimension(21,-1:nx/8,3,3) :: Rie_1D_1
    real(8), dimension(21,-1:nx/4,3,3) :: Rie_1D_2
    real(8), dimension(21,-1:nx/2,3,3) :: Rie_1D_3
    real(8), dimension(21,-1:nx,3,3) :: Rie_1D_4
    real(8), dimension(5,nx) :: Rie_1D_exact
    integer :: l, m, n
    real(8), dimension(nx) :: x
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    real(8) :: max_wave_speed, dt, dx
    real(8), dimension(4) :: dxes, rmserrors
    real(8), dimension(2) :: fitted_poly

    GodConvergenceTester1D = 0
    dx = 1d0
    do n = 1, nx
       x(n) = ((n-1)*dx-.5d0*nx)
    end do
!    call riemann_solve(test_1_left,test_1_right,1,1,[0d0],out,max_wave_speed,&
!         riemann_middle_states,riemann_wave_speeds)
    do n = 1, 3
       do m = 1, 3
          dx = 1d0/(nx/8-1)
          do l = -1, nx/8/2-1
             Rie_1D_1(1:5,l,m,n) = left
             Rie_1D_1(18,l,m,n) = dx*l
          end do
          do l = nx/8/2, nx/8
             Rie_1D_1(1:5,l,m,n) = right
             Rie_1D_1(18,l,m,n) = dx*l
          end do
          Rie_1D_1( 6,:,:,:) = 8d0/nx
          Rie_1D_1(10,:,:,:) = 8d0/nx
          Rie_1D_1(14,:,:,:) = 8d0/nx
          Rie_1D_1(21,:,:,:) = (8d0/nx)**3
          dx = 1d0/(nx/4-1)
          do l = -1, nx/4/2-1
             Rie_1D_2(1:5,l,m,n) = left
             Rie_1D_2(18,l,m,n) = dx*l
          end do
          do l = nx/4/2, nx/4
             Rie_1D_2(1:5,l,m,n) = right
             Rie_1D_2(18,l,m,n) = dx*l
          end do
          Rie_1D_2( 6,:,:,:) = 4d0/nx
          Rie_1D_2(10,:,:,:) = 4d0/nx
          Rie_1D_2(14,:,:,:) = 4d0/nx
          Rie_1D_2(21,:,:,:) = (4d0/nx)**3
          dx = 1d0/(nx/2-1)
          do l = -1, nx/2/2-1
             Rie_1D_3(1:5,l,m,n) = left
             Rie_1D_3(18,l,m,n) = dx*l
          end do
          do l = nx/2/2, nx/2
             Rie_1D_3(1:5,l,m,n) = right
             Rie_1D_3(18,l,m,n) = dx*l
          end do
          Rie_1D_3( 6,:,:,:) = 2d0/nx
          Rie_1D_3(10,:,:,:) = 2d0/nx
          Rie_1D_3(14,:,:,:) = 2d0/nx
          Rie_1D_3(21,:,:,:) = (2d0/nx)**3
          dx = 1d0/(nx-1)
          do l = -1, nx/2-1
             Rie_1D_4(1:5,l,m,n) = left
             Rie_1D_4(18,l,m,n) = dx*l
          end do
          do l = nx/2, nx
             Rie_1D_4(1:5,l,m,n) = right
             Rie_1D_4(18,l,m,n) = dx*l
          end do
          Rie_1D_4( 6,:,:,:) = 1d0/nx
          Rie_1D_4(10,:,:,:) = 1d0/nx
          Rie_1D_4(14,:,:,:) = 1d0/nx
          Rie_1D_4(21,:,:,:) = (1d0/nx)**3
       end do
    end do
    Rie_1D_exact = 0d0
    call riemann_solve(Rie_1D_4(:,-1,2,2),Rie_1D_4(:,nx,2,2),1,nx,&
         ((Rie_1D_4(18,0:nx-1,2,2)-.5d0)/Rie_1D_4(6,0:nx-1,2,2))/t_out,&
         Rie_1D_exact,max_wave_speed,riemann_middle_states,riemann_wave_speeds)
    write(*,*) "Middle States ",riemann_middle_states
!!$    Rie_1D_1 = right; Rie_1D_1(:,1:size(Rie_1D_1,2)/2,:,:) = left
!!$    Rie_1D_2 = right; Rie_1D_2(:,1:size(Rie_1D_2,2)/2,:,:) = left
!!$    Rie_1D_3 = right; Rie_1D_3(:,1:size(Rie_1D_3,2)/2,:,:) = left
!!$    Rie_1D_4 = right; Rie_1D_4(:,1:size(Rie_1D_4,2)/2,:,:) = left
    dt = 1d-4
    do n = 1, 100000
!!$       write(*,*) "t = ",n*dt
       if(n*dt.ge.t_out) exit
       call prim_update(Rie_1D_1,1,dt,.7d0,nx/8,1,1)
       call prim_update(Rie_1D_2,1,dt,.7d0,nx/4,1,1)
       call prim_update(Rie_1D_3,1,dt,.7d0,nx/2,1,1)
       call prim_update(Rie_1D_4,1,dt,.7d0,nx  ,1,1)
    end do
    dxes = [8d0/(nx-1),4d0/(nx-1),2d0/(nx-1),1d0/(nx-1)]
    rmserrors = [&
         sqrt(sum((Rie_1D_1(2,0:nx/8-1,2,2)-Rie_1D_exact(2,1:nx:8))**2)/(nx/8)),&
         sqrt(sum((Rie_1D_2(2,0:nx/4-1,2,2)-Rie_1D_exact(2,1:nx:4))**2)/(nx/4)),&
         sqrt(sum((Rie_1D_3(2,0:nx/2-1,2,2)-Rie_1D_exact(2,1:nx:2))**2)/(nx/2)),&
         sqrt(sum((Rie_1D_4(2,0:nx  -1,2,2)-Rie_1D_exact(2,:     ))**2)/(nx  ))]
    
    fitted_poly = polyfit(log(dxes),log(rmserrors),1)
    write(*,*) "Polyfit results: ",fitted_poly
    open(unit=92922,file="exact.dat")
    open(unit=92923,file="Rie1.dat")
    open(unit=92924,file="Rie2.dat")
    open(unit=92925,file="Rie3.dat")
    open(unit=92926,file="Rie4.dat")
    write(92922,*) Rie_1D_4(18,0:nx-1,2,2)
    write(92922,*) Rie_1D_exact(2,:)
    write(92923,*) Rie_1D_1(18,0:nx/8-1,2,2)
    write(92923,*) Rie_1D_1( 2,0:nx/8-1,2,2)
    write(92924,*) Rie_1D_2(18,0:nx/4-1,2,2)
    write(92924,*) Rie_1D_2( 2,0:nx/4-1,2,2)
    write(92925,*) Rie_1D_3(18,0:nx/2-1,2,2)
    write(92925,*) Rie_1D_3( 2,0:nx/2-1,2,2)
    write(92926,*) Rie_1D_4(18,0:nx  -1,2,2)
    write(92926,*) Rie_1D_4( 2,0:nx  -1,2,2)
    close(92922); close(92923); close(92924); close(92925); close(92926)

    if(fitted_poly(2) < 1d0) GodConvergenceTester1D = 2
    !       call prim_update(riemann_test_array_1d,1,dt,.7d0,nx-2,ny-2,1)

  end function GodConvergenceTester1D

  integer function GodTester()
    use TimeAdvancementStuff
    implicit none
    real(8), dimension(21,3,4,3) :: main
    integer :: i, j, nx, ny
    real(8), dimension(5) :: riemann_test
    real(8), dimension(21) :: geom_test
    real(8), dimension(21,3,3,3) :: constoprimtest1, constoprimtest2
    real(8), dimension(5) :: left_test, right_test
    real(8) :: max_wave_speed
    real(8), dimension(3) :: grad, norm, tan1, tan2
    real(8) :: test1, test2, test3, dt
    real(8), dimension(21), parameter :: flux_seed = [0.270991,0.0613936,&
         0.748627,0.24336,0.374409,-0.423353,0.323678,-0.333788,-0.53982,&
         -0.991525,0.75888,-0.921888,-0.631486,0.801117,0.634729,-0.750043,&
         0.412858,-0.0734701,-0.0563167,0.97121,1.]
    real(8), dimension(5), parameter :: flux_test_x = [-0.0171431,-0.0982248,&
         -0.0765652,-0.161747,-0.335386]
    real(8), dimension(5), parameter :: flux_test_y = [-0.0384556,-0.0419378,&
         -0.184655,-0.167707,-0.548871]
    real(8), dimension(5), parameter :: flux_test_z = [0.0285833,-0.00172451,&
         0.142847,0.171804,0.402354]
    real(8), dimension(21,62,102,3) :: riemann_test_array
    real(8), dimension(21,102,3,3) :: riemann_test_array_1d
    real(8), dimension(5) :: left, right
    integer :: n
    
    GodTester = 0

    ! Check that primtocons and constoprim are mutual inverses
    ! Are there invertibility problems from non-physical variables,
    ! such as negative pressures, etc? Random_number only returns
    ! positive numbers, I suppose.
    call random_number(constoprimtest1)
    constoprimtest2 = constoprimtest1
    call constoprim(constoprimtest1)
    call primtocons(constoprimtest1)
    if(sqrt(sum((constoprimtest1-constoprimtest2)**2)&
         /size(constoprimtest1))>EPS) GodTester = 1

    ! Ensure that the coordinate system returned from grid_coords
    ! is orthonormal, and that norm is parallel to grad.
    
    call random_number(grad)
    grad = grad - .5d0
    call grid_coords(grad, norm, tan1, tan2)
    if(maxval(abs([norm2(norm), norm2(tan1), norm2(tan2)]-1.d0))>EPS)then
       GodTester = 2
    end if
    if(dot_product(grad/norm2(grad),norm)-1.d0>EPS) GodTester = 2
    if(abs(dot_product(norm,tan1))>EPS) GodTester = 2
    if(abs(dot_product(norm,tan2))>EPS) GodTester = 2
    if(abs(dot_product(tan1,tan2))>EPS) GodTester = 2
       
    ! Test flux routine against Mathematica-generated test problem & solution
    ! Mathematica only returns to a set precision.
    if(sqrt(sum((flux(flux_seed(1:5),flux_seed,1) - flux_test_x)**2)/5.)>EPSs)then
       GodTester = 3
    end if
    if(sqrt(sum((flux(flux_seed(1:5),flux_seed,2) - flux_test_y)**2)/5.)>EPSs)then
       GodTester = 3
    end if
    if(sqrt(sum((flux(flux_seed(1:5),flux_seed,3) - flux_test_z)**2)/5.)>EPSs)then
       GodTester = 3
    end if

    ! It is important to also test prim_update, though the only known way 
    ! to do this is via specific convergence testing. For now, that testing
    ! must be approved manually.
    ! subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)

!!$    left = [0d0,0d0,0d0,0d0,0d0,&
!!$         1d-2,0d0,0d0,&
!!$         0d0,1d-2,0d0,&
!!$         0d0,0d0,1d-2,&
!!$         0d0,0d0,0d0,&
!!$         0d0,0d0,0d0,1d-6]
!!$    right = left
    left(1:5)  = [1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    right(1:5) = [.1d0, .125d0, 0.d0, 0.d0, 0.d0]
    GodTester = GodConvergenceTester1D(left,right,.25d0)
    stop
    left(1:5)  = [.4d0, 1.d0, -2.d0, 0.d0, 0.d0]
    right(1:5) = [.4d0, 1.d0, 2.d0, 0.d0, 0.d0]
    GodTester = GodConvergenceTester1D(left,right,.15d0)
    left(1:5)  = [1.d3, 1.d0, 0.d0, 0.d0, 0.d0]
    right(1:5) = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    GodTester = GodConvergenceTester1D(left,right,.012d0)
    left(1:5)  = [.01d0, 1.d0, 0.d0, 0.d0, 0.d0]
    right(1:5) = [1.d2, 1.d0, 0.d0, 0.d0, 0.d0]
    GodTester = GodConvergenceTester1D(left,right,.035d0)
    left(1:5)  = [460.894d0, 5.99924d0, 19.5975d0, 0d0, 0d0]
    right(1:5) = [46.095d0, 5.99242d0, -6.19633d0, 0d0, 0d0]
    GodTester = GodConvergenceTester1D(left,right,.035d0)
    stop

    call RieInit1D(riemann_test_array_1d)
    nx = size(riemann_test_array_1d,2)
    ny = size(riemann_test_array_1d,3)
    call write_files_matlab(riemann_test_array_1d(:,:,:,2),1d-4*n,&
         nx,ny,.true.)!(x,y,E,h,t,dt)
    dt = 1d-3
    do n = 1, 50000
       if(n*dt>.25d0) exit
       write(*,*) "time step = ", n
       riemann_test_array_1d(:,:,1,:) = riemann_test_array_1d(:,:,2,:)
       riemann_test_array_1d(:,:,3,:) = riemann_test_array_1d(:,:,2,:)
       riemann_test_array_1d(:,:,:,1) = riemann_test_array_1d(:,:,:,2)
       riemann_test_array_1d(:,:,:,3) = riemann_test_array_1d(:,:,:,2)
       call prim_update(riemann_test_array_1d,1,dt,.7d0,nx-2,ny-2,1)
       call write_files_matlab(riemann_test_array_1d(:,:,:,2),dt*n,&
            nx,ny,.false.)!(x,y,E,h,t,dt)
!!$       write(*,*) riemann_test_array_1d(1,20:60,1,1) 
!!$       read(*,*)
    end do
       
         ! subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)

    write(*,*) "Got this far!!"
    write(*,*) riemann_test_array_1d(1,:,2,2)
  end function GodTester

  subroutine RieInit1D(main)
    implicit none
    real(8), dimension(21,-1:100,-1:1,-1:1) :: main
    integer :: i
    main = 0d0
    main(1,-1:49,:,:) = 1d0 ; main(1,50:100,:,:) = 1d-1
    main(2,-1:49,:,:) = 1d0 ; main(2,50:100,:,:) = 125d-3
    main(3,-1:49,:,:) = 0d0 ; main(3,50:100,:,:) = 0d0
    main(6,:,:,:) = 1d-2 ; main(10,:,:,:) = 1d-2 ; main(14,:,:,:) = 1d-2
!    main(6,:,:,:) = 1d0; main(10,:,:,:) = 1d0; main(14,:,:,:) = 1d0
    do i = -1, 100
       main(18,i,:,:) = (real(i,8)+.5d0)*1d-2
    end do
    main(21,:,:,:) = 1d-6
!    main(21,:,:,:) = 1d0
  end subroutine RieInit1D

  subroutine GodRieInit(main)
    implicit none
    real(8), dimension(21,62,102,3), intent(inout) :: main
    integer :: i, j

    main( 1,:,:,:) = 1.d0
    main( 2,:,:,:) = 1.d0
    main( 3,:,:,:) = 2.4d0*sqrt(1.4d0*main(1,:,:,:)/main(2,:,:,:))
    main( 4,:,:,:) = 0.d0
    main( 5,:,:,:) = 0.d0
    main( 6,:,:,:) = 1.d0/99.d0
    main( 7,:,:,:) = 0.d0
    main( 8,:,:,:) = 0.d0
    main( 9,:,:,:) = 0.d0
    main(10,:,:,:) = 1.d0/99.d0
    main(11,:,:,:) = 0.d0
    main(12,:,:,:) = 0.d0
    main(13,:,:,:) = 0.d0
    main(14,:,:,:) = 1.d0
    main(15,:,:,:) = .25d0*main(3,:,:,:)
    main(16,:,:,:) = 0.d0
    main(17,:,:,:) = 0.d0
    main(18,:,:,:) = 0.d0
    do j = 2, size(main,3)-1
       do i = 2, size(main,2)-1
          main(19,i,j,2) = 1.d0/100.d0*(i-1.5d0)
          main(20,i,j,2) = 0.d0
          main(21,i,j,2) = Jacobian(main(6:14,i,j,2))
       end do
    end do
    main(1,:,size(main,3)/2+1:size(main,3),:) = .25d0
    main(2,:,size(main,3)/2+1:size(main,3),:) = .5d0
    main(3,:,size(main,3)/2+1:size(main,3),:) = 7.d0*&
         sqrt(1.4d0*.25d0/.5d0)
    main(15,:,:,:) = .25d0*main(3,:,:,:)
  end subroutine GodRieInit
end module Godunov_tester

program Godunov_tester_program
  use Godunov_tester
  use Riemann_tester
  integer :: result, junk
  result = RieTester()
  write(*,*) "Riemann_tester_program: "
  junk = RieErrorReader(result)
  result = GodTester()
  write(*,*) "Godunov_tester_program: "
  junk = GodErrorReader(result)

end program Godunov_tester_program

!!$module grid_velocity
!!$contains
!!$  subroutine grid_velocity_update(main)
!!$    implicit none
!!$    
!!$
!!$  end subroutine grid_velocity_update
!!$end module grid_velocity
