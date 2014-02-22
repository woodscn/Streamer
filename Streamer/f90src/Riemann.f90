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
    DL =  left(2); PL =  left(1);
    DR = right(2); PR = right(1);
!!$    select case(dir)
!!$    case(1)
       UL =  left(3); VL =  left(4); WL =  left(5)
       UR = right(3); VR = right(4); WR = right(5)
       Uavg = .5d0*(left(15)+right(15))
!!$    case(2)
!!$       UL =  left(4); VL =  left(5); WL =  left(3)
!!$       UR = right(4); VR = right(5); WR = right(3)
!!$       Uavg = .5d0*(left(16)+right(16))
!!$    case(3)
!!$       UL =  left(5); VL =  left(3); WL =  left(4)
!!$       UR = right(5); VR = right(3); WR = right(4)
!!$       Uavg = .5d0*(left(17)+right(17))
!!$    end select
    AL = sqrt(gamma_const*PL/DL); AR = sqrt(gamma_const*PR/DR)

    met_inv = MetricInverse(.5d0*(left(6:14)+right(6:14)));
    J = Jacobian(.5d0*(left(6:14)+right(6:14)))
!!$    if(abs(J-1d-6)>1d-7)then
!!$       write(*,*) "J = ",J,"  Metric = ",.5d0*(left(6:14)+right(6:14))
!!$       read(*,*)
!!$    end if
    S = sqrt(sum((J*met_inv(dir:9:3))**2))
!    write(*,*) "S/J = ", S/J
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

  function wave_speeds(left,right,dir,riemann_middle_states,Uavg,J,S)
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
            *sqrt(&
            (gamma_const+1d0)/(2d0*gamma_const)*(Pstar/right(1)-1d0)+1d0))
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
       f = 2.d0*a0/(gamma_const-1.d0)*(psi**((gamma_const-1.d0)/&
            (2.d0*gamma_const))-1.d0)
       df= 1.d0/(in(2)*a0)*psi**((gamma_const+1.d0)/(-2.d0*gamma_const))
    end if
  end subroutine u_fun

  pure function beta(psi)
    implicit none
    real(8), intent(in) :: psi
    real(8) :: beta
    if( psi .gt. 1 )then
       beta = ((gamma_const+1.)*psi+gamma_const-1.)/&
            (gamma_const+1.+(gamma_const-1.)*psi)
    else
       beta = psi**(1.d0/gamma_const)
    end if
  end function beta

  pure subroutine sample(x,left,right,riemann_middle_states,&
       riemann_wave_speeds,Uavg,J,S,out)    
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
             out(1) = left(1)*(2d0*gamma6+gamma5/aL*&
                  (left(3)-Uavg-J/S*x))**gamma7
             out(3) = left(3) - 2.d0*aL/(gamma_const-1.d0)*((out(1)/left(1))&
                  **((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = left(2)*(out(1)/left(1))**(1.d0/gamma_const)
          end if
       end if
    else
!!$       write(*,*) " The boundary lies to the right of the contact wave"
       out(4) = right(4); out(5) = right(5)
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
             out(1) = right(1)*(2d0*gamma6-gamma5/&
                  aR*(right(3)-Uavg-J/S*x))**gamma7
             out(3) = right(3) + 2.d0*aR/(gamma_const-1.d0)*((out(1)/right(1))&
                  **((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = right(2)*(out(1)/right(1))**(1.d0/gamma_const)
          end if
       end if
    end if
  end subroutine sample
end module Riemann

