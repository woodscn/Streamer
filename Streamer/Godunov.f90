module Godunov
  use GeneralUtilities
  implicit none
!!$  integer, parameter :: SGL = 4
!!$  integer, parameter :: DBL = 8
   real(8), parameter :: EPS = 5.d-15
   real(8), parameter :: max_dt = 1.d0
!!$  integer, parameter :: eqn_of_state = 1 ! 1 = Ideal Gas Law
!!$  integer, parameter :: prim_length = 21 ! Number of elements in primitive variables vector
!!$  integer, parameter :: cons_length =  5 ! Number of elements in conservative variables vector
  real(8), parameter :: gamma_const = 1.4d0
  real(8), parameter :: gamma1 = 1.d0/(gamma_const-1.d0)
  real(8), parameter :: gamma2 = (gamma_const-1.d0)
  real(8), parameter :: gamma3 = (gamma_const - 1.d0)/(2.d0*gamma_const)
  real(8), parameter :: gamma4 = 1.d0/gamma3
  real(8), parameter :: gamma5 = (gamma_const-1.d0)/(gamma_const+1.d0)
  real(8), parameter :: gamma6 = 1.d0/(gamma_const+1.d0)
  real(8), parameter :: gamma7 = 1.d0/gamma3
!!$  integer :: nx, ny, nz
  real(8), parameter :: dxi   = 1.d0
  real(8), parameter :: deta  = 1.d0
  real(8), parameter :: dzeta = 1.d0
  real(8), parameter :: dxi_inv   = 1.d0/dxi
  real(8), parameter :: deta_inv  = 1.d0/deta
  real(8), parameter :: dzeta_inv = 1.d0/dzeta
  real(8), parameter :: dV_inv = dxi_inv*deta_inv*dzeta_inv
contains
!  integer elemental function gt0(x)
!    implicit none
!    real(8), intent(in) :: x
!    gt0 = ishft( int(sign(1.d0,x) + 1) , -1 )
!  end function gt0

  function riemann_solve( left, right, max_wave_speed, verbose_flag, t_out )
    implicit none
    real(8), dimension(:), intent(in) :: left, right
    real(8), dimension(5) :: riemann_solve
    real(8), intent(out) :: max_wave_speed
    logical, intent(in), optional :: verbose_flag
    real(8), intent(in), optional :: t_out
    real(8) :: tout = 1.d0
    logical :: verbose = .false.
    integer :: n
    integer, parameter :: nx = 1000

    real(8) :: DL, PL, UL, VL, AL
    real(8) :: DR, PR, UR, VR, AR
    real(8) :: Pstar, Ustar, DstarL, DstarR
    real(8) :: tol = 1.d-10, x
    real(8) :: PsiL, PsiR, temp=1.d0, fL, fR, dfL, dfR
    real(8), dimension(nx,5) :: data

    if(present(verbose_flag))verbose=verbose_flag
    if(present(t_out)) tout = t_out
!!$    left = tempL ; right = tempR
    DL =  left(2) ; PL =  left(1) ; UL =  left(3) ; VL =  left(4)
    DR = right(2) ; PR = right(1) ; UR = right(3) ; VR = right(4)

    AL = sqrt(gamma_const*PL/DL) ; AR = sqrt(gamma_const*PR/DR)
    Pstar = guessp(left,right)
    if(verbose)write(*,*) "Initial guess P = " , Pstar
    temp = 1.d0 ; n = 0
    do while(abs(temp) .gt. tol)
       if(verbose)write(*,*) n , "Pstar = " , Pstar; 
       n = n + 1
       if(n .gt. 10)then ;  write(*,*) "Failed Convergence" ; stop ; end if
       PsiL = Pstar/PL
       call u_fun( left,Pstar,fL,dfL)
       PsiR = Pstar/PR
       call u_fun(right,Pstar,fR,dfR)
       temp = ( UR - UL + fR + fL )/( dfL + dfR )
       Pstar = max( Pstar - temp , tol )
       if(verbose)write(*,*) "fL , fR , Pstar = " , fL , fR , Pstar
    end do
    Ustar = .5*(UR+fR+UL-fL)
    DstarL = beta(Pstar/PL)*DL
    DstarR = beta(Pstar/PR)*DR
    if(PsiL<1.d0)then
       max_wave_speed = abs(left(3) - .5d0*(UL+UR) - aL)
    else
       max_wave_speed = abs(left(3)-.5d0*(UL+UR)-aL*sqrt((gamma_const+1.d0)/&
            (2.d0*gamma_const)*(PsiL-1.d0)+1.d0))
    end if
    if(PsiR<1.d0)then
       max_wave_speed = max(abs(right(3) - .5d0*(UL+UR) + aR),max_wave_speed)
    else
       max_wave_speed = max(abs(right(3)-.5*(UL+UR)+aR*sqrt((gamma_const+1.d0)/&
            (2.d0*gamma_const)*(PsiR-1.d0)+1.d0)),max_wave_speed)
    end if
    if(verbose)write(*,*) Ustar , DstarL , DstarR
    call sample(0.D0,left,right,Pstar,Ustar,DstarL,DstarR,riemann_solve,verbose)
    if(verbose)write(*,*)"RiemannSolve = ", riemann_solve
  end function riemann_solve

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
    if(.not.( guessp .gt. min(left(1),right(1)) .and. guessp .lt. max(left(1),right(1)) &
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
          gL=sqrt( 2./((gamma_const+1.)* left(2))/(guessp+ left(1)*(gamma_const-1.)/(gamma_const+1.)))
          gR=sqrt( 2./((gamma_const+1.)*right(2))/(guessp+right(1)*(gamma_const-1.)/(gamma_const+1.)))
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

  subroutine sample(x,left,right,Pstar,Ustar,DstarL,DstarR,out,verbose)
    implicit none
    real(8), dimension(:), intent(in) :: left, right
    real(8), intent(in) :: Pstar, Ustar, DstarL, DstarR, x
    real(8), dimension(:) , intent(out) :: out
    logical, intent(in) :: verbose
    real(8) :: PsiL , SL , SR , DeltaL, DeltaR
    real(8) :: PsiR , aL , aR , betaL , betaR
    real(8) :: cL , cR , cLT , cLH , cRT , cRH , h
    logical :: test_flag = .false.
    real(8) :: Uavg
    if(verbose) test_flag = .true.
    Uavg = .5d0*(left(15)+right(15)) ! The average of the normal grid velocity
    PsiL = Pstar/left(1)
    PsiR = Pstar/right(1)
    aL   = sqrt(gamma_const* left(1)/ left(2))
    aR   = sqrt(gamma_const*right(1)/right(2))
    betaL= DstarL/ left(2)
    betaR= DstarR/right(2)
    if( Ustar .gt. x )then
       if(test_flag) write(*,*) "The boundary lies to the left of the contact wave"
       out(4) = left(4) ; out(5) = left(5)
       if( PsiL .gt. 1.d0 )then
          if(test_flag) write(*,*) " Left shock"
          cL = left(3)-Uavg-aL*sqrt((gamma_const+1.d0)/(2.d0*gamma_const)*(PsiL-1.d0)+1.d0)
          if(test_flag)cL = cL*SL/DeltaL
          if(test_flag) write(*,*) "Left shock speed = " , cL
          if( cL .gt. x )then
             if(test_flag) write(*,*) " The boundary lies to the left of the shock"
             out(1) = left(1)
             out(3) = left(3)
             out(2) = left(2)
          else
             if(test_flag) write(*,*) " The boundary lies in the left central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarL
          end if
       else
          if(test_flag) write(*,*) "Left rarefaction wave"
          cLT = left(3) - Uavg - aL
          if(test_flag) cLT = cLT*SL/DeltaL
          if(test_flag) write(*,*) "Left rarefaction tail speed = " , cLT
          cLH = Ustar-Uavg - sqrt(gamma_const*Pstar/DstarL)
          if(test_flag) write(*,*) "Left rarefaction head speed = " , cLH
          if( cLT .gt. x )then
             if(test_flag) write(*,*) " The boundary lies to the left of the wave"
             out(1) = left(1)
             out(3) = left(3)
             out(2) = left(2)
          elseif( cLH .lt. x )then
             if(test_flag) write(*,*) "The boundary lies in the left central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarL
          else
             if(test_flag) write(*,*) "The boundary lies within the left expansion wave"
             out(1) = left(1)*( 2.d0*gamma6 + gamma5/aL*(left(3)-Uavg) )**gamma7
             out(3) = left(3) - 2.d0*aL/(gamma_const-1.d0)*((out(1)/left(1))&
                  **((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = left(2)*(out(1)/left(1))**(1.d0/gamma_const)
          end if
       end if
    else
       if(test_flag) write(*,*) " The boundary lies to the right of the contact wave"
       out(4) = right(4) ; out(5) = right(5)
       if( PsiR .gt. 1. )then
          if(test_flag) write(*,*) " Right shock"
          cR = right(3)-Uavg+aR*sqrt((gamma_const+1.d0)/(2.d0*gamma_const)*(PsiR-1.d0)+1.d0)
          if(test_flag) cR = cR*SR/DeltaR
          if(test_flag) write(*,*) " Right shock speed = " , cR
          if( cR .lt. x )then
             if(test_flag) write(*,*) " The boundary lies to the right of the shock"
             out(1) = right(1)
             out(3) = right(3)
             out(2) = right(2)
          else
             if(test_flag) write(*,*) " The boundary lies in the right central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarR
          end if
       else
          if(test_flag) write(*,*) " Right rarefaction wave"
          cRT = right(3) - Uavg + aR
          if(test_flag) cRT = cRT*SR/DeltaR
          if(test_flag) write(*,*) "Right rarefaction tail speed = " , cRT
          cRH = Ustar - Uavg + sqrt(gamma_const*Pstar/DstarR)
          if(test_flag) write(*,*) "Right rarefaction head speed = " , cRH
          if( cRT .lt. x )then
             if(test_flag) write(*,*) " The boundary lies to the right of the wave"
             out(1) = right(1)
             out(3) = right(3)
             out(2) = right(2)
          elseif( cRH .gt. x )then
             if(test_flag) write(*,*) "The boundary lies in the right central region"
             out(1) = Pstar
             out(3) = Ustar
             out(2) = DstarR
          else
             if(test_flag) write(*,*) " The boundary lies within the right expansion wave"
             out(1) = right(1)*( 2.d0*gamma6 - gamma5/aR*(right(3)-Uavg) )**gamma7
             out(3) = right(3) + 2.d0*aR/(gamma_const-1.d0)*((out(1)/right(1))**&
                  ((gamma_const-1.d0)/(2.d0*gamma_const))-1.d0)
             out(2) = right(2)*(out(1)/right(1))**(1.d0/gamma_const)
          end if
       end if
    end if
  end subroutine sample

  ! Computes specific energy, given the array [ p, rho, u, v, w ]
  real(8) function energy_func(in)
    implicit none
    real(8), dimension(:), intent(in) :: in
    energy_func = 0.5d0*(in(3)**2+in(4)**2+in(5)**2) + in(1)/(in(2)*gamma2)
  end function energy_func

  subroutine primtocons(main,nx,ny,nz)
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
    integer :: nx, ny, nz
    real(8), dimension(nx,ny,nz) :: energy
    energy = main(21,:,:,:)*( main(1,:,:,:)*gamma1 + 0.5d0*main(2,:,:,:)&
         *( main(3,:,:,:)**2 + main(4,:,:,:)**2 + main(5,:,:,:)**2 ) )
    main(1,:,:,:) = main(2,:,:,:)*main(21,:,:,:)
    main(2,:,:,:) = main(1,:,:,:)*main(3,:,:,:)
    main(3,:,:,:) = main(1,:,:,:)*main(4,:,:,:)
    main(4,:,:,:) = main(1,:,:,:)*main(5,:,:,:)
    main(5,:,:,:) = energy
  end subroutine primtocons
  
  subroutine constoprim(main,nx,ny,nz)
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
    integer, intent(in) :: nx, ny, nz
    real(8), dimension(nx,ny,nz) :: temp1, temp2, p
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
  end subroutine constoprim

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
  
  subroutine vels_transform(in, matrix)
    real(8), dimension(21), intent(inout) :: in
    real(8), dimension(3,3), intent(in) :: matrix
    in(3:5) = matmul(matrix,in(3:5))
    in(15:17) = matmul(matrix,in(15:17))
  end subroutine vels_transform
    
  subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)
! Advance the solution using the integral form of the equations
! This subroutine assumes that main is the full array of primitive variables. main 
! must also include the boundary cell values. That is, main contains a 0-index and 
! an n + 1 index containing the contents prescribed by the boundary conditions.
    implicit none
    real(8), dimension(21,nx+2*bcextent,ny+2*bcextent,nz+2*bcextent), intent(inout) :: main
!f2py intent(in,out) :: main
    real(8), intent(in), optional :: dt_in
    real(8), intent(in), optional :: CFL
    integer, intent(in) :: nx,ny,nz,bcextent
    integer :: i, j, k
    real(8), dimension(14,3,nx+1,ny+1,nz+1) :: fluxes
    real(8), dimension(14,nx,ny,nz) :: temp
    real(8) :: max_wave_speed, dt, max_dt_grid
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8), dimension(5) :: junk
    fluxes = 0.d0
    do k = 1, nz + 1
       do j = 1, ny + 1
          do i = 1, nx + 1
             fluxes( 1: 5,1,i,j,k) = flux(riemann_solve(main(:,i,j+1,k+1),&
                  main(:,i+1,j+1,k+1),max_wave_speed),&
                  .5d0*(main(:,i,j+1,k+1)+main(:,i+1,j+1,k+1)),1)
             fluxes( 1: 5,2,i,j,k) = flux(riemann_solve(main(:,i+1,j,k+1),&
                  main(:,i+1,j+1,k+1),max_wave_speed),&
                  .5d0*(main(:,i+1,j,k+1)+main(:,i+1,j+1,k+1)),2)
             fluxes( 1: 5,3,i,j,k) = flux(riemann_solve(main(:,i+1,j+1,k),&
                  main(:,i+1,j+1,k+1),max_wave_speed),&
                  .5d0*(main(:,i+1,j+1,k)+main(:,i+1,j+1,k+1)),3)
             fluxes( 6: 8,1,i,j,k) = 0.5d0*(main(15:17,i+1,j+1,k+1)+&
                  main(15:17,i,j+1,k+1))
             fluxes( 9:11,2,i,j,k) = 0.5d0*(main(15:17,i+1,j+1,k+1)+&
                  main(15:17,i+1,j,k+1))
             fluxes(12:14,3,i,j,k) = 0.5d0*(main(15:17,i+1,j+1,k+1)+&
                  main(15:17,i+1,j+1,k))
          end do
       end do
    end do
    max_wave_speed = max(max_wave_speed,EPS)
    if(dt_in > 0.d0)then
       dt = dt_in
    elseif(CFL>0.d0)then
       dt = min(CFL/max_wave_speed,max_dt)
       dt = min(maxval((main(6,2,2:ny+1,2:nz+1)-main(18,2,2:ny+1,2:nz+1))&
            /main(15,2,2:ny+1,2:nz+1)),dt)
    else
       write(*,*) "Error in Godunov prim_update, no time step information given"
       read(*,*)
       stop
    end if
!!$    write(*,*) "Dt = ", dt

!!$    write(*,*) "I need to implement a system to keep the grid from advancing too far in a single step!"
!!$    stop
    call primtocons(main(:,2:nx+1,2:ny+1,2:nz+1),nx,ny,nz)
    main(1:14,2:nx+1,2:ny+1,2:nz+1) = main(1:14,2:nx+1,2:ny+1,2:nz+1) - (&
         (fluxes(:,1,2:nx+1,:,:)-fluxes(:,1,1:nx,:,:))*deta*dzeta + &
         (fluxes(:,2,:,2:ny+1,:)-fluxes(:,2,:,1:ny,:))*dxi*dzeta + &
         (fluxes(:,3,:,:,2:nz+1)-fluxes(:,3,:,:,1:nz))*dxi*deta &
         )*dt*dV_inv
    call constoprim(main(:,2:nx+1,2:ny+1,2:nz+1),nx,ny,nz)
! Update grid velocity
!!$    call grid_velocity_update(main)
!!$    main(15:17,:,:,:) = main(3:5,:,:,:)*.25d0
! Update grid position
    main(18:20,2:nx+1,2:ny+1,2:nz+1) = main(18:20,2:nx+1,2:ny+1,2:nz+1) + &
         main(15:17,2:nx+1,2:ny+1,2:nz+1)*dt
! Update extra variables
  do k = 2, nz+1
     do j = 2, ny+1
        do i = 2, nx+1
           main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
        end do
     end do
  end do
  end subroutine prim_update

  function flux(in,geom_avg,case_no)
    implicit none
    real(8), dimension(:), intent(in) :: in, geom_avg
    integer, intent(in) :: case_no
    real(8), dimension(5) :: flux
    real(8) :: D, grad(3), J, Jinv
!!$    write(*,*) "Got this far"
    J = geom_avg(21)!Jacobian(geom_avg(6:14))
    Jinv = 1.d0/geom_avg(21)
    select case(case_no)
    case(1)
       grad = [&
            geom_avg(10)*geom_avg(14) - geom_avg(11)*geom_avg(13) , &
            geom_avg(11)*geom_avg(12) - geom_avg( 9)*geom_avg(14) , &
            geom_avg( 9)*geom_avg(13) - geom_avg(10)*geom_avg(12) ]*Jinv
    case(2)
       grad = [&
            geom_avg( 8)*geom_avg(13) - geom_avg( 7)*geom_avg(14) , &
            geom_avg( 6)*geom_avg(14) - geom_avg( 8)*geom_avg(12) , &
            geom_avg( 7)*geom_avg(12) - geom_avg( 6)*geom_avg(13) ]*Jinv
    case(3)
       grad = [&
            geom_avg( 7)*geom_avg(11) - geom_avg( 8)*geom_avg(10) , &
            geom_avg( 8)*geom_avg( 9) - geom_avg( 6)*geom_avg(11) , &
            geom_avg( 6)*geom_avg(10) - geom_avg( 7)*geom_avg( 9) ]*Jinv
!!$    case default
!!$       write(*,*) "Error in flux -- invalid case_no: ",case_no
    end select
    D = sum(([ in(3), in(4), in(5) ] - [ geom_avg(15), geom_avg(16), geom_avg(17) ])*grad)
    flux = [&
         in(2)*J*D, &
         in(2)*J*D*in(3) + J*grad(1)*in(1), &
         in(2)*J*D*in(4) + J*grad(2)*in(1), &
         in(2)*J*D*in(5) + J*grad(3)*in(1), &
         in(2)*J*D*energy_func(in(1:5)) + J*sum(grad*in(3:5))*in(1) ]
  end function flux

end module Godunov
!!$module grid_velocity
!!$contains
!!$  subroutine grid_velocity_update(main)
!!$    implicit none
!!$    
!!$
!!$  end subroutine grid_velocity_update
!!$end module grid_velocity
