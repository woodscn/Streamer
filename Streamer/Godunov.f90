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

  function riemann_solve( left, right, geom_avg, max_wave_speed, verbose_flag, t_out )
    implicit none
    real(8), dimension(5), intent(in) :: left, right
    real(8), dimension(21), intent(in) :: geom_avg
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
    call sample(0.D0,left,right,geom_avg,Pstar,Ustar,DstarL,DstarR,riemann_solve,verbose)
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

  subroutine sample(x,left,right,geom_avg,Pstar,Ustar,DstarL,DstarR,out,verbose)    
    implicit none
    real(8), dimension(:), intent(in) :: left, right, geom_avg
    real(8), intent(in) :: Pstar, Ustar, DstarL, DstarR, x
    real(8), dimension(:) , intent(out) :: out
    logical, intent(in) :: verbose
    real(8) :: PsiL , SL , SR , DeltaL, DeltaR
    real(8) :: PsiR , aL , aR , betaL , betaR
    real(8) :: cL , cR , cLT , cLH , cRT , cRH , h
    logical :: test_flag = .false.
    real(8) :: Uavg
    if(verbose) test_flag = .true.
    Uavg = geom_avg(15) ! The average of the normal grid velocity
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
  
  subroutine prim_update(main,bcextent,dt_in,CFL,nx,ny,nz)
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
             call compute_fluxes(main(1:5,i-1,j,k), main(1:5,i,j,k),&
                  .5d0*(main(:,i-1,j,k)+main(:,i,j,k)),fluxx(:,i,j,k),&
                  1,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
          end do
       end do
    end do
    do k = 0, nz-1
       do j = 0, ny
          do i = 0, nx-1
             call compute_fluxes(main(1:5,i,j-1,k), main(1:5,i,j,k),&
                  .5d0*(main(:,i,j-1,k)+main(:,i,j,k)),fluxy(:,i,j,k),&
                  2,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta],debug_flag=.true.)
          end do
       end do
    end do
    do k = 0, nz
       do j = 0, ny-1
          do i = 0, nx-1
             call compute_fluxes(main(1:5,i,j,k-1), main(1:5,i,j,k),&
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
!!$    write(*,*) "Dt = ", dt

!!$    write(*,*) "I need to implement a system to keep the grid from advancing too far in a single step!"
!!$    stop
!!$    write(*,*) "Compute fluxes inputs:"
!!$    write(*,*) "Left state = "
!!$    write(*,*) main(1:5,0,0,0)
!!$    write(*,*) "Center state = "
!!$    write(*,*) main(1:5,0,1,0)
!!$    write(*,*) "Right state = "
!!$    write(*,*) main(1:5,0,2,0)
!!$    write(*,*) "Left interface = "
!!$  subroutine compute_fluxes(inL, inR, geom_avg, flux_vec, case_no,&
!!$       max_wave_speed,dt,dV_in,debug_flag)
!!$    write(*,*) 
    call compute_fluxes(main(1:5,0,0,0),main(1:5,0,1,0),&
         .5d0*(main(:,0,0,0)+main(:,0,1,0)), junk, 2, max_wave_speed,&
         dt,[dxi,deta,dzeta],debug_flag=.false.)
!!$    write(*,*) "Right interface = "
    call compute_fluxes(main(1:5,0,1,0),main(1:5,0,2,0),&
         .5d0*(main(:,0,1,0)+main(:,0,2,0)), junk, 2, max_wave_speed,&
         dt,[dxi,deta,dzeta],debug_flag=.false.)
!!$    
!!$    write(*,*) "Geometry = "
!!$    write(*,*) .5d0*(main(6:14,1,2,1)+main(6:14,1,1,1))
!!$    write(*,*) "dV = ", [dxi,deta,dzeta]
!!$    write(*,*) "Grid motion = "
!!$    write(*,*) main(15:17,0,0,0)
!!$    write(*,*) main(15:17,0,1,0)
!!$    write(*,*) main(15:17,0,2,0)
!!$    write(*,*) "Left flux = "
!!$    write(*,*) fluxy(1:5,0,1,0)
!!$    write(*,*) "Right flux = "
!!$    write(*,*) fluxy(1:5,0,2,0)
!!$    write(*,*) "Total change in conservative variables (y):"
!!$    write(*,*) "( Mass, Momentum(x), Momentum(y), Momentum(z), Energy)"
!!$    write(*,*) - ( &
!!$         (fluxy(1:5,0,2,0)-fluxy(1:5,0,1,0)) &
!!$         )
!!$    read(*,*)
!!!$             call compute_fluxes(main(1:5,i,j,k-1), main(1:5,i,j,k),&
!!!$                  .5d0*(main(:,i,j,k-1)+main(:,i,j,k)),fluxz(:,i,j,k),&
!!!$                  3,max_wave_speed,dt=dt,dV_in=[dxi,deta,dzeta])
    call primtocons(main(:,0:nx-1,0:ny-1,0:nz-1),nx,ny,nz)
    main(1:14,0:nx-1,0:ny-1,0:nz-1) = main(1:14,0:nx-1,0:ny-1,0:nz-1) - (&
         (fluxx(:,1:nx,:,:)-fluxx(:,0:nx-1,:,:))*deta*dzeta + &
         (fluxy(:,:,1:ny,:)-fluxy(:,:,0:ny-1,:))*dxi*dzeta + &
         (fluxz(:,:,:,1:nz)-fluxz(:,:,:,0:nz-1))*dxi*deta &
         )*dt*dV_inv

    call constoprim(main(:,0:nx-1,0:ny-1,0:nz-1),nx,ny,nz)
! Update grid velocity
!!$    call grid_velocity_update(main)
    main(15:17,:,:,:) = main(3:5,:,:,:)*.25d0
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
  end subroutine prim_update

  subroutine compute_fluxes(inL, inR, geom_avg, flux_vec, case_no,&
       max_wave_speed,dt,dV_in,debug_flag)
    implicit none
    real(8), dimension(:), intent(in) :: inL, inR
    real(8), dimension(5) :: StateL, StateR
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
    select case(case_no)
    case(1)
       row_ops_mat = reshape([&
            1,0,0,&
            0,1,0,&
            0,0,1],&
            [3,3])
       flux_vec(6:8) = -geom_avg(15:17)
       dA = dV_in(2)*dV_in(3)
    case(2)
       row_ops_mat = reshape([&
            0,1,0,&
            0,0,1,&
            1,0,0],&
            [3,3])
       flux_vec(9:11) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(3)
    case(3)
       row_ops_mat = reshape([&
            0,0,1,&
            1,0,0,&
            0,1,0],&
            [3,3])
       flux_vec(12:14) = -geom_avg(15:17)
       dA = dV_in(1)*dV_in(2)
    case default
       write(*,*) "Invalid case_no in compute_fluxes -- case_no = ",case_no
       stop
    end select
    dXi_dX_u = matmul(dXi_dX_u,transpose(row_ops_mat))
    dX_dXi_u = matmul(row_ops_mat,dX_dXi_u)
    if(maxval(matmul(dXi_dX_u,dX_dXi_u)-reshape([1,0,0,0,1,0,0,0,1],[3,3]))>1d-14)then
       write(*,*) "Inverses not working!"
       write(*,*) 
       write(*,*) dXi_dX_u
       write(*,*) dX_dXi_u
       read(*,*)
    end if
    StateL(3:5) = matmul(dXi_dX_u,StateL(3:5))
    StateR(3:5) = matmul(dXi_dX_u,StateR(3:5))
    geom_avg(15:17) = matmul(dXi_dX_u,geom_avg(15:17))
    interface_vars = riemann_solve(StateL,StateR,geom_avg,temp_wave_speed)
    interface_vars(3:5) = matmul(dX_dXi_u,interface_vars(3:5))
!!$    write(*,*) "Interface = "
!!$    write(*,*) interface_vars
    geom_avg(15:17) = matmul(dX_dXi_u,geom_avg(15:17))
!!$    geom_avg(6:14) = geom_avg(6:14)-flux_vec(6:14)*dt*dV_inv
!!$    flux_vec(6:14) = 0.d0
    flux_vec(1:5) = flux(interface_vars,geom_avg,case_no)
    max_wave_speed = max(max_wave_speed,temp_wave_speed)
  end subroutine compute_fluxes

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

end module Godunov

program Godunov_tester
  use Godunov
  implicit none
  real(8), dimension(21,3,4,3) :: main
  real(8), dimension(3) :: gradxi, gradeta, gradzeta
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
  call prim_update(main,1,.001d0,.25d0,size(main,2)-2,size(main,3)-2,&
       size(main,4)-2)

end program Godunov_tester
!!$module grid_velocity
!!$contains
!!$  subroutine grid_velocity_update(main)
!!$    implicit none
!!$    
!!$
!!$  end subroutine grid_velocity_update
!!$end module grid_velocity
