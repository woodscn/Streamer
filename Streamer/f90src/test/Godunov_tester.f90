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
  use GeneralUtilitiesTester
  use GodunovDriver
  use TimeAdvancementStuff
contains
  real(8) function norm2(in)
    implicit none
    real(8), dimension(:), intent(in) :: in
    norm2 = sqrt(sum(in**2))
  end function norm2
  integer function GodErrorReader(in)
    integer, intent(in) :: in
    write(*,*) 
    select case(in)
    case(0) 
       write(*,*) "   All tests passed"
    case(1)
       write(*,*) "   Primitive-Conservative mutual inverse test failed"
    case(2)
       write(*,*) "   Grid_coords failed to return an orthonormal system"
    case(3)
       write(*,*) "   Flux fails provided test problem"
    case(4)
       write(*,*) "   GodConvergenceTester1D returned unexpectedly low convergence rate"
    case default
       write(*,*) "   Unexpected error code"
    end select
    GodErrorReader = 0
  end function GodErrorReader

  integer function GodConvergenceTester1D(left,right,x0,t_out,dt,nmax,prim_update_options,base_filename)
    implicit none
    real(8), dimension(5), intent(in) :: left, right
    real(8), dimension(5) :: left_init, right_init
    real(8), intent(in) :: x0, t_out, dt
    integer, intent(in) :: nmax
    integer, dimension(:), intent(in) :: prim_update_options
    character(len=*), intent(in), optional :: base_filename
    integer, parameter :: nx = 100
    real(8), dimension(0:10) :: dxes, rmserrors, fvec
    real(8), dimension(0:10,3) :: fjac
    integer :: nxmax
    real(8), dimension(:,:,:,:), allocatable :: Rie_1D_dir
    real(8), dimension(:,:), allocatable :: Rie_1D
    real(8), dimension(:,:), allocatable :: Rie_1D_exact
    integer :: filenum
    integer :: dir, m, n, i, j, k
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    real(8) :: max_wave_speed, t
    real(8), dimension(3) :: fitted_poly, fitted_poly2
    integer :: info
    logical, parameter :: verbose = .true.
    real(8) :: maxvalue
    real(8), dimension(:), allocatable :: exact_x
    real(8), dimension(21) :: left_full, right_full
    integer, dimension(3) :: nxes
    real(8), dimension(3,3) :: row_ops_mat
    real(8) :: dt_out

    GodConvergenceTester1D = 0
    if(nmax>10)then
       write(*,*) "Error in GodConvergenceTester1D: Too many refinements!"
       stop
    end if
    maxvalue = 0d0
    do dir = 1, 1
       row_ops_mat = row_ops_mat_func(dir)
       left_init = left; left_init(3:5) = matmul(row_ops_mat,left(3:5))
       right_init = right; right_init(3:5) = matmul(row_ops_mat,right(3:5))
       do n = 0, nmax
          nxmax = nx*2**n
          allocate(Rie_1D_exact(5,0:nxmax-1),exact_x(0:nxmax-1),&
               Rie_1D(21,-1:nxmax))
          select case(dir)
          case(1)
             allocate(Rie_1D_dir(21,-1:nxmax,-1:1,-1:1))
             nxes = [nx*2**n,1,1]
          case(2)
             allocate(Rie_1D_dir(21,-1:1,-1:nxmax,-1:1))
             nxes = [1,nx*2**n,1]
          case(3)
             allocate(Rie_1D_dir(21,-1:1,-1:1,-1:nxmax))
             nxes = [1,1,nx*2**n]
          end select
          call RieInit1D(left_init,right_init,x0,nxmax,dir,[0d0,1d0],Rie_1D_dir)
          t = 0d0
          select case(dir)
          case(1) 
             exact_x = (Rie_1D_dir(18,0:nxmax-1,0,0)-x0)&
                  /Rie_1D_dir(6,0:nxmax-1,0,0)/t_out
          case(2)
             exact_x = (Rie_1D_dir(19,0,0:nxmax-1,0)-x0)&
                  /Rie_1D_dir(10,0,0:nxmax-1,0)/t_out
          case(3)
             exact_x = (Rie_1D_dir(20,0,0,0:nxmax-1)-x0)&
                  /Rie_1D_dir(14,0,0,0:nxmax-1)/t_out
          end select
!!$          call write_files_matlab(Rie_1D_dir,t,nxes(1),nxes(2),.true.)
          do 
             call FreeExitConditions(Rie_1D_dir,size(Rie_1D_dir,2),&
                  size(Rie_1D_dir,3),size(Rie_1D_dir,4))
             call prim_update(Rie_1D_dir,dt_out,dt,&
                  .7d0,nxes(1),nxes(2),nxes(3),prim_update_options)
             t = t + dt_out
!!$             call write_files_matlab(Rie_1D_dir(:,:,:,1),t,nxes(1),nxes(2),.false.)
             if(t .ge. t_out - 1d-13) exit
          end do
          select case(dir)
          case(1)
             Rie_1D = Rie_1D_dir(:,:,0,0)
             dxes(n) = Rie_1D(18,2) - Rie_1D(18,1)
             exact_x = ((Rie_1D(18,0:nxmax-1)-x0)/t_out&
                  -.5d0*(Rie_1D(15,0)+Rie_1D(15,nxmax-1)))&
                  /(.5d0*(Rie_1D(6,0)+Rie_1D(6,nxmax-1)))
          case(2)             
             Rie_1D(:,:) = Rie_1D_dir(:,0,:,0)
             Rie_1D(3,:) = Rie_1D(4,:); Rie_1D(4,:) = 0d0;
             dxes(n) = Rie_1D(19,2) - Rie_1D(19,1)
             exact_x = ((Rie_1D(19,0:nxmax-1)-x0)/t_out&
                  -.5d0*(Rie_1D(16,0)+Rie_1D(16,nxmax-1)))/Rie_1D(10,0:nxmax-1)
          case(3)
             Rie_1D(:,:) = Rie_1D_dir(:,0,0,:)
             Rie_1D(3,:) = Rie_1D(5,:); Rie_1D(5,:) = 0d0;
             dxes(n) = Rie_1D(20,2) - Rie_1D(20,1)
             exact_x = ((Rie_1D(20,0:nxmax-1)-x0)/t_out&
                  -.5d0*(Rie_1D(17,0)+Rie_1D(17,nxmax-1)))/Rie_1D(14,0:nxmax-1)
          end select
          left_full = Rie_1D(:,-1)
          right_full = Rie_1D(:,nxmax)
          Rie_1D_exact = 0d0
!!$          write(*,*) " Left_full = ",left_full
!!$          write(*,*) "Right_full = ",right_full
          call riemann_solve(left_full,right_full,dir,nxmax,exact_x,&
               Rie_1D_exact,max_wave_speed,riemann_middle_states,riemann_wave_speeds)
!!$          write(*,*) left_full(1:3)
!!$          write(*,*) right_full(1:3)
!!$          write(*,*) riemann_middle_states
!!$          write(*,*) riemann_wave_speeds
!!$          read(*,*)
          rmserrors(n) = sqrt(sum(((Rie_1D(2,0:nxmax-1)-Rie_1D_exact(2,:))/&
               maxval(abs(Rie_1D_exact(2,:))))**2)/nxmax)
          maxvalue = max(maxvalue,maxval(abs(Rie_1D(2,:))))
          if(present(base_filename))then
             filenum = 92920 + n
             open(unit=filenum,file=base_filename//achar(48+n)//"_"//achar(87+dir)//".dat")
             do m = 1, 20
                write(filenum,*) Rie_1D(m,0:nxmax-1)
             end do
             close(filenum)
             open(unit=92919,file=base_filename//achar(48+n)//"_"//achar(87+dir)//"_exact.dat")
             do m = 1, 5
                write(92919,*) Rie_1D_exact(m,:)
             end do
             close(92919)
          end if
          deallocate(Rie_1D,Rie_1D_exact,exact_x,Rie_1D_dir)
       end do
    fitted_poly(1:2) = polyfit(log(dxes(0:nmax)),log(rmserrors(0:nmax)),1)
    fitted_poly(3) = 0d0
!!$    call minpack_function_fitting(dxes(0:nmax),rmserrors(0:nmax),&
!!$         exponential_with_y_offset,fitted_poly,fvec(0:nmax),fjac(0:nmax,:),1d-3,info)
    call ToroRiemannConvergence(fitted_poly2)
    info = 1
    if(info .eq. 1 .or. info .eq. 2 .or. info .eq. 3)then
       if(verbose)then
          write(*,*) "Simulation converges with order = ",fitted_poly(2)
          write(*,*) "Simulation converges to within ",&
               fitted_poly(3)/maxvalue*100,&
               "% of exact solution."
          write(*,*) "As a comparison, Numerica's hyper_eul converges with ",&
               "order = ", fitted_poly2(2)
          write(*,*) "CLAWPACK converges with order = 0.608111"
          write(*,*) "CLAWPACK converges to within 0.49% of exact solution"
       end if
    else
       write(*,*) "Function fitting returned unexpected value"
       write(*,*) "info = ", info
       write(*,*) "The meanings of various info values are found in lmder1.f"
       stop
    end if
    if(fitted_poly(2) < .5d0 .or. &
         fitted_poly(3)/maxvalue > .5d0)&
         GodConvergenceTester1D = 4
    end do
  end function GodConvergenceTester1D

  function RotateCoords(in,phi,theta)
    implicit none
    real(8), dimension(3), intent(in) :: in
    real(8), intent(in) :: phi, theta
    optional theta
    real(8), dimension(3) :: RotateCoords
    
    RotateCoords = matmul(reshape(&
         [cos(phi),sin(phi),0d0,&
         -sin(phi),cos(phi),0d0,&
         0d0,0d0,1d0],[3,3]),in)
    if(present(theta))then
       write(*,*) "Error in RotateCoords, 3-D rotations not yet supported!!!"
       stop
    end if
  end function RotateCoords
  
!!$  integer function GodConvergenceTester2D(left,right,t_out,dt,nmax,base_filename)
!!$    use GeneralUtilities, only : PI
!!$    use TimeAdvancementStuff
!!$    implicit none
!!$    real(8), dimension(5), intent(in) :: left, right
!!$    real(8), intent(in) :: t_out, dt
!!$    integer, intent(in) :: nmax
!!$    character(len=*), intent(in), optional :: base_filename
!!$    integer, parameter :: nx = 50
!!$    real(8), dimension(0:10) :: dxes, rmserrors, fvec
!!$    real(8), dimension(0:10,3) :: fjac
!!$    integer :: nxmax
!!$    real(8), dimension(:,:,:,:), allocatable :: Rie_2D
!!$    real(8), dimension(:,:,:), allocatable :: Rie_2D_exact
!!$    real(8), parameter :: phi = 0d0!PI*.5d0!25d0/180d0
!!$    real(8), dimension(2) :: xy
!!$    integer :: filenum
!!$    integer :: m, n, i, j, k
!!$    real(8), dimension(4) :: riemann_middle_states
!!$    real(8), dimension(5) :: riemann_wave_speeds
!!$    real(8) :: max_wave_speed, t
!!$    real(8), dimension(3) :: fitted_poly
!!$    integer :: info
!!$    real(8), dimension(3) :: rotated_coords, coords_shift
!!$    real(8) :: dt_out
!!$
!!$    GodConvergenceTester2D = 0
!!$    if(nmax>10)then
!!$       write(*,*) "Error in GodConvergenceTester2D: Too many refinements!"
!!$       stop
!!$    end if
!!$    do n = 0, nmax
!!$       nxmax = nx*2**n
!!$       allocate(Rie_2D(21,-1:nxmax,-1:nxmax,-1:1),&
!!$            Rie_2D_exact(5,0:nxmax-1,0:nxmax-1))
!!$       call RieInit2D(left,right,phi,nxmax,[0d0,1d0],Rie_2D)
!!$       t = 0d0
!!$!       call write_files_matlab(Rie_2D(:,0:nx-1,0:nx-1,0),t,nx,nx,.true.)
!!$       do 
!!$          call prim_update(Rie_2D,dt_out,FreeExitConditions,1,dt,.7d0,nxmax,nxmax,1,[2,0])
!!$          t = t + dt
!!$          write(*,*) "t = ",t
!!$          if(t .ge. t_out) exit
!!$       end do
!!$       Rie_2D_exact = 0d0
!!$       do i = 0, nxmax-1
!!$          do j = 0, nxmax-1
!!$             coords_shift = .5d0*(maxval(maxval(Rie_2D(18:20,:,:,0),3),2)-&
!!$                  minval(minval(Rie_2D(18:20,:,:,0),3),2))
!!$             rotated_coords = RotateCoords((Rie_2D(18:20,i,j,0)-coords_shift)&
!!$                  /Rie_2D(6:14:4,i,j,0),phi)/t_out
!!$             call riemann_solve(Rie_2D(:,-1,-1,-1),Rie_2D(:,nxmax,nxmax,-1),1,1,&
!!$                  rotated_coords,Rie_2D_exact(:,i,j),max_wave_speed)
!!$          end do
!!$       end do
!!$       dxes(n) = Rie_2D(18,1,0,0)-Rie_2D(18,0,0,0)
!!$       rmserrors(n) = sqrt(sum(((Rie_2D(2,0:nxmax-1,0:nxmax-1,0)&
!!$            -Rie_2D_exact(2,:,:))/&
!!$            maxval(abs(Rie_2D_exact(2,:,:))))**2)/nxmax**2)
!!$       if(present(base_filename))then
!!$          filenum = 93020 + n
!!$          open(unit=filenum,file=base_filename//achar(48+n)//".dat")
!!$          do m = 1, 20
!!$             write(filenum,*) Rie_2D(m,:,:,0)
!!$          end do
!!$          close(filenum)
!!$          if(n .eq. nmax)then
!!$             open(unit=93019,file=base_filename//achar(48+n)//"_exact.dat")
!!$             do m = 1, 5
!!$                write(93019,*) Rie_2D_exact(m,:,:)
!!$             end do
!!$             close(93019)
!!$          end if
!!$       end if
!!$       deallocate(Rie_2D,Rie_2D_exact)
!!$    end do
!!$    write(*,*) "dxes = ",dxes(0:nmax)
!!$    write(*,*) "rmserrors = ",rmserrors(0:nmax)
!!$    fitted_poly(1:2) = polyfit(log(dxes(0:nmax)),log(rmserrors(0:nmax)),1)
!!$    fitted_poly(3) = 0d0
!!$    call minpack_function_fitting(dxes(0:nmax),rmserrors(0:nmax),&
!!$         exponential_with_y_offset,fitted_poly,fvec(0:nmax),fjac(0:nmax,:),1d-3,info)
!!$    write(*,*) "2-D data"
!!$    if(info .eq. 1)then
!!$       write(*,*) "Simulation converges with order = ",fitted_poly(2)
!!$       write(*,*) "Simulation converges to within ",fitted_poly(3)*100,&
!!$            "% of exact solution."
!!$       write(*,*) "CLAWPACK converges with order = 0.608111"
!!$       write(*,*) "CLAWPACK converges to within 0.49% of exact solution"
!!$       stop
!!$    else
!!$       write(*,*) "Function fitting returned unexpected value"
!!$       write(*,*) "info = ", info
!!$       write(*,*) "The meanings of various info values are found in lmder1.f"
!!$       stop
!!$    end if
!!$    if(fitted_poly(2) < 1d0) GodConvergenceTester2D = 3
!!$    
!!$  end function GodConvergenceTester2D

  integer function GodTester(output_dir)
    use TimeAdvancementStuff
    implicit none
    character(len=*), intent(in) :: output_dir
    real(8), dimension(21,3,4,3) :: main
    integer :: i, j, nx, ny
    real(8), dimension(5) :: riemann_test
    real(8), dimension(21) :: geom_test
    real(8), dimension(21,3,3,3) :: constoprimtest1, constoprimtest2
    real(8), dimension(21) :: constoprimtest3
    real(8), dimension(5) :: left_test, right_test
    real(8) :: max_wave_speed
    real(8), dimension(3) :: grad, norm, tan1, tan2
    integer, dimension(1000) :: prim_update_options
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
    ! Also check against a Matlab solution.
    constoprimtest3 = [&
         0.081125768865785,0.929385970968730,0.775712678608402,&
         0.486791632403172,0.435858588580919,0.446783749429806,&
         0.306349472016557,0.508508655381127,0.510771564172110,&
         0.817627708322262,0.794831416883453,0.644318130193692,&
         0.378609382660268,0.811580458282477,0.532825588799455,&
         0.350727103576883,0.939001561999887,0.875942811492984,&
         0.550156342898422,0.622475086001227,0.022367194643555]
    call primtocons(constoprimtest3)
    constoprimtest3=[&
         0.020787756911647,0.016125326596194,0.010119306121021,&
         0.009060522387274,0.015228249823238,0.446783749429806,&
         0.306349472016557,0.508508655381127,0.510771564172110,&
         0.817627708322262,0.794831416883453,0.644318130193692,&
         0.378609382660268,0.811580458282477,0.532825588799455,&
         0.350727103576883,0.939001561999887,0.875942811492984,&
         0.550156342898422,0.622475086001227,0.022367194643555]&
         -constoprimtest3
    if(any(constoprimtest3**2>1d-15)) GodTester = 1
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

    prim_update_options = 0
    prim_update_options(1) = 1
    prim_update_options(101) = 1
    prim_update_options(102) = 1
    prim_update_options(103) = 0
    prim_update_options(104) = 1
    left(1:5)  = [1d0, 1d0, .75d0, 0d0, 0d0]
    right(1:5) = [.1d0, .125d0, 0d0, 0d0, 0d0]
    GodTester = GodConvergenceTester1D(left,right,.3d0,.18d0,1d-4,3,&
         prim_update_options,base_filename=trim(output_dir)//"1D_Rie_Test_1")
!!$    left(1:5)  = [.4d0, 1d0,-2d0, 0d0, 0d0]
!!$    right(1:5) = [.4d0, 1d0, 2d0, 0d0, 0d0]
!!$    GodTester = GodConvergenceTester1D(left,right,.5d0,.12d0,1d-4,3,&
!!$         prim_update_options,base_filename=trim(output_dir)//"1D_Rie_Test_2")
!!$    left(1:5)  = [1000d0, 1d0, 0d0, 0d0, 0d0]
!!$    right(1:5) = [.01d0, 1d0, 0d0, 0d0, 0d0]
!!$    GodTester = GodConvergenceTester1D(left,right,.5d0,.01d0,1d-4,3,&
!!$         prim_update_options,base_filename=trim(output_dir)//"1D_Rie_Test_3")
!!$    left(1:5)  = [460.894d0, 5.99924d0, 19.5975d0, 0d0, 0d0]
!!$    right(1:5) = [46.095d0, 5.99242d0,-6.19633d0, 0d0, 0d0]
!!$    GodTester = GodConvergenceTester1D(left,right,.4d0,.03d0,1d-4,3,&
!!$         prim_update_options,base_filename=trim(output_dir)//"1D_Rie_Test_4")
!!$    left(1:5)  = [1000d0, 1d0,-19.59745d0, 0d0, 0d0]
!!$    right(1:5) = [.01d0, 1d0,-19.59745d0, 0d0, 0d0]
!!$    GodTester = GodConvergenceTester1D(left,right,.8d0,.01d0,1d-4,3,&
!!$         prim_update_options,base_filename=trim(output_dir)//"1D_Rie_Test_5")
!    GodTester = GodConvergenceTester2D(left,right,.18d0,1d-4,2,"2D_Rie_Test_1")
  end function GodTester

  subroutine NormShockInit(nx,main)
    implicit none
    integer, intent(in) :: nx
    real(8), dimension(21,-1:nx,-1:1,-1:1), intent(out) :: main
    real(8), dimension(3) :: upstream, downstream
    real(8), dimension(2) :: xrange
    integer :: half_nx, i,j,k

    xrange = [main(18,0,0,0),main(18,nx-1,0,0)]
    upstream = [1d0,1d0,2d0*sqrt(1.4d0)]
    call NormalShockRelations(upstream,1.4d0,downstream)
    
    if(mod(nx,2).eq.0)then
       half_nx = nx/2
    else
       write(*,*) "Warning: NormShockInit called with an odd value for nx."
       write(*,*) "         Press any key to continue, but the normal shock"
       write(*,*) "         will be slightly offset from center."
       read(*,*)
       half_nx = nx/2
       write(*,*) "         The true halfway point is ",.5*(xrange(2)-xrange(1))
       write(*,*) "         The initial shock location is ",&
            .5*(xrange(2)-xrange(1))-.5*(xrange(2)-xrange(1))/(nx-1)
    end if
    main = 0d0
    forall (i=-1:half_nx-1, j=-1:1, k=-1:1)
       main(1:3,i,j,k) = upstream
    end forall
    forall (i=half_nx:nx, j=-1:1, k=-1:1)
       main(1:3,i,j,k) = downstream
    end forall
    main(6,:,:,:) = 1d0
    main(10,:,:,:) = main(6,:,:,:)
    main(14,:,:,:) = main(6,:,:,:)
    forall (i=-1:nx,j=-1:1,k=-1:1)
       main(18,i,j,k) = real(i,8)
       main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
    end forall
    main(18,:,:,:) = main(18,:,:,:) - .5d0*maxval(main(18,:,:,:))
  end subroutine NormShockInit

  subroutine NormalShockRelations(in,gamma,out)
    implicit none
    real(8), dimension(3), intent(in) :: in
    real(8), intent(in) :: gamma
    real(8), dimension(3), intent(out) :: out
    real(8) :: a1, M1, p, d, u, a2, M2
    
    p = in(1); d = in(2); u = in(3)
    a1 = sqrt(gamma*p/d); M1 = u/a1
    M2 = sqrt((1+((gamma-1d0)*.5d0)*M1**2)/(gamma*M1**2-.5d0*(gamma-1d0)))
    out(1) = p*(1d0 + 2d0*gamma/(gamma-1d0)*(M1**2-1d0))
    out(2) = d*((gamma+1d0)*M1**2)/(2d0 + (gamma-1d0)*M1**2)
    a2 = sqrt(gamma*out(1)/out(2))
    out(3) = M2*a2
  end subroutine NormalShockRelations

  subroutine RieInit1D(left,right,x0,nx,dir,xrange,main)
    implicit none
    real(8), dimension(5), intent(in) :: left, right
    real(8), intent(in) :: x0
    integer, intent(in) :: nx, dir
    real(8), dimension(2), intent(in) :: xrange
    real(8), dimension(:,:,:,:), intent(out) :: main
    real(8), dimension(21,-1:nx,-1:1,-1:1) :: main_temp
    integer :: i,j,k,half_nx 
    real(8) :: dx
    if(mod(nx*x0,1d0).eq.0)then
       half_nx = int(floor(nx*x0))
    else
       write(*,*) "Warning: RieInit1D called with an odd value for nx."
       write(*,*) "         Press any key to continue, but the Riemann"
       write(*,*) "         problem will be slightly offset."
       read(*,*)
       half_nx = int(nx*x0)
!!$       write(*,*) "         The true halfway point is ",.5*(xrange(2)-xrange(1))
!!$       write(*,*) "         The Riemann problem is centered at ",&
!!$            .5*(xrange(2)-xrange(1))-.5*(xrange(2)-xrange(1))/(nx-1)
    end if
    main = 0d0
    main_temp = 0d0
    forall (i=-1:half_nx-1, j=-1:1, k=-1:1)
       main_temp(1:5,i,j,k) = left
    end forall
    forall (i=half_nx:nx, j=-1:1, k=-1:1)
       main_temp(1:5,i,j,k) = right
    end forall
    dx = (xrange(2)-xrange(1))/(nx)
    main_temp(6,:,:,:) = dx
    main_temp(10,:,:,:) = dx
    main_temp(14,:,:,:) = dx
!!$    main_temp(6:14:4,:,:,:) = 1d0
    forall (i=-1:nx, j=-1:1, k=-1:1)
       main_temp(18,i,j,k) = dx*(i+.5) + xrange(1)
       main_temp(21,i,j,k) = Jacobian(main_temp(6:14,i,j,k))
    end forall
    select case(dir)
    case(1)
       main = main_temp
    case(2)
       forall (i=-1:nx,j=-1:1,k=-1:1)
          main(:,k+2,i+2,j+2) = main_temp(:,i,j,k)
       end forall
       main(19,:,:,:) = main(18,:,:,:)
       main(18,:,:,:) = 0d0
    case(3)
       forall (i=-1:nx,j=-1:1,k=-1:1)
          main(:,j+2,k+2,i+2) = main_temp(:,i,j,k)
       end forall
       main(20,:,:,:) = main(18,:,:,:)
       main(18,:,:,:) = 0d0
    case default
       write(*,*) "Error in RieInit1D, invalid dir value!!!"
       stop
    end select
  end subroutine RieInit1D

  subroutine RieInit2D(left,right,phi,nx,xrange,main)
    use GeneralUtilities, only : PI
    implicit none
    real(8), dimension(5) :: left, right
    real(8), intent(in) :: phi
    integer, intent(in) :: nx
    real(8), dimension(2), intent(in) :: xrange
    real(8), dimension(21,-1:nx,-1:nx,-1:1), intent(out) :: main
    integer :: i,j,k,half_nx
    real(8) :: dx
    real(8), dimension(3,3) :: vels_transform

    vels_transform = reshape(&
         [cos(phi),-sin(phi),0d0,sin(phi),cos(phi),0d0,0d0,0d0,1d0],[3,3])
    main = 0d0
    dx = (xrange(2)-xrange(1))/(nx-1)
    main(6,:,:,:) = dx
    main(10,:,:,:) = dx
    main(14,:,:,:) = main(6,:,:,:)
    left(3:5) = matmul(vels_transform,left(3:5))
    right(3:5) = matmul(vels_transform,right(3:5))
    do i = -1, nx
       do j = -1, nx
          do k = -1, 1
             main(18,i,j,k) = dx*i + xrange(1)
             main(19,i,j,k) = dx*j + xrange(1)
             main(21,i,j,k) = Jacobian(main(6:14,i,j,k))
             if( (main(19,i,j,k)-.5d0*(xrange(2)-xrange(1)))&
                  - tan(phi-PI*.5d0)*( main(18,i,j,k)-.5d0&
                  *(xrange(2)-xrange(1))) < 0d0)then
                main(1:5,i,j,k) = left
             else
                main(1:5,i,j,k) = right
             end if
          end do
       end do
    end do
        
  end subroutine RieInit2D

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
  subroutine ToroRiemannConvergence(out)
    implicit none
    real(8), dimension(3), intent(out) :: out
    real(8), dimension(:,:), allocatable :: num
    real(8), dimension(:,:), allocatable :: exact
    real(8), dimension(:), allocatable :: x
    real(8), dimension(21) :: left, right
    real(8) :: max_wave_speed
    real(8), dimension(4) :: riemann_middle_states
    real(8), dimension(5) :: riemann_wave_speeds
    integer :: n, nx, ierr
    real(8), dimension(0:4) :: rmserrors, dxes, x0, t_out
    real(8), dimension(3) :: fitted_poly
    real(8) :: max_value
    max_value = 0d0
    x0 = [.3d0,.5d0,.5d0,.4d0,.8d0]
    t_out = [.18d0,.12d0,.01d0,.03d0,.01d0]
    x0 = .3d0
    t_out = .2d0
    do n = 0, 4
       nx = 100*2**n
       allocate(exact(5,0:nx-1),x(0:nx-1),num(3,0:nx-1))
       num = 0d0
       call ReadToroDat("datToro/1D_Rie_Test_1"//char(48+n)//".dat",num,x,ierr)
       left = 0d0                           ; right = 0d0
       left(1:3) = num(:,0)                 ; right(1:3) = num(:,nx-1)
       left(6) = (x(2)-x(1))                ; right(6) = left(6)
       left(10) = left(6)                   ; right(10) = left(10)
       left(14) = left(6)                   ; right(14) = left(14)
       left(21) = left(6)*left(10)*left(14) ; right(21) = left(21)
       call riemann_solve(left,right,1,nx,(x-x0(n))/left(6)/t_out(n),exact,max_wave_speed,&
            riemann_middle_states,riemann_wave_speeds)
       open(file="datToro/exact"//char(n+48)//".dat",unit=100)
       write(100,*) exact(2,:)
       close(100)
!!$       write(*,*) exact(2,:) - num(2,:)
!!$       read(*,*)
       rmserrors(n) = sqrt(sum(&
            ((num(2,:)-exact(2,:))/maxval(abs(exact(2,:))))**2)/nx)
       max_value = max(max_value,maxval(abs(num(2,:))))
       dxes(n) = x(2)-x(1)
       deallocate(exact,x,num)
    end do
    fitted_poly(1:2) = polyfit(log(dxes),log(rmserrors),1)
    fitted_poly(3) = 0d0
!!$    call minpack_function_fitting(dxes,rmserrors,exponential_with_y_offset,&
!!$         fitted_poly,fvec,fjac,1d-3,info)
    out = fitted_poly
  end subroutine ToroRiemannConvergence
  subroutine ReadToroDat(file,out,x,ierr)
    implicit none
    integer, parameter :: maxnx = 10000
    character(len=*), intent(in) :: file
    real(8), intent(out), dimension(:,:) :: out
    real(8), intent(out), dimension(:) :: x
    integer, intent(out) :: ierr
    real(8), dimension(5) :: read_buffer
    integer :: iostatus, n
    open(unit=837347,file=file,status="old")
    do n = 1, maxnx
       read(837347,*,IOSTAT=iostatus) read_buffer
       if(iostatus.eq.0)then
          out(1,n) = read_buffer(4)
          out(2,n) = read_buffer(2)
          out(3,n) = read_buffer(3)
          x(n) = read_buffer(1)
       else if(iostatus<0)then
          exit
       else
          write(*,*) "Error reading ",trim(file)," in ReadToroDat!!"
          write(*,*) "n = ", n, "iostatus = ", iostatus
          stop
       end if
    end do
  end subroutine ReadToroDat
end module Godunov_tester
