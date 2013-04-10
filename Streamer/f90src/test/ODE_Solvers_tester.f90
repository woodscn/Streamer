module ODE_Solvers_tester
  ! Test the convergence of ODE_Solvers routines
  ! subroutine Forward_Euler(main,dt,options)
  use ODE_Solvers
  use GeneralUtilitiesTester
contains
  function ODETester()
  use ODE_Solvers
  implicit none
!  character(len=*), intent(in) :: output_dir
  integer :: ODETester 
  real(8), dimension(21,1,1,1) :: init, main
  integer, parameter :: nt = 20
  real(8), dimension(21,1,1,1) :: exact
  real(8), dimension(21,1,1,1,nt,3) :: error
  real(8), dimension(21,1,1,1,nt,2) :: fitted_poly
  integer :: n, m, i, j, k
  real(8) :: dt, t, dt_base
  integer, parameter :: n_iter = 3
  real(8), dimension(n_iter) :: dts

  dt_base = .0001d0
  call random_number(init)
  do m = 1, n_iter
     t=0d0 
     main = init
     dt = dt_base*.01d0**(m-1)
     dts(m) = dt
     do n = 1, nt*2**(m-1)
        call Forward_Euler(main,dt)
        call exact_sol_exp(init,dt*n,exact)
        t = t+dt
        if(mod(n,2**(m-1)).eq.0) error(:,:,:,:,n/(2**(m-1)),m) = main-exact
     end do
  end do
  open(unit = 919191,file="odejunk"//achar(48+m)//".dat")
  write(919191,*) error(1,1,1,1,:,:)
  close(919191)
  do m = 1, size(error,1)
     do i = 1, size(error,2)
        do j = 1, size(error,3)
           do k = 1, size(error,4)
              do n = 1, size(error,5)
                 fitted_poly(m,i,j,k,n,:) = polyfit(log(dts),&
                      log(abs(error(m,i,j,k,n,:)/error(m,i,j,k,nt,1))),1)
              end do
           end do
        end do
     end do
  end do
  write(*,*) fitted_poly(1,1,1,1,1,:)


end function ODETester
subroutine exact_sol_exp(in,t,out)
  real(8), dimension(:,:,:,:), intent(in) :: in
  real(8), intent(in) :: t
  real(8), dimension(:,:,:,:), intent(out) :: out
  out = (in - 1d0) + exp(t)
end subroutine exact_sol_exp
   
!!$  subroutine ODE_Convergence(init,exact,out,options)
!!$    real(8), dimension(:,:,:,:), intent(in) :: init
!!$    integer, dimension(:), intent(in) :: options
!!$    real(8), dimension(:,:,:,:), intent(out) :: out
!!$    integer :: nx, ny, nz, nt, n, m, number_of_refinements
!!$    real(8) :: t_out, dt_base, dt, t
!!$    real(8), dimension(:,:,:,:), allocatable :: main
!!$    real(8), dimension(:,:,:,:,:), allocatable :: convergence_rates
!!$    real(8), dimension(:,:,:,:,:,:), allocatable :: error
!!$    
!!$    nx = size(init,2)
!!$    ny = size(init,3)
!!$    nz = size(init,4)
!!$    nt = 20
!!$    do n = 1, 
!!$    select case(options(1))
!!$    case(1) ! Forward_Euler
!!$       dt_base = .01d0
!!$    case default
!!$       write(*,*) "Invalid algorithm specification in ODE_Convergence!"
!!$       stop
!!$    end select
!!$    t_out = dt_base*nt
!!$    allocate(main(size(init,1),nx,ny,nz),&
!!$         convergence_rates(size(init,1),nx,ny,nz,nt),&
!!$         error(size(init,1),nx,ny,nz,nt,options(2)))
!!$    main = init
!!$
!!$    do m = 0, options(2)
!!$       dt = dt_base*.5d0**m
!!$       do n = 1, nt*2**m
!!$          select case(options(1))
!!$          case(1)
!!$             call Forward_Euler(main,dt)
!!$          end select
!!$          if(mod(n,2**m)==0)then
!!$             error(:,:,:,:,n/(2**m),m+1) = main - exact
!!$          end if
!!$       end do
!!$    end do
!!$    
!!$    deallocate(main,convergence_rates)
!!$  end function ODE_Convergence_1D
end module ODE_Solvers_tester

