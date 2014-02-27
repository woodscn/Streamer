module TimeAdvancementStuff
use GeneralUtilities
contains
  logical function CheckCreateColumn(main_data, bc_state, nx, ny)
    implicit none
    real(8), intent(in), dimension(21,nx,ny) :: main_data, bc_state
    integer, intent(in) :: nx, ny
    integer :: i, j
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    
    CheckCreateColumn = .false.
    do i = 1, size(main_data,2)
       do j = 1, size(main_data,3)
          call ComputationalGrads(main_data(6:14,i,j),main_data(21,i,j),&
               gradXi, gradEta, gradZeta)
          if(dot_product(gradXi,main_data(18:20,i,j)-bc_state(18:20,i,j))&
               >= 1.d0)then
!!$          write(*,*) main_data(6:14,i,j)
!!$          write(*,*) gradXi,gradEta,gradZeta
!!$          write(*,*) main_data(18:20,i,j)-bc_state(18:20,i,j)
!!$          write(*,*) i,j
!!$          read(*,*)
             CheckCreateColumn = .true.
             return
          end if
       end do
    end do
  end function CheckCreateColumn

  subroutine CreateColumn(main_data,bc_state,ny,nz,out)
    implicit none
    integer, intent(in) :: ny, nz
    real(8), intent(in), dimension(21,ny,nz) :: main_data, bc_state
    real(8), intent(out), dimension(21,ny,nz) :: out
    integer :: i,j
    logical :: fixed_grid = .true.
    real(8), parameter :: dx = 1.d0, dy = 1.d0, dz = 1.d0
    
    out(:,:,:) = bc_state
    if(fixed_grid)then
       ! This is for computing the metric based on fixed grid position
       out(6,:,:) = (main_data(18,:,:)-bc_state(18,:,:))/dx
       out(7,:,:) = (main_data(19,:,:)-bc_state(19,:,:))/dx
       out(8,:,:) = (main_data(20,:,:)-bc_state(20,:,:))/dx
       call TwoDGradient(bc_state(18,:,:),dy,dz,ny,nz,&
            out(9,:,:),out(12,:,:))
       call TwoDGradient(bc_state(19,:,:),dy,dz,ny,nz,&
            out(10,:,:),out(13,:,:))
       call TwoDGradient(bc_state(20,:,:),dy,dz,ny,nz,&
            out(11,:,:),out(14,:,:))
       if(nz == 1) out(14,:,:) = 1.d0
       do i = 1, ny
          do j = 1, nz
             out(21,i,j) = Jacobian(out(6:14,i,j))
          end do
       end do
    else
       ! This is for computing the points based on a fixed metric
       do i = 1, ny
          do j = 1, nz
             out(18:20,i,j) = main_data(18:20,i,j)+&
                  matmul(reshape(out(6:14,i,j),[3,3]),[-1.d0,0.d0,0.d0])
          end do
       end do
    end if
    out(:,:,:) = bc_state
  end subroutine CreateColumn

  subroutine write_files_matlab(main_data,time,nx,ny,first_flag)!(x,y,E,h,t,dt)
    ! Inputs: x, y, E (array of conserved variables), t, dt
    ! Outputs: two_D.dat, an ASCII file that contains the primitive variables together
    !          with their coordinate values and the associated time. Data is designed
    !          to be read by the matlab file geom_data_reader.m. It should be noted
    !          that the file is opened and headers are written in the main program.
!!$  use global_data, only: dimensions, output_file_name
!!$  use types, only: node_data
!!$  use node_array, only: get_node

    implicit none
    real(8), intent(in), dimension(21,nx,ny) :: main_data
    real(8), intent(in) :: time
    logical, intent(in), optional :: first_flag
    integer, intent(in) :: nx, ny

    real(8), dimension(21) :: current
    logical :: first
    integer :: i,j
    character(len=128) :: error_message
    integer :: error

    first = .false.
    if(present(first_flag))first = first_flag
    if(first)then
       open( unit = 3141, file = 'output_mat.dat', iostat = error, iomsg = error_message, action = 'write' )
       write(3141,*) ' nt= ' , 1000
       !        write(3141,*) ' '
       if( error /= 0 )then
          write(*,*) 'Error opening file for output'
          write(*,*) error_message
          stop
       end if
    end if
    ! We output primitive variables
    write(3141,*) ' '
    write(3141,*) 't=',time,'nx=',nx,'ny=',ny,'dt=',0.
    write(3141,*) ' '
    write(3141,*) 'x= ','y= ','u= ','v= ','rho= ','P= '
    do i = 1 , nx
       do j = 1 , ny
          current = main_data(:,i,j)
          write(3141,*) ' ', current(18) , current(19) , current(3) , current(4) , current(2) , current(1) &
               , current(6) , current(7) , current(9) , current(10) , current(21)
       end do
    end do
  end subroutine write_files_matlab

end module TimeAdvancementStuff
