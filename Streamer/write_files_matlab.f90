
subroutine write_files_matlab(main_data,time,first_flag,nx,ny)!(x,y,E,h,t,dt)
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
