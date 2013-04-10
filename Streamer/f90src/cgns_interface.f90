module cgns_interface
  include 'cgnslib_f.h'
  character(len=32) :: CGNS_FILE_NAME_UCS
! So what do I need to accomplish with this interface? 
!
! What I have been used to having, from my Matlab days,
! is a printout of my primary flow variables at certain
! time steps. So, my most-commonly-used routine will 
! open my file, write those variables, and then close 
! the file again. 
!
! However, CGNS is capable of much more, and I would do
! well to use that capability. So, I'm also going to 
! include my boundary condition specifications. 
! Eventually, I will want to handle multiple zones,
! and CGNS is well equipped for that as well.
contains
  integer function write_initial_data(main_data,filename_in,nx,ny,nz)
    character(len=*), intent(in) :: filename_in
    real(8), intent(in), dimension(21,nx,ny,nz) :: main_data

    character(len=32) :: basename
    character(len=32) :: zonename
    character(len=32) :: solname
    integer, dimension(3,3) :: isize
    real, allocatable, dimension(:,:,:) :: grid_x, grid_y, grid_z
    
    ni = size(main_data,2)
    nj = size(main_data,3)
    nk = size(main_data,4)

    allocate(grid_x(ni,nj,nk), grid_y(ni,nj,nk), grid_z(ni,nj,nk))

    grid_x = main_data(18,:,:,:)
    grid_y = main_data(19,:,:,:)
    grid_z = main_data(20,:,:,:)

    write_initial_data = 0 ! 0 corresponds to success
    CGNS_FILE_NAME_UCS = filename_in
    call cg_open_f(CGNS_FILE_NAME_UCS,CG_MODE_WRITE,index_file,ier)  
    basename = "Base"
    icelldim = 3 ! We work only in volume cells
    iphysdim = 3 ! For 3-D flows
    call cg_base_write_f(index_file,basename,icelldim,iphysdim,&
         index_base, ier)
    zonename = "Zone 1"
    isize(:,1)=[ni,nj,nk] ! Vertex size
    isize(:,2) = isize(:,1)-1 ! Cell size
    isize(:,3) = 0 ! Boundary vertex size (0 for structured grids)
    call cg_zone_write_f(index_file,index_base,zonename,isize,&
         Structured, index_zone, ier)
    ! Write grid coordinates (user must use SIDS-standard names here)
!!$    call cg_grid_write_f(index_file,index_base,index_zone,&
!!$         'GridCoordinates',index_grid)
    call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
         'CoordinateX',grid_x,index_coord,ier)
    call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
         'CoordinateY',grid_y,index_coord,ier)
    call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
         'CoordinateZ',grid_z,index_coord,ier)

    ! Close CGNS file
    call cg_close_f(index_file,ier)
    deallocate(grid_x, grid_y, grid_z)
  end function write_initial_data
end module cgns_interface

!!$program cgns_interface_test
!!$!
!!$! Example taken from CGNS User Guide
!!$!
!!$! Notice that implicit variables are being used.
!!$!
!!$! Example compilation for 64-bit operation is
!!$!
!!$! gfortran -fdefault-real-8 -fdefault-integer-8 -I /usr/local/include
!!$! -L /usr/local/lib -l cgns ../cgns_interface.f90
!!$!
!!$! the cgns_interface module contains the interface functions as well
!!$! as the shared data contained in the cgnslib_f.h file.
!!$  use cgns_interface
!!$  real(8), dimension(21,17,9) :: x, y, z
!!$  integer, dimension(3,3) :: isize
!!$  character(len=32) :: basename, zonename
!!$  character(len=512) :: error_message
!!$  integer, dimension(3) :: irmax, irmin
!!$! Create gridpoints for a simple example:
!!$  ni = 21
!!$  nj = 17
!!$  nk = 9
!!$  do k = 1,nk
!!$     do j = 1, nj
!!$        do i = 1, ni
!!$           x(i,j,k) = real(i-1,kind=8)
!!$           y(i,j,k) = real(j-1,kind=8)
!!$           z(i,j,k) = real(k-1,kind=8)
!!$        end do
!!$     end do
!!$  end do
!!$  write(*,*) "Simple 3-D grid points created"
!!$!
!!$! Write x, y, z to CGNS file
!!$! Open CGNS file for write
!!$  call cg_open_f('grid.cgns',CG_MODE_WRITE,index_file,ier)  
!!$  basename = "Base"
!!$  icelldim = 3
!!$  iphysdim = 3
!!$  call cg_base_write_f(index_file,basename,icelldim,iphysdim,&
!!$       index_base, ier)
!!$  zonename = 'Zone 1'
!!$  isize(:,1)=[21,17,9] ! Vertex size
!!$  isize(:,2) = [21,17,9]-1 ! Cell size
!!$  isize(:,3) = 0 ! Boundary vertex size (0 for structured grids)
!!$  call cg_zone_write_f(index_file,index_base,zonename,isize,&
!!$       Structured, index_zone, ier)
!!$! Write grid coordinates (user must use SIDS-standard names here)
!!$  call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
!!$       'CoordinateX',x,index_coord,ier)
!!$  call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
!!$       'CoordinateY',y,index_coord,ier)
!!$  call cg_coord_write_f(index_file, index_base, index_zone, RealDouble,&
!!$       'CoordinateZ',z,index_coord,ier)
!!$! Close CGNS file
!!$  call cg_close_f(index_file,ier)
!!$  write(*,*) "Successfully wrote grid to file grid.cgns"
!!$  x=0.
!!$  call cg_open_f('grid.cgns',CG_MODE_READ,index_file,ier)
!!$    if (ier .ne. CG_OK)then
!!$     call cg_get_error_f(error_message)
!!$     call cg_error_exit_f()
!!$  end if
!!$  index_base = 1
!!$  index_zone = 1
!!$  call cg_zone_read_f(index_file,index_base,index_zone,zonename,&
!!$       isize,ier)
!!$  irmin = [1,1,1]
!!$  irmax = isize(:,1)
!!$  call cg_coord_read_f(index_file,index_base,index_zone,&
!!$       'CoordinateX',RealDouble,irmin,irmax,x,ier)
!!$  call cg_coord_read_f(index_file,index_base,index_zone,&
!!$       'CoordinateY',RealDouble,irmin,irmax,y,ier)
!!$  call cg_coord_read_f(index_file,index_base,index_zone,&
!!$       'CoordinateZ',RealDouble,irmin,irmax,z,ier)
!!$  call cg_close_f(index_file,ier)
!!$  write(*,*) x(:,1,1)
!!$
!!$end program cgns_interface_test
