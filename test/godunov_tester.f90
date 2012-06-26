program godunov_tester
  use Streamer
  implicit none
  integer, parameter :: nx=4,ny=1,nz=1
  integer :: n
  real(8) :: sample_main_avg(21), sample_interface(5), dt, temp(21,nx+2,ny+2,nz+2)
  sample_main_avg = [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1 ]
  sample_interface= [ 1, 1, 1, 1, 1 ]
  write(*,*) "Given sample geometric variables: "
  write(*,*) "[ A B C L M N P Q R U V W ] = :"
  write(*,*) real(sample_main_avg( 6: 8))
  write(*,*) real(sample_main_avg( 9:11))
  write(*,*) real(sample_main_avg(12:14))
  write(*,*) real(sample_main_avg(15:17))
  write(*,*) "And sample interface variables: "
  write(*,*) "[ p rho u v w ]"
  write(*,*) real(sample_interface)
  write(*,*) "The flux vectors are given by: "
  write(*,*) size(sample_interface), size(sample_main_avg)
  write(*,*) real(flux(sample_interface,sample_main_avg,1))
  write(*,*) real(flux(sample_interface,sample_main_avg,2))
  write(*,*) real(flux(sample_interface,sample_main_avg,3))

dt = .0001
temp = 0.d0
temp( 1,:,:,:) = 460.894
temp( 2,:,:,:) = 5.99924
temp( 3,:,:,:) = 19.5975
temp( 4,:,:,:) = 0
temp( 5,:,:,:) = 0
temp( 6,:,:,:) = .01
temp( 7,:,:,:) = 0
temp( 8,:,:,:) = 0
temp( 9,:,:,:) = 0
temp(10,:,:,:) = .01
temp(11,:,:,:) = 0
temp(12,:,:,:) = 0
temp(13,:,:,:) = 0
temp(14,:,:,:) = .01
temp(15,:,:,:) = 0
temp(16,:,:,:) = 0
temp(17,:,:,:) = 0
temp(18,:,:,:) = 0
temp(19,:,:,:) = 0
temp(20,:,:,:) = 0
temp(21,:,:,:) = .000001

temp( 1,nx/2+2:,:,:) = 46.0950
temp( 2,nx/2+2:,:,:) = 5.99242
temp( 3,nx/2+2:,:,:) =-6.19633

do n = 1,350
   call prim_update(temp,.0001d0,1,nx,ny,nz)
end do

! flux(in,geom_avg,case_no)
end program godunov_tester
