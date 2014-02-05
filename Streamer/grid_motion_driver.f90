module grid_motion_driver
  use h_update_mod
contains
  subroutine h_update_driver(main,options)
    ! Only supports 2-D in xi-eta. 
    ! Assumes that main is size (21 x nxi+2 x neta+2 x 3)
    real(8), dimension(:,:,:,:) :: main
    integer, dimension(:), intent(in) :: options
    !f2py intent(in,out) :: main
    real(8) :: dxi, deta, h0
    real(8), dimension(:,:), allocatable :: h
    integer nxi, neta, grid_motion
    real(8), dimension(3), parameter :: h0_array = [.25, .5, .999]
    real(8), dimension(7), parameter :: dx_array = [1.,.5,.25,.2,2.,4.,5.]
    grid_motion = options(6)
    h0 = h0_array(options(7))
    nxi = size(main,2)
    neta = size(main,3)
    nzeta = size(main,4)
    dxi = dx_array(options(3))
    deta = dx_array(options(4))
    dzeta = dx_array(options(5))

    select case(grid_motion)
    case(0)
       main(15:17,:,:,:) = 0d0
    case(1)
       main(15:17,:,:,:) = h0*main(3:5,:,:,:)
    case(2)
       main(15:17,:,:,:) = h0
    case(3)
       allocate(h(neta-2,nxi-2))
       call h_update(&
            transpose(main( 6,2:nxi-1,2:neta-1,2)),&
            transpose(main( 7,2:nxi-1,2:neta-1,2)),&
            transpose(main( 9,2:nxi-1,2:neta-1,2)),&
            transpose(main(10,2:nxi-1,2:neta-1,2)),&
            h,& 
            transpose(main(3,2:nxi-1,2:neta-1,2)),&
            transpose(main(4,2:nxi-1,2:neta-1,2)),& 
            dxi,deta,h0)
       main(15,2:nxi-1,2:neta-1,2) = transpose(h)*main(3,2:nxi-1,2:neta-1,2)
       main(16,2:nxi-1,2:neta-1,2) = transpose(h)*main(4,2:nxi-1,2:neta-1,2)
       deallocate(h)
    case default
       write(*,*) "Invalid grid motion specification!"
       stop
    end select

  end subroutine h_update_driver
end module grid_motion_driver
