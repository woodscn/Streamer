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
               > 1.d0)then
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
  end subroutine CreateColumn

end module TimeAdvancementStuff
