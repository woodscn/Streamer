module TimeAdvancementStuff
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

  subroutine ComputationalGrads(Metric,Jac,gradXi,gradEta,gradZeta)
    implicit none
    ! Metric has the form:
    !   1  2  3  4  5  6  7  8  9  
    ! [ A, B, C, L, M, N, P, Q, R ]
    ! Gradients are taken with respect to global Cartesian system
    real(8), intent(in), dimension(9) :: Metric
    real(8), intent(in) :: Jac
    real(8), intent(out), dimension(3) :: gradXi
    real(8), intent(out), dimension(3) :: gradEta
    real(8), intent(out), dimension(3) :: gradZeta
    real(8) :: J

!!$    J = Metric(1)*Metric(5)*Metric(9)&
!!$         + Metric(3)*Metric(4)*Metric(8)&
!!$         + Metric(2)*Metric(6)*Metric(7)&
!!$         - Metric(3)*Metric(5)*Metric(7)&
!!$         - Metric(1)*Metric(6)*Metric(8)&
!!$         - Metric(2)*Metric(4)*Metric(9)
!!$    write(*,*) "Jac = ",Jac
!!$    write(*,*) "J = ",J

    J = Jac
    
    gradXi   = [ Metric(5)*Metric(9) - Metric(6)*Metric(8),&
                 Metric(6)*Metric(7) - Metric(4)*Metric(9),&
                 Metric(4)*Metric(8) - Metric(5)*Metric(7) ]/J

    gradEta  = [ Metric(3)*Metric(8) - Metric(2)*Metric(9),&
                 Metric(1)*Metric(9) - Metric(3)*Metric(7),&
                 Metric(2)*Metric(7) - Metric(1)*Metric(8) ]/J

    gradZeta = [ Metric(2)*Metric(6) - Metric(3)*Metric(5),&
                 Metric(3)*Metric(4) - Metric(1)*Metric(6),&
                 Metric(1)*Metric(5) - Metric(2)*Metric(4) ]/J
    if(sum((matmul(transpose(reshape(metric,[3,3])),&
         reshape([gradXi,gradEta,gradZeta],[3,3]))&
         - reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2) > 1.d-10)then
       write(*,*) "Failed matrix inverse test!"
       write(*,*) (matmul(transpose(reshape(metric,[3,3])),&
            reshape([gradXi,gradEta,gradZeta],[3,3]))&
            - reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2
       read(*,*)
    end if
    
  end subroutine ComputationalGrads

  subroutine CreateColumn(main_data,bc_state,nx,ny,out)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in), dimension(21,nx,ny) :: main_data, bc_state
    real(8), intent(out), dimension(21,nx,ny) :: out
    integer :: i,j
    logical :: fixed_grid = .false.

    out(:,:,:) = bc_state
    if(fixed_grid)then
       ! This is for computing the metric based on fixed grid position
       out(6,:,:) = (main_data(18,:,:)-bc_state(18,:,:))
       out(7,:,:) = (main_data(19,:,:)-bc_state(19,:,:))
       out(8,:,:) = (main_data(20,:,:)-bc_state(20,:,:))
       call TwoDGradient(bc_state(18,:,:),1.d0,1.d0,nx,ny,&
            out(9,:,:),out(12,:,:))
       call TwoDGradient(bc_state(19,:,:),1.d0,1.d0,nx,ny,&
            out(10,:,:),out(13,:,:))
       call TwoDGradient(bc_state(20,:,:),1.d0,1.d0,nx,ny,&
            out(11,:,:),out(14,:,:))
       do i = 1, nx
          do j = 1, ny
             out(21,i,j) = Jacobian(out(6:14,i,j))
          end do
       end do
    else
       ! This is for computing the points based on a fixed metric
       do i = 1, nx
          do j = 1, ny
             out(18:20,i,j) = main_data(18:20,i,j)+&
                  matmul(reshape(out(6:14,i,j),[3,3]),[-1.d0,0.d0,0.d0])
          end do
       end do
    end if
  end subroutine CreateColumn

  real(8) function Jacobian(in)
    implicit none
    real(8), dimension(9), intent(in) :: in
    Jacobian = &
         in(1)*in(5)*in(9) - in(1)*in(6)*in(8) + & ! A*M*R - A*N*Q
         in(2)*in(6)*in(7) - in(2)*in(4)*in(9) + & ! B*N*P - B*L*R
         in(3)*in(4)*in(8) - in(3)*in(5)*in(7)     ! C*L*Q - C*M*P
  end function Jacobian

  subroutine TwoDGradient(in,dx,dy,nx,ny,gradx,grady)
    implicit none
    real(8), intent(in), dimension(nx,ny) :: in
    real(8), intent(in) :: dx, dy
    integer, intent(in) :: nx, ny
    real(8), intent(out), dimension(nx,ny) :: gradx, grady
    integer :: i, j

    ! Central differencing where possible
    gradx(2:nx-1,:) = .5d0*(in(3:nx,:)-in(1:nx-2,:))/dx
    grady(:,2:ny-1) = .5d0*(in(:,3:ny)-in(:,1:ny-2))/dy
    ! Forward & backward differencing elsewhere
    gradx(1,:)  = (in( 2,:)-in(   1,:))/dx
    gradx(nx,:) = (in(nx,:)-in(nx-1,:))/dx
    grady(:,1)  = (in(:, 2)-in(:,   1))/dy
    grady(:,ny) = (in(:,ny)-in(:,ny-1))/dy
  end subroutine TwoDGradient

end module TimeAdvancementStuff
