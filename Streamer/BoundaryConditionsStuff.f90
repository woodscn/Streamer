module BoundaryConditionsStuff
  use GeneralUtilities
  implicit none
contains
! This function does not work in f2py. Perhaps it has something to do with the 
! use of too many callback functions (>1). This means that I will basically have
! a more complicated calling sequence when I use normal vectors, but that should
! be acceptable. (18 May 2012)
  function WallReflect(point, normal, vertices,dim)
    ! I need to work out a way to include the effects of metric components
    ! and the actual location of the .stl boundary. I'm not currently sure
    ! how to do that, but it will probably involve modifying the metric 
    ! components of both WallReflect and point (bad practice...).
    ! 
    ! In fact, I need to figure out how to specify a number of other things:
    ! X,Y,Z; J; and so on, as well as the requirement that the metric yield 
    ! a zero flux across the boundary.
    ! point and WallReflect have the form: 
    !   1  2   3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    ! [ p rho  u  v  w  A  B  C  L  M  N  P  Q  R  U  V  W  x  y  z  J ]
    implicit none
    real(8), intent(in), dimension(21) :: point
    real(8), intent(in), dimension(3) :: normal
    real(8), intent(in), dimension(3,3) :: vertices
    integer, intent(in) :: dim
    real(8), dimension(21) :: WallReflect
    real(8), dimension(3,3) :: matrix_inverse
    
    ! Okay, so I need solid wall boundary conditions to do a few things. 
    ! First, I need to ensure that a condition of no wall penetration is 
    ! maintained by reflecting the physical variables.
    WallReflect = 0.d0
    WallReflect(1:2) = point(1:2)
    WallReflect(3:5) = point(3:5) - 2.d0*vectorProjection(point(3:5),normal)
    ! This actually isn't enough by itself. The Riemann solver requires the
    ! velocity vector be expressed in normal and tangential components. 
    ! In order to do that, I need to project the vector using the metric. 
    ! I have to use an average between the two points. So, in order to 
    ! correctly enforce the no-penetration condition, I have to ensure that
    ! the average normal vector is in fact normal to the wall. That part is 
    ! simple: reflect the normal vector across the wall and negate it.
    ! Unfortunately, the normal vector (Grad(xi) e.g.) isn't a simple 
    ! function of the metric components, but rather the components of the 
    ! metric inverse. So, in order to reflect the vector, I need to reflect
    ! the whole thing. I also have to negate the normal vector, which 
    ! depends on which one that is, xi, eta, zeta.
    matrix_inverse = MatrixInverse(reshape(point(6:14),[3,3]))
    matrix_inverse(1,:) = ReflectionOperator(matrix_inverse(1,:),normal)
    matrix_inverse(2,:) = ReflectionOperator(matrix_inverse(2,:),normal)
    matrix_inverse(3,:) = ReflectionOperator(matrix_inverse(3,:),normal)
    matrix_inverse(dim,:) = -matrix_inverse(dim,:)
    WallReflect(6:14) = reshape(MatrixInverse(matrix_inverse),[9])
    ! I reflect the grid velocity as well. Not sure about this bit.
    WallReflect(15:17) = point(15:17) &
         - 2.d0*vectorProjection(point(15:17),normal)
    ! I also have a hunch that I can help control the grid by using the 
    ! actual location of the grid boundaries. I'm not sure exactly how
    ! to implement that, but it's on the menu. For now, I simply copy 
    ! the Jacobian.
    WallReflect(21) = point(21)
    ! I haven't used the vertices because I don't know yet what to do with
    ! them. Therefore, throw a warning if the vertices given are .ne. 0.
    if(sum(vertices**2) > 1e-14)then
       write(*,*) "Warning, vertex functionality not yet implemented in wallreflect"
       stop
    end if
  end function WallReflect
  
  function ReflectionOperator(in,normal)
    ! Though this is written specifically for three-dimensional arrays, the 
    ! only change that should be required to make it n-dimensional is to 
    ! change the dimension of the arrays. Using (:) notation breaks f2py, 
    ! so a generalized version needs to pass the dimensionality as an 
    ! integer argument.
    implicit none
    real(8), dimension(3), intent(in) :: in, normal
    real(8), dimension(3) :: ReflectionOperator
    
    ReflectionOperator = in-2.d0*normal*&
         (dot_product(in,normal))/(dot_product(normal,normal))
  end function ReflectionOperator
  
  function MatrixInverse(in)
    ! Compute the inverse of a 3x3 matrix
    implicit none
    real(8), dimension(3,3), intent(in) :: in
    real(8), dimension(3,3) :: MatrixInverse
    real(8) :: J
    
    J = ( &
         in(1,1)*in(2,2)*in(3,3) + &
         in(1,2)*in(2,3)*in(3,1) + &
         in(1,3)*in(2,1)*in(3,2) - &
         in(1,1)*in(2,3)*in(3,2) - &
         in(1,2)*in(2,1)*in(3,3) - &
         in(1,3)*in(2,2)*in(3,1) )

    MatrixInverse = transpose(reshape( [ & 
         in(3,3)*in(2,2)-in(3,2)*in(2,3) , &
         in(3,2)*in(1,3)-in(3,3)*in(1,2) , &
         in(2,3)*in(1,2)-in(2,2)*in(1,3) , &
         in(3,1)*in(2,3)-in(3,3)*in(2,1) , &
         in(3,3)*in(1,1)-in(3,1)*in(1,3) , &
         in(2,1)*in(1,3)-in(2,3)*in(1,1) , &
         in(3,2)*in(2,1)-in(3,1)*in(2,2) , &
         in(3,1)*in(1,2)-in(3,2)*in(1,1) , &
         in(2,2)*in(1,1)-in(2,1)*in(1,2) ] &
         ,[3,3])/J)
    if(.true. .and. maxval(matmul(in,MatrixInverse)&
            -reshape([1,0,0,0,1,0,0,0,1],[3,3]))**2>1.d-15)then
          write(*,*) "MatrixInverse failed!!"
          write(*,*) matmul(in,MatrixInverse)
          stop
    end if
  end function MatrixInverse

  function vectorProjection(in,normal)
    implicit none
    real(8), intent(in), dimension(3) :: in
    real(8), intent(in), dimension(3) :: normal
    real(8), dimension(3) :: vectorProjection
    vectorProjection = normal*&
         (dot_product(in,normal)/dot_product(normal,normal))
  end function vectorProjection

  function FindClosestPoint(starting_point, grid_points, nx, ny, nz)
    implicit none
    real(8), intent(in), dimension(3) :: starting_point
    real(8), intent(in), dimension(3,nx,ny,nz) :: grid_points
    integer, intent(in) :: nx, ny, nz
    integer :: i, j, k
    real(8) :: displacement, min_displ
    integer, dimension(3) :: min_index, FindClosestPoint

    min_displ = 5.d0**300
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             displacement = sum((starting_point - grid_points(:,i,j,k))**2)
             if(displacement < min_displ)then
                min_displ = displacement
                min_index = [i,j,k]
             end if
          end do
       end do
    end do    
    FindClosestPoint = min_index
  end function FindClosestPoint

  function ComputationalDisplacement(displacement,metric,jacobian)
    implicit none
    real(8), intent(in), dimension(3) :: displacement
    real(8), intent(in), dimension(9) :: metric
    real(8), intent(in) :: jacobian
    real(8), dimension(3) :: ComputationalDisplacement
    real(8), dimension(3) :: gradXi, gradEta, gradZeta

    call ComputationalGrads(metric, jacobian, gradXi, gradEta, gradZeta)
    ComputationalDisplacement = [dot_product(gradXi,displacement),&
         dot_product(gradEta,displacement),dot_product(gradZeta,displacement)]
  end function ComputationalDisplacement

  integer function FortranXiOffset(xi_offset)
    implicit none
    integer, intent(in) :: xi_offset
    FortranXiOffset = xi_offset - 1
  end function FortranXiOffset

  subroutine ApplyInflowConditions(main_data,bc_state,out,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in), dimension(21,nx,ny) :: main_data, bc_state
    real(8), intent(out), dimension(21,nx,ny) :: out
    integer :: i, j
    out = bc_state
    do i = 1, size(main_data,2)
       do j = 1, size(main_data,3)
          if(.not.CheckSupersonic(main_data(:,i,j)))then
             out(1,i,j) = main_data(1,i,j)
          end if
       end do
    end do
  end subroutine ApplyInflowConditions
  
  subroutine ApplyOutflowConditions(main_data,bc_state,out,nx,ny)
    implicit none
    integer, intent(in) :: nx, ny
    real(8), intent(in), dimension(21,nx,ny) :: main_data, bc_state
    optional :: bc_state
    real(8), intent(out), dimension(21,nx,ny) :: out
    integer :: i, j
    out = main_data
    do i = 1, size(main_data,2)
       do j = 1, size(main_data,3)
          if(.not.CheckSupersonic(main_data(:,i,j)))then
             out(1,i,j) = main_data(1,i,j)
          end if
       end do
    end do
  end subroutine ApplyOutflowConditions

!!$  subroutine ApplyReflectionConditions(main_data,patch_id,out,nx,ny,dim,t)
!!$    implicit none
!!$    integer, intent(in) :: nx, ny, dim
!!$    real(8), intent(out), dimension(21,nx,ny) :: out
!!$    real(8), intent(in), dimension(21,nx,ny) :: main_data
!!$    integer, intent(in) :: patch_id
!!$    integer :: i, j
!!$    real(8) :: x, y, z, t
!!$    intent(in) :: t
!!$    real(8), dimension(3) :: normal
!!$    !f2py(
!!$    do i = 1, size(main_data,2)
!!$       do j = 1, size(main_data,3)
!!$          !  function WallReflect(point, normal, vertices,dim)
!!$          x = main_data(18,i,j)
!!$          y = main_data(19,i,j)
!!$          z = main_data(20,i,j)
!!$!          call normal_func(x,y,z,t,normal)
!!$          normal = normal_func(x,y,z,t)
!!$          out(:,i,j) = WallReflect(main_data(:,i,j), normal,&
!!$               reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],[3,3]),&
!!$               dim)
!!$       end do
!!$    end do
!!$  end subroutine ApplyReflectionConditions
!!$
  logical function CheckSupersonic(point,direction_in)
    implicit none
    real(8), intent(in), dimension(21) :: point
    integer, intent(in), optional :: direction_in ! 1: xi, 2: eta, 3: zeta
    integer :: direction
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8) :: normal_velocity

    direction = 1
    if(present(direction_in)) direction = direction_in
    call ComputationalGrads(point(6:15),point(21),gradXi,gradEta,gradZeta)
    select case(direction)
    case(1)
       normal_velocity = dot_product(  gradXi,point(3:5))
    case(2)
       normal_velocity = dot_product( gradEta,point(3:5))
    case(3)
       normal_velocity = dot_product(gradZeta,point(3:5))
    end select
    CheckSupersonic = normal_velocity > SoundSpeed(point)
  end function CheckSupersonic

  real(8) function SoundSpeed(point)
    implicit none
    real(8), intent(in), dimension(21) :: point
    SoundSpeed = sqrt(1.4d0*point(1)/point(2))
  end function SoundSpeed

  subroutine LeadingEdgePointSearch(&
       nodes_x,faces_x,face_nodes_x,point_in,num_nodes,num_faces,inda)
    implicit none
    real(8), intent(in), dimension(3,num_nodes) :: nodes_x
    real(8), intent(in), dimension(3,num_faces) :: faces_x
    integer, intent(in), dimension(3,num_faces) :: face_nodes_x
    real(8), intent(in), dimension(21) :: point_in
    integer, intent(in) :: num_nodes, num_faces
    integer, intent(out) :: inda
    real(8), dimension(3) :: normal, xi_normal, point
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8), dimension(3) :: nodea, nodeb, nodec
    real(8), dimension(9) :: metric
    real(8), dimension(3) :: diff, unit_normal
    integer :: i
    real(8), dimension(3,3) :: nodes
    real(8) :: test2, v3, u, v, test2min
    real(8), dimension(3) :: v0, v1, v2
    real(8) :: jacobian
    logical :: debug = .false.
    
    inda = 0
    metric = point_in(6:14)
    point = point_in(18:20)
    jacobian = point_in(21)
    call ComputationalGrads(metric, jacobian, gradXi, gradEta, gradZeta)
    do i = 1, num_faces
       if(.true.)then
          nodea = nodes_x(:,face_nodes_x(1,i))
          nodeb = nodes_x(:,face_nodes_x(2,i))
          nodec = nodes_x(:,face_nodes_x(3,i))
          nodes = nodes_x(:,face_nodes_x(:,i))
          normal = faces_x(:,i)
          unit_normal = normal/sqrt(sum(normal**2))
          diff = point - nodea
!          xi_normal = [dot_product(gradXi,normal),&
!               dot_product(gradEta,normal),dot_product(gradZeta,normal)]
!          dx = abs(dot_product(dxin,normal))
          v3 = dot_product(diff,normal)/dot_product(normal,normal)
          xi_normal = matmul(reshape([gradXi,gradEta,gradZeta],[3,3]),v3*normal)
          if(sum(xi_normal**2)>1.)then
             ! The triangle is too far from the point in the normal direction,
             ! or the point is on the wrong side of the triangle.
             !write(*,*) "i = ", i, v3*xi_normal
          else
             v2 = diff - v3*normal
             v0 = nodec - nodea
             v1 = nodeb - nodea
             u = (dot_product(v1,v1)*dot_product(v2,v0) - &
                  dot_product(v1,v0)*dot_product(v2,v1))/ &
                  (dot_product(v0,v0)*dot_product(v1,v1) - &
                  dot_product(v0,v1)**2)
             v = (dot_product(v0,v0)*dot_product(v2,v1) - &
                  dot_product(v1,v0)*dot_product(v2,v0))/ &
                  (dot_product(v0,v0)*dot_product(v1,v1) - &
                  dot_product(v0,v1)**2)
             test2 = sum(((point-v3*normal)-matmul(nodes,&
                  [.3333333333333333,.3333333333333333,.3333333333333333]))**2)
!             if(abs(test2) < test2min)then
             if( u>=0. .and. v>=0. .and. u+v<=1.)then
!                test2min = min(abs(test2), test2min)
                inda = i
                if(debug)then
                   write(*,*) 'Point = [',point(1),',',point(2),',',point(3),'];'
                   write(*,*) 'V2 = [',v2(1),',',v2(2),',',v2(3),'];'
                   write(*,*) 'NodeA = [',nodea(1),',',nodea(2),',',nodea(3),'];'
                   write(*,*) 'NodeB = [',nodeb(1),',',nodeb(2),',',nodeb(3),'];'
                   write(*,*) 'NodeC = [',nodec(1),',',nodec(2),',',nodec(3),'];'
                   write(*,*) 'Normal = [',normal(1),',',normal(2),',',normal(3),'];'
                   write(*,*) 'XiNormal = [',xi_normal(1),',',xi_normal(2),','&
                        ,xi_normal(3),'];'
                   read(*,*)
                end if
             end if
          end if
       end if
    end do
  end subroutine LeadingEdgePointSearch

  function ComputeCompuCoordsDelta(point, metric, jacobian, base_point)
    real(8), intent(in), dimension(3) :: point
    real(8), intent(in), dimension(9) :: metric
    real(8), intent(in) :: jacobian
    real(8), intent(in), dimension(3) :: base_point

    integer, dimension(3) :: ComputeCompuCoordsDelta
    
    real(8), dimension(3) :: gradXi, gradEta, gradZeta
    real(8), dimension(3) :: diff

    call ComputationalGrads(metric, jacobian, gradXi, gradEta, gradZeta)
    diff = point-base_point
    ComputeCompuCoordsDelta = [&
         sum(gradXi*diff), sum(gradEta*diff), sum(gradZeta*diff) ]
  end function ComputeCompuCoordsDelta

  logical function pnpoly(npol,xp,yp,x,y)
    ! Check to see if a point lies on the interior of a polygon.
    ! The polygon is given by  x,y points in xp & yp. The point
    ! to test is given by x, y.
    ! Returns true if within the polygon, false otherwise.

    ! Uses the method of counting the number of times a ray from 
    ! the point crosses the polygon. Original C code given on
    ! http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/
    ! and written by Randolph Franklin.

    ! int pnpoly(int npol, float *xp, float *yp, float x, float y)
    !     {
    !       int i, j, c = 0;
    !       for (i = 0, j = npol-1; i < npol; j = i++) {
    !         if ((((yp[i] <= y) && (y < yp[j])) ||
    !              ((yp[j] <= y) && (y < yp[i]))) &&
    !             (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
    !           c = !c;
    !       }
    !       return c;
    !     }
    implicit none
    integer :: npol
    real(8) :: x, y
    real(8), dimension(npol) :: xp, yp
    integer :: i, j, n

    pnpoly = .false.
    n = size(xp)
    do i = 1, size(xp)
       if(i==1)then
          j = n
       else
          j = i-1
       end if
       if((((yp(i)<=y).and.(y<yp(j))).or.&
            ((yp(j)<=y).and.(y<yp(i)))).and.&
            (x<(xp(j)-xp(i))*(y-yp(i))/(yp(j)-yp(i))+xp(i)))&
            pnpoly = .not. pnpoly
       write(*,*) "i = ",i,"j = ",j       
    end do

  end function pnpoly
  
end module BoundaryConditionsStuff
