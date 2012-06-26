module TriangulatedBoundaryConditions
  implicit none
contains
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
    call ComputationalGrads(metric, jacobian, &
         gradXi, gradEta, gradZeta)
!!$    write(*,*) "GradEta = ",gradEta
!!$    read(*,*)
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
  
end module TriangulatedBoundaryConditions
