!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_generalutilities_jacobian (jacobianf2pywrap, i&
     &n)
      use generalutilities, only : jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      jacobianf2pywrap = jacobian(in)
      end subroutine f2pywrap_generalutilities_jacobian
      
      subroutine f2pyinitgeneralutilities(f2pysetupfunc)
      use generalutilities, only : computationalgrads
      use generalutilities, only : twodgradient
      interface 
      subroutine f2pywrap_generalutilities_jacobian (jacobianf2pywrap, j&
     &acobian, in)
      real(kind=8) jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      end subroutine f2pywrap_generalutilities_jacobian
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(computationalgrads,f2pywrap_generalutilities_ja&
     &cobian,twodgradient)
      end subroutine f2pyinitgeneralutilities

      subroutine f2pywrap_boundaryconditionsstuff_wallreflect (wallrefle&
     &ctf2pywrap, point, normal, vertices, dim)
      use boundaryconditionsstuff, only : wallreflect
      integer dim
      real(kind=8) point(21)
      real(kind=8) normal(3)
      real(kind=8) vertices(3,3)
      real(kind=8) wallreflectf2pywrap(21)
      wallreflectf2pywrap = wallreflect(point, normal, vertices, dim)
      end subroutine f2pywrap_boundaryconditionsstuff_wallreflect
      subroutine f2pywrap_boundaryconditionsstuff_reflectionoperator (re&
     &flectionoperatorf2pywrap, in, normal)
      use boundaryconditionsstuff, only : reflectionoperator
      real(kind=8) in(3)
      real(kind=8) normal(3)
      real(kind=8) reflectionoperatorf2pywrap(3)
      reflectionoperatorf2pywrap = reflectionoperator(in, normal)
      end subroutine f2pywrap_boundaryconditionsstuff_reflectionoperator
      subroutine f2pywrap_boundaryconditionsstuff_matrixinverse (matrixi&
     &nversef2pywrap, in)
      use boundaryconditionsstuff, only : matrixinverse
      real(kind=8) in(3,3)
      real(kind=8) matrixinversef2pywrap(3,3)
      matrixinversef2pywrap = matrixinverse(in)
      end subroutine f2pywrap_boundaryconditionsstuff_matrixinverse
      subroutine f2pywrap_boundaryconditionsstuff_vectorprojection (vect&
     &orprojectionf2pywrap, in, normal)
      use boundaryconditionsstuff, only : vectorprojection
      real(kind=8) in(3)
      real(kind=8) normal(3)
      real(kind=8) vectorprojectionf2pywrap(3)
      vectorprojectionf2pywrap = vectorprojection(in, normal)
      end subroutine f2pywrap_boundaryconditionsstuff_vectorprojection
      subroutine f2pywrap_boundaryconditionsstuff_findclosestpoint (find&
     &closestpointf2pywrap, starting_point, grid_points, nx, ny, nz)
      use boundaryconditionsstuff, only : findclosestpoint
      integer nx
      integer ny
      integer nz
      real(kind=8) starting_point(3)
      real(kind=8) grid_points(3,nx,ny,nz)
      integer findclosestpointf2pywrap(3)
      findclosestpointf2pywrap = findclosestpoint(starting_point, grid_p&
     &oints, nx, ny, nz)
      end subroutine f2pywrap_boundaryconditionsstuff_findclosestpoint
      subroutine f2pywrap_boundaryconditionsstuff_computationaldisplacem&
     &ent (computationaldisplacementf2pywrap, displacement, metric, jaco&
     &bian)
      use boundaryconditionsstuff, only : computationaldisplacement
      real(kind=8) jacobian
      real(kind=8) displacement(3)
      real(kind=8) metric(9)
      real(kind=8) computationaldisplacementf2pywrap(3)
      computationaldisplacementf2pywrap = computationaldisplacement(disp&
     &lacement, metric, jacobian)
      end subroutine f2pywrap_boundaryconditionsstuff_computationaldispl&
     &acement
      subroutine f2pywrap_boundaryconditionsstuff_fortranxioffset (fortr&
     &anxioffsetf2pywrap, xi_offset)
      use boundaryconditionsstuff, only : fortranxioffset
      integer xi_offset
      integer fortranxioffsetf2pywrap
      fortranxioffsetf2pywrap = fortranxioffset(xi_offset)
      end subroutine f2pywrap_boundaryconditionsstuff_fortranxioffset
      subroutine f2pywrap_boundaryconditionsstuff_checksupersonic (check&
     &supersonicf2pywrap, point, direction_in)
      use boundaryconditionsstuff, only : checksupersonic
      integer direction_in
      real(kind=8) point(21)
      logical checksupersonicf2pywrap
      checksupersonicf2pywrap = .not.(.not.checksupersonic(point, direct&
     &ion_in))
      end subroutine f2pywrap_boundaryconditionsstuff_checksupersonic
      subroutine f2pywrap_boundaryconditionsstuff_soundspeed (soundspeed&
     &f2pywrap, point)
      use boundaryconditionsstuff, only : soundspeed
      real(kind=8) point(21)
      real(kind=8) soundspeedf2pywrap
      soundspeedf2pywrap = soundspeed(point)
      end subroutine f2pywrap_boundaryconditionsstuff_soundspeed
      subroutine f2pywrap_boundaryconditionsstuff_computecompucoordsdelt&
     &a (computecompucoordsdeltaf2pywrap, point, metric, jacobian, base_&
     &point)
      use boundaryconditionsstuff, only : computecompucoordsdelta
      real(kind=8) jacobian
      real(kind=8) point(3)
      real(kind=8) metric(9)
      real(kind=8) base_point(3)
      integer computecompucoordsdeltaf2pywrap(3)
      computecompucoordsdeltaf2pywrap = computecompucoordsdelta(point, m&
     &etric, jacobian, base_point)
      end subroutine f2pywrap_boundaryconditionsstuff_computecompucoords&
     &delta
      subroutine f2pywrap_boundaryconditionsstuff_pnpoly (pnpolyf2pywrap&
     &, npol, xp, yp, x, y)
      use boundaryconditionsstuff, only : pnpoly
      integer npol
      real(kind=8) x
      real(kind=8) y
      real(kind=8) xp(npol)
      real(kind=8) yp(npol)
      logical pnpolyf2pywrap
      pnpolyf2pywrap = .not.(.not.pnpoly(npol, xp, yp, x, y))
      end subroutine f2pywrap_boundaryconditionsstuff_pnpoly
      
      subroutine f2pyinitboundaryconditionsstuff(f2pysetupfunc)
      use boundaryconditionsstuff, only : applyinflowconditions
      use boundaryconditionsstuff, only : applyoutflowconditions
      use boundaryconditionsstuff, only : leadingedgepointsearch
      interface 
      subroutine f2pywrap_boundaryconditionsstuff_wallreflect (wallrefle&
     &ctf2pywrap, wallreflect, point, normal, vertices, dim)
      integer dim
      real(kind=8) point(21)
      real(kind=8) normal(3)
      real(kind=8) vertices(3,3)
      real(kind=8) wallreflect(21)
      real(kind=8) wallreflectf2pywrap(21)
      end subroutine f2pywrap_boundaryconditionsstuff_wallreflect 
      subroutine f2pywrap_boundaryconditionsstuff_reflectionoperator (re&
     &flectionoperatorf2pywrap, reflectionoperator, in, normal)
      real(kind=8) in(3)
      real(kind=8) normal(3)
      real(kind=8) reflectionoperator(3)
      real(kind=8) reflectionoperatorf2pywrap(3)
      end subroutine f2pywrap_boundaryconditionsstuff_reflectionoperator&
     & 
      subroutine f2pywrap_boundaryconditionsstuff_matrixinverse (matrixi&
     &nversef2pywrap, matrixinverse, in)
      real(kind=8) in(3,3)
      real(kind=8) matrixinverse(3,3)
      real(kind=8) matrixinversef2pywrap(3,3)
      end subroutine f2pywrap_boundaryconditionsstuff_matrixinverse 
      subroutine f2pywrap_boundaryconditionsstuff_vectorprojection (vect&
     &orprojectionf2pywrap, vectorprojection, in, normal)
      real(kind=8) in(3)
      real(kind=8) normal(3)
      real(kind=8) vectorprojection(3)
      real(kind=8) vectorprojectionf2pywrap(3)
      end subroutine f2pywrap_boundaryconditionsstuff_vectorprojection 
      subroutine f2pywrap_boundaryconditionsstuff_findclosestpoint (find&
     &closestpointf2pywrap, findclosestpoint, starting_point, grid_point&
     &s, nx, ny, nz)
      integer nx
      integer ny
      integer nz
      real(kind=8) starting_point(3)
      real(kind=8) grid_points(3,nx,ny,nz)
      integer findclosestpoint(3)
      integer findclosestpointf2pywrap(3)
      end subroutine f2pywrap_boundaryconditionsstuff_findclosestpoint 
      subroutine f2pywrap_boundaryconditionsstuff_computationaldisplacem&
     &ent (computationaldisplacementf2pywrap, computationaldisplacement,&
     & displacement, metric, jacobian)
      real(kind=8) jacobian
      real(kind=8) displacement(3)
      real(kind=8) metric(9)
      real(kind=8) computationaldisplacement(3)
      real(kind=8) computationaldisplacementf2pywrap(3)
      end subroutine f2pywrap_boundaryconditionsstuff_computationaldispl&
     &acement 
      subroutine f2pywrap_boundaryconditionsstuff_fortranxioffset (fortr&
     &anxioffsetf2pywrap, fortranxioffset, xi_offset)
      integer fortranxioffset
      integer xi_offset
      integer fortranxioffsetf2pywrap
      end subroutine f2pywrap_boundaryconditionsstuff_fortranxioffset 
      subroutine f2pywrap_boundaryconditionsstuff_checksupersonic (check&
     &supersonicf2pywrap, checksupersonic, point, direction_in)
      logical checksupersonic
      integer direction_in
      real(kind=8) point(21)
      logical checksupersonicf2pywrap
      end subroutine f2pywrap_boundaryconditionsstuff_checksupersonic 
      subroutine f2pywrap_boundaryconditionsstuff_soundspeed (soundspeed&
     &f2pywrap, soundspeed, point)
      real(kind=8) soundspeed
      real(kind=8) point(21)
      real(kind=8) soundspeedf2pywrap
      end subroutine f2pywrap_boundaryconditionsstuff_soundspeed 
      subroutine f2pywrap_boundaryconditionsstuff_computecompucoordsdelt&
     &a (computecompucoordsdeltaf2pywrap, computecompucoordsdelta, point&
     &, metric, jacobian, base_point)
      real(kind=8) jacobian
      real(kind=8) point(3)
      real(kind=8) metric(9)
      real(kind=8) base_point(3)
      integer computecompucoordsdelta(3)
      integer computecompucoordsdeltaf2pywrap(3)
      end subroutine f2pywrap_boundaryconditionsstuff_computecompucoords&
     &delta 
      subroutine f2pywrap_boundaryconditionsstuff_pnpoly (pnpolyf2pywrap&
     &, pnpoly, npol, xp, yp, x, y)
      logical pnpoly
      integer npol
      real(kind=8) x
      real(kind=8) y
      real(kind=8) xp(npol)
      real(kind=8) yp(npol)
      logical pnpolyf2pywrap
      end subroutine f2pywrap_boundaryconditionsstuff_pnpoly
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_boundaryconditionsstuff_wallreflect,f2&
     &pywrap_boundaryconditionsstuff_reflectionoperator,f2pywrap_boundar&
     &yconditionsstuff_matrixinverse,f2pywrap_boundaryconditionsstuff_ve&
     &ctorprojection,f2pywrap_boundaryconditionsstuff_findclosestpoint,f&
     &2pywrap_boundaryconditionsstuff_computationaldisplacement,f2pywrap&
     &_boundaryconditionsstuff_fortranxioffset,applyinflowconditions,app&
     &lyoutflowconditions,f2pywrap_boundaryconditionsstuff_checksuperson&
     &ic,f2pywrap_boundaryconditionsstuff_soundspeed,leadingedgepointsea&
     &rch,f2pywrap_boundaryconditionsstuff_computecompucoordsdelta,f2pyw&
     &rap_boundaryconditionsstuff_pnpoly)
      end subroutine f2pyinitboundaryconditionsstuff


