!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_timeadvancementstuff_checkcreatecolumn (checkc&
     &reatecolumnf2pywrap, main_data, bc_state, nx, ny)
      use timeadvancementstuff, only : checkcreatecolumn
      integer nx
      integer ny
      real(kind=8) main_data(21,nx,ny)
      real(kind=8) bc_state(21,nx,ny)
      logical checkcreatecolumnf2pywrap
      checkcreatecolumnf2pywrap = .not.(.not.checkcreatecolumn(main_data&
     &, bc_state, nx, ny))
      end subroutine f2pywrap_timeadvancementstuff_checkcreatecolumn
      subroutine f2pywrap_timeadvancementstuff_jacobian (jacobianf2pywra&
     &p, in)
      use timeadvancementstuff, only : jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      jacobianf2pywrap = jacobian(in)
      end subroutine f2pywrap_timeadvancementstuff_jacobian
      
      subroutine f2pyinittimeadvancementstuff(f2pysetupfunc)
      use timeadvancementstuff, only : computationalgrads
      use timeadvancementstuff, only : createcolumn
      use timeadvancementstuff, only : twodgradient
      interface 
      subroutine f2pywrap_timeadvancementstuff_checkcreatecolumn (checkc&
     &reatecolumnf2pywrap, checkcreatecolumn, main_data, bc_state, nx, n&
     &y)
      logical checkcreatecolumn
      integer nx
      integer ny
      real(kind=8) main_data(21,nx,ny)
      real(kind=8) bc_state(21,nx,ny)
      logical checkcreatecolumnf2pywrap
      end subroutine f2pywrap_timeadvancementstuff_checkcreatecolumn 
      subroutine f2pywrap_timeadvancementstuff_jacobian (jacobianf2pywra&
     &p, jacobian, in)
      real(kind=8) jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      end subroutine f2pywrap_timeadvancementstuff_jacobian
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_timeadvancementstuff_checkcreatecolumn&
     &,computationalgrads,createcolumn,f2pywrap_timeadvancementstuff_jac&
     &obian,twodgradient)
      end subroutine f2pyinittimeadvancementstuff

