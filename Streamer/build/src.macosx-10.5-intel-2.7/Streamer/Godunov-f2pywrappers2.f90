!     -*- f90 -*-
!     This file is autogenerated with f2py (version:2)
!     It contains Fortran 90 wrappers to fortran functions.

      subroutine f2pywrap_generalutilities_metricinverse (metricinversef&
     &2pywrap, metric)
      use generalutilities, only : metricinverse
      real(kind=8) metric(9)
      real(kind=8) metricinversef2pywrap(9)
      metricinversef2pywrap = metricinverse(metric)
      end subroutine f2pywrap_generalutilities_metricinverse
      subroutine f2pywrap_generalutilities_metrictomatrix (metrictomatri&
     &xf2pywrap, metric)
      use generalutilities, only : metrictomatrix
      real(kind=8) metric(9)
      real(kind=8) metrictomatrixf2pywrap(3,3)
      metrictomatrixf2pywrap = metrictomatrix(metric)
      end subroutine f2pywrap_generalutilities_metrictomatrix
      subroutine f2pywrap_generalutilities_jacobian (jacobianf2pywrap, i&
     &n)
      use generalutilities, only : jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      jacobianf2pywrap = jacobian(in)
      end subroutine f2pywrap_generalutilities_jacobian
      subroutine f2pywrap_generalutilities_gradstomatrix (gradstomatrixf&
     &2pywrap, grad1, grad2, grad3)
      use generalutilities, only : gradstomatrix
      real(kind=8) grad1(3)
      real(kind=8) grad2(3)
      real(kind=8) grad3(3)
      real(kind=8) gradstomatrixf2pywrap(3,3)
      gradstomatrixf2pywrap = gradstomatrix(grad1, grad2, grad3)
      end subroutine f2pywrap_generalutilities_gradstomatrix
      
      subroutine f2pyinitgeneralutilities(f2pysetupfunc)
      use generalutilities, only : gamma7
      use generalutilities, only : gamma6
      use generalutilities, only : gamma5
      use generalutilities, only : gamma4
      use generalutilities, only : gamma3
      use generalutilities, only : gamma2
      use generalutilities, only : gamma1
      use generalutilities, only : eps
      use generalutilities, only : gamma_const
      use generalutilities, only : epss
      use generalutilities, only : computationalgrads
      use generalutilities, only : twodgradient
      interface 
      subroutine f2pywrap_generalutilities_metricinverse (metricinversef&
     &2pywrap, metricinverse, metric)
      real(kind=8) metric(9)
      real(kind=8) metricinverse(9)
      real(kind=8) metricinversef2pywrap(9)
      end subroutine f2pywrap_generalutilities_metricinverse 
      subroutine f2pywrap_generalutilities_metrictomatrix (metrictomatri&
     &xf2pywrap, metrictomatrix, metric)
      real(kind=8) metric(9)
      real(kind=8) metrictomatrix(3,3)
      real(kind=8) metrictomatrixf2pywrap(3,3)
      end subroutine f2pywrap_generalutilities_metrictomatrix 
      subroutine f2pywrap_generalutilities_jacobian (jacobianf2pywrap, j&
     &acobian, in)
      real(kind=8) jacobian
      real(kind=8) in(9)
      real(kind=8) jacobianf2pywrap
      end subroutine f2pywrap_generalutilities_jacobian 
      subroutine f2pywrap_generalutilities_gradstomatrix (gradstomatrixf&
     &2pywrap, gradstomatrix, grad1, grad2, grad3)
      real(kind=8) grad1(3)
      real(kind=8) grad2(3)
      real(kind=8) grad3(3)
      real(kind=8) gradstomatrix(3,3)
      real(kind=8) gradstomatrixf2pywrap(3,3)
      end subroutine f2pywrap_generalutilities_gradstomatrix
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(gamma7,gamma6,gamma5,gamma4,gamma3,gamma2,gamma&
     &1,eps,gamma_const,epss,f2pywrap_generalutilities_metricinverse,f2p&
     &ywrap_generalutilities_metrictomatrix,computationalgrads,f2pywrap_&
     &generalutilities_jacobian,f2pywrap_generalutilities_gradstomatrix,&
     &twodgradient)
      end subroutine f2pyinitgeneralutilities

      subroutine f2pywrap_generalutilitiestest_guerrorreader (guerrorrea&
     &derf2pywrap, in)
      use generalutilitiestest, only : guerrorreader
      integer in
      integer guerrorreaderf2pywrap
      guerrorreaderf2pywrap = guerrorreader(in)
      end subroutine f2pywrap_generalutilitiestest_guerrorreader
      subroutine f2pywrap_generalutilitiestest_gutest (gutestf2pywrap)
      use generalutilitiestest, only : gutest
      integer gutestf2pywrap
      gutestf2pywrap = gutest()
      end subroutine f2pywrap_generalutilitiestest_gutest
      subroutine f2pywrap_generalutilitiestest_polyfit (polyfitf2pywrap,&
     & vx, vy, d, f2py_vx_d0, f2py_vy_d0)
      use generalutilitiestest, only : polyfit
      integer d
      integer f2py_vx_d0
      integer f2py_vy_d0
      real(kind=8) vx(f2py_vx_d0)
      real(kind=8) vy(f2py_vy_d0)
      real(kind=8) polyfitf2pywrap(d + 1)
      polyfitf2pywrap = polyfit(vx, vy, d)
      end subroutine f2pywrap_generalutilitiestest_polyfit
      
      subroutine f2pyinitgeneralutilitiestest(f2pysetupfunc)
      interface 
      subroutine f2pywrap_generalutilitiestest_guerrorreader (guerrorrea&
     &derf2pywrap, guerrorreader, in)
      integer guerrorreader
      integer in
      integer guerrorreaderf2pywrap
      end subroutine f2pywrap_generalutilitiestest_guerrorreader 
      subroutine f2pywrap_generalutilitiestest_gutest (gutestf2pywrap, g&
     &utest)
      integer gutest
      integer gutestf2pywrap
      end subroutine f2pywrap_generalutilitiestest_gutest 
      subroutine f2pywrap_generalutilitiestest_polyfit (polyfitf2pywrap,&
     & polyfit, vx, vy, d, f2py_vx_d0, f2py_vy_d0)
      integer d
      integer f2py_vx_d0
      integer f2py_vy_d0
      real(kind=8) vx(f2py_vx_d0)
      real(kind=8) vy(f2py_vy_d0)
      real(kind=8) polyfit(d + 1)
      real(kind=8) polyfitf2pywrap(d + 1)
      end subroutine f2pywrap_generalutilitiestest_polyfit
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_generalutilitiestest_guerrorreader,f2p&
     &ywrap_generalutilitiestest_gutest,f2pywrap_generalutilitiestest_po&
     &lyfit)
      end subroutine f2pyinitgeneralutilitiestest

      subroutine f2pywrap_riemann_riemann_solve (riemann_solvef2pywrap, &
     &left, right, geom_avg, max_wave_speed, t_out)
      use riemann, only : riemann_solve
      real(kind=8) max_wave_speed
      real(kind=8) t_out
      real(kind=8) left(5)
      real(kind=8) right(5)
      real(kind=8) geom_avg(21)
      real(kind=8) riemann_solvef2pywrap(5)
      riemann_solvef2pywrap = riemann_solve(left, right, geom_avg, max_w&
     &ave_speed, t_out)
      end subroutine f2pywrap_riemann_riemann_solve
      subroutine f2pywrap_riemann_guessp (guesspf2pywrap, left, right, f&
     &2py_left_d0, f2py_right_d0)
      use riemann, only : guessp
      integer f2py_left_d0
      integer f2py_right_d0
      real(kind=8) left(f2py_left_d0)
      real(kind=8) right(f2py_right_d0)
      real(kind=8) guesspf2pywrap
      guesspf2pywrap = guessp(left, right)
      end subroutine f2pywrap_riemann_guessp
      subroutine f2pywrap_riemann_u_fun (in, pstar, f, df, f2py_in_d0)
      use riemann, only : u_fun
      real(kind=8) pstar
      real(kind=8) f
      real(kind=8) df
      integer f2py_in_d0
      real(kind=8) in(f2py_in_d0)
      call u_fun(in, pstar, f, df)
      end subroutine f2pywrap_riemann_u_fun
      subroutine f2pywrap_riemann_beta (betaf2pywrap, psi)
      use riemann, only : beta
      real(kind=8) psi
      real(kind=8) betaf2pywrap
      betaf2pywrap = beta(psi)
      end subroutine f2pywrap_riemann_beta
      subroutine f2pywrap_riemann_sample (x, left, right, geom_avg, psta&
     &r, ustar, dstarl, dstarr, out, f2py_left_d0, f2py_right_d0, f2py_g&
     &eom_avg_d0, f2py_out_d0)
      use riemann, only : sample
      real(kind=8) x
      real(kind=8) pstar
      real(kind=8) ustar
      real(kind=8) dstarl
      real(kind=8) dstarr
      integer f2py_left_d0
      integer f2py_right_d0
      integer f2py_geom_avg_d0
      integer f2py_out_d0
      real(kind=8) left(f2py_left_d0)
      real(kind=8) right(f2py_right_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) out(f2py_out_d0)
      call sample(x, left, right, geom_avg, pstar, ustar, dstarl, dstarr&
     &, out)
      end subroutine f2pywrap_riemann_sample
      
      subroutine f2pyinitriemann(f2pysetupfunc)
      use riemann, only : test_flag
      use riemann, only : riemann_test_flag
      use riemann, only : verbose
      use riemann, only : test_sol
      interface 
      subroutine f2pywrap_riemann_riemann_solve (riemann_solvef2pywrap, &
     &riemann_solve, left, right, geom_avg, max_wave_speed, t_out)
      real(kind=8) max_wave_speed
      real(kind=8) t_out
      real(kind=8) left(5)
      real(kind=8) right(5)
      real(kind=8) geom_avg(21)
      real(kind=8) riemann_solve(5)
      real(kind=8) riemann_solvef2pywrap(5)
      end subroutine f2pywrap_riemann_riemann_solve 
      subroutine f2pywrap_riemann_guessp (guesspf2pywrap, guessp, left, &
     &right, f2py_left_d0, f2py_right_d0)
      real(kind=8) guessp
      integer f2py_left_d0
      integer f2py_right_d0
      real(kind=8) left(f2py_left_d0)
      real(kind=8) right(f2py_right_d0)
      real(kind=8) guesspf2pywrap
      end subroutine f2pywrap_riemann_guessp 
      subroutine f2pywrap_riemann_u_fun (in, pstar, f, df, f2py_in_d0)
      real(kind=8) pstar
      real(kind=8) f
      real(kind=8) df
      integer f2py_in_d0
      real(kind=8) in(f2py_in_d0)
      end subroutine f2pywrap_riemann_u_fun 
      subroutine f2pywrap_riemann_beta (betaf2pywrap, beta, psi)
      real(kind=8) beta
      real(kind=8) psi
      real(kind=8) betaf2pywrap
      end subroutine f2pywrap_riemann_beta 
      subroutine f2pywrap_riemann_sample (x, left, right, geom_avg, psta&
     &r, ustar, dstarl, dstarr, out, f2py_left_d0, f2py_right_d0, f2py_g&
     &eom_avg_d0, f2py_out_d0)
      real(kind=8) x
      real(kind=8) pstar
      real(kind=8) ustar
      real(kind=8) dstarl
      real(kind=8) dstarr
      integer f2py_left_d0
      integer f2py_right_d0
      integer f2py_geom_avg_d0
      integer f2py_out_d0
      real(kind=8) left(f2py_left_d0)
      real(kind=8) right(f2py_right_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) out(f2py_out_d0)
      end subroutine f2pywrap_riemann_sample
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(test_flag,riemann_test_flag,verbose,test_sol,f2&
     &pywrap_riemann_riemann_solve,f2pywrap_riemann_guessp,f2pywrap_riem&
     &ann_u_fun,f2pywrap_riemann_beta,f2pywrap_riemann_sample)
      end subroutine f2pyinitriemann

      subroutine f2pywrap_riemann_tester_rieerrorreader (rieerrorreaderf&
     &2pywrap, in)
      use riemann_tester, only : rieerrorreader
      integer in
      integer rieerrorreaderf2pywrap
      rieerrorreaderf2pywrap = rieerrorreader(in)
      end subroutine f2pywrap_riemann_tester_rieerrorreader
      subroutine f2pywrap_riemann_tester_rietester (rietesterf2pywrap)
      use riemann_tester, only : rietester
      integer rietesterf2pywrap
      rietesterf2pywrap = rietester()
      end subroutine f2pywrap_riemann_tester_rietester
      
      subroutine f2pyinitriemann_tester(f2pysetupfunc)
      use riemann_tester, only : riemann_sol
      use riemann_tester, only : test_geom
      use riemann_tester, only : mws_test
      interface 
      subroutine f2pywrap_riemann_tester_rieerrorreader (rieerrorreaderf&
     &2pywrap, rieerrorreader, in)
      integer rieerrorreader
      integer in
      integer rieerrorreaderf2pywrap
      end subroutine f2pywrap_riemann_tester_rieerrorreader 
      subroutine f2pywrap_riemann_tester_rietester (rietesterf2pywrap, r&
     &ietester)
      integer rietester
      integer rietesterf2pywrap
      end subroutine f2pywrap_riemann_tester_rietester
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(riemann_sol,test_geom,mws_test,f2pywrap_riemann&
     &_tester_rieerrorreader,f2pywrap_riemann_tester_rietester)
      end subroutine f2pyinitriemann_tester

      subroutine f2pywrap_godunov_energy_func (energy_funcf2pywrap, in, &
     &f2py_in_d0)
      use godunov, only : energy_func
      integer f2py_in_d0
      real(kind=8) in(f2py_in_d0)
      real(kind=8) energy_funcf2pywrap
      energy_funcf2pywrap = energy_func(in)
      end subroutine f2pywrap_godunov_energy_func
      subroutine f2pywrap_godunov_primtocons (main, nx, ny, nz, f2py_mai&
     &n_d0, f2py_main_d1, f2py_main_d2, f2py_main_d3)
      use godunov, only : primtocons
      integer nx
      integer ny
      integer nz
      integer f2py_main_d0
      integer f2py_main_d1
      integer f2py_main_d2
      integer f2py_main_d3
      real(kind=8) main(f2py_main_d0,f2py_main_d1,f2py_main_d2,f2py_main&
     &_d3)
      call primtocons(main, nx, ny, nz)
      end subroutine f2pywrap_godunov_primtocons
      subroutine f2pywrap_godunov_constoprim (main, nx, ny, nz, f2py_mai&
     &n_d0, f2py_main_d1, f2py_main_d2, f2py_main_d3)
      use godunov, only : constoprim
      integer nx
      integer ny
      integer nz
      integer f2py_main_d0
      integer f2py_main_d1
      integer f2py_main_d2
      integer f2py_main_d3
      real(kind=8) main(f2py_main_d0,f2py_main_d1,f2py_main_d2,f2py_main&
     &_d3)
      call constoprim(main, nx, ny, nz)
      end subroutine f2pywrap_godunov_constoprim
      subroutine f2pywrap_godunov_invnorm3 (invnorm3f2pywrap, in)
      use godunov, only : invnorm3
      real(kind=8) in(3)
      real(kind=8) invnorm3f2pywrap
      invnorm3f2pywrap = invnorm3(in)
      end subroutine f2pywrap_godunov_invnorm3
      subroutine f2pywrap_godunov_grid_coords (grad, normal, tangential1&
     &, tangential2, f2py_grad_d0)
      use godunov, only : grid_coords
      integer f2py_grad_d0
      real(kind=8) grad(f2py_grad_d0)
      real(kind=8) normal(3)
      real(kind=8) tangential1(3)
      real(kind=8) tangential2(3)
      call grid_coords(grad, normal, tangential1, tangential2)
      end subroutine f2pywrap_godunov_grid_coords
      subroutine f2pywrap_godunov_compute_fluxes (inl, inr, geom_avg, fl&
     &ux_vec, case_no, max_wave_speed, dt, dv_in, debug_flag, f2py_inl_d&
     &0, f2py_inr_d0, f2py_geom_avg_d0, f2py_flux_vec_d0)
      use godunov, only : compute_fluxes
      integer case_no
      real(kind=8) max_wave_speed
      real(kind=8) dt
      logical debug_flag
      integer f2py_inl_d0
      integer f2py_inr_d0
      integer f2py_geom_avg_d0
      integer f2py_flux_vec_d0
      real(kind=8) inl(f2py_inl_d0)
      real(kind=8) inr(f2py_inr_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) flux_vec(f2py_flux_vec_d0)
      real(kind=8) dv_in(3)
      call compute_fluxes(inl, inr, geom_avg, flux_vec, case_no, max_wav&
     &e_speed, dt, dv_in, debug_flag)
      end subroutine f2pywrap_godunov_compute_fluxes
      subroutine f2pywrap_godunov_flux (fluxf2pywrap, in, geom_avg, case&
     &_no, f2py_in_d0, f2py_geom_avg_d0)
      use godunov, only : flux
      integer case_no
      integer f2py_in_d0
      integer f2py_geom_avg_d0
      real(kind=8) in(f2py_in_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) fluxf2pywrap(5)
      fluxf2pywrap = flux(in, geom_avg, case_no)
      end subroutine f2pywrap_godunov_flux
      
      subroutine f2pyinitgodunov(f2pysetupfunc)
      use godunov, only : dzeta
      use godunov, only : deta
      use godunov, only : max_dt
      use godunov, only : deta_inv
      use godunov, only : dxi_inv
      use godunov, only : dv_inv
      use godunov, only : dxi
      use godunov, only : dzeta_inv
      use godunov, only : prim_update
      interface 
      subroutine f2pywrap_godunov_energy_func (energy_funcf2pywrap, ener&
     &gy_func, in, f2py_in_d0)
      real(kind=8) energy_func
      integer f2py_in_d0
      real(kind=8) in(f2py_in_d0)
      real(kind=8) energy_funcf2pywrap
      end subroutine f2pywrap_godunov_energy_func 
      subroutine f2pywrap_godunov_primtocons (main, nx, ny, nz, f2py_mai&
     &n_d0, f2py_main_d1, f2py_main_d2, f2py_main_d3)
      integer nx
      integer ny
      integer nz
      integer f2py_main_d0
      integer f2py_main_d1
      integer f2py_main_d2
      integer f2py_main_d3
      real(kind=8) main(f2py_main_d0,f2py_main_d1,f2py_main_d2,f2py_main&
     &_d3)
      end subroutine f2pywrap_godunov_primtocons 
      subroutine f2pywrap_godunov_constoprim (main, nx, ny, nz, f2py_mai&
     &n_d0, f2py_main_d1, f2py_main_d2, f2py_main_d3)
      integer nx
      integer ny
      integer nz
      integer f2py_main_d0
      integer f2py_main_d1
      integer f2py_main_d2
      integer f2py_main_d3
      real(kind=8) main(f2py_main_d0,f2py_main_d1,f2py_main_d2,f2py_main&
     &_d3)
      end subroutine f2pywrap_godunov_constoprim 
      subroutine f2pywrap_godunov_invnorm3 (invnorm3f2pywrap, invnorm3, &
     &in)
      real(kind=8) invnorm3
      real(kind=8) in(3)
      real(kind=8) invnorm3f2pywrap
      end subroutine f2pywrap_godunov_invnorm3 
      subroutine f2pywrap_godunov_grid_coords (grad, normal, tangential1&
     &, tangential2, f2py_grad_d0)
      integer f2py_grad_d0
      real(kind=8) grad(f2py_grad_d0)
      real(kind=8) normal(3)
      real(kind=8) tangential1(3)
      real(kind=8) tangential2(3)
      end subroutine f2pywrap_godunov_grid_coords 
      subroutine f2pywrap_godunov_compute_fluxes (inl, inr, geom_avg, fl&
     &ux_vec, case_no, max_wave_speed, dt, dv_in, debug_flag, f2py_inl_d&
     &0, f2py_inr_d0, f2py_geom_avg_d0, f2py_flux_vec_d0)
      integer case_no
      real(kind=8) max_wave_speed
      real(kind=8) dt
      logical debug_flag
      integer f2py_inl_d0
      integer f2py_inr_d0
      integer f2py_geom_avg_d0
      integer f2py_flux_vec_d0
      real(kind=8) inl(f2py_inl_d0)
      real(kind=8) inr(f2py_inr_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) flux_vec(f2py_flux_vec_d0)
      real(kind=8) dv_in(3)
      end subroutine f2pywrap_godunov_compute_fluxes 
      subroutine f2pywrap_godunov_flux (fluxf2pywrap, flux, in, geom_avg&
     &, case_no, f2py_in_d0, f2py_geom_avg_d0)
      integer case_no
      integer f2py_in_d0
      integer f2py_geom_avg_d0
      real(kind=8) in(f2py_in_d0)
      real(kind=8) geom_avg(f2py_geom_avg_d0)
      real(kind=8) flux(5)
      real(kind=8) fluxf2pywrap(5)
      end subroutine f2pywrap_godunov_flux
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(dzeta,deta,max_dt,deta_inv,dxi_inv,dv_inv,dxi,d&
     &zeta_inv,f2pywrap_godunov_energy_func,f2pywrap_godunov_primtocons,&
     &f2pywrap_godunov_constoprim,f2pywrap_godunov_invnorm3,f2pywrap_god&
     &unov_grid_coords,prim_update,f2pywrap_godunov_compute_fluxes,f2pyw&
     &rap_godunov_flux)
      end subroutine f2pyinitgodunov

      subroutine f2pywrap_godunov_tester_goderrorreader (goderrorreaderf&
     &2pywrap, in)
      use godunov_tester, only : goderrorreader
      integer in
      integer goderrorreaderf2pywrap
      goderrorreaderf2pywrap = goderrorreader(in)
      end subroutine f2pywrap_godunov_tester_goderrorreader
      subroutine f2pywrap_godunov_tester_godtester (godtesterf2pywrap)
      use godunov_tester, only : godtester
      integer godtesterf2pywrap
      godtesterf2pywrap = godtester()
      end subroutine f2pywrap_godunov_tester_godtester
      
      subroutine f2pyinitgodunov_tester(f2pysetupfunc)
      use godunov_tester, only : godinit
      interface 
      subroutine f2pywrap_godunov_tester_goderrorreader (goderrorreaderf&
     &2pywrap, goderrorreader, in)
      integer goderrorreader
      integer in
      integer goderrorreaderf2pywrap
      end subroutine f2pywrap_godunov_tester_goderrorreader 
      subroutine f2pywrap_godunov_tester_godtester (godtesterf2pywrap, g&
     &odtester)
      integer godtester
      integer godtesterf2pywrap
      end subroutine f2pywrap_godunov_tester_godtester
      end interface
      external f2pysetupfunc
      call f2pysetupfunc(f2pywrap_godunov_tester_goderrorreader,f2pywrap&
     &_godunov_tester_godtester,godinit)
      end subroutine f2pyinitgodunov_tester


