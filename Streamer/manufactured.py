import numpy
import sympy
import scipy.integrate
import Euler_UCS_manufactured
import multiprocessing
from sympy.core.cache import * # Needed for clear_cache
#from pygsl import monte
#import pygsl.rng
import gc
t=sympy.Symbol('t')
xi=sympy.Symbol('xi')
eta=sympy.Symbol('eta')
zeta=sympy.Symbol('zeta')
def mc_integrate(f,ranges,calls,args=()):
#    pool = multiprocessing.Pool()
    f_sum = 0
    for n in xrange(calls):
        coords_lst = []
        for rng in ranges:
            coords_lst.append(rng[0] + rng[1]*numpy.random.random())
        f_sum += f(coords_lst,args)
#        clear_cache()
    vol = numpy.prod([float(rng[1])-float(rng[0]) for rng in ranges])
    return vol/calls*f_sum    
def sample_point(eqns,t_in,xi_in,eta_in,zeta_in):
    out = []
    for eqn in eqns:
        out.append(eqn.evalf(subs={t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}))
    return out
def smooth_sols_1():
    return [1+t,  # pressure
     1+t,  # density
     xi,    # x-velocity
     eta,    # y-velocity
     zeta,    # z-velocity
     1,    # dx_dxi
     0,    # dy_dxi
     0,    # dz_dxi
     0,    # dx_deta
     1,    # dy_deta
     0,    # dz_deta
     0,    # dx_dzeta
     0,    # dy_dzeta
     1,    # dz_dzeta
     0,    # dx_dt
     0,    # dy_dt
     0,    # dz_dt
     0,    # x
     0,    # y
     0]    # z
def jacobian(a,b,c,l,m,n,p,q,r):
    return a*(m*r-n*q)+b*(n*p-l*r)+c*(l*q-m*p)
def gradxi(a,b,c,l,m,n,p,q,r):
    return sympy.Matrix([(m*r-n*q)/jacobian(a,b,c,l,m,n,p,q,r),
            (n*p-l*r)/jacobian(a,b,c,l,m,n,p,q,r),
            (l*q-m*p)/jacobian(a,b,c,l,m,n,p,q,r)])
def shock_speed(pressure,density,velocity,a,b,c,l,m,n,p,q,r,u,
                gamma,pressure_hi,pm):
    return(sympy.sqrt(sum(gradxi(a,b,c,l,m,n,p,q,r).applyfunc(lambda x:x**2)))*(
        velocity-u+pm*sympy.sqrt(gamma*pressure/density)*sympy.sqrt(
            (gamma+1)/(2*gamma)*(pressure_hi/pressure-1)+1)))
def shocked_sol_1():
    xmin, xmax = 0, 1
    nx = 100
    dx = (xmax-xmin)/(nx)
    gamma = 1.4
    shock_pos_0 = .1
    pressure_lo, density_lo = 1, 1
    velocity_lo = 1.5*numpy.sqrt(1.4*pressure_lo/density_lo)
    a_lo, b_lo, c_lo, l_lo, m_lo, n_lo, p_lo, q_lo, r_lo = (1,0,0,0,1,0,0,0,1)
    u_lo, v_lo, w_lo = 0, 0, 0
    pressure_hi = 2#2.458333333333333*pressure_lo
    u_hi, v_hi, w_hi = 0, 0, 0
    density_hi = density_lo*(pressure_hi/pressure_lo*(gamma+1)+(gamma-1)
                             )/(pressure_hi/pressure_lo*(gamma-1)+(gamma+1))
    velocity_hi = velocity_lo-(pressure_hi/pressure_lo-1)*numpy.sqrt(
        pressure_lo*gamma/density_lo)/numpy.sqrt(.5*gamma*(
            (gamma+1)*pressure_hi/pressure_lo+(gamma-1)))
    shock_function = sympy.functions.Heaviside((xi-shock_pos_0)-t*shock_speed(
            pressure_lo,density_lo,velocity_lo,a_lo,b_lo,c_lo,l_lo,m_lo,n_lo,
            p_lo,q_lo,r_lo,u_lo,gamma,pressure_hi,pm=-1))
#    print "pressure ratio = "+str(pressure_hi/pressure_lo)
#    print "density ratio = "+str(density_hi/density_lo)
#    print "mach hi = "+str(velocity_hi/numpy.sqrt(1.4*pressure_hi/density_hi))
#    print "mach lo = "+str(velocity_lo/numpy.sqrt(1.4*pressure_lo/density_lo))
    out = []
    out.append(pressure_lo+(pressure_hi-pressure_lo)*shock_function)
    out.append(density_lo+(density_hi-density_lo)*shock_function)
    out.append(velocity_lo+(velocity_hi-velocity_lo)*shock_function)
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(xi)
    out.append(eta)
    out.append(zeta)
    return out
def sym_sol():
    return (sympy.Symbol('p'),
            sympy.Symbol('rho'),
            sympy.Symbol('u'),
            sympy.Symbol('v'),
            sympy.Symbol('w'),
            sympy.Symbol('A'),
            sympy.Symbol('B'),
            sympy.Symbol('C'),
            sympy.Symbol('L'),
            sympy.Symbol('M'),
            sympy.Symbol('N'),
            sympy.Symbol('P'),
            sympy.Symbol('Q'),
            sympy.Symbol('R'),
            sympy.Symbol('U'),
            sympy.Symbol('V'),
            sympy.Symbol('W'),
            sympy.Symbol('x'),
            sympy.Symbol('y'),
            sympy.Symbol('z'))
def cons_sample((xi_in,eta_in,zeta_in),(t_in,func)):
    out = func.subs({t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}).evalf()
    clear_cache()
    return out
def flx1_sample((t_in,eta_in,zeta_in),(xi_in,func)):
    out = func.subs({t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}).evalf()
    clear_cache()
    return out
def flx2_sample((t_in,xi_in,zeta_in),(eta_in,func)):
    out = func.subs({t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}).evalf()
    clear_cache()
    return out
def flx3_sample((t_in,xi_in,eta_in),(zeta_in,func)):
    out = func.subs({t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}).evalf()
    clear_cache()
    return out
def quad_sample(x,y,z,t,func,mc_sample_func):
    return mc_sample_func((x,y,z),(t,func))
def int_eqn_sum(eqn_obj,sol,t_range,xi_range,eta_range,zeta_range):
    cons,flux1,flux2,flux3,source = eqn_obj(sol)
    calls = 10000
    cons_int_lst = [];flux1_int_lst = [];flux2_int_lst = [];flux3_int_lst = []
    for elem in cons:
#        left = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                            calls,args=(t_range[0],elem))
#        right = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                             calls,args=(t_range[1],elem))
        left = scipy.integrate.tplquad(
            quad_sample,xi_range[0],xi_range[1],lambda x:eta_range[0],
            lambda x:eta_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(t_range[0],elem,cons_sample))
        right = scipy.integrate.tplquad(
            quad_sample,xi_range[0],xi_range[1],lambda x:eta_range[0],
            lambda x:eta_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(t_range[1],elem,cons_sample))
        cons_int_lst.append(right[0]-left[0])
        print left,right
    print "Done with cons"
    for elem in flux1:
#        left = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                            calls,args=(xi_range[0],elem))
#        right = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                             calls,args=(xi_range[1],elem))
        left = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:eta_range[0],
            lambda x:eta_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(xi_range[0],elem,flx1_sample))
        right = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:eta_range[0],
            lambda x:eta_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(xi_range[1],elem,flx1_sample))
        flux1_int_lst.append(right[0]-left[0])
        print left,right
    print "Done with flx1"
    for elem in flux2:
#        left = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                            calls,(eta_range[0],elem))
#        right = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                             calls,(eta_range[1],elem))
        left = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:xi_range[0],
            lambda x:xi_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(eta_range[0],elem,flx2_sample))
        right = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:xi_range[0],
            lambda x:xi_range[1],lambda x,y:zeta_range[0],
            lambda x,y:zeta_range[1],args=(eta_range[1],elem,flx2_sample))
        flux2_int_lst.append(right[0]-left[0])
        print left,right
    print "Done with flx2"
    for elem in flux3:
#        left = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                            calls,(zeta_range[0],elem))
#        right = mc_integrate(flx3_sample,(t_range,xi_range,eta_range),
#                             calls,(zeta_range[1],elem))
        left = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:xi_range[0],
            lambda x:xi_range[1],lambda x,y:eta_range[0],
            lambda x,y:eta_range[1],args=(zeta_range[0],elem,flx3_sample))
        right = scipy.integrate.tplquad(
            quad_sample,t_range[0],t_range[1],lambda x:xi_range[0],
            lambda x:xi_range[1],lambda x,y:eta_range[0],
            lambda x,y:eta_range[1],args=(zeta_range[1],elem,flx3_sample))
        flux3_int_lst.append(right[0]-left[0])
        print left,right
    print "Done with flx3"
    source = [0 for elem in cons]
#    import pdb;pdb.set_trace()
    return (numpy.array(cons_int_lst)+numpy.array(flux1_int_lst)+
            numpy.array(flux2_int_lst)+numpy.array(flux3_int_lst)+source)
if __name__=="__main__":
    junk = Euler_UCS_manufactured.Euler_UCS()
    print int_eqn_sum(junk,shocked_sol_1(),(0,1),(0,1),(0,1),(0,1))
