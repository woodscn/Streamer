import numpy
import sympy
import scipy.integrate
import Euler_UCS_manufactured
import multiprocessing
from sympy.core.cache import * # Needed for clear_cache

import myquad
from shocked_solutions import shocked_sol_1
from smooth_solutions import smooth_sols_1
from symbolic_solution import sym_sol

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
def int_eqn_sum(eqn_obj,sol,t_range,xi_range,eta_range,zeta_range,shocks=()):
    cons,flux1,flux2,flux3,source = eqn_obj(sol)
    calls = 10000
    cons_int_lst = [];flux1_int_lst = [];flux2_int_lst = [];flux3_int_lst = []
#    for elem in cons:
#        left = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                            calls,args=(t_range[0],elem))
#        right = mc_integrate(cons_sample,(xi_range,eta_range,zeta_range),
#                             calls,args=(t_range[1],elem))
#        left = scipy.integrate.tplquad(
#            quad_sample,xi_range[0],xi_range[1],lambda x:eta_range[0],
#            lambda x:eta_range[1],lambda x,y:zeta_range[0],
#            lambda x,y:zeta_range[1],args=(t_range[0],elem,cons_sample))
#        right = scipy.integrate.tplquad(
#            quad_sample,xi_range[0],xi_range[1],lambda x:eta_range[0],
#            lambda x:eta_range[1],lambda x,y:zeta_range[0],
#            lambda x,y:zeta_range[1],args=(t_range[1],elem,cons_sample))
#        cons_int_lst.append(right[0]-left[0])
#        print left,right
#    print "Done with cons"
#    for elem in flux1:
#        left = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                            calls,args=(xi_range[0],elem))
#        right = mc_integrate(flx1_sample,(t_range,eta_range,zeta_range),
#                             calls,args=(xi_range[1],elem))
#        left = scipy.integrate.tplquad(
#            quad_sample,t_range[0],t_range[1],lambda x:eta_range[0],
#            lambda x:eta_range[1],lambda x,y:zeta_range[0],
#            lambda x,y:zeta_range[1],args=(xi_range[0],elem,flx1_sample))
#        right = scipy.integrate.tplquad(
#            quad_sample,t_range[0],t_range[1],lambda x:eta_range[0],
#            lambda x:eta_range[1],lambda x,y:zeta_range[0],
#            lambda x,y:zeta_range[1],args=(xi_range[1],elem,flx1_sample))
#        flux1_int_lst.append(right[0]-left[0])
#        print left,right
#    print "Done with flx1"
    for elem in flux2:
#        left = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                            calls,(eta_range[0],elem))
#        right = mc_integrate(flx2_sample,(t_range,xi_range,zeta_range),
#                             calls,(eta_range[1],elem))
        points=[]
# This doesn't work if the shock point is outside the range of integration.
        for shock in shocks:
            points.append(
                float(
                    sympy.solve(shocks[0].subs({xi:.5}),sympy.Symbol('t'))[0]))
        import pdb;pdb.set_trace()
        test = scipy.integrate.quad(quad_sample,t_range[0],t_range[1],
                                    args=(.5,.5,.5,elem,flx2_sample),
                                    points=points)
        pdb.set_trace()
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
    sol=shocked_sol_1()
    print int_eqn_sum(junk,sol[0],(0,1),(0,1),(0,1),(0,1),shocks=sol[1])
