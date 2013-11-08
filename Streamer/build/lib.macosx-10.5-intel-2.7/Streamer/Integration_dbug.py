import numpy
import sympy
from sympy.core.cache import *
t=sympy.Symbol('t')
xi=sympy.Symbol('xi')
eta=sympy.Symbol('eta')
zeta=sympy.Symbol('zeta')
def mc_integrate(f,ranges,calls,args=()):
    f_sum = 0
    for n in xrange(calls):
        coords_lst = []
        for rng in ranges:
            coords_lst.append(rng[0] + rng[1]*numpy.random.random())
        f_sum += f(coords_lst,args)
        del coords_lst
    vol = numpy.prod([float(rng[1])-float(rng[0]) for rng in ranges])
    return vol/calls*f_sum    
def test_func():
    return xi**2-4*eta*sympy.sqrt(xi)+zeta*xi*eta**2+t*zeta
def test_func_sample((xi_in,eta_in,zeta_in),(t_in,func)):
    return func().subs({t:t_in,xi:xi_in,eta:eta_in,zeta:zeta_in}).evalf()
if __name__=="__main__":
    out=[]
    for n in range(100):
        calls = 10000
        t_range = (0,1)
        xi_range = (0,1)
        eta_range = (0,1)
        zeta_range = (0,1)
        left = mc_integrate(test_func_sample,(xi_range,eta_range,zeta_range),
                            calls,args=(t_range[0],test_func))
        right = mc_integrate(test_func_sample,(xi_range,eta_range,zeta_range),
                             calls,args=(t_range[1],test_func))
        out.append(right-left)
        clear_cache()
        print "n = ",n,right-left
