import numpy
import sympy

t=sympy.Symbol('t')
xi=sympy.Symbol('xi')
eta=sympy.Symbol('eta')
zeta=sympy.Symbol('zeta')

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
