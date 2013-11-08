import numpy
import sympy

t=sympy.Symbol('t')
xi=sympy.Symbol('xi')
eta=sympy.Symbol('eta')
zeta=sympy.Symbol('zeta')

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
