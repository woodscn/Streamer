from sympy import *

from sympy.utilities.codegen import codegen

from numpy import array
def normalvec(functionStr,x_in,y_in,z_in,out):
    x=Symbol('x')
    y=Symbol('y')
    z=Symbol('z')
    out[0]=diff(functionStr,x).evalf(subs={x:x_in,y:y_in,z:z_in})
    out[1]=diff(functionStr,y).evalf(subs={x:x_in,y:y_in,z:z_in})
    out[2]=diff(functionStr,z).evalf(subs={x:x_in,y:y_in,z:z_in})
    return 0#array([diff(functionStr,x),diff(functionStr,y),diff(functionStr,z)])

if __name__ == "__main__":
    str = 'x**2'
    out = array([0,0,0])
    args = Symbol('x'),Symbol('y'),Symbol('z'),out
    print codegen(('tester',normalvec(str,,out)),'F95','junk',argument_sequence=args
)
