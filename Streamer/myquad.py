"""
Recursive integration using SciPy's quad function.

Contains the following:

class Error - module wrapper for Exception. For future expansion only.

function mul_quad - Evaluate multiple integrals using recursive calls to quad.
"""
from scipy.integrate import quad
import numpy
class Error(Exception):
    pass
def mul_quad(func,ranges,args=([],[]),opts=(),depth=0):
    """
    Evaluate multiple integrals through recursive calls to scipy.integrate.quad.
    
    mul_quad(func,ranges,args=([],[]),opts=(),depth=0)

    Inputs:
      func - callable, acceptable to SciPy's quad, returning a number. 
        Should accept a float, followed by the contents of args[0] and 
        args[1], e.g. if args=([1,1,3],[1.4,string]), then func should 
        be of the form func(x,1,1,3,1.4,string).
      ranges - sequence of sequences describing the ranges of integration. 
        Integrals are performed in order, so ranges[0] corresponds to the 
        first argument of func, ranges[1] to the second, and so on.
      args - optional sequences of arguments. args[0] should be left empty, 
        to be used in the recursive construction of new functions, while 
        args[1] contains any additional arguments required by func 
        beyond those over which integration is being performed, similarly
        to how args functions in quad.
      opts - optional sequence of dictionaries. Each dictionary contains any 
        optional arguments to be passed to quad. Options are passed in order.
        If omitted, default options are used:
        - full_output = 0
        - epsabs      = 1.49e-08
        - epsrel      = 1.49e-08
        - limit       = 50
        - points      = None
        - weight      = None
        - wvar        = None
        - wopts       = None
        Currently, only points is fully supported. The user must specify
        appropriate functions for points, rather than specific values.
      depth - used to determine level of integration. Should be omitted
        by the user, except for debugging purposes.
        
    Returns:
      out - value of multiple integral

    mul_quad takes care of the programming required to perform nested
    integration using the 1-dimensional integration routines provided
    by SciPy's adaptive quadrature function. It extends the capabilities
    of dblquad and tplquad by allowing for deeper integration, up to
    six nested integrals. It also allows the user to specify the full
    range of options allowed by quad, independently, for each level of
    integration. Support is not provided for functional dependence of 
    limits of integration as in dblquad and tplquad, so limits must be
    specified as constants. Users are cautioned that nested Gaussian 
    integration of this kind is computationally intensive, and may be 
    unsuitable for many nested integrals. 

    TODO: 
      - implement better support for options other than points.
      - implement limits of integration as functions, as in dblquad.
      - simplify user input and make interface more robust.

    """
    if len(ranges)>6:
        raise Error(
            'greater than six-dimensional integration not supported!')
    ind = -depth-1
    # Set to default options if not given
    opts_default={'full_output':0,'epsabs':1.49e-08,'epsrel':1.49e-08,
                  'limit':50,'points':None,'weight':None,'wvar':None,
                  'wopts':None}
    for k in (opts_default.keys()):
        try:
            exec(k+'=opts[ind][k]')
        except(KeyError,IndexError):
            exec(k+'=opts_default[k]')
    current_range = ranges[ind]
    try:
        points_sample=[point(args[0]) for point in points]
#        import pdb;pdb.set_trace()
    except(TypeError):
        points_sample=None
    if ranges[0] is ranges[ind]:# Check to see if down to 1-D integral
        out = quad(func,current_range[0],current_range[1],
                   args=tuple(sum(args,[])),full_output=full_output,
                   epsabs=epsabs,epsrel=epsrel,limit=limit,points=points_sample,
                   weight=weight,wvar=wvar,wopts=wopts)
    else:
# Generate a new integrand. This is a recursive operation.
        vars_base=['x0','x1','x2','x3','x4','x5','x6']
        vars_str = [vars_base.pop(0) for n in range(depth+1)]
# Create a new function, written as a Python string.
        func_new_str = (
            'lambda '+','.join(vars_str)+
            ',func=func,ranges=ranges,args=args,opts=opts,depth=depth'+
            ':mul_quad(func,ranges,args=(['+','.join(vars_str)+'],'+
            str(args[1])+'),opts=opts,depth=depth+1)')
        func_new = eval(func_new_str)
# Integrate the new function.
        out = quad(func_new,current_range[0],current_range[1],
                   args=tuple(sum(args,[])),full_output=full_output,
                   epsabs=epsabs,epsrel=epsrel,limit=limit,points=points_sample,
                   weight=weight,wvar=wvar,wopts=wopts)
    return out[0]

if __name__=='__main__':
    func = lambda x0,x1,x2,x3 : x0**2+x1*x2-x3**3+numpy.sin(x0)+(
        1 if (x0-.2*x3-.5-.25*x1>0) else 0)
    points=[[lambda (x1,x2,x3) : .2*x3+.5+.25*x1],[],[],[]]
    print mul_quad(func,[[0,1],[-1,1],[.13,.8],[-.15,1]],
                   opts=[{'points':points[0]},{},{},{}])
