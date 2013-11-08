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
def mul_quad(func,ranges,args=(),opts=(),_depth=0,_int_vars=()):
    """
    Evaluate multiple integrals through recursive calls to scipy.integrate.quad.
    
    mul_quad takes care of the programming required to perform nested
    integration using the 1-dimensional integration routines provided
    by SciPy's adaptive quadrature function. It extends the capabilities
    of dblquad and tplquad by allowing for more levels of integration, 
    and allowing the user to specify nearly the full range of options 
    allowed by quad, for each level of integration. Users are cautioned 
    that nested Gaussian integration of this kind is computationally 
    intensive, and may be unsuitable for many nested integrals. 

    Usage: mul_quad(func,ranges,args=(),opts=(),_depth=0,_int_vars=())

    Inputs:
      func - callable, acceptable to SciPy's quad, returning a number. 
        Should accept a float, followed by the contents of _int_vars and 
        args, e.g. if x is a float, args=(1.4,string), and _int_vars = 
        (1,1,3), then func should be of the form 
        func(x,1,1,3,1.4,string).
      ranges - sequence describing the ranges of integration. Integrals
        are performed in order, so ranges[0] corresponds to the first 
        argument of func, ranges[1] to the second, and so on. Each 
        element of ranges may be either a constant sequence of length 2 
        or else a function that returns such a sequence. If a function, 
        then it will be called with all of the integration arguments 
        available to that point. e.g. for func = f(x0,x1,x2,x3), the 
        range of integration for x0 may be defined as either a constant 
        such as (0,1) or as a function range0(x1,x2,x3). The functional 
        range of integration for x1 will be range1(x2,x3), x2 will be 
        range2(x3), and so on.
      args - optional sequence of arguments. Contains only arguments of 
        func beyond those over which the integration is being performed.
      opts - optional sequence of options for SciPy's quad. Each element 
        of opts may be specified as either a dictionary or as a function 
        that returns a dictionary similarly to ranges. opts must either 
        be left empty (), or it must be the same length as ranges. 
        Options are passed in the same order as ranges, so opts[0] 
        corresponds to integration over the first argument of func, and 
        so on. The full_output option from quad is not available, due to 
        the difficulty of consolidating the large number of additional 
        outputs. For reference, the default options from quad are:
        - epsabs      = 1.49e-08
        - epsrel      = 1.49e-08
        - limit       = 50
        - points      = None
        - weight      = None
        - wvar        = None
        - wopts       = None
        (As of Apr 2013)
      _depth - used to determine level of integration. Should be omitted
        by the user, except for debugging purposes.
      _int_vars - contains values of integration variables in inner 
        integration loops. Should not be used manually except for
        debugging.
        
    Returns:
      out - value of multiple integral in the specified range.
      abserr - estimate of the absolute error in the result. The 
        maximum value of abserr among all the SciPy quad evaluations.

    """
    global abserr
    if _depth == 0:
        abserr = None
    if not (len(opts) in [0,len(ranges)]):
        raise Error('opts must be given for all integration levels or none!')
    total_args = _int_vars+args
    # Select the range and opts for the given depth of integration.
    ind = -_depth-1
    if callable(ranges[ind]):
        current_range = ranges[ind](*total_args)
    else:
        current_range = ranges[ind]
    if len(opts) != 0:
        if callable(opts[ind]):
            current_opts = opts[ind](*total_args)
        else:
            current_opts = opts[ind]
    else:
        current_opts = {}
    try:
        if current_opts["full_output"] != 0:
            raise Error('full_output option is disabled!')
    except(KeyError):
        pass
    if current_range is ranges[0]:# Check to see if down to 1-D integral
        func_new = func
    else:
        # Define a new integrand.
        def func_new(*_int_vars):
            return mul_quad(func,ranges,args=args,opts=opts,
                            _depth=_depth+1,_int_vars=_int_vars)
    out = quad(func_new,*current_range,args=_int_vars+args,**current_opts)
    if abserr is None:
        abserr = out[1]
    if out[1] > abserr:
        abserr = out[1]
    if _depth == 0:
        return out[0],abserr
    else:
        return out[0]

if __name__=='__main__':
    func = lambda x0,x1,x2,x3 : x0**2+x1*x2-x3**3+numpy.sin(x0)+(
        1 if (x0-.2*x3-.5-.25*x1>0) else 0)
    points=[[lambda (x1,x2,x3) : .2*x3+.5+.25*x1],[],[],[]]
    def opts0(*args,**kwargs):
        return {'points':[.2*args[2]+.5+.25*args[0]]}
    print mul_quad(func,[[0,1],[-1,1],[.13,.8],[-.15,1]],
                   opts=[opts0,{},{},{}])
