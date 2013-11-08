from scipy.integrate import quad
import numpy
class NQuadError(Exception):
    pass
def nquad(func,ranges,args=(),opts=[]):
    """
Integration over multiple variables.

Wraps scipy.integrate.quad to enable integration over multiple variables.
Various options allow improved integration of discontinuous functions, as 
well as the use of weighted integration, and generally finer control of the
integration process.

Parameters
----------
func : callable object
    The function to be integrated. Has arguments of x0, ... xn, t0, ... tn,
    where integration is carried out over x0, ... xn, which must be floats.
    Function signature should be func(x0,x1,...xn,t0,t1,...tn). Integration
    is performed in order: x0, x1, ... xn.
ranges : iterable object
    Each element of ranges may be either a sequence 2 numbers, or else a 
    callable that returns such a sequence. ranges[0] corresponds to 
    integration over x0, and so on. If an element of ranges is a callable,
    then it will be called with all of the integration arguments available.
    e.g. if func = f(x0,x1,x2), then ranges[0] may be defined as either
    (a,b) or else as (a,b) = range0(x1,x2). 
args : iterable object, optional
    Additional arguments t0, ... tn, required by func.
opts : iterable object or dict, optional
    Options to be passed to scipy.integrate.quad. May be empty, a dict, or
    a sequence of dicts or functions that return a dict. If empty, the 
    default options from scipy.integrate.quadare used. If a dict, the same 
    options are used for all levels of integraion. If a sequence, then each
    element of the sequence corresponds to a particular integration. e.g. 
    opts[0] corresponds to integration over x0, and so on. The available 
    options together with their default values (as of Apr 2013) are:
    - epsabs = 1.49e-08
    - epsrel = 1.49e-08
    - limit  = 50
    - points = None
    - weight = None
    - wvar   = None
    - wopts  = None
    The full_output option from quad is unavailable, due to the complexity
    of handling the large amount of data such an option would return for 
    this kind of nested integration. For more information on these options,
    consult the documentation for scipy.integrate.quad

Returns
-------
result : float
    The result of the integration
abserr : float
    The maximum of the estimates of the absolute error in the various 
    integration results

See Also
--------
scipy.integrate.quad : 1-dimensional Gaussian integration.

Examples
--------
>>> func = lambda x0,x1,x2,x3 : x0**2+x1*x2-x3**3+numpy.sin(x0)+(
        1 if (x0-.2*x3-.5-.25*x1>0) else 0)
>>> points=[[lambda (x1,x2,x3) : .2*x3+.5+.25*x1],[],[],[]]
>>> def opts0(*args,**kwargs):
        return {'points':[.2*args[2]+.5+.25*args[0]]}
>>> nquad(func,[[0,1],[-1,1],[.13,.8],[-.15,1]],opts=[opts0,{},{},{}])
(1.5267454070738635, 2.943736000140233e-14)

>>> scale = .1
>>> def func(x0,x1,x2,x3,t0,t1):
        return x0*x1*x3**2+numpy.sin(x2)+1+(1 if x0+t1*x1-t0>0 else 0)
>>> def lim0(x1,x2,x3,t0,t1):
        return [scale*(x1**2+x2+numpy.cos(x3)*t0*t1+1)-1, 
                scale*(x1**2+x2+numpy.cos(x3)*t0*t1+1)+1]
>>> def lim1(x2,x3,t0,t1):
        return [scale*(t0*x2+t1*x3)-1,
                scale*(t0*x2+t1*x3)+1]
>>> def lim2(x3,t0,t1):
        return [scale*(x3+t0**2*t1**3)-1, 
                scale*(x3+t0**2*t1**3)+1]
>>> def lim3(t0,t1):
        return [scale*(t0+t1)-1,
                scale*(t0+t1)+1]
>>> def opts0(x1,x2,x3,t0,t1):
        return {'points':[t0-t1*x1]}
>>> def opts1(x2,x3,t0,t1):
        return {}
>>> def opts2(x3,t0,t1):
        return {}
>>> def opts3(t0,t1):
        return {}
>>> nquad(func,[lim0,lim1,lim2,lim3],args=(0,0),opts=[opts0,opts1,opts2,opts3])

"""
    new_ranges = [
        range_ if callable(range_) else RangeFunc(range_) for range_ in ranges]
    if isinstance(opts,dict):
        new_opts = [opts for ind in range(len(ranges))]
    else:
        if len(opts) == 0:
            new_opts = opts
        else:
            new_opts = [
                opt if callable(opt) else OptFunc(opt) for opt in opts]
    return _NQuad(func,new_ranges,args,new_opts).integrate()
class RangeFunc(object):
    def __init__(self,range_):
        self.range_ = range_
    def __call__(self,*args):
        return self.range_
class OptFunc(object):
    def __init__(self,opt):
        self.opt = opt
    def __call__(self,*args):
        return self.opt
class _NQuad(object):
    def __init__(self,func,ranges,args,opts):
        self.abserr = 0
        self.func = func
        self.ranges = ranges
        self.args = args
        self.opts = opts
    def integrate(self):
        args_and_depth = self.args+(0,)
        return (self._int(*args_and_depth),self.abserr)
    def _int(self,*args):
        depth = args[-1]
        ind = -depth-1
        range_ = self.ranges[ind](*args[0:-1])
        opt = self.opts[ind](*args[0:-1])
        try:
            for point in opt["points"]:
                if point < range_[0] or point > range_[1]:
                    opt["points"].remove(point)
        except(KeyError):
            pass
        if self.ranges[ind] is self.ranges[0]:
            newfunc = self.func
            newargs = args[0:-1]
        else:
            def newfunc(*args):
                return self._int(*args)
            newargs = list(args)
            newargs[-1] = newargs[-1] + 1
        out = quad(newfunc,*range_,args=tuple(newargs),**opt)
        self.abserr = max(self.abserr,out[1])
        return out[0]

if __name__=='__main__':
    func = lambda x0,x1,x2,x3 : x0**2+x1*x2-x3**3+numpy.sin(x0)+(
        1 if (x0-.2*x3-.5-.25*x1>0) else 0)
    points=[[lambda (x1,x2,x3) : .2*x3+.5+.25*x1],[],[],[]]
    def opts0(*args,**kwargs):
        return {'points':[.2*args[2]+.5+.25*args[0]]}
    print nquad(func,[[0,1],[-1,1],[.13,.8],[-.15,1]],opts=[opts0,{},{},{}])

    scale = .1
    def func(x0,x1,x2,x3,t0,t1):
        return x0*x1*x3**2+numpy.sin(x2)+1+(1 if x0+t1*x1-t0>0 else 0)
    def lim0(x1,x2,x3,t0,t1):
        return [scale*(x1**2+x2+numpy.cos(x3)*t0*t1+1)-1, 
                scale*(x1**2+x2+numpy.cos(x3)*t0*t1+1)+1]
    def lim1(x2,x3,t0,t1):
        return [scale*(t0*x2+t1*x3)-1,
                scale*(t0*x2+t1*x3)+1]
    def lim2(x3,t0,t1):
        return [scale*(x3+t0**2*t1**3)-1, 
                scale*(x3+t0**2*t1**3)+1]
    def lim3(t0,t1):
        return [scale*(t0+t1)-1,
                scale*(t0+t1)+1]
    def opts0(x1,x2,x3,t0,t1):
        return {'points':[t0-t1*x1]}
    def opts1(x2,x3,t0,t1):
        return {}
    def opts2(x3,t0,t1):
        return {}
    def opts3(t0,t1):
        return {}
    print nquad(func,[lim0,lim1,lim2,lim3],args=(0,0),
                opts=[opts0,opts1,opts2,opts3])
