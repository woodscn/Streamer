from scipy.integrate import quad
import numpy
class Error(Exception):
    pass
def mul_quad(func,ranges,args=([],[]),opts=(),depth=0):
    if len(ranges)>6:
        raise Error('myquad does not currently support more than six-dimensional integrals!')
    opts_default={'full_output':0,'epsabs':1.49e-08,'epsrel':1.49e-08,
                  'limit':50,'points':None,'weight':None,'wvar':None,
                  'wopts':None}
# Set to default options if not given
    for k in (opts_default.keys()):
        try:
            exec(k+'=opts[0][k]')
        except(KeyError,IndexError):
            exec(k+'=opts_default[k]')

    current_range = ranges[depth]
    if ranges[-1] is ranges[depth]:# Check to see if down to 1-D integral
        out = quad(func,current_range[0],current_range[1],
                   args=tuple(sum(args,[])),full_output=full_output,
                   epsabs=epsabs,epsrel=epsrel,limit=limit,points=points,
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
                   epsabs=epsabs,epsrel=epsrel,limit=limit,points=points,
                   weight=weight,wvar=wvar,wopts=wopts)
    return out[0]
if __name__=='__main__':
    func = lambda x0,x1,x2,x3 : x0**2+x1*x2-x3**3+numpy.sin(x0)+1 if x0-.0*x3-.5>0 else 0
    points=[[],[],[],[]]
    print mul_quad(func,[[0,1],[0,1],[0,1],[0,1]])
