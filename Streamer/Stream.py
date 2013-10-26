import os, sys
import numpy
import scipy.integrate
from functools import partial
#import Euler_UCS
libpath = os.path.abspath('/Users/woodscn/')
sys.path.insert(0,libpath)
from manufactured import Euler_UCS
Euler_UCS = Euler_UCS.Euler_UCS(Euler_UCS.MASA_solution_E())
class Stream(object):
    '''
    A self-contained, structured-grid flow. Includes flow data, initial 
    conditions, boundary conditions, equation sets, and any other necessary 
    options.
    '''
    import Godunov_driver
    Godunov=Godunov_driver.godunovdriver
    def __init__(self,boundary_conditions,initial_conditions,options):
        self.bounds = None
        self.main_data = None
        self.nt = 0
    def advance(self,dt,options):
        # Advance the conservation equations
        [self.main_data,dt_out] = self.Godunov.prim_update(
            numpy.asfortranarray(self.main_data),dt,cfl=.25,
            nx=self.main_data.shape[1]-2,ny=self.main_data.shape[2]-2,
            nz=self.main_data.shape[3]-2,options=options)
        # Advance the source ODEs
        t_array = numpy.zeros(2)
        t_array[1]=dt
        source=Euler_UCS
        import pdb;pdb.set_trace()
        source.balance_diff()
        for i in range(self.main_data.shape[1]):
            for j in range(self.main_data.shape[2]):
                for k in range(self.main_data.shape[3]):
                    ODE_func = partial(source.sample_point_diff,i,j,k)
                    out = scipy.integrate.odeint(
                        # What are the source terms?
                        source.sample_point_diff,
#                        manufactured.sample_point(
#                            source(smooth_sols_1()),nt,i,j,k),
                        self.main_data[:,i,j,k],t_array,
                        args=((0,-1*pow(10,-4),0),21))
                    self.main_data[:,i,j,k] = out[1,:]

        self.nt = self.nt + 1
        return dt_out
