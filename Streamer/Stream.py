import numpy
import scipy.integrate
import manufactured
import Euler_UCS_manufactured
#    gravity = Source_functions.source_functions.gravity
#import Godunov_driver
#Godunov = Godunov_driver.godunovdriver
class Stream(object):
    '''
    A self-contained, structured-grid flow. Includes flow data, initial 
    conditions, boundary conditions, equation sets, and any other necessary 
    options.
    '''
    import Godunov_driver
#    stream_dict = {Godunov:Godunov_driver.godunovdriver}
    Godunov=Godunov_driver.godunovdriver
#    Sources={gravity:
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
        source=Euler_UCS_manufactured.Euler_UCS(
            Euler_UCS_manufactured.MASA_solution_E())
        source.balance_diff()
        for i in range(self.main_data.shape[1]):
            for j in range(self.main_data.shape[2]):
                for k in range(self.main_data.shape[3]):
                    out = scipy.integrate.odeint(
                        # What are the source terms?
                        source.sample_point_diff(nt,i,j,k),
#                        manufactured.sample_point(
#                            source(smooth_sols_1()),nt,i,j,k),
                        self.main_data[:,i,j,k],t_array,
                        args=((0,-1*pow(10,-4),0),21))
                    self.main_data[:,i,j,k] = out[1,:]
        self.nt = self.nt + 1
        return dt_out
#def Riemann_init_1D_1:
#    nx = 100
#    dx = 1/(nx-1.)
#    temp = numpy.zeros([21,nx],dtype=numpy.float64,order='F')
#    temp[0] = 1
#    temp[1] = 1
#    temp[2] = .75
#    temp[0,30:] = .1
#    temp[1,30:] = .125
#    temp[2,30:] = 0
#    temp[5:14:4] = dx
#    for n in range(nx):
#        temp[17,n] = dx*n
#    temp[20] = dx*dx*dx
#    return temp
#
#if __name__=='__main__':
