import numpy
import scipy.integrate
import Source_functions
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
    def advance(self,dt,options):
        [self.main_data,dt_out] = self.Godunov.prim_update(
            numpy.asfortranarray(self.main_data),dt,cfl=.25,
            nx=self.main_data.shape[1]-2,ny=self.main_data.shape[2]-2,
            nz=self.main_data.shape[3]-2,options=options)
        for i in range(self.main_data.shape[1]):
            for j in range(self.main_data.shape[2]):
                for k in range(self.main_data.shape[3]):
                    out = scipy.integrate.odeint(Source_functions.source_functions.gravity,self.main_data[:,i,j,k],(0,dt),args=((0,1,0),3))
                    self.main_data[:,i,j,k] = out[1,:]
                    import pdb;pdb.set_trace()
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
