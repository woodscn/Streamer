import numpy
import scipy.integrate

sys.path.insert(0,'/Users/woodscn/')
from manufactured import Euler_UCS
Euler_UCS = Euler_UCS.Euler_UCS(Euler_UCS.MASA_solution_E())

def advance_source_ODE():
    t_array = numpy.zeros(2)
    t_array[1] = dt
    source_sym = Euler_UCS.balance_diff()
    for i in range(self.main_data.shape[1]):
        for j in range(self.main_data.shape[2]):
            for k in range(self.main_data.shape[3]):
                ODE_func = partial(source.sample_point_diff,i,j,k)
                out = scipy.integrate.odeint(
                    # What are the source terms?
                    source.sample_point_diff,
                    self.main_data[:,i,j,k],t_array,
                    args=((0,-1*pow(10,-4),0),21))
                self.main_data[:,i,j,k] = out[1,:]
