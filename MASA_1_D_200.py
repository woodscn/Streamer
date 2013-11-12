import os, sys

import numpy

libpath = os.path.abspath('/Users/woodscn/')
sys.path.insert(0,libpath)
from manufactured import Euler_UCS
Euler_UCS = Euler_UCS.Euler_UCS(
    Euler_UCS.MASA_with_pinned_bounds([[0,100],[0,1],[0,1]])
    )
manufactured_source_function = Euler_UCS.balance_lambda_init()

class PatchInit:
    def __init__(self,t,bp,dim,bs=None,fs=None):
        self.type = t
        self.bounding_points = bp
        self.boundary_surface = bs
        self.flow_state = fs
        self.dim = dim
class BoundsInit:
    def __init__(self,lf,rf,bof,tf,baf,ff):
        self.left_face = lf
        self.right_face = rf
        self.bottom_face = bof
        self.top_face = tf
        self.back_face = baf
        self.front_face = ff
def init():
    left_face_init = []
    right_face_init = []
    top_face_init = []
    bottom_face_init = []
    back_face_init = []
    front_face_init = []
    nx = 200
    ny = 1
    nz = 1
    initial_conds = numpy.zeros((21,nx,ny,nz),order="F")
    source_funcs = numpy.zeros((nx,ny,nz),dtype=object)
    Dirichlet_pin = 0.01
    x_bcs = numpy.zeros((21,ny,nz))
    y_bcs = numpy.zeros((21,nx,nz))
    z_bcs = numpy.zeros((21,nx,ny))
    x_bcs[:5,:,:] += Dirichlet_pin
    x_bcs[5:14:4,:,:] = 1.
    y_bcs[:5,:,:] += Dirichlet_pin
    y_bcs[5:14:4,:,:] = 1.
    z_bcs[:5,:,:] += Dirichlet_pin
    z_bcs[5:14:4,:,:] = 1.
    x_bcs[-1],y_bcs[-1],z_bcs[-1] = 1.,1.,1.

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
#                out = Euler_UCS.sol.subs(
#                    dict(zip(Euler_UCS.vars_,(0,i,j,k)))).evalf()
                initial_conds[:-1,i,j,k] = Euler_UCS.sol.subs(dict(zip(
                            Euler_UCS.vars_,
                            (0.,float(i),float(j),float(k))))).evalf()[:]
                source_funcs[i,j,k] = (lambda y, t, i=i, j=j, k=k : numpy.array(
                        [item[0,0] for item in 
                         manufactured_source_function(t,i,j,k)]+[0]
                        ,dtype=numpy.float64))
    initial_conds[-1,:,:,:] = (initial_conds[5,:,:,:]*initial_conds[9,:,:,:]*
                               initial_conds[13,:,:,:])
    exact_solution = lambda t,x,y,z : (
        numpy.array(
            Euler_UCS.sol.subs(dict(zip(Euler_UCS.vars_,(t,x,y,z)))).evalf()[:]
            )       
        )
    solver_options = numpy.zeros(300)
#Options meanings
# [1]: controls which prim_update algorithm to use
# [101]: reports how many boundary ghost points are present
# [102]: controls spatial order of accuracy
# [103]: controls grid motion
# [104]: Controls type of time step (constant or CFL)
# [201-203]: same as [101-103]
    solver_options[0] = 1
    solver_options[100] = 1
    solver_options[101] = 1
    solver_options[102] = 0
    solver_options[103] = 0
    stream_options = {
        'solver_type':'euler',
        'boundary_layers':False,
        'multistream':False,
        'solver_options':solver_options,
        'manufactured':True,
        'source_funcs' : source_funcs,
        'exact_sol_func' : exact_solution
        }
    left_face_init.append(
        PatchInit('Dirichlet',
                  ((0.0,0.0,-.1),(0.0,0.0,0.1),(0.0,1.0,0.1),(0.0,1.0,-.1)),
                  1,'f = x',x_bcs))
    right_face_init.append(
        PatchInit('Dirichlet',
                  ((1.0,0.0,-.1),(0.0,1.0,-.1),(0.0,1.0,0.1),(0.0,0.0,0.1)),
                  1,'f = -x+1.0',x_bcs))
    bottom_face_init.append(
        PatchInit('Transmissive',#'Dirichlet',
                  ((0.0,0.0,-.1),(100.1,0.0,-.1),(100.1,0.0,0.1),(0.0,0.0,0.1)),
                  2,"f = y",y_bcs))
    top_face_init.append(
        PatchInit('Transmissive',#'Dirichlet',
#        PatchInit('Dirichlet',
                  ((0.0,1.0,-.1),(100.1,1.0,-.1),(100.1,1.0,0.1),(0.0,1.0,0.1)),
                  2,"f = 1.-y",y_bcs))
    back_face_init.append(
        PatchInit('Transmissive',#'Dirichlet',
#        PatchInit('Dirichlet',
                  ((0.0,0.0,-.1),(100.1,0.0,-.1),(100.1,1.0,-.1),(0.0,1.0,-.1)),
                  3,'f = z+.1',z_bcs))
    front_face_init.append(
        PatchInit('Transmissive',#'Dirichlet',
#        PatchInit('Dirichlet',
                  ((0.0,0.0,0.1),(100.1,0.0,0.1),(100.1,1.0,0.1),(0.0,1.0,0.1)),
                  3,'f = z-.1',z_bcs))
    bounds_init = BoundsInit(lf=left_face_init,rf=right_face_init,
                             bof=bottom_face_init,tf=top_face_init,
                             baf=back_face_init,ff=front_face_init)
    return bounds_init, initial_conds, stream_options

def MASA_generator():
    import numpy as np
    inputs = np.zeros((21,25,1))
    inputs[0,:,:] = 1.
    inputs[1,:,:] = 1.
    inputs[2,:,:] = 1.8*np.sqrt(1.4*inputs[0,:,:]/inputs[1,:,:])
    inputs[3,:,:] = 0.
    inputs[4,:,:] = 0.
    inputs[ 5,:,:] = 1./24.
    inputs[ 6,:,:] = 0.
    inputs[ 7,:,:] = 0.
    inputs[ 8,:,:] = 0.
    inputs[ 9,:,:] = 1./24.
    inputs[10,:,:] = 0.
    inputs[11,:,:] = 0.
    inputs[12,:,:] = 0.
    inputs[13,:,:] = 1.
    inputs[14,:,:] = .25*inputs[2,:,:]
    inputs[15,:,:] = 0.
    inputs[16,:,:] = 0.
    inputs[17,:,:] = 0.
    for inda in range(inputs.shape[1]):
        for indb in range(inputs.shape[2]):
            inputs[18,inda,indb] = 1./25.*(inda+.5)
            inputs[19,inda,indb] = 0.
    inputs[20,:,:] = 1./24./24.
    # np.savetxt("transonic_duct_bounds_x.txt",np.reshape(inputs[17,:,:],-1))
    # np.savetxt("transonic_duct_bounds_y.txt",np.reshape(inputs[18,:,:],-1))
    # np.savetxt("transonic_duct_bounds_z.txt",np.reshape(inputs[19,:,:],-1))
    # np.save("transonic_duct_bounds",inputs)
    return(inputs)

if __name__=='__main__':
    test = init()
