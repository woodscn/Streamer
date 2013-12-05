import os, sys

import numpy

libpath = os.path.abspath('/Users/woodscn/')
sys.path.insert(0,libpath)
from manufactured import Euler_UCS
xmin,xmax,ymin,ymax,zmin,zmax = 0,100,0,1,0,1
nx = 24
ny = 2
nz = 2
nxis = [nx,ny,nz]
dxis = [4,1,1]
#Euler_UCS = Euler_UCS.Euler_UCS(
#    Euler_UCS.MASA_with_pinned_bounds(
#        ranges=[[xmin,xmax],[ymin,ymax],[zmin,zmax]],nxes=(nx,ny,nz),dxis=dxis))
MMS = Euler_UCS.MMS_solution([[xmin,xmax],[ymin,ymax],[zmin,zmax]],nxis,dxis)
index_to_xi = MMS['MMS_obj'].index_to_xi
Euler_UCS = Euler_UCS.Euler_UCS(MMS)
    
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
    initial_conds = numpy.zeros((21,nx,ny,nz),order="F")
    source_funcs = numpy.zeros((nx,ny,nz),dtype=object)
    Dirichlet_pin = 0.01
    x_min_bcs = numpy.zeros((21,ny,nz))
    x_max_bcs = numpy.zeros((21,ny,nz))
    for inda in range(ny):
        for indb in range(nz):
            xi_in = index_to_xi(0,inda,indb)
            x_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]
            xi_in = index_to_xi(nx+1,inda,indb)
            x_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]
    y_min_bcs = numpy.zeros((21,nx,nz))
    y_max_bcs = numpy.zeros((21,nx,nz))
    for inda in range(nx):
        for indb in range(nz):
            xi_in = index_to_xi(inda,0,indb)
            y_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]
            xi_in = index_to_xi(inda,ny+1,indb)
            y_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]
    z_min_bcs = numpy.zeros((21,nx,ny))
    z_max_bcs = numpy.zeros((21,nx,ny))
    for inda in range(nx):
        for indb in range(ny):
            xi_in = index_to_xi(inda,indb,0)
            z_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]
            xi_in = index_to_xi(inda,indb,nz+1)
            z_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xi_in[0],xi_in[1],xi_in[2])))
                                                          ).evalf()[:]

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xi_in = index_to_xi(i+1,j+1,k+1)
                initial_conds[:-1,i,j,k] = Euler_UCS.sol.subs(dict(zip(
                            Euler_UCS.vars_,
                            (0.,xi_in[0],xi_in[1],xi_in[2])))).evalf()[:]
                source_funcs[i,j,k] = (lambda y, t, i=(xi_in[0]), 
                                       j=(xi_in[1]), 
                                       k=(xi_in[2]) : numpy.array(
                        [item[0,0] for item in 
                         manufactured_source_function(t,i,j,k)]+[0]
                        ,dtype=numpy.float64))
    initial_conds[-1,:,:,:] = (initial_conds[5,:,:,:]*initial_conds[9,:,:,:]*
                               initial_conds[13,:,:,:])
    def exact_solution(t,i,j,k):
        xi_in = index_to_xi(i,j,k)
        return numpy.array(
            Euler_UCS.sol.subs(dict(
                    zip(Euler_UCS.vars_,(t,xi_in[0],xi_in[1],xi_in[2])))
            ).evalf()[:])
#    exact_solution = lambda t,x,y,z : (
#        numpy.array(
#            Euler_UCS.sol.subs(dict(zip(Euler_UCS.vars_,(t,dxis[0]*x,
#                                                         dxis[1]*y,dxis[2]*z
#                                                         )))).evalf()[:]
#            )       
#        )
    solver_options = numpy.zeros(300)
#Options meanings
# [1]: controls which prim_update algorithm to use
# [101]: reports how many boundary ghost points are present
# [102]: controls spatial order of accuracy
# [103]: controls grid motion
# [104]: Controls type of time step (constant or CFL)
# [201-203]: same as [101-103]
    solver_options[0] = 1
    solver_options[2] = 6
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
                  ((0.0,-.1,-.1),(0.0,-.1,1.1),(0.0,1.1,1.1),(0.0,1.1,-.1)),
                  1,'f = x',x_min_bcs))
    right_face_init.append(
        PatchInit('Dirichlet',
                  ((100.0,-.1,-.1),(100.0,1.1,-.1),
                   (100.0,1.1,1.1),(100.0,-.1,1.1)),
                  1,'f = -x+100.0',x_max_bcs))
    bottom_face_init.append(
        PatchInit('Transmissive',
                  ((-.1,0.0,-.1),(100.1,0.0,-.1),(100.1,0.0,1.1),(-0.1,0.0,1.1)),
                  2,"f = y",y_min_bcs))
    top_face_init.append(
        PatchInit('Transmissive',
                  ((-.1,1.0,-.1),(100.1,1.0,-.1),(100.1,1.0,1.1),(-.1,1.0,1.1)),
                  2,"f = 1.-y",y_max_bcs))
    back_face_init.append(
        PatchInit('Transmissive',
                  ((-.1,-.1,0.),(100.1,-0.1,0.),(100.1,1.1,0.),(-.1,1.1,0.)),
                  3,'f = z+1',z_min_bcs))
    front_face_init.append(
        PatchInit('Transmissive',
                  ((-.1,-.1,1.),(100.1,-.1,1.),(100.1,1.1,1.),(-0.1,1.1,1.)),
                  3,'f = z-1',z_max_bcs))
    bounds_init = BoundsInit(lf=left_face_init,rf=right_face_init,
                             bof=bottom_face_init,tf=top_face_init,
                             baf=back_face_init,ff=front_face_init)
    import pdb;pdb.set_trace()
    return bounds_init, initial_conds, stream_options

if __name__=='__main__':
    test = init()
