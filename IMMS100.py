import os, sys

import numpy

libpath = os.path.abspath('/Users/woodscn/')
sys.path.insert(0,libpath)
libpath = os.path.abspath('/home/woodscn/')
sys.path.insert(0,libpath)
from manufactured import Euler_UCS
xmin,xmax,ymin,ymax,zmin,zmax = 0,100,0,1,0,1
nx = 101
ny = 1
nz = 1
dxis = [1,1,1]
Euler_UCS = Euler_UCS.Euler_UCS(
    Euler_UCS.MASA_solution_full(
        ranges=[[xmin,xmax],[ymin,ymax],[zmin,zmax]],nxes=(nx,ny,nz),dxis=dxis,
        disc=True))
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
            x_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xmin,
                                         dxis[1]*float(inda),
                                         dxis[2]*float(indb))))
                                                          ).evalf()[:]
            x_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,xmax,
                                         dxis[1]*float(inda),
                                         dxis[2]*float(indb))))
                                                          ).evalf()[:]
    y_min_bcs = numpy.zeros((21,nx,nz))
    y_max_bcs = numpy.zeros((21,nx,nz))
    for inda in range(nx):
        for indb in range(nz):
            y_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,dxis[0]*float(inda),
                                         ymin,dxis[2]*float(indb))))
                                                          ).evalf()[:]
            y_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,dxis[0]*float(inda),
                                         ymax,dxis[1]*float(indb))))
                                                          ).evalf()[:]
    z_min_bcs = numpy.zeros((21,nx,ny))
    z_max_bcs = numpy.zeros((21,nx,ny))
    for inda in range(nx):
        for indb in range(ny):
            z_min_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,dxis[0]*float(inda),
                                         dxis[1]*float(indb),zmin)))
                                                          ).evalf()[:]
            z_max_bcs[:-1,inda,indb] = Euler_UCS.sol.subs(dict(zip(
                        Euler_UCS.vars_,(0.,dxis[0]*float(inda),
                                         dxis[1]*float(indb),zmax)))
                                                          ).evalf()[:]
#    x_bcs[:5,:,:] += Dirichlet_pin
#    x_bcs[5:14:4,:,:] = 1.
#    y_bcs[:5,:,:] += Dirichlet_pin
#    y_bcs[5:14:4,:,:] = 1.
#    z_bcs[:5,:,:] += Dirichlet_pin
#    z_bcs[5:14:4,:,:] = 1.
#    x_bcs[-1],y_bcs[-1],z_bcs[-1] = 1.,1.,1.

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
#                out = Euler_UCS.sol.subs(
#                    dict(zip(Euler_UCS.vars_,(0,i,j,k)))).evalf()
                initial_conds[:-1,i,j,k] = Euler_UCS.sol.subs(dict(zip(
                            Euler_UCS.vars_,
                            (0.,dxis[0]*float(i),dxis[1]*float(j),
                             dxis[2]*float(k))))).evalf()[:]
                source_funcs[i,j,k] = (lambda y, t, i=(dxis[0]*i), 
                                       j=(dxis[1]*j), 
                                       k=(dxis[2]*k) : numpy.array(
                        [item[0,0] for item in 
                         manufactured_source_function(t,i,j,k)]+[0]
                        ,dtype=numpy.float64))
    initial_conds[-1,:,:,:] = (initial_conds[5,:,:,:]*initial_conds[9,:,:,:]*
                               initial_conds[13,:,:,:])
    exact_solution = lambda t,x,y,z : (
        numpy.array(
            Euler_UCS.sol.subs(dict(zip(Euler_UCS.vars_,(t,dxis[0]*x,
                                                         dxis[1]*y,dxis[2]*z
                                                         )))).evalf()[:]
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
    solver_options[2] = 5
    solver_options[100] = 1
    solver_options[101] = 1
    solver_options[102] = 0
    solver_options[103] = 0
    stream_options = {
        'solver_type':'euler',
        'boundary_layers':False,
        'multistream':False,
        'solver_options':solver_options,
        'manufactured':'IMMS',
        'source_funcs' : source_funcs,
        'exact_sol_func' : exact_solution,
        'manufactured_object' : Euler_UCS,
        'dxis' : dxis,
        'discs' : [50.0]
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
    return bounds_init, initial_conds, stream_options

if __name__=='__main__':
    test = init()
