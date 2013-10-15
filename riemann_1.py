import numpy
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
    initial_conds = initial_conditions()
    solver_options = numpy.zeros(300)
    solver_options[0] = 1
    solver_options[101] = 1
    solver_options[102] = 1
    solver_options[103] = 0
    solver_options[104] = 1
    stream_options = {
        'solver_type':'euler',
        'boundary_layers':False,
        'multistream':False,
        'solver_options':solver_options
        }
    left_face_init.append(
        PatchInit('Transmissive',
                  ((0.0,-.1,-.1),(0.0,-.1,.1),(0.0,.1,0.1),(0.0,.1,-.1)),
                  1,'f = x'))
    right_face_init.append(
        PatchInit('Transmissive',
                  ((1.0,.1,-.1),(1.0,-.1,-.1),(1.0,-.1,.1),(1.0,.1,.1)),
                  1,'f = -x+1', 
                  None))
    bottom_face_init.append(
        PatchInit('Transmissive',
                  ((0.0,-.1,-.1),(1.0,-.1,-.1),(1.0,0.1,0.1),(0.0,0.1,0.1)),
            2,"f = y"))
    top_face_init.append(
        PatchInit('SolidWall',
                  ((-.1,.1,-.1),(1.0,.1,-.1),(1.0,.1,0.1),(-.1,.1,0.1)),
            2,"f = -y"))
    back_face_init.append(
        PatchInit('Transmissive',
                  ((0.0,-.1,-.1),(1.0,-.1,-.1),(1.0,.1,-.1),(0.0,.1,-.1)),3))
    front_face_init.append(
        PatchInit('Transmissive',
                  ((0.0,-.1,0.1),(1.0,-.1,0.1),(1.0,.1,0.1),(0.0,.1,0.1)),3))
    bounds_init = BoundsInit(lf=left_face_init,rf=right_face_init,
                             bof=bottom_face_init,tf=top_face_init,
                             baf=back_face_init,ff=front_face_init)
    return bounds_init, initial_conds, stream_options

def initial_conditions():
    import numpy as np
    inputs = np.zeros((21,50,1,1))
    inputs[0,:,:] = 1.
    inputs[1,:,:] = 1.
    inputs[2,:,:] = .75
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
            inputs[18,inda,indb] = 1./50.*(inda+.5)
            inputs[19,inda,indb] = 0.
    inputs[20,:,:] = 1./24./24.
    inputs[0,25:,:] = .1
    inputs[1,25:,:] = .125
    inputs[2,25:,:]
    # np.savetxt("transonic_duct_bounds_x.txt",np.reshape(inputs[17,:,:],-1))
    # np.savetxt("transonic_duct_bounds_y.txt",np.reshape(inputs[18,:,:],-1))
    # np.savetxt("transonic_duct_bounds_z.txt",np.reshape(inputs[19,:,:],-1))
    # np.save("transonic_duct_bounds",inputs)
    return(inputs)

if __name__=='__main__':
    inputs = transonic_duct_inflow_generator()
#    for inda in range(inputs.shape[1]):
#        print inputs[5:15,inda,0]
    import numpy as np
    import matplotlib
#    print matplotlib.__version__
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    def randrange(n, vmin, vmax):
        return (vmax-vmin)*np.random.rand(n) + vmin
    
    fig = plt.figure()
    ax = Axes3D(fig)
    n = 100
    for i in range(inputs.shape[1]):
        ax.scatter(inputs[17,i,:],inputs[18,i,:],inputs[19,i,:])            
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()
    print init.__doc__