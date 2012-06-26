def round_duct_init():
    import numpy as np
    inputs = np.zeros((21,25,25))
    inputs[0,:,:] = 1.
    inputs[1,:,:] = 1.
    inputs[2,:,:] = 1.
    inputs[3,:,:] = 0.
    inputs[4,:,:] = 0.
    for inda in range(25):
        for indb in range(25):
            inputs[ 9,inda,indb] = 1./12.*np.sqrt(1.-(indb/12.-1)**2/2.)
            inputs[12,inda,indb] = -(indb/12.-1.)*(inda/12.-1.)/(
                24.*np.sqrt(1.-(inda/12.-1.)**2/2.))
            inputs[10,inda,indb] = -(indb/12.-1.)*(inda/12.-1.)/(
                24.*np.sqrt(1.-(indb/12.-1.)**2/2.))
            inputs[13,inda,indb] = 1./12.*np.sqrt(1.-(inda/12.-1)**2/2.)
    inputs[ 5,:,:] = 1.
    inputs[ 6,:,:] = 0.
    inputs[ 7,:,:] = 0.
    inputs[ 8,:,:] = 0.
    inputs[11,:,:] = 0.
    inputs[14,:,:] = 0.
    inputs[15,:,:] = 0.
    inputs[16,:,:] = 0.
    inputs[17,:,:] = 0.
    for inda in range(25):
        for indb in range(25):
            inputs[18,inda,indb] = (inda/12.-1.)*np.sqrt(1.-(indb/12.-1.)**2/2.)
            inputs[19,inda,indb] = (indb/12.-1.)*np.sqrt(1.-(inda/12.-1.)**2/2.)
#            print inputs[18:20,inda,indb]
    inputs[20,:,:] = inputs[9,:,:]*inputs[13,:,:]-inputs[10,:,:]*inputs[12,:,:]
    np.savetxt("round_duct_bounds_x.txt",np.reshape(inputs[17,:,:],-1))
    np.savetxt("round_duct_bounds_y.txt",np.reshape(inputs[18,:,:],-1))
    np.savetxt("round_duct_bounds_z.txt",np.reshape(inputs[19,:,:],-1))
    print '''
    Warning: The corner nodes ([0,0],[0,-1],[-1,0],[-1,-1]) may be
        ill-conditioned, because the computed Jacobian for those
        nodes is very small.
    '''
    return(inputs)

if __name__=='__main__':
#    from pylab import *
    inputs = round_duct_init()
#    N = 30
#    x = 0.9*rand(N)
#    y = 0.9*rand(N)
#    area = pi*(10 * rand(N))**2 # 0 to 10 point radiuses
    for inda in range(inputs.shape[1]):
        print inputs[5:15,inda,0]
#    print inputs[12,:2,:2]
#    print inputs[12,:2,-2:]
#    print inputs[12,-2:,:2]
#    print inputs[12,-2:,-2:]
#    contourf(inputs[18,:,:],inputs[19,:,:],inputs[10,:,:])
#    xlim((-1.01,1.01))
#    ylim((-1.01,1.01))
#    show()
    import numpy as np
    import matplotlib
    print matplotlib.__version__
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    def randrange(n, vmin, vmax):
        return (vmax-vmin)*np.random.rand(n) + vmin
    
    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    n = 100
#    for c, m, zl, zh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
    for i in range(inputs.shape[1]):
#        for j in range(inputs.shape[2]):
        ax.scatter(inputs[17,i,:],inputs[18,i,:],inputs[19,i,:])            
#        xs = randrange(n, 23, 32)
#        ys = randrange(n, 0, 100)
#        zs = randrange(n, zl, zh)
#        ax.scatter(xs, ys, zs, c=c, marker=m)
#        print 'xs = ',xs
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()
