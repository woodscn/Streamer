import numpy as np

def steady_riemann(main_data):
    nx,ny,nz = 60,100,1#main_data.shape[1:]
    xmin,xmax = 0.,.6
    ymin,ymax = -.5,.5
    zmin,zmax = 0.,0.
    dx,dy,dz = (xmax-xmin)/(nx),(ymax-ymin)/(ny),1
    main_data[:,0,:,:] = 0.
    main_data[0,0,:,:] = 1.
    main_data[0,0,ny/2:,:] = .25
    main_data[1,0,:,:] = 1.
    main_data[1,0,ny/2:,:] = .5
    main_data[2,0,:,:] = 2.4*(1.4)**.5
    main_data[2,0,ny/2:,:] = 7*(.25/.5*1.4)**.5
    main_data[5,0,:,:] = dx
    main_data[9,0,:,:] = dy
    main_data[13,0,:,:] = dz
    for ind in range(ny):
        main_data[18,0,ind+1,:] = (ind+.5)*dy
    main_data[20,0,:,:] = dx*dy*dz
    main_data[:,:,0,:] = main_data[:,:,1,:]
    main_data[:,:,-1,:] = main_data[:,:,-2,:]
    main_data[:,:,:,0] = main_data[:,:,:,1]
    main_data[:,:,:,-1] = main_data[:,:,:,-2]
    main_data[:,-1,:,:] = main_data[:,-2,:,:]
    return main_data

def transonic_duct(main_data):
    nx,ny,nz = main_data.shape[1:]
    xes = 0.,.5,1.,3.6
    ymin,ymax = 0.,1.
    zmin,zmax = 0.,0.
    dx,dy,dz = .02,.02,1.
    nxstars = []
    xiter = iter(xes[1:-1])
    xstar = xiter.next()
    main_data[:,0,:,:] = 0.
    for ind, x in enumerate(main_data[17,:,1,1]):
        if x > xstar:
            nxstars.append(ind)
            try:
                xstar = xiter.next()
            except(StopIteration):
                break
    # Front and back walls
    main_data[:,1:-1,1:-1,0] = main_data[:,1:-1,1:-1,1]
    main_data[:,1:-1,1:-1,-1] = main_data[:,1:-1,1:-1,-2]
    # Bottom wall
    if len(nxstars)>0:
        main_data[:,1:nxstars[0],0,1] = wall_reflection(
            main_data[:,1:nxstars[0],1,1],0)
        if len(nxstars)>1:
            main_data[:,nxstars[0]:nxstars[1],0,1] = wall_reflection(
                main_data[:,nxstars[0]:nxstars[1],1,1],15./180.*np.pi)
            main_data[:,nxstars[1]:-1,0,1] = wall_reflection(
                main_data[:,nxstars[1]:-1,1,1],0)
        else:
            main_data[:,nxstars[0]:-1,0,1] = wall_reflection(
                main_data[:,nxstars[0]:-1,1,1],15./180.*np.pi)
    else:
        main_data[:,1:-1,0,1] = wall_reflection(main_data[:,1:-1,1,1],0)
    # Top wall
    main_data[:,1:-1,-1,1] = wall_reflection(
        main_data[:,1:-1,-2,1],0)
    # Outflow boundary
    main_data[:,-1,:,:] = main_data[:,-2,:,:]
    # Inflow boundary
    main_data[:,0,:,:] = 0.
    main_data[0,0,:,:] = 1.
    main_data[1,0,:,:] = 1.
    main_data[2,0,:,:] = 1.8*(1.4)**.5
    main_data[5,0,:,:] = dx
    main_data[9,0,:,:] = dy
    main_data[13,0,:,:] = dz
    main_data[20,0,:,:] = dx*dy*dz
    for ind in range(len(main_data[18,0,1:-1,1])):
        main_data[18,0,ind+1,1] = (float(ind)+.5)*dy
    # Inflow grid motion
    main_data[14:17,0,:,:] = .25*main_data[2:5,0,:,:]
    return main_data

def wall_reflection(interior,angle):
    import numpy as np
    cos = np.cos
    sin = np.sin
    tan = np.tan
    bound = 0*interior
    bound[0,:] = interior[0,:]
    bound[1,:] = interior[1,:]
    bound[2,:] = np.cos(2.*angle)*interior[2,:]+np.sin(2.*angle)*interior[3,:]
    bound[3,:] = np.sin(2.*angle)*interior[2,:]-np.cos(2.*angle)*interior[3,:]
    if(np.abs(angle)<=.5*np.pi):
        bound[5,:] = interior[5,:]
        bound[9,:] = (
            interior[9,:]*(interior[5,:]+bound[5,:]*sin(angle)**2)-sin(angle)
            *cos(angle)*(interior[5,:]*interior[8,:]+interior[6,:]
                         *interior[9,:]+bound[5,:]*interior[8,:])
            )/(bound[5,:]*(cos(angle)**2-sin(angle)**2)-interior[5,:]
               *sin(angle)**2+interior[6,:]*sin(angle)*cos(angle))
        bound[8,:] = tan(angle)*(interior[9,:]+bound[9,:])-interior[8,:]
        bound[6,:] = tan(angle)*(interior[5,:]+bound[5,:])-interior[6,:]
    else:
        bound[6,:] = interior[6,:]
        bound[8,:] = (
            interior[8,:]*(interior[6,:]+bound[6,:]*cos(angle)**2)-sin(angle)
            *cos(angle)*(interior[5,:]*interior[8,:]+interior[6,:]
                         *interior[9,:]+bound[6,:]*interior[9,:])
            )/(bound[6,:]*(sin(angle)**2-cos(angle)**2)-interior[6,:]
               *cos(angle)**2+interior[5,:]*sin(angle)*cos(angle))
        bound[5,:] = 1./tan(angle)*(interior[6,:]+bound[6,:])-interior[5,:]
        bound[9,:] = 1./tan(angle)*(interior[8,:]+bound[8,:])-interior[9,:]
#    bound[5:14,:] = interior[5:14,:]
#
#    if(np.abs(angle)<=.5*np.pi):
#        bound[5,:] = interior[5,:]
#        numer = (
#            interior[9,:]*(interior[5,:]+bound[5,:]*sin(angle)**2)
#            -sin(angle)*cos(angle)*
#            (interior[5,:]*interior[8,:]+interior[6,:]*interior[9,:]+
#             bound[5,:]*interior[8,:]))
#        denom = (bound[5,:]*(cos(angle)**2-sin(angle)**2)-interior[5,:]*
#                 sin(angle)**2+interior[6,:]*sin(angle)*cos(angle))
#        bound[9,:] = numer/denom
#        bound[8,:] = tan(angle)*(interior[9,:]+bound[9,:])-interior[8,:]
#        bound[6,:] = tan(angle)*(interior[5,:]+bound[5,:])-interior[6,:]
#    else:
#        bound[6,:] = interior[6,:]
#        numer = (
#            interior[8,:]*(interior[6,:]+bound[6,:]*cos(angle)**2)
#            -sin(angle)*cos(angle)*
#            (interior[5,:]*interior[8,:]+interior[6,:]*interior[9,:]+
#             bound[6,:]*interior[9,:]))
#        denom = (bound[6,:]*(sin(angle)**2-cos(angle)**2)-interior[6,:]
#                 *cos(angle)**2+interior[5,:]*sin(angle)*cos(angle))
#        bound[8,:]  = numer/denom
#        bound[5,:] = 1./tan(angle)*(interior[6,:]+bound[6,:])-interior[5,:]
#        bound[9,:] = 1./tan(angle)*(interior[8,:]+bound[8,:])-interior[9,:]
    bound[13,:] = interior[13,:]
    bound[14:17] = interior[14:17]
    bound[20] = interior[20]
    return bound
        
