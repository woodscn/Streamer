import Stream
import numpy
import Godunov_driver
main = Stream.Stream(None,None,None)
ny = 100
dy = 1/(ny-1.)
temp = numpy.zeros([21,ny],dtype=numpy.float64,order='F')
temp[0] = 1
temp[1] = 1
temp[2] = 0
temp[5:14:4] = dy
for n in range(ny):
    temp[18,n] = dy*n
temp[20] = dy*dy*dy
main.main_data = numpy.zeros([21,3,ny+2,3],dtype=numpy.float64,order='F')
main.main_data[:,1,1:-1,1] = temp
advance_options = numpy.zeros(300)
advance_options[0]=1
advance_options[100]=1
advance_options[101]=1
advance_options[102]=0
advance_options[103]=1
#opt
dt = .0001
for n in range(1800):
    main.main_data[:,0,:,:] = main.main_data[:,1,:,:]
    main.main_data[:,-1,:,:] = main.main_data[:,-2,:,:]
    main.main_data[:,:,0,:] = main.main_data[:,:,1,:]
    main.main_data[:,:,-1,:] = main.main_data[:,:,-2,:]
    main.main_data[:,:,:,0] = main.main_data[:,:,:,1]
    main.main_data[:,:,:,-1] = main.main_data[:,:,:,-2]
#    [main.main_data,dt_out] = Godunov_driver.godunovdriver.prim_update(
#        numpy.asfortranarray(main.main_data),dt,.7,100,1,1,advance_options)
    dt_out = main.advance(dt,advance_options)
print main.main_data[0,1,1:-1,1]
