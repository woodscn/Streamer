import Stream
import numpy
import Godunov_driver
main = Stream.Stream(None,None,None)
nx = 100
dx = 1/(nx-1.)
temp = numpy.zeros([21,nx],dtype=numpy.float64,order='F')
main=[Stream.Stream(None,None,None),Stream.Stream(None,None,None)]
left = numpy.zeros(21)
left[5:14:4] = dx
left[20] = dx*dx*dx
right = left.copy()
left[0:3] = [1,1,.75]
right[0:3] = [.1,.125,0]
main[0].main_data = numpy.zeros([21,32,3,3],dtype=numpy.float64,order='F')
for n in range(32):
    main[0].main_data[:,n,1,1] = left
    main[0].main_data[17,n,1,1] = (n-1)*dx
main[1].main_data = numpy.zeros([21,72,3,3],dtype=numpy.float64,order='F')
for n in range(72):
    main[1].main_data[:,n,1,1] = right
    main[1].main_data[17,n,1,1] = (n+30)*dx
advance_options = numpy.zeros(300)
advance_options[0]=1
advance_options[100]=1
advance_options[101]=1
advance_options[102]=0
advance_options[103]=0
dt = .0001
for n in range(1800):
    for stream in main:
        stream.main_data[:,:,0,1] = stream.main_data[:,:,1,1]
        stream.main_data[:,:,2,1] = stream.main_data[:,:,1,1]
        stream.main_data[:,:,1,0] = stream.main_data[:,:,1,1]
        stream.main_data[:,:,1,2] = stream.main_data[:,:,1,1]
    main[0].main_data[:,-1,:,:] = main[1].main_data[:,1,:,:]
    main[1].main_data[:,0 ,:,:] = main[0].main_data[:,-2,:,:]
    [main[0].main_data,dt_out] = Godunov_driver.godunovdriver.prim_update(
        numpy.asfortranarray(main[0].main_data),dt,.7,30,1,1,advance_options)
    [main[1].main_data,dt_out] = Godunov_driver.godunovdriver.prim_update(
        numpy.asfortranarray(main[1].main_data),dt,.7,70,1,1,advance_options)
temp = numpy.zeros(100)
temp[0:30] = main[0].main_data[1,1:-1,1,1]
temp[30:] = main[1].main_data[1,1:-1,1,1]
print temp
