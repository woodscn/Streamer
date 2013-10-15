import Stream
import numpy
import Godunov_driver
main = Stream.Stream(None,None,None)
nx = 100
dx = 1/(nx-1.)
temp = numpy.zeros([21,nx],dtype=numpy.float64,order='F')
temp[0] = 1
temp[1] = 1
temp[2] = .75
temp[0,30:] = .1
temp[1,30:] = .125
temp[2,30:] = 0
temp[5:14:4] = dx
for n in range(nx):
    temp[17,n] = dx*n
temp[20] = dx*dx*dx
main.main_data = numpy.zeros([21,nx+2,3,3],dtype=numpy.float64,order='F')
main.main_data[:,1:nx+1,1,1] = temp
main.main_data[:,0] = main.main_data[:,1]
main.main_data[:,nx+1] = main.main_data[:,nx]
advance_options = numpy.zeros(300)
advance_options[0]=1
advance_options[100]=1
advance_options[101]=1
advance_options[102]=0
advance_options[103]=1
#opt
dt = .0001
for n in range(1800):
    main.main_data[:,:,0,1] = main.main_data[:,:,1,1]
    main.main_data[:,:,2,1] = main.main_data[:,:,1,1]
    main.main_data[:,:,1,0] = main.main_data[:,:,1,1]
    main.main_data[:,:,1,2] = main.main_data[:,:,1,1]
#    [main.main_data,dt_out] = Godunov_driver.godunovdriver.prim_update(
#        numpy.asfortranarray(main.main_data),dt,.7,100,1,1,advance_options)
    dt_out = main.advance(dt,advance_options)
print main.main_data[1,1:-1,1,1]