import sys
import numpy
import BoundaryConditions
import TimeAdvancementStuff as TAS
TAS = TAS.timeadvancementstuff
import Godunov
Godunov = Godunov.godunov
import CGNS_Interface
cgns = CGNS_Interface.cgns_interface
import imp, os
#import getopt
global interactive_flag
#prim_update = TAS.timeadvancementstuff.prim_update
#create_column = TAS.timeadvancementstuff.create_column
from sys import exit as sysexit

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Options(object):
    def __init__(self):
        return None
    
class Stream(object):
    '''
    A complete stream of flow. The Stream class contains the primary variable
    array, boundary condition specifications, simulation parameters, and
    anything else required for a fluid simulation of a single, contiguous,
    stream of fluid. More complex flows involving embedded surfaces, bleed
    flows, etc. must be handled by multiple Stream objects, communicating via
    boundary conditions.

    Stream objects are currently initialized using files, one for boundary
    conditions and one for everything else. This is unwieldy for dynamic stream
    creation, as would be required for complex flows, so other functionality
    will be implemented in the future.
    '''
    def __init__(self,bounds_init,initial_conds,stream_options):
        # Initialize boundary conditions
        self.bounds = BoundaryConditions.Bounds(bounds_init)
        bcextent = 1
        try:
            test = bool(initial_conds)
        except ValueError:
            test = initial_conds.any()
        if test:
            print "Using prescribed initial conditions."
            temp = initial_conds.shape
            self.main_data = numpy.zeros((temp[0],temp[1]+2*bcextent,
                                          temp[2]+2*bcextent,temp[3]+2*bcextent))
            self.main_data[:,0,1:-1,1:-1] = initial_conds
        else:
            print "No initial conditions given; using boundary conditions."
            temp2 = self._set_inflow_conditions_from_bounds(self.bounds)
            temp = temp2.shape
            # Assume that the inflow boundary is in the xi-direction
            self.main_data = numpy.zeros((temp[0],1+2*bcextent,
                                          temp[1]+2*bcextent,temp[2]+2*bcextent))
            self.main_data[:,1,1:-1,1:-1] = temp2
        self.xi_offset = 0
    def _set_inflow_conditions_from_bounds(self,bounds):
        '''
        Set initial conditions based on Inflow-type boundary conditions.

        Search a Bounds object for Inflow boundary conditions, and generate
        an initial condition for the simulation just inside those boundaries.
        At present, only Bounds objects with only one Inflow boundary are 
        supported.
        
        Args:
          bounds: Initialized Bounds object

        Returns:
          out: Initialized array corresponding the points just bordering the
                 Inflow boundary condition.
        '''
        inflows = []
        for face in bounds:
            for patch in face:
                if patch.type is 'Inflow':
                    if not (patch.which_face == 'left' 
                            or patch.which_face == 'right'):
                        print "Inflow condition detected on eta or zeta boundary!"
                        sysexit()
                    inflows.append(patch)
            sorted(inflows, key=lambda inflow: min(inflow.bounding_points[:][2]))
#        try:
#            if len(inflows)>1:
#                raise IndexError('More than 1 Inflow condition!')
#        except IndexError:
#            print "Multiple Inflow conditions not supported!"
#            sysexit()
        initial_condition = numpy.concatenate(
            [inflow.flow_state.copy() for inflow in inflows],axis=1)
        return initial_condition
#    def _find_computational_coords(self,
        
    def advance(self,dt,opts=Options()):
        opts.xi_offset = self.xi_offset
#        import pdb;pdb.set_trace()
        self.bounds(self.main_data,opts)
        for ind,element in numpy.ndenumerate(self.main_data):
            if numpy.isnan(element):
                print ind
                import pdb;pdb.set_trace()
#        print self.main_data[1,:,:,1]
#        import pdb;pdb.set_trace()
        self.main_data = Godunov.prim_update(
            numpy.asfortranarray(self.main_data),dt_in=.001,cfl=.25,bcextent=1,
            nx=self.main_data.shape[1]-2,
            ny=self.main_data.shape[2]-2,
            nz=self.main_data.shape[3]-2)
        #! Check to see if a new column needs to be created
        check_create_column = TAS.checkcreatecolumn(
            self.main_data[:,1,1:-1,1:-1],self.main_data[:,0,1:-1,1:-1])
        if check_create_column:
            # print "Creating a column"
            pass
            new_column = TAS.createcolumn(self.main_data[:,1,1:-1,1:-1],
                                          self.main_data[:,0,1:-1,1:-1],
                self.main_data.shape[2]-2,self.main_data.shape[3]-2)
            dims = self.main_data.shape
            self.main_data = numpy.concatenate((
                numpy.empty((dims[0],1,dims[2],dims[3])),self.main_data),axis=1)
            self.main_data[:,1,1:-1,1:-1] = new_column
            self.xi_offset = self.xi_offset + 1
        TAS.write_files_matlab(self.main_data[:,1:-1,1:-1,1],0.)
        return None

def run(input_file,interactive=False):
    global interactive_flag
    interactive_flag = interactive
    usage = ""
    fullpath = os.path.abspath(input_file)
    init_module = os.path.splitext(os.path.basename(fullpath))[0]
    init_ = imp.load_source(init_module,fullpath)
    try:
        bounds_init, initial_conds, stream_options = init_.init()
    except AttributeError:
        print "Input file does not contain an init() function"
        sys.exit()
    streams = [Stream(bounds_init, initial_conds,stream_options)]
    dt = .0001 
    TAS.write_files_matlab(streams[0].main_data[:,1:-1,1:-1,1],0.,first_flag=True)
    for step in range(50000):
        print "Time step = ",step
        for stream in streams:
#            if numpy.amax(stream.main_data[17,1:-1,1:-1,1:-1]) >= .5:
#                import pdb;pdb.set_trace()
            stream.advance(dt)
    cgns.write_initial_data(stream.main_data,'test.cgns')
    if interactive_flag:
        import pdb; pdb.set_trace()

if __name__=='__main__':
    import getopt
    filename = sys.argv[1]
    usage = ''
    try:
        if sys.argv[1].startswith('-'):
            raise IndexError()
    except IndexError:
        print "Error: First argument must be an input file!"
        sys.exit()
    options, args = getopt.getopt(sys.argv[2:],'ihp',
                                  ['interactive','help','profile'])
    if len(options)>1:
        print "No support yet for multiple options"
        sys.exit()
    elif len(options) == 0:
        run(filename)
    for option, value in options:
        if option in ('-i', '--interactive'):
            run(filename,interactive=True)
        elif option in ('-h', '--help'):
            print usage; sysexit(0)
        elif option in ('-p','--profile'):
            import cProfile
            cProfile.run('run(filename)','cProfout')
            import pstats
            p=pstats.Stats('cProfout')
            p.sort_stats('time').print_stats(10)        
#    usage = ""
#    try:
#        if sys.argv[1].startswith('-'):
#            raise IndexError()
#        fullpath = os.path.abspath(sys.argv[1])
#    except IndexError:
#        print "Error: First argument must be an input file!"
#        sys.exit()
#    options, args = getopt.getopt(sys.argv[2:],'ih', ['interactive','help'])
#    interactive_flag = False
#    for option, value in options:
#        if option in ('-i', '--interactive'):
#            interactive_flag = True
#        elif option in ('-h', '--help'):
#            print usage; sysexit(0)
#    init_module = os.path.splitext(os.path.basename(fullpath))[0]
#    init_=imp.load_source(init_module,fullpath)
#    try:
#        bounds_init, initial_conds, stream_options = init_.init()
#    except AttributeError:
#        print "Input file does not contain an init() function"
#        sys.exit()
#    streams = [Stream(bounds_init, initial_conds,stream_options)]
#    dt = .0001
#    for step in range(2500):
#        for stream in streams:
#            stream.advance(dt)
#    cgns.write_initial_data(stream.main_data,'test.cgns')
#    #    import pdb; pdb.set_trace()
        

#class Grid(object):
#    def __init__():
#        
#    def marchStep(dt):
#        
#    def applyBounds():    
            
