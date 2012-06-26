import sys
import numpy
import BoundaryConditions
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
        try:
            test = bool(initial_conds)
        except ValueError:
            test = initial_conds.any()
        if test:
            print "Using prescribed initial conditions."
            self.main_data = initial_conds
        else:
            print "No initial conditions given; using boundary conditions."
            self.main_data = self._set_inflow_conditions_from_bounds(
                self.bounds)
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
        print bounds
        sys.exit()
if __name__=='__main__':
    try:
        init_module=sys.argv[1].rstrip('.py')
    except IndexError:
        print "Error: No input file given!"
        sys.exit()
    init_=__import__(init_module)
    try:
        bounds_init, initial_conds, stream_options = init_.init()
    except AttributeError:
        print "Input file does not contain an init() function"
        sys.exit()
    stream = Stream(bounds_init, initial_conds,stream_options)
    print "success!"

#class Grid(object):
#    def __init__():
#        
#    def marchStep(dt):
#        
#    def applyBounds():    
