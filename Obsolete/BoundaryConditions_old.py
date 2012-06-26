import numpy as np
import scipy.io as sio
import stla_io as stl
import BoundaryConditions as bc
from sys import exit as sysexit
bctri = bc.triangulatedboundaryconditions

class Error(Exception):
    pass

class STLError(Error):
    def __init__(self,msg):
        self.msg = msg
        
class InputFormatError(Error):
    def __init__(self,msg):
        self.msg = msg

class Bounds(object):
    '''
    Boundary Conditions object for a given Stream.
    
    A complete description of the boundary conditions associated with a given
    Stream. Since all Streams are structured, the boundary conditions can be
    divided into 6 Sides, one for each face of the Stream.
    '''
    def __init__(self,default_bounds,stl_bounds_file=None,num_faces=0,
                 init=dict()):
        '''
        Initialize Bounds object.

        Args:
          default_bounds: Default_bounds is a description of a boundary surface.
            It may be superceded by other boundary conditions such as solid
            walls defined in an STL file. It may be either a filename string or
            or a list of strings (the results of a file.readlines() command).
          stl_bounds_file: Filename string of an ASCII .stl file describing any
            solid wall geometries present. Will eventually support the results 
            of an stl_read command (a list of arrays of nodes and faces).
          num_faces: Integer number of individual faces in the .stl file. Must
            be present if an STL file is given.
        Returns:
          self.left_face
          self.right_face
          self.top_face
          self.bottom_face
          self.back_face
          self.front_face
        Raises:
          None as yet.
        '''
        # Read STL file, returning lists of triangulation vertices, indices of
        # the vertices associated with triangles, and the normal vectors of
        # those triangles.
        if stl_bounds_file:
            num_nodes = num_faces*3
            [self.stl_nodes,self.stl_face_nodes,self.stl_face_normal,
             error]=stl.stla_read(stl_bounds_file,num_nodes,num_faces)
        if error:
            try:
                str="STLError: stla_read failed in BoundaryConditions.__init__"
                raise STLError(str)
            except STLError as e:
                print e.msg
                sysexit()

        # Isolate the parts of the boundary specification pertaining to each
        # Side and pass them on.
        try:
            infile = open(default_bounds,'r')
            lines = infile.readlines()
        except TypeError:
            print '''
                  Note: Bounds.__init__ received default boundaries
                  as lines instead of as a file
                  '''
            lines = default_bounds
        for inda in range(6):
            substr = self._pick_sides(lines)
            cmd = 'self.'+substr[0].split()[0].strip().lower()
            cmd = cmd+'_face = Face(substr)'
#            print cmd
            exec(cmd)
    def __call__(self,main,t=0):
        '''
        Apply boundary conditions to the Stream.

        Args:
          self: Bounds object containing the initialized boundary conditions.
          main: Array object containing the state information for the Stream.
          t: Simulation time; used for time-dependent boundary conditions.

        Returns:
          main: Updated Array of state variables.
        '''
#        print main[0,0,1:-1,1:-1]
        self.left_face(main[:, 0,1:-1,1:-1],t)
        self.right_face(main[:,-1,1:-1,1:-1],t)
        self.bottom_face(main[:,1:-1, 0,1:-1],t)
        self.top_face(main[:,1:-1,-1,1:-1],t)
        self.back_face(main[:,1:-1,1:-1, 0],t)
        self.front_face(main[:,1:-1,1:-1,-1],t)
    def _pick_sides(self,lines):
        '''
        Extract the specification of a boundary Face from a list of strings
        
        A boundary Face specification begins with "____ boundary:", perhaps
        with leading or trailing whitespace. _pick_sides extracts and returns
        the substring containing Face specification, also as a list of strings.
        
        Args:
          lines: List of strings (as in the result of file.readlines()). First
            line must begin a boundary Face specification.

        Returns:
          side_lines: List of strings, containing exactly one complete Face
            specification      
        '''
        side_lines=[]
        side_lines.append(lines.pop(0).rstrip())
        if not side_lines[-1].endswith('boundary:'):
            print "Error in pick_sides: Invalid first entry"
        try:
            while not lines[0].rstrip().endswith('boundary:'):
                side_lines.append(lines.pop(0).rstrip())
        except IndexError:
            # This simply means that we've encountered the end of the list.
            # print 'End of list'
            pass
        return(side_lines)
class Face(object):
    '''
    Boundary Conditions object for an individual face of a Stream.
    '''
    types = {'Inflow':'inflow_init',
             'Outflow':'outflow_init',
             'Wall':'wall_init',
             'Pressure':'pressure_init',
             'Transmissive':'transmissive_init'}
    def __init__(self,lines):
        '''
        Initialize Face object.
        '''
        # Look for field descriptions in the input file.

        # This identifies all of the different types of boundary condition that
        # are associated with a particular Face. It returns a list of the 
        # indices of those types, so that the boundary condition can be read in
        # distinct types.
        inds = []
        for type in self.types.keys():
            inds.append(index_substring(lines,type))
        flat_inds = [item for indlist in inds for item in indlist]
        flat_inds.sort()
        flat_inds.append(-1)
        junk = Patch(lines[flat_inds[0]:flat_inds[1]])
        sysexit()
#        for ind in range(len(flat_inds)):
#            types[lines[flat_inds[ind]].strip()](
#                lines[flat_inds[ind]:flat_inds[ind+1]])
    def __call__(self,inputs,t):
        '''
        '''
class Patch:
    '''
    An individual boundary condition; a piece of a Face.

    Each Face contains a single Patch. A Patch is a particular type of boundary
    condition, such as a pressure condition or an inflow condition. Each
    simulation can have, at present, a single inflow condition with a specified
    grid.

    Args:
      lines: Description of a Patch.
    Returns:
      self.bounding_points: Nested tuple, containing three dimensional points in
        space. At present, only sets of four points are supported.
      self.boundary_surface: String describing the surface of the Patch in the
        form "f = f(x,y,z)". The surface is taken to correspond to f = 0.
      self.flow_state: Array of flow variables. May be either a 1-dimensional
        array representing a constant flow state throughout the Patch, or else a
        three-dimensional array containing a general boundary description with
        an associated grid.
    '''
    # What to do to make this a parameter-based function?
    def __init__(self,lines):
        self.type = lines[0].strip()
        self.bounding_points = None
        self.boundary_surface = None
        self.flow_state = None
        self.fields = {'Bounding Points:':'junk',
                       'Boundary Surface:':'junk',
                       'Flow State:':'junk'}
        inds = []
        for field in self.fields.keys():
            inds.append(index_substring(lines,field))
            inds = index_substring(lines,field)
            if len(inds)>1:
                msg = "Duplicate field entries detected in Patch.__init__!"
                try:
                    raise InputFormatError(msg)
                except InputFormatError as e:
                    print e.msg
                    print 'Inds = ',inds
                    print 'Field = ',field
                    print 'Lines = ',lines
                    sysexit()
            elif not inds:
                self.fields[field] = None
            for ind in inds:
                self.read_field(lines[ind:ind+3])
        sysexit()
    def wall_init(self,lines):
        pass
    def read_field(self,lines):
        if lines[0].rstrip().endswith('Bounding Points:'):
            self.bounding_points = eval(lines[1])
        elif lines[0].rstrip().endswith('Boundary Surface:'):
            self.boundary_surface_string = lines[1]
        elif lines[0].rstrip().endswith('Flow State:'):
            if lines[1].strip() == 'Initialized:':
                self.flow_state = np.load(lines[2].strip())
            else:
                self.flow_state = eval(lines[1])
        else:
            print "Bad field value in read_field!"
            print "Value = ", lines[0]
            sysexit()
    
def index_substring(strings, substring):
    '''
    Find all elements of a list containing a given substring.
    
    Args:
      strings: List of strings so be searched
      substring: Substring being sought
    Returns:
      Indices for elements of strings containing substring.
      Returns None if substring is not found in strings.
    '''
    return list(i for i, string in enumerate(strings) if substring in string)
            
def leadingEdgeSearchSTL(leading_face_array, nodes_x, faces_x, face_nodes_x):
    '''LeadingEdgeSearchSTL searches for boundary conditions that may apply to
    a two-dimensional array of grid points representing the leading edge of 
    the simulation. The boundaries are assumed to come from a .stl file, and 
    represent solid walls.
    Inputs: 
       leading_face_array - Array of leading-edge grid points with their
                              associated flow and metric values
       nodes_x            - Array of boundary vertices read from a .stl file
       faces_x            - Array of normal vectors read from a .stl file
       face_nodes_x       - Integer indexing array connecting nodes to faces 
                              (Fortran-style indexing, starting from 1)
     Outputs:
        indx              - Array of indices connecting grid points with 
                              individual boundary faces (Fortran-style)
     '''
# Each point has the form:
#   0   1   2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
# [ p, rho, u, v, w, A, B, C, L, M,  N, P, Q, R, U, V, W, x, y, z, J ]

    lfa = leading_face_array
    indx = np.zeros((leading_face_array.shape[1],leading_face_array.shape[2])
                    ,'int')
    for i in range(leading_face_array.shape[1]):
        for j in range(leading_face_array.shape[2]):
            indx[i,j] = bctri.leadingedgepointsearch(nodes_x,faces_x,
                               face_nodes_x,leading_face_array[:,i,j])
    return(indx)

if __name__=='__main__':
    from transonic_duct_bounds import transonic_duct_init as tdinit
    inputs = tdinit()
    sio.savemat('transonic_duct',{'inputs':inputs})
    infile = open('transonic_duct_bounds.in','r')
    test = Bounds(default_bounds='transonic_duct_bounds.in',
                  stl_bounds_file='transonic_duct.stl',num_faces=8)
#    num_faces = 8
#    num_nodes = num_faces*3
#    [nodes,face_nodes,face_normal,error]=stl.stla_read(
#        'transonic_duct.stl',num_nodes,num_faces)
#    out= leadingEdgeSearchSTL(inputs,nodes,face_normal,face_nodes)
