import numpy as np
from numpy import f2py
import scipy.io as sio
import STLA_IO as stl
import BoundaryConditionsStuff
from sys import exit as sysexit
from main import Stream
bcs = BoundaryConditionsStuff.boundaryconditionsstuff

class Error(Exception):
    pass

class STLError(Error):
    def __init__(self,msg):
        self.msg = msg
        
class InputFormatError(Error):
    def __init__(self,msg):
        self.msg = msg

class UpperBoundsError(Error):
    pass

class Bounds(object):
    '''
    Boundary Conditions object for a given Stream.
    
    A complete description of the boundary conditions associated with a given
    Stream. Since all Streams are structured, the boundary conditions can be
    divided into 6 Sides, one for each face of the Stream.
    '''
    def __init__(self,default_bounds,stl_bounds_file=None,num_faces=0):
        '''
        Initialize Bounds object.

        Args:
          default_bounds: Default_bounds is a description of a boundary surface.
            It may be superceded by other boundary conditions such as solid
            walls defined in an STL file. It is an object containing one 
            element for each Face.
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
          STLError: There was an error reading the .stl file.
        '''
        # Read STL file, returning lists of triangulation vertices, indices of
        # the vertices associated with triangles, and the normal vectors of
        # those triangles.
        if stl_bounds_file:
            print " Warning: STL boundaries are not yet implemented."
            sysexit()
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
        self.left_face = Face('left',default_bounds.left_face)
        self.right_face = Face('right',default_bounds.right_face)
        self.bottom_face = Face('bottom',default_bounds.bottom_face)
        self.top_face = Face('top',default_bounds.top_face)
        self.back_face = Face('back',default_bounds.back_face)
        self.front_face = Face('front',default_bounds.front_face)
        fortran_normal_src = "! -*- f90 -*-\n"
        for face in self:
            for patch in face:
                try:
                    fortran_normal_src += patch.gradsrc
                except AttributeError: # Some patches don't have normal vectors
                    pass
        f2py.compile(fortran_normal_src,modulename='FortranNormalVectors',
                     verbose=False,source_fn='FortranNormalVectors.f90',
                     extra_args='--quiet')
    def __iter__(self):
        faces = (self.left_face,self.right_face,self.bottom_face,self.top_face,
                 self.back_face,self.front_face)
        for face in faces:
            yield face
    def __call__(self,main,opts,t=0):
        '''
        Update boundary computational coordinates and fill in ghost cells.
        
        Args:
          self: Bounds object containing the initialized boundary conditions.
          main: Array object containing the state information for the Stream.
          t: Simulation time; used for time-dependent boundary conditions.

        Returns:
          self: Updated boundary computational coordinates.
          main: Updated Array of state variables.
        '''
        main[:, 0,1:-1,1:-1] =   self.left_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:, 0,1:-1,1:-1],opts)
        main[:,-1,1:-1,1:-1] =  self.right_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:,-1,1:-1,1:-1],opts)
        main[:,1:-1, 0,1:-1] = self.bottom_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:,1:-1, 0,1:-1],opts)
        main[:,1:-1,-1,1:-1] =    self.top_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:,1:-1,-1,1:-1],opts)
        main[:,1:-1,1:-1, 0] =   self.back_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:,1:-1,1:-1, 0],opts)
        main[:,1:-1,1:-1,-1] =  self.front_face(main[:,1:-1,1:-1,1:-1], t,
                                                main[:,1:-1,1:-1,-1],opts)
class Face(object):
    '''
    Boundary Conditions object for an individual face of a Stream.
    '''
    def __init__(self,which_face,face_init):
        '''
        Initialize Face object.
        '''
        # Look for field descriptions in the input file.
        self.patches = []
        self.which_face = which_face
        for patch in face_init:
            self.patches.append(Patch(patch))
        for patch in self.patches:
            patch.which_face = self.which_face
        if self.which_face == "left" or self.which_face == "right":
            # Sort in ascending order, first by y and then by z coordinates.
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,2]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,1]))
        elif self.which_face == "bottom" or self.which_face == "top":
            # Sort first by z, then by x.
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,0]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,2]))
        elif self.which_face == "back" or self.which_face == "front":
            # Sort first by x, then by y
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,1]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,0]))
    def __iter__(self):
        for patch in self.patches:
            yield patch
    def __call__(self,main_data,t,out,opts):
        '''
        '''
        if self.which_face == "left":
            main_pass = main_data[:,0,eta_min:eta_max,zeta_min:zeta_max]
        elif self.which_face == "right":
            main_pass = main_data[:,-1,eta_min:eta_max,zeta_min:zeta_max]
        elif self.which_face == "bottom":
            main_pass = main_data[:,xi_min:xi_max,0,zeta_min:zeta_max]
        elif self.which_face == "top":
            main_pass = main_data[:,xi_min:xi_max,-1,zeta_min:zeta_max]
        elif self.which_face == "back":
            main_pass = main_data[:,xi_min:xi_max,eta_min:eta_max,0]
        elif self.which_face == "front":
            main_pass = main_data[:,xi_min:xi_max,eta_min:eta_max,-1]
        else:
            raise(Error)
        for patch in self.patches:
            # First, I need to advance the positions of boundary patches.
            try:
                test = bool(patch.bounding_points_xi)
            except ValueError:
                test = patch.bounding_points_xi.any()
            if not test:
                # Compute the appropriate computational coordinates
                # self.bounding_points_xi will be None unless it has been
                # initialized previously.
                ccc = self._compute_computational_coordinates
                patch.bounding_points_xi = ccc(patch,main_data,
                                               opts)
                # Seriously, how should I do this? Which point should I pick to
                # extrapolate my grid metric? Grr...! No wonder I'm job hunting!
            else:
                # Eventually, just time-advance the points. For now, do the same
                # as above
                ccc = self._compute_computational_coordinates
                patch.bounding_points_xi = ccc(patch,main_data,
                                               opts)
            # Next, I need to pass the appropriate grid points to the patch.
            xi_min   = np.floor(min(patch.bounding_points_xi[:,0]))
            xi_max   =  np.ceil(max(patch.bounding_points_xi[:,0]))
            eta_min  = np.floor(min(patch.bounding_points_xi[:,1]))
            eta_max  =  np.ceil(max(patch.bounding_points_xi[:,1]))
            zeta_min = np.floor(min(patch.bounding_points_xi[:,2]))
            zeta_max =  np.ceil(max(patch.bounding_points_xi[:,2]))

            if self.which_face == "left":
                out = patch(main_data[:,0,eta_min:eta_max,
                                        zeta_min:zeta_max],opts,out)
            elif self.which_face == "right":
                out = patch(main_data[:,-1,eta_min:eta_max,
                                        zeta_min:zeta_max],opts,out)
            elif self.which_face == "bottom":
                out = patch(main_data[:,xi_min:xi_max,0,
                                        zeta_min:zeta_max],opts,out)
            elif self.which_face == "top":
                out = patch(main_data[:,xi_min:xi_max,-1,
                                        zeta_min:zeta_max],opts,out)
            elif self.which_face == "back":
                out = patch(main_data[:,xi_min:xi_max,eta_min:eta_max,
                                        0],opts,out)
            elif self.which_face == "front":
                out = patch(main_data[:,xi_min:xi_max,eta_min:eta_max,
                                        -1],opts,out)
            else:
                raise(Error)
        # Can a patch be duplicated on multiple faces? I guess so.
        # Do I need to have the capability to dynamically add patches to a face?
        # That would be one way to organize .stl surfaces, I suppose.
        # Second, I need to call the patches.
        return out
    def _compute_computational_coordinates(self,patch,main_data,opts):
        '''
        Compute the computational coordinates of a particular boundary patch.
        
        The location of a boundary patch in computational space is useful in 
        order to reduce the computational load of dividing computational points
        among various boundary patches. However, accuracy is not particularly
        important, except to further improve efficiency. A simple extrapolation
        from the nearest point is therefore used to compute the computational
        coordinates.
                
        Args: 
           patch: the Patch object for which coordinates are being computed. In
              particular, must contain the field self.bounding_points. 
           input_points: the NumPy array containing all simulation data 
              associated with the current stream. The metric components
              A,B,C,L,M,N,P,Q,R are stored in input_points[5:14,:,:,:] (Python
              numbering), and the global coordinates X,Y,Z are stored in 
              input_points[17:20,:,:,:].
        Returns:
           patch.bounding_points_xi: computational coordinates computed using
              1st-order extrapolation of the metric components from the nearest
              computational point.
        '''
        out = []
        for bounding_point in patch.bounding_points:
            main_point_inds = bcs.findclosestpoint(bounding_point,
                                                   main_data[17:20,1:-1,
                                                             1:-1,1:-1])
            # Since main_point_inds is found using only a subset of the
            # main_data array, you don't have to do the usual index+1 stuff
            # to convert to python indexing.
            main_point = main_data[:,main_point_inds[0],main_point_inds[1],
                                   main_point_inds[2]]
            # main_point_inds happens to also closely correspond with the 
            # computational coordinates of the points. Therefore:
            out.append((main_point_inds + [opts.xi_offset,0,0])+
                       bcs.computationaldisplacement(bounding_point-
                                                     main_point[17:20],
                                                     main_point[5:14],
                                                     main_point[20]))
            if min(point[0] for point in out) > 1 and len(out)>=4 and patch.which_face is 'bottom':
                import pdb; pdb.set_trace()
        return np.array(out)
        
class Patch(object):
    '''
    An individual boundary condition; a piece of a Face.

    Each Face contains a single Patch. A Patch is a particular type of boundary
    condition, such as a pressure condition or an inflow condition.

    Args:
      init: Description of a Patch.
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
    def __init__(self,init):
        from numpy import zeros, f2py
        from sympy.utilities.codegen import codegen
        import sympy as sp
        self.type = init.type
        self.bounding_points = np.array(init.bounding_points)
        self.bounding_points_xi = None
        self.boundary_surface = init.boundary_surface
        self.flow_state = init.flow_state
        if self.boundary_surface:
            gradxsrc = codegen(
                ('gradx'+str(id(self)),
                 sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),
                                    sp.Symbol('x'))),'F95','junk')[0]
            gradysrc = codegen(
                ('grady'+str(id(self)),
                 sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),
                                    sp.Symbol('y'))),'F95','junk')[0]
            gradzsrc = codegen(
                ('gradz'+str(id(self)),
                 sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),
                                    sp.Symbol('z'))),'F95','junk')[0]
            self.gradsrc = gradxsrc[1]+'\n'+gradysrc[1]+'\n'+gradzsrc[1]+'\n'
        try:
            iftest = bool(self.flow_state)
        except ValueError: # This means that the truth value is ambiguous (array).
            iftest = bool(self.flow_state.ndim == 1)
        if iftest:
            temp = zeros((points.shape[1],points.shape[2]))
            for i in range(points.shape[1]):
                for j in range(points.shape[2]):
                    temp[:,i,j] = self.flow_state
            self.flow_state = temp

    def __call__(self,main_data,opts,out):
        '''
        Apply boundary conditions to the main data points bounding the patch
        '''
        # Sort out what points belong to the patch. Remember that main_data
        # contains ghost points. Therefore, main_data[:,i,j,k] (python indexing)
        # has xi-coordinates (i+xi_offset, j, k). So, I'll just grab any point
        # that falls within my bounding box. 
        # Identify the maximum values of xi, eta, zeta among the bounding
        # points. The minimum allowable value is 1.
#        min_xi   = max(np.ceil(min(self.bounding_points_xi[:,0])),1)
#        max_xi   = np.ceil(max(self.bounding_points_xi[:,0]))
#        min_eta  = max(np.ceil(min(self.bounding_points_xi[:,1])),1)
#        max_eta  = np.ceil(max(self.bounding_points_xi[:,1]))
#        min_zeta = max(np.ceil(min(self.bounding_points_xi[:,2])),1)
#        max_zeta = np.ceil(max(self.bounding_points_xi[:,2]))
#        # Check to make sure that your max values are within the array bounds.
#        max_xi, max_eta, max_zeta = map(lambda x,shape,minx: (x if x<=shape
#                                    else shape-1) if x > minx else minx+1,
#                                    (max_xi, max_eta, max_zeta),
#                                    (main_data.shape[1:]),(
#                                        min_xi,min_eta,min_zeta))
#        if self.which_face == 'left':
#            out = self._apply_patch(main_data[:,1,min_eta:max_eta,
#                                                min_zeta:max_zeta],opts,out)
#        elif self.which_face == 'right':
#            out = self._apply_patch(main_data[:,-2,min_eta:max_eta,
#                                                 min_zeta:max_zeta],opts,out)
#        elif self.which_face == 'bottom':
#            out = self._apply_patch(main_data[:,min_xi:max_xi,1,
#                                                min_zeta:max_zeta],opts,out)
#        elif self.which_face == 'top':
#            out = self._apply_patch(main_data[:,min_xi:max_xi,-2,
#                                                min_zeta:max_zeta],opts,out)
#        elif self.which_face == 'back':
#            out = self._apply_patch(main_data[:,min_xi:max_xi,min_eta:max_eta,
#                                                1],opts,out)
#        elif self.which_face == 'front':
#            out = self._apply_patch(main_data[:,min_xi:max_xi,min_eta:max_eta,
#                                                -2],opts,out)
        out = self._apply_patch(main_data,opts,out)
        return out
    def _apply_patch(self,main_data,opts,out):
        '''
        '''
        if self.type == 'Inflow':
            out=bcs.applyinflowconditions(main_data=main_data,
                                          bc_state=self.flow_state)
        elif self.type == 'Outflow':
            if self.flow_state is None:
                out=bcs.applyoutflowconditions(main_data=main_data)
            else:
                out=bcs.applyoutflowconditions(main_data=main_data,
                                           bc_state=self.flow_state)
        
        elif self.type == 'Transmissive':
            out = main_data
        return out
                
#if __name__ == "__main__":
#    str = 'x**2'
#    callable = autowrap(normalvec(sympify(str),'x','y','z')[0],args=('x','y','z','t'))
#    print callable(1),callable(0)

            
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
            indx[i,j] = tbctri.leadingedgepointsearch(nodes_x,faces_x,
                               face_nodes_x,leading_face_array[:,i,j])
    return(indx)

if __name__=='__main__':
    import transonic_duct
    bounds = Bounds(transonic_duct.init())
    inputs = transonic_duct.transonic_duct_inflow_generator()
    initial_array = np.empty((inputs.shape[0],1+2,
                              inputs.shape[1]+2,inputs.shape[2]+2))
    initial_array[:,1,1:-1,1:-1] = inputs
    bounds(initial_array)
#    sio.savemat('transonic_duct',{'inputs':inputs})
#    infile = open('transonic_duct_bounds.in','r')
#    test = Bounds(default_bounds='transonic_duct_bounds.in',
#                  stl_bounds_file='transonic_duct.stl',num_faces=8)
#    num_faces = 8
#    num_nodes = num_faces*3
#    [nodes,face_nodes,face_normal,error]=stl.stla_read(
#        'transonic_duct.stl',num_nodes,num_faces)
#    out= leadingEdgeSearchSTL(inputs,nodes,face_normal,face_nodes)
