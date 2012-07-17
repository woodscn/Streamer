import numpy as np
from numpy import f2py
import scipy.io as sio
import STLA_IO as stl
import BoundaryConditionsStuff
from sys import exit as sysexit
from main import Stream
bcs = BoundaryConditionsStuff.boundaryconditionsstuff
import os

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
        fortran_normal_src +='''
module FortranNormalVectors
contains
subroutine ApplyReflectionConditions(main_data,patch_id,out,nx,ny,dim,t)
use BoundaryConditionsStuff, only: WallReflect
implicit none
integer, intent(in) :: nx, ny, dim
real(8), intent(out), dimension(21,nx,ny) :: out
real(8), intent(in), dimension(21,nx,ny) :: main_data
integer(8), intent(in) :: patch_id
integer :: i, j
real(8) :: x, y, z, t
intent(in) :: t
real(8), dimension(3) :: normal
do i = 1, size(main_data,2)
do j = 1, size(main_data,3)
x = main_data(18,i,j)
y = main_data(19,i,j)
z = main_data(20,i,j)
select case (patch_id)'''
        for face in self:
            for patch in face:
                try:
                    fortran_normal_src += patch.gradsrc
                except AttributeError: # Some patches don't have normal vectors
                    pass
        fortran_normal_src +='''
end select
out(:,i,j) = WallReflect(main_data(:,i,j), normal,&
reshape([0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],[3,3]),&
int(dim,4))
end do
end do
end subroutine ApplyReflectionConditions
end module FortranNormalVectors
'''
        if os.access("FortranNormalVectors.f90",os.F_OK):
            os.remove("FortranNormalVectors.f90")
            os.remove("FortranNormalVectors.so")
        f2py.compile(fortran_normal_src,modulename='FortranNormalVectors',
                     verbose=False,source_fn='FortranNormalVectors.f90',
                     extra_args='--quiet --f90flags=-Wno-unused-dummy-argument')
    def __iter__(self):
        faces = (self.left_face,self.right_face,self.bottom_face,self.top_face,
                 self.back_face,self.front_face)
        for face in faces:
            yield face
    def __call__(self,main,opts,t=0):
        '''
Call each Face object with the appropriate subset of your return array.
The main_data array is passed in it's entirety, to simplify any future
higher-order boundary condition implementation that may require more
than just the first interior point.

Ex:
main[:,0,1:-1,1:-1] = left_face(main_data = main[:,1:-1,1:-1,1:-1],
                                time = t,
                                out = main[:, 0,1:-1,1:-1],
                                options = opts)
        
Args:
  self: Bounds object containing the initialized boundary conditions.
  main: Array object containing the state information for the Stream.
  t: Simulation time; used for time-dependent boundary conditions.

Returns:
  self: Updated boundary computational coordinates.
  main: Updated Array of state variables.'''
        opts.t = t
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
        dim_dict = {'left':1,'right':1,'bottom':2,'top':2,
                    'back':3,'front':3}
        for patch in face_init:
            self.patches.append(Patch(patch))
        for patch in self.patches:
            patch.dim = dim_dict[self.which_face]
            patch.which_face = self.which_face
        if self.which_face == "left" or self.which_face == "right":
            # Sort in ascending order, first by y and then by z coordinates.
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,2]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,1]))
        elif self.which_face == "bottom" or self.which_face == "top":
            # Sort first by x, then by z.
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,2]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,0]))
        elif self.which_face == "back" or self.which_face == "front":
            # Sort first by x, then by y
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,1]))
            self.patches.sort(key=lambda patch: min(patch.bounding_points[:,0]))
    def __iter__(self):
        for patch in self.patches:
            yield patch
    def __call__(self,main_data,t,out,opts):
        '''
Call each Patch object with the appropriate subset of main and
the appropriate subset of your return array.

Since the computational coordinates of a given point are intimately
connected with the array indices of that point, this is done by first
computing the computational coordinates of each Patch, and then using
those coordinates to pass along the appropriate array subset.

The current implementation is only valid for 2-dimensional simulations.'''
        prev_pri_max = 0
        prev_sec_max = 0
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
            else:
                # Eventually, just time-advance the points. For now, do the same
                # as above
                ccc = self._compute_computational_coordinates
                patch.bounding_points_xi = ccc(patch,main_data,opts)
            # Next, I need to pass the appropriate grid points to the patch.
            # Only pass the points that lie completely within the patch.
            # I believe that using only computational coordinates may obviate
            # the problems with overlapping boundaries and such.
            xi_min   =  np.ceil(min(patch.bounding_points_xi[:,0]))+opts.xi_offset
            xi_max   = np.floor(max(patch.bounding_points_xi[:,0]))+opts.xi_offset
            eta_min  =  np.ceil(min(patch.bounding_points_xi[:,1]))
            eta_max  = np.floor(max(patch.bounding_points_xi[:,1]))
            zeta_min =  np.ceil(min(patch.bounding_points_xi[:,2]))
            zeta_max = np.floor(max(patch.bounding_points_xi[:,2]))
            # Because of the way Python indexes arrays, we must add 1 to each
            # max value.
            if self.which_face == "left":
                main_pass = main_data[:,0,eta_min:eta_max+1,zeta_min:zeta_max+1]
                out_pass = out[:,eta_min:eta_max+1,zeta_min:zeta_max+1]
                key_ind, key_ind_max = eta_min, eta_max
            elif self.which_face == "right":
                if len(self.patches) > 1:
                    print "Multiple outflow patches are currently unsupported."
                    raise(Error)
                main_pass = main_data[:,-1,:,:]
                out_pass = out[:,:,:]
                key_ind, key_ind_max = eta_min, eta_max
            elif self.which_face == "bottom":
                main_pass = main_data[:,xi_min:xi_max+1,0,zeta_min:zeta_max+1]
                out_pass = out[:,xi_min:xi_max+1,zeta_min:zeta_max+1]
                key_ind, key_ind_max = xi_min, xi_max
            elif self.which_face == "top":
                main_pass = main_data[:,xi_min:xi_max+1,-1,zeta_min:zeta_max+1]
                out_pass = out[:,xi_min:xi_max+1,zeta_min:zeta_max+1]
                key_ind, key_ind_max = xi_min, xi_max
            elif self.which_face == "back":
                main_pass = main_data[:,xi_min:xi_max+1,eta_min:eta_max+1,0]
                out_pass = out[:,xi_min:xi_max+1,eta_min:eta_max+1]
                key_ind, key_ind_max = xi_min, xi_max
            elif self.which_face == "front":
                main_pass = main_data[:,xi_min:xi_max+1,eta_min:eta_max+1,-1]
                out_pass = out[:,xi_min:xi_max+1,eta_min:eta_max+1]
                key_ind, key_ind_max = xi_min, xi_max
            else:
                raise(Error)
            if self.which_face == "left" or self.which_face == "right":
                pass
            elif key_ind > prev_pri_max+1:
                # Average of some kind. For now, just raise an error.
                print "Side boundary patches do not match up perfectly."
                import pdb; pdb.set_trace()
                raise Error()
                sysexit()
            if main_pass.size>0:
		out_pass[:,:,:] = patch(main_pass,opts,out_pass)
                # What happens here? Do I get all of the patches?
                # import pdb;pdb.set_trace()
            prev_pri_max = key_ind_max
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
              associated with the dcurrent stream. The metric components
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
            # main_point_inds returns Fortran-style indices that must be
            # adjusted for use in Python. 
            main_point_inds = bcs.findclosestpoint(bounding_point,
                                                   main_data[17:20,:,:,:]) - 1
            try:
                main_point = main_data[:,main_point_inds[0],main_point_inds[1],
                                       main_point_inds[2]]
            except:
                import pdb;pdb.set_trace()
            # main_point_inds happens to also closely correspond with the 
            # computational coordinates of the points. Therefore:
            main_point_xi_coords = main_point_inds - [opts.xi_offset,0,0]
            displacement = bcs.computationaldisplacement(
                bounding_point-main_point[17:20],main_point[5:14],main_point[20]
                )
            out.append(main_point_xi_coords + displacement)
            # It is also possible here to do some kind of higher-order
            # extrapolation in order to compute the computational
            # coordinates of the patch bounding points.
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
        # from sympy.utilities.codegen import codegen
        from sympy.printing.fcode import fcode
        import sympy as sp
        self.type = init.type
        self.bounding_points = np.array(init.bounding_points)
        self.bounding_points_xi = None
        self.boundary_surface = init.boundary_surface
        self.flow_state = init.flow_state
        if self.boundary_surface:
            # gradxsrc = codegen(
            #     ('gradx'+str(id(self)),
            #      sp.diff(sp.sympify(
            #          self.boundary_surface.split('=')[1].strip()),
            #          sp.Symbol('x'))),'F95','junk',
            #          argument_sequence=(sp.Symbol('x'),sp.Symbol('y'),
            #                             sp.Symbol('z'),sp.Symbol('t')))[0]
            # gradysrc = codegen(
            #     ('grady'+str(id(self)),
            #      sp.diff(sp.sympify(
            #          self.boundary_surface.split('=')[1].strip()),
            #          sp.Symbol('y'))),'F95','junk',
            #          argument_sequence=(sp.Symbol('x'),sp.Symbol('y'),
            #                             sp.Symbol('z'),sp.Symbol('t')))[0]
            # gradzsrc = codegen(
            #     ('gradz'+str(id(self)),
            #      sp.diff(sp.sympify(
            #          self.boundary_surface.split('=')[1].strip()),
            #          sp.Symbol('z'))),'F95','junk',
            #          argument_sequence=(sp.Symbol('x'),sp.Symbol('y'),
            #                             sp.Symbol('z'),sp.Symbol('t')))[0]
            # self.gradsrc = gradxsrc[1]+'\n'+gradysrc[1]+'\n'+gradzsrc[1]+'\n'
            self.gradsrc = '''
case('''+str(id(self))+'''_8)
normal=[real('''+fcode(sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),sp.Symbol('x')),source_format='free')+',8),real('+fcode(sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),sp.Symbol('y')),source_format='free')+',8),real('+fcode(sp.diff(sp.sympify(self.boundary_surface.split('=')[1].strip()),sp.Symbol('z')),source_format='free')+''',8)]'''
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
        out[:,:,:] = self._apply_patch(main_data,opts,out)
        return out
    def _apply_patch(self,main_data,opts,out):
        '''
        '''
        # The [:,:,:] is necessary to keep the array from being reassigned
        # to a different object in memory.
        if self.type == 'Inflow':
            out[:,:,:]=bcs.applyinflowconditions(main_data=main_data,
                                          bc_state=self.flow_state)
        elif self.type == 'Outflow':
            if self.flow_state is None:
                out[:,:,:]=bcs.applyoutflowconditions(main_data=main_data)
            else:
                out[:,:,:]=bcs.applyoutflowconditions(main_data=main_data,
                                               bc_state=self.flow_state)
        elif self.type == 'Transmissive':
            out[:,:,:] = main_data
        elif self.type == 'SolidWall':
            import FortranNormalVectors
            fnv = FortranNormalVectors.fortrannormalvectors
            # ApplyReflectionConditions(main_data,patch_id,out,nx,ny,dim,t)
            try:
                out[:,:,:] = fnv.applyreflectionconditions(main_data,str(id(self)),
                                                           self.dim,opts.t)
            except ValueError:
                import pdb;pdb.set_trace()
        else:
            print "Undefined patch type!"
            raise(Error)
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
