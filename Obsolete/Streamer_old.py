import numpy as np
import sympy as sp
from sympy.utilities.autowrap import autowrap
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
    def __init__(self,bounds_file,init_file):
        self.bounds = Bounds(bounds_file)# Initialize boundary conditions
        self.xi_offset = 0

#class Grid(object):
#    def __init__():
#        
#    def marchStep(dt):
#        
#    def applyBounds():    

class Bounds(object):
    '''
    A complete description of the boundary conditions associated with a given
    Stream. Since all Streams are structured, the boundary conditions can be
    divided into 6 Sides, one for each face of the Stream.
    '''
    boundsList = ["Left","Right","Bottom","Top","Back","Front"]
    def __init__(self,bounds_file):
        '''
        Isolate the parts of the boundary specification pertaining to each Side
        and pass them on.
        '''
        self.boundsDict = {}
        infile = open(bounds_file,'r')
        for inda in range(6):
            infile.seek(0)
            infile.seek(infile.read().find(self.boundsList[inda]+" boundary"))
            location = infile.tell()
#            print boundsList[inda]+" boundary:"
            self.boundsDict[self.boundsList[inda]] = Side(self.readBoundary(infile,location))
#        print self.boundsDict["Left"]
    def readBoundary(self,file,location):
        file.seek(location)
        lines = [file.readline().strip(),file.readline().strip()]
        while not lines[-1].strip().endswith('boundary:'):
            lines.append(file.readline().strip())
            if not lines[-1]: break
        lines=lines[0:-1]
        return lines
    def __call__(self,main,t=0):
        print main[0,0,1:-1,1:-1]
        self.boundsDict["Left"  ](main[:, 0,1:-1,1:-1],t)
        self.boundsDict["Right" ](main[:,-1,1:-1,1:-1],t)
        self.boundsDict["Bottom"](main[:,1:-1, 0,1:-1],t)
        self.boundsDict["Top"   ](main[:,1:-1,-1,1:-1],t)
        self.boundsDict["Back"  ](main[:,1:-1,1:-1, 0],t)
        self.boundsDict["Front" ](main[:,1:-1,1:-1,-1],t)
class Side(object):
    '''
    A collection of Segment objects specific to a given face of a Stream.
    Initialized with a list of strings corresponding to the section from the
    input file that applies to the given face.
    '''
    def __init__(self,input):
        '''
        Isolate the specifications of each Segment and pass them on.
        '''
        segmentsInd = []
        for inda in range(len(input)):
            if input[inda].strip().startswith("Bounding Points:"):
                segmentsInd.append(inda)
        segmentsInd.append(len(input))
#        print len(input), len(segmentsInd)
#        print segmentsInd[0],segmentsInd[-1],input[segmentsInd[0]:segmentsInd[1]]
        self.Segments = [Segment(input[segmentsInd[inda]:segmentsInd[inda+1]])
                         for inda in range(len(segmentsInd)-1)]
    def __call__(self,main,t):
        print "Got here!\n", main[0,:,:]
        
def is_numbers(line):
    '''
    Check to see whether a string represents a number.
    '''
    for digit in line.split():
        try:
            float(digit)
            return True
        except ValueError:
            return False
def index_substring(strings, substring):
    '''
    Find all elements of a list containing a given substring.
    '''
    return next(i for i, string in enumerate(strings) if substring in string)
class Segment(object):
    '''
    Contains: 
       Points (List of lists)
       Type (string)
       Surface Specification 0 = f(x,y,z,t) (function)
       Surface Normal [dx,dy,dz] = grad(f)
       Flow State (ndarray)
    '''
    typeDict = {"Inflow":("Boundary Surface","Flow State"),
                "Outflow":("Boundary Surface","Flow State"),
                "Transmissive":(),
                "Pressure":("Flow State"),
                "Reflective":("Boundary Surface")}
    masks = {}
    masks["SupersonicIn"]  = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    masks["SubsonicIn"]    = [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    masks["SupersonicOut"] = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    masks["SubsonicOut"]   = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    masks["Pressure"] = masks["SubsonicOut"]
    masks["Communication"] = masks["SupersonicIn"]
#    funcDict = {"Inflow":inflow,
#                "Outflow":outflow,
#                "Transmissive":outflow,
#                "Pressure":pressure,
#                "Reflective":reflection}
    
    def __init__(self,input):
#        attrDict={"Boundary Surface":boundarySurface,"Flow State":Flow}
        self.Type = input.pop(0).split(":")[1].strip()
# Read in boundary surface points
        self.Points = []
        while input and is_numbers(input[0]):
            line = input.pop(0)
            self.Points.append([float(digit) for digit in line.split()])
        for attr in self.typeDict[self.Type]:
            if attr == "Boundary Surface":
                try:
                    segmentInd = index_substring(input,attr)
                    input.pop(segmentInd)
                    self.surfaceStr = input.pop(segmentInd)
                    print "surfaceStr = ",self.surfaceStr.split('=')[1].strip()
                    Xes = [sp.Symbol('x'),sp.Symbol('y'),sp.Symbol('z')]
                    symGrad = sympyGradient(
                        self.surfaceStr.split('=')[1].strip(),Xes)
                    print "sSymGrad = ",str(symGrad)
#                    Xes.append(sp.Symbol('t'))
#                    self.Gradx = str(symGrad[0])#autowrap(symGrad[0],args=Xes)
#                    self.Grady = str(symGrad[1])#autowrap(symGrad[1],args=Xes)
#                    self.Gradz = str(symGrad[2])#autowrap(symGrad[2],args=Xes)
                    self.Grad = str(symGrad)
#                    print self.normalVec(0,0,0,0)
                except StopIteration:
                    print attr+" was missing in input"
            elif attr == "Flow State":
                try:
                    segmentInd = index_substring(input,attr)
                except:
                    print "Error"
                input.pop(segmentInd)
                # This signifies a definite flow state
                if not (input[segmentInd].strip() == 'None'):
                    exec(input[segmentInd])
                    print "Hi"
                # This is for something without a state, like supersonic outflow.
                else:
                    self.state = None
                    print "None"
            else:
                print "Error, attribute "+attr+" not found!"
    def normalVec(self,x,y,z,t):
        return(eval(self.Grad))
#        return([self.Gradx(Xes[0],Xes[1],Xes[2],Xes[3]),
#                 self.Grady(Xes[0],Xes[1],Xes[2],Xes[3]),
#                 self.Gradz(Xes[0],Xes[1],Xes[2],Xes[3])])
    def __call__(self,point):
        out = numpy.ndarray(
            [var in self.state if masks[self.Type] else var in point]
            ,order='F')
        return out
#        def boundarySurface(self,input):
#            self.surfaceFunc(
#class Surface(object):
#    def __init__(self,input):
#        self.function = input
#        print "Surface function: 0 = "+self.function

#class functionXYZT(str,x,y,z,t):
#    exec(str)
#    return f

def sympyGradient(str,Xes):
#    Xes=(sp.Symbol(x_in),sp.Symbol(y_in),sp.Symbol(z_in))
    [dx,dy,dz] = [sp.diff(sp.sympify(str),X) for X in Xes]
    return([dx,dy,dz])

def allindex(value,list):
    indices = []
    indx = -1
    while 1:
        try:
            indx = list.index(value,indx+1)
            indices.append(indx)
        except ValueError:
            break
    return indices

if __name__ == "__main__":
    test = Stream('riemann1d_bounds.in','riemann1d_init.in')
    main = np.zeros((4,6,6,6),order='F')
    main[:,:,:,0] = 1
    main[:,:,:,-1] = 1
    main[:,:,0,:] = 1
    main[:,:,-1,:] = 1
#    main[:,0,:,:] = 1
#    main[:,-1,:,:] = 1
    
    print main[0,:,:,:]
    test.bounds(main,0)



