import numpy as np
import stla_io as stl

def find_boundary(P):
    testx = 10
    CUTOFF = .01
    for inda in range(num_faces):
        [nodea,nodeb,nodec] = np.transpose(nodes[:,face_nodes[:,inda]])
        normal = face_normal[:,inda]
        v3 = np.dot((P - nodea),normal)
        v2 = (P - nodea) - v3*normal
        v0 = (nodec - nodea)
        v1 = (nodeb - nodea)
        u = (
            np.dot(v1,v1)*np.dot(v2,v0)-np.dot(v1,v0)*np.dot(v2,v1)
            )/(
                np.dot(v0,v0)*np.dot(v1,v1)-np.dot(v0,v1)**2
                )
        v = (
            np.dot(v0,v0)*np.dot(v2,v1)-np.dot(v1,v0)*np.dot(v2,v0)
            )/(
                np.dot(v0,v0)*np.dot(v1,v1)-np.dot(v0,v1)**2
                )
        test1 = v3**2
        if test1 < CUTOFF:
            test2 = (u-.3333333333333333)**2+(v-.3333333333333333)**2
            if test2 < testx:
                testx = test2
                indx = inda

    return(indx)



def computationalGrads(p):
    '''
    Inputs: point = single element of primitive variables array:
     0   1   2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20
    [p, rho, u, v, w, A, B, C, L, M, N, P, Q, R, U, V, W, X, Y, Z, J]
    '''
    gradXi   = (( p[ 9]*p[13] - p[10]*p[12] )/p[20],
                ( p[10]*p[11] - p[ 8]*p[13] )/p[20],
                ( p[ 8]*p[12] - p[ 9]*p[11] )/p[20])
    gradEta  = (( p[ 7]*p[12] - p[ 6]*p[13] )/p[20],
                ( p[ 5]*p[13] - p[ 7]*p[11] )/p[20],
                ( p[ 6]*p[11] - p[ 5]*p[12] )/p[20])
    gradZeta = (( p[ 6]*p[10] - p[ 7]*p[ 9] )/p[20],
                ( p[ 7]*p[ 8] - p[ 5]*p[10] )/p[20],
                ( p[ 5]*p[ 9] - p[ 6]*p[ 8] )/p[20])
    return(gradXi, gradEta, gradZeta)

def compute_computational_coordinates(base,point):
    gradXi,gradEta,gradZeta = computationalGrads(base)
    diff = point-base[17:20]
    return( sum(gradXi*diff), sum(gradEta*diff), sum(gradZeta*diff) )

#def build_boundary_mask(in):
#    mask = np.ones(
#    mask[:,1:-1,1:-1,1:-1] = 

if __name__=='__main__':
    from round_duct_bounds import round_duct_init as rdinit
    inputs = rdinit()
    num_faces = 506
    num_nodes = num_faces*3
    [nodes,face_nodes,face_normal,error]=stl.stla_read(
        'round_duct.stl',num_nodes,num_faces)
    face_nodes = face_nodes - 1
    primaries = np.zeros((21,1,25,25))
    primaries[:,0,:,:] = inputs
    primariesplus = np.zeros((21,3,27,27))
#    primariesplus[:,1,0,:]
    for inda in range(len(primaries[0,0,0,:])):
        primariesplus[:,1,0,inda] = primaries[:,0,0,inda].copy()
#        print compute_computational_coordinates(primaries[:,0,0,inda],
#                                                [1.,1.,0.])
    indx = find_boundary([.01,.99,-.99])
    print indx
    print 'Nodes = '
    print nodes[:,face_nodes[:,indx]]
    print 'Normal vector = '
    print face_normal[:,indx]
    print 'Computational coordinates (xi, eta, zeta) = '
    print compute_computational_coordinates(inputs[:,0,0],[1.0,1.0,0.0])
    import BoundaryConditions as bc
    print bc.triangulatedboundaryconditions.leadingedgesearch(nodes,face_normal,face_nodes,[10.,.99,-.99],[1.])
