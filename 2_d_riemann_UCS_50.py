import numpy as np
import scipy.optimize as opt

class PatchInit:
    def __init__(self,t,bp,bs=None,fs=None):
        self.type = t
        self.bounding_points = bp
        self.boundary_surface = bs
        self.flow_state = fs
class BoundsInit:
    def __init__(self,lf,rf,bof,tf,baf,ff):
        self.left_face = lf
        self.right_face = rf
        self.bottom_face = bof
        self.top_face = tf
        self.back_face = baf
        self.front_face = ff
def init():
    left_face_init = []
    right_face_init = []
    top_face_init = []
    bottom_face_init = []
    back_face_init = []
    front_face_init = []
#    inflow_condition = inflow_generator()
    initial_conds = initial_condition()
    left_face_init.append(
        PatchInit('Inflow',
                  ((0.0,0.0,0.0),(0.0,0.0,1.0),(0.0,1.0,1.0),(0.0,1.0,0.0)),
                  'f = x',
                  initial_conds[:,0,:,:]))
    right_face_init.append(
        PatchInit('Outflow',
                  ((.6,-1.0,-1.0),(.6,2.0,-1.0),(.6,2.0,2.0),(.6,-1.0,2.0)),
                  'f = -x+.6', 
                  None))
    bottom_face_init.append(
        PatchInit('Transmissive',
                  ((-1.0,0.0,0.0),(1.0,0.0,0.0),(1.0,0.0,1.0),(-1.0,0.0,1.0))))
    top_face_init.append(
        PatchInit('Transmissive',
                  ((-1.0,1.0,0.0),(0.8,1.0,0.0),(0.8,1.0,1.0),(-1.0,1.0,1.0))))
    back_face_init.append(
        PatchInit('Transmissive',
                  ((-1.0,0.0,0.0),(0.8,0.0,0.0),(0.8,1.0,0.0),(-1.0,1.0,0.0))))
    front_face_init.append(
        PatchInit('Transmissive',
                  ((-1.0,0.0,1.0),(0.8,0.0,1.0),(0.8,1.0,1.0),(-1.0,1.0,1.0))))
    bounds_init = BoundsInit(lf=left_face_init,rf=right_face_init,
                             bof=bottom_face_init,tf=top_face_init,
                             baf=back_face_init,ff=front_face_init)

    exact_sol_kwargs = {'pR':.25,'dR':.5,'MR':7,'alphaR':np.pi*.5,
                        'pL':1.,'dL':1.,'ML':2.4,'alphaL':np.pi*.5,
                        'gamma':1.4}
    exact_solution_obj = SteadyRiemannSolution(**exact_sol_kwargs)
    exact_solution = lambda t,indxi,indeta,indzeta :(
        np.array(
            list(exact_solution_obj(
                        initial_conds[17,indxi,indeta,indzeta],
                        initial_conds[18,indxi,indeta,indzeta]-.5))
            +list(initial_conds[5:,indxi,indeta,indzeta])
            )
        )

    solver_options = np.zeros(300)
    solver_options[0] = 1
    solver_options[100] =  1
    solver_options[101] = 1
    solver_options[102] = 1
    solver_options[103] = 0
    stream_options = {
        'solver_type':'euler',
        'boundary_layers':False,
        'multistream':False,
        'solver_options':solver_options,
        'manufactured':False,
        'source_funcs':None,
        'exact_sol_func':exact_solution}
    return bounds_init, initial_conds, stream_options

def initial_condition():
    nx,ny,nz = 30,50,1
    xmin,xmax = 0.,.6
    ymin,ymax = -.5,.5
    zmin,zmax = 0.,0.
    dx,dy,dz = (xmax-xmin)/(nx),(ymax-ymin)/(ny),1
    inputs = np.zeros((21,nx,ny,nz))
    inputs[0,:,:,:] = 1.
    inputs[0,:,25:,:] = .25
    inputs[1,:,:,:] = 1.
    inputs[1,:,25:,:] = .5
    inputs[2,:,:,:] = 2.4*np.sqrt(1.4)
    inputs[2,:,25:,:] = 7*np.sqrt(.25/.5*1.4)
    inputs[5,:,:,:] = dx
    inputs[9,:,:,:] = dy
    inputs[13,:,:,:] = dz
    inputs[14:17,:,:,:] = 0.5*inputs[2:5,:,:,:]
    for inda in range(ny):
        for indb in range(nx):
            inputs[17,indb,inda,0] = (indb+.5)*dx
            inputs[18,indb,inda,0] = (inda+.5)*dy 
    inputs[20,:,:,:] = dx*dy*dz
    return inputs
    
class SteadyRiemannSolution(object):
    def __init__(self,pL,dL,ML,alphaL,pR,dR,MR,alphaR,gamma):
        self.pL,self.dL,self.ML,self.alphaL = pL,dL,ML,alphaL
        self.pR,self.dR,self.MR,self.alphaR = pR,dR,MR,alphaR
        self.gamma = gamma
        (self.pstar,self.alphastar,self.dstarL,self.dstarR,
         self.MstarL,self.MstarR) = self.star_state()
        if self.pstar/self.pL <= 1.:
            self.left_wave = {
                'type':'fan',
                'head_angle':self.left_fan_angle(self.MstarL,self.alphastar),
                'tail_angle':self.left_fan_angle(self.ML,self.alphaL)}
        else:
            self.left_wave = {
                'type':'shock',
                'shock_angle':self.left_shock_angle(
                    self.pstar,self.pL,self.alphaL,self.ML)}
        if self.pstar/self.pR <=1.:
            self.right_wave = {
                'type':'fan',
                'head_angle':self.right_fan_angle(self.MstarR,self.alphastar),
                'tail_angle':self.right_fan_angle(self.MR,self.alphaR)}
        else:
            self.right_wave = {
                'type':'shock',
                'shock_angle':self.right_shock_angle(
                    self.pstar,self.pR,self.alphaR,self.MR)}

    def __call__(self,x,y):
        theta = np.arctan(y/x)
        p,d,M,alpha = self.sample(theta)
        a = (p/d*self.gamma)**.5
        umag = M*a
        u = umag*np.sin(alpha)
        v = umag*np.cos(alpha)
        return p,d,u,v

    def sample(self,theta):
        if theta <= np.pi*.5-self.alphastar: #Left of slip line
            if self.left_wave['type'] == 'fan':
                if theta <= self.left_wave['tail_angle']:
                    return self.pL,self.dL,self.ML,self.alphaL
                else:
                    if theta >= self.left_wave['head_angle']: #Left star region
                        return self.pstar,self.dstarL,self.MstarL,self.alphastar
                    else: #Inside left fan
                        p_fan,alpha_fan,d_fan,M_fan = self.left_fan_state(
                            theta,self.pL,self.dL,self.ML,self.alphaL)
                        return p_fan,d_fan,M_fan,alpha_fan
            else: #Left shock 
                if theta <= self.left_wave['shock_angle']: 
                    return self.pL,self.dL,self.ML,self.alphaL
                else:
                    return self.pstar,self.dstarL,self.MstarL,self.alphastar
        else: #Right of slip line
            if self.right_wave['type'] == 'fan':
                if theta >= self.right_wave['tail_angle']:
                    return self.pR,self.dR,self.MR,self.alphaR
                else:
                    if theta <= self.right_wave['head_angle']:
                        return self.pstar,self.dstarR,self.MstarR,self.alphastar
                    else:
                        p_fan,alpha_fan,d_fan,M_fan = self.right_fan_state(
                            theta,self.pR,self.dR,self.MR,self.alphaR)
                        return p_fan,d_fan,M_fan,alpha_fan
            else: #Right shock
                if theta >= self.right_wave['shock_angle']:
                    return self.pR,self.dR,self.MR,self.alphaR
                else:
                    return self.pstar,self.dstarR,self.MstarR,self.alphastar

    def Prandtl_Meyer(self,M):
        out = (((self.gamma+1)/(self.gamma-1))**.5*
               np.arctan(((self.gamma-1)/(self.gamma+1)*(M**2-1))**.5)
               -np.arctan((M**2-1)**.5))
        return out
    
    def h_func(self,eta):
        if eta <= 1.:
            out = eta**(1./self.gamma)
        else:
            out = ((1+.5*(self.gamma+1)/self.gamma*(eta-1))/
                   (1+.5*(self.gamma-1)/self.gamma*(eta-1)))
        return out

    def g_func(self,eta,M0):
        out = (2/(self.gamma-1)*(self.h_func(eta)/
                                 eta*(1+.5*(self.gamma-1)*M0**2)-1))**.5
        return out

    def D_func(self,eta,M0):
        out = (np.arcsin(1./self.g_func(eta,M0)*
                         (1+(self.gamma+1)/self.gamma*.5*(1./eta-1))**.5)
               - np.arcsin(1./M0*(1+(self.gamma+1)/self.gamma*.5*(eta-1))**.5))
        return out

    def f_func(self,eta,M0):
        if eta<=0.:
            out = self.Prandtl_Meyer(self.gfunc(eta,M0))-self.Prandtl_Meyer(M0)
        else:
            out = self.D_func(eta,M0)
        return out

    def pressure_func(self,p):
        return (self.f_func(p/self.pL,self.ML)+
                self.f_func(p/self.pR,self.MR)+
                self.alphaR-self.alphaL)

    def pressure_guess(self):
        return .5*(self.pL+self.pR)
    
    def pressure_solve(self):
        guess = self.pressure_guess()
        return opt.newton(self.pressure_func,guess)

    def star_state(self):
        pstar = self.pressure_solve()
        alphastar = self.f_func(pstar/self.pR,self.MR)+self.alphaR
        dstarL = self.h_func(pstar/self.pL)*self.dL
        dstarR = self.h_func(pstar/self.pR)*self.dR
        MstarL = self.g_func(pstar/self.pL,self.ML)
        MstarR = self.g_func(pstar/self.pR,self.MR)
        return (pstar,alphastar,dstarL,dstarR,MstarL,MstarR)

    def left_fan_state(self,theta,p,d,M,alpha):
        return self.fan_state(theta,p,d,M,alpha,'left')

    def right_fan_state(self,theta,p,d,M,alpha):
        return self.fan_state(theta,p,d,M,alpha,'right')

    def fan_state(self,theta,p,d,M,alpha,leftright_in):
        leftright = {'left':-1,'right':1}
        M_fan = opt.newton(self.fan_mach_func,M,args=(theta,M,alpha,leftright_in))
        p_fan = p*opt.newton(self.fan_eta_func,p,args=(M_fan,M))
        d_fan = d*(p_fan/p)**(1./self.gamma)
        alpha_fan = alpha + leftright[leftright_in]*(
            self.Prandtl_Meyer(M_fan)-self.Prandtl_Meyer(M))
        return p_fan,alpha_fan,d_fan,M_fan
    
    def fan_mach_func(self,M,theta,M0,alpha0,leftright_in):
        leftright = {'left':-1,'right':1}
        return (np.pi*.5-alpha0)-theta-leftright[leftright_in]*(
            self.Prandtl_Meyer(M)-self.Prandtl_Meyer(M0)-np.arcsin(1./M))

    def fan_eta_func(self,eta,M,M0):
        return (2./(self.gamma-1)*(eta**((1-self.gamma)/self.gamma)*
                                   (1+.5*(self.gamma-1)*M0**2)-1))**.5-M
    
    def left_shock_angle(self,p,p0,alpha0,M0):
        return np.pi*.5-alpha0-np.arcsin(((p/p0-1)*.5*(self.gamma+1)
                                          /self.gamma+1)**.5/M0)

    def right_shock_angle(self,p,p0,alpha0,M0):
        return np.pi*.5-alpha0+np.arcsin(((p/p0-1)*.5*(self.gamma+1)/
                                          self.gamma+1)**.5/M0)

    def left_fan_angle(self,M,alpha):
        return np.pi*.5-alpha-np.arcsin(1./M)

    def right_fan_angle(self,M,alpha):
        return np.pi*.5-alpha+np.arcsin(1./M)

if __name__=='__main__':
    kwargs = {'pR':.25,'dR':.5,'MR':7,'alphaR':np.pi*.5,
              'pL':1.,'dL':1.,'ML':2.4,'alphaL':np.pi*.5,
              'gamma':1.4}
    sol = SteadyRiemannSolution(**kwargs)
#    import pdb;pdb.set_trace()
    inputs = transonic_duct_inflow_generator()
#    for inda in range(inputs.shape[1]):
#        print inputs[5:15,inda,0]
    import numpy as np
    import matplotlib
#    print matplotlib.__version__
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    
    def randrange(n, vmin, vmax):
        return (vmax-vmin)*np.random.rand(n) + vmin
    
    fig = plt.figure()
    ax = Axes3D(fig)
    n = 100
    for i in range(inputs.shape[1]):
        ax.scatter(inputs[17,i,:],inputs[18,i,:],inputs[19,i,:])            
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()
    print init.__doc__
