import numpy
np = numpy
import scipy.optimize as opt
class PatchInit:
    def __init__(self,t,bp,dim,bs=None,fs=None):
        self.type = t
        self.bounding_points = bp
        self.boundary_surface = bs
        self.flow_state = fs
        self.dim = dim
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
    initial_conds = initial_conditions()
    solver_options = numpy.zeros(300)
    solver_options[0] = 1
    solver_options[100] = 1
    solver_options[101] = 1
    solver_options[102] = 1
    solver_options[103] = 0
    solver_options[104] = 1
    exact_solution_obj = UnsteadyRiemannSolution()
    exact_solution = lambda t,x,y,z:(
        np.array(
            list(exact_solution_obj(x,t))+[0 for ind in range(17)]))
                    
    stream_options = {
        'solver_type':'euler',
        'boundary_layers':False,
        'multistream':False,
        'solver_options':solver_options,
        'manufactured':False,
        'source_funcs':None,
        'exact_sol_func':exact_solution
        }
    left_face_init.append(
        PatchInit('Transmissive',
                  ((-.5,-.1,-.1),(-.5,-.1,.1),(-.5,.1,0.1),(-.5,.1,-.1)),1))
    right_face_init.append(
        PatchInit('Transmissive',
                  ((0.5,.1,-.1),(0.5,-.1,-.1),(0.5,-.1,0.1),(0.5,.1,.1)),1))
    bottom_face_init.append(
        PatchInit('Transmissive',
                  ((-.5,-.1,-.1),(0.5,-.1,-.1),(0.5,0.1,0.1),(-.5,0.1,0.1)),2))
    top_face_init.append(
        PatchInit('Transmissive',
                  ((-.5,0.1,-.1),(0.5,0.1,-.1),(0.5,0.1,0.1),(-.5,0.1,0.1)),2))
    back_face_init.append(
        PatchInit('Transmissive',
                  ((-.5,-.1,-.1),(0.5,-.1,-.1),(0.5,.1,-.1),(-.5,.1,-.1)),3))
    front_face_init.append(
        PatchInit('Transmissive',
                  ((-.5,-.1,0.1),(0.5,-.1,0.1),(0.5,.1,0.1),(-.5,.1,0.1)),3))
    bounds_init = BoundsInit(lf=left_face_init,rf=right_face_init,
                             bof=bottom_face_init,tf=top_face_init,
                             baf=back_face_init,ff=front_face_init)
    return bounds_init, initial_conds, stream_options

def initial_conditions():
    import numpy as np
    inputs = np.zeros((21,101,2,2))
    inputs[0,:,:] = 1.
    inputs[1,:,:] = 1.
    inputs[2,:,:] = .75
    inputs[3,:,:] = 0.
    inputs[4,:,:] = 0.
    inputs[ 5,:,:] = 0.01
    inputs[ 6,:,:] = 0.
    inputs[ 7,:,:] = 0.
    inputs[ 8,:,:] = 0.
    inputs[ 9,:,:] = 1.
    inputs[10,:,:] = 0.
    inputs[11,:,:] = 0.
    inputs[12,:,:] = 0.
    inputs[13,:,:] = 1.
    inputs[14,:,:] = .25*inputs[2,:,:]
    inputs[15,:,:] = 0.
    inputs[16,:,:] = 0.
    inputs[17,:,0,0] = [-.5 + 0.01*ind for ind in range(101)]
    inputs[18,:,:] = 0.
    inputs[19,:,:] = 0.
    inputs[20,:,:] = 0.01
    inputs[0,25:,:] = .1
    inputs[1,25:,:] = .125
    inputs[2,25:,:] = 0
    # np.savetxt("transonic_duct_bounds_x.txt",np.reshape(inputs[17,:,:],-1))
    # np.savetxt("transonic_duct_bounds_y.txt",np.reshape(inputs[18,:,:],-1))
    # np.savetxt("transonic_duct_bounds_z.txt",np.reshape(inputs[19,:,:],-1))
    # np.save("transonic_duct_bounds",inputs)
    return(inputs)

class UnsteadyRiemannSolution(object):
    def __init__(self,pL=1,dL=1,uL=0,pR=.1,dR=.125,uR=0,gamma=1.4):
        self.pL,self.dL,self.uL = pL,dL,uL
        self.pR,self.dR,self.uR = pR,dR,uR
        self.gamma = gamma
        self.pstar,self.ustar,self.dstarL,self.dstarR = self.star_state()
        self.wave_speeds()
#        if self.pstar/self.pL <= 1.:
#            self.left_wave = {
#                'type':'fan',
#                'head_speed':self.left_fan_speed(),
#                'tail_speed':self.left_fan_speed()}
#        else:
#            self.left_wave = {
#                'type':'shock',
#                'shock_speed':self.left_shock_speed()}
#        if self.pstar/self.pR <= 1.:
#            self.right_wave = {
#                'type':'fan',
#                'head_speed':self.right_fan_speed(),
#                'tail_speed':self.right_fan_speed()}
#        else:
#            self.right_wave = {
#                'type':'shock'
#                'shock_speed':self.right_shock_speed()}

    def __call__(self,x,t):
        if t > 0:
            s = x/t
        else:
            s = 100000000000
        p,d,u = self.sample(s)
        return p,d,u

    def sample(self,s):
        if s <= self.ustar: #Left of slip line
            if self.left_wave['type'] == 'fan':
                if s <= self.left_wave['tail_speed']:
                    return self.pL,self.dL,self.uL
                else:
                    if s >= self.left_wave['head_speed']:
                        return self.pstar,self.dstarL,self.ustar
                    else:
                        p_fan,d_fan,u_fan = self.left_fan_state(s)
                        return p_fan,d_fan,u_fan
            else:
                if s <= self.left_wave['shock_speed']:
                    return self.pL,dL,uL
                else:
                    return self.pstar,self.dstarL,self.ustar
        else:
            if self.right_wave['type'] == 'fan':
                if s >= self.right_wave['tail_speed']:
                    return self.pR,self.dR,self.uR
                else:
                    if s <= self.right_wave['head_speed']:
                        return self.pstar,self.dstarR,self.ustar
                    else:
                        p_fan,d_fan,u_fan = self.right_fan_state(s)
                        return p_fan,d_fan,u_fan
            else:
                if s >= self.right_wave['shock_speed']:
                    return self.pR,self.dR,self.uR
                else:
                    return self.pstar,self.dstarR,self.ustar

    def star_state(self):
        pstar = self.pressure_solve()
        ustar = .5*(self.uR+self.u_fun_right(pstar)+self.uL-self.u_fun_left(pstar))
        dstarL = self.dL*self.beta(pstar/self.pL)
        dstarR = self.dR*self.beta(pstar/self.pR)
        return pstar, ustar, dstarL, dstarR

    def pressure_guess(self):
        return (self.pL+self.pR)*.5
    
    def pressure_func(self,p):
        return (self.uR-self.uL+self.u_fun_right(p)+self.u_fun_left(p))
    
    def pressure_solve(self):
        guess = self.pressure_guess()
        return opt.newton(self.pressure_func,guess)
    
    def u_fun_left(self,pstar):
        return self.u_fun(pstar,-1)

    def u_fun_right(self,pstar):
        return self.u_fun(pstar,1)

    def u_fun(self,pstar,leftright):
        if leftright == -1:
            p0,d0,u0 = self.pL,self.dL,self.uL
        else:
            p0,d0,u0 = self.pR,self.dR,self.uR
        psi = pstar/p0
        a0 = np.sqrt(self.gamma*p0/d0)
        if psi > 1:
            A = 2/((self.gamma+1)*d0)
            B = (self.gamma-1)/(self.gamma+1)*p0
            f = p0*(psi-1)*np.sqrt(A/(p0*psi+B))
        else:
            f = 2*a0/(self.gamma-1)*(psi**((self.gamma-1)/(2*self.gamma))-1)
        return f

    def beta(self,psi):
        if(psi>1):
            return ((self.gamma+1)*psi+self.gamma-1)/(self.gamma+1+(self.gamma-1)*psi)
        else:
            return psi**(1/self.gamma)

    def wave_speeds(self):
        self.left_wave = {}
        self.right_wave = {}
        if self.pstar/self.pL>1:
            self.left_wave['type'] = 'shock'
            self.left_wave['shock_speed'] = (self.uL-np.sqrt(self.pL*self.gamma/self.dL)
                                             *np.sqrt((self.gamma+1)/(2*self.gamma)
                                                   *(self.pstar/self.pL-1)+1))
        else:
            self.left_wave['type'] = 'fan'
            self.left_wave['head_speed'] = self.ustar-np.sqrt(
                self.pstar*self.gamma/self.dstarL)
            self.left_wave['tail_speed'] = self.uL-np.sqrt(
                self.pL*self.gamma/self.dL)

        if self.pstar/self.pR>1:
            self.right_wave['type'] = 'shock'
            self.right_wave['shock_speed'] = (self.uR+np.sqrt(self.pR*self.gamma/self.dR)
                                              *np.sqrt((self.gamma+1)/(2*self.gamma)
                                                    *(self.pstar/self.pR-1)+1))
        else:
            self.right_wave['type'] = 'fan'
            self.right_wave['head_speed'] = self.ustar+np.sqrt(
                self.pstar*self.gamma/self.dstarR)
            self.right_wave['tail_speed'] = self.uR+np.sqrt(
                self.pR*self.gamma/self.dR)
        return None

    def left_fan_state(self,x):
        return self.fan_state(x,-1)
    def right_fan_state(self,x):
        return self.fan_state(x,1)

    def fan_state(self,s,leftright):
        if leftright == -1:
            p0,d0,u0 = self.pL,self.dL,self.uL
        else:
            p0,d0,u0 = self.pR,self.dR,self.uR
        a0 = np.sqrt(self.gamma*p0/d0)
        p = p0*(2/(self.gamma+1)-leftright*(self.gamma-1)/(self.gamma+1)/a0*(u0-s)
                )**((2*self.gamma)/(self.gamma-1))
        u = u0+leftright*2*a0/(self.gamma-1)*((p/p0)*((self.gamma-1)/(2*self.gamma))-1)
        d = d0*(p/p0)**(1/self.gamma)
        return p,d,u
    

if __name__=='__main__':
    pass
