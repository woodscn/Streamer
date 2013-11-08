import numpy
import sympy

t=sympy.Symbol('t')
xi=sympy.Symbol('xi')
eta=sympy.Symbol('eta')
zeta=sympy.Symbol('zeta')

def shocked_sol_1():
    xmin, xmax, nx = 0, 1, 100
    dx = (xmax-xmin)/(nx)
    gamma = 1.4
    shock_pos_0 = .1
    pressure_lo, density_lo = 1, 1
    velocity_lo = 1.5*numpy.sqrt(1.4*pressure_lo/density_lo)
    a_lo, b_lo, c_lo, l_lo, m_lo, n_lo, p_lo, q_lo, r_lo = (1,0,0,0,1,0,0,0,1)
    u_lo, v_lo, w_lo = 0, 0, 0
    pressure_hi = 2#2.458333333333333*pressure_lo
    u_hi, v_hi, w_hi = 0, 0, 0
    density_hi = density_lo*(pressure_hi/pressure_lo*(gamma+1)+(gamma-1)
                             )/(pressure_hi/pressure_lo*(gamma-1)+(gamma+1))
    velocity_hi = velocity_lo-(pressure_hi/pressure_lo-1)*numpy.sqrt(
        pressure_lo*gamma/density_lo)/numpy.sqrt(
        .5*gamma*((gamma+1)*pressure_hi/pressure_lo+(gamma-1))) 
    def jacobian(a,b,c,l,m,n,p,q,r):
        return a*(m*r-n*q)+b*(n*p-l*r)+c*(l*q-m*p)
    def gradxi(a,b,c,l,m,n,p,q,r):
        return sympy.Matrix([(m*r-n*q)/jacobian(a,b,c,l,m,n,p,q,r),
                             (n*p-l*r)/jacobian(a,b,c,l,m,n,p,q,r),
                             (l*q-m*p)/jacobian(a,b,c,l,m,n,p,q,r)])
    def shock_speed(pressure,density,velocity,a,b,c,l,m,n,p,q,r,u,
                    gamma,pressure_hi,pm):
        return(sympy.sqrt(
                sum(gradxi(a,b,c,l,m,n,p,q,r).applyfunc(lambda x:x**2)))*(
                velocity-u+pm*sympy.sqrt(gamma*pressure/density)*sympy.sqrt(
                    (gamma+1)/(2*gamma)*(pressure_hi/pressure-1)+1)))
    shock_function_arg = ((xi-shock_pos_0)-t*shock_speed(
            pressure_lo,density_lo,velocity_lo,a_lo,b_lo,c_lo,l_lo,m_lo,n_lo,
            p_lo,q_lo,r_lo,u_lo,gamma,pressure_hi,pm=-1))
    shock_function = sympy.functions.Heaviside(shock_function_arg)
    out = []
    out.append(pressure_lo+(pressure_hi-pressure_lo)*shock_function)
    out.append(density_lo+(density_hi-density_lo)*shock_function)
    out.append(velocity_lo+(velocity_hi-velocity_lo)*shock_function)
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(1))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(sympy.Integer(0))
    out.append(xi)
    out.append(eta)
    out.append(zeta)
    out2=[out]
    out2.append([shock_function_arg])
    return out2

if __name__=="__main__":
    print shocked_sol_1()
