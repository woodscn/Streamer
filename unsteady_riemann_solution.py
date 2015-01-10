import scipy.optimize as opt
from numpy.testing import assert_almost_equal


class UnsteadyRiemannSolution(object):
    def __init__(self, pL, dL, uL, pR, dR, uR, gamma):
        self.pL, self.dL, self.uL = pL, dL, uL
        self.pR, self.dR, self.uR = pR, dR, uR
        self.gamma = gamma
        self.pstar, self.ustar, self.dstarL, self.dstarR = self.star_state()
        if self.pstar / self.pL <= 1.:
            self.left_wave = {
                'type': 'fan',
                'head_speed': self.left_fan_speed(
                    self.pL, self.dL, self.uL),
                'tail_speed': self.left_fan_speed(
                    self.pstar, self.dstarL, self.ustar)}
        else:
            self.left_wave = {
                'type': 'shock',
                'shock_speed': self.left_shock_speed()}
        if self.pstar / self.pR <= 1.:
            self.right_wave = {
                'type': 'fan',
                'head_speed': self.right_fan_speed(
                    self.pR, self.dR, self.uR),
                'tail_speed': self.right_fan_speed(
                    self.pstar, self.dstarR, self.ustar)}
        else:
            self.right_wave = {
                'type': 'shock',
                'shock_speed': self.right_shock_speed()}

    def __call__(self, x, t):
        speed = x / t
        p, d, u = self.sample(speed)
        return p, d, u

    def sample(self, speed):
        if speed < self.ustar: # Left of slip line
            if self.left_wave['type'] == 'fan': # Left region
                if speed <= self.left_wave['tail_speed']:
                    out = self.pL, self.dL, self.uL
                else:
                    if speed >= self.left_wave['head_angle']:# Left star region
                        out = self.pstar, self.dstarL, self.ustar
                    else: # Inside left fan
                        p, d, u = self.left_fan_state(
                            speed, self.pL, self.dL, self.uL)
                        out = p, d, u
            else: # Left shock
                if speed <= self.left_wave['shock_speed']:
                    out = self.pL, self.dL, self.uL
                else:
                    out = self.pstar, self.dstarL, self.ustar
        else: # Right of slip line
            if self.right_wave['type'] == 'fan':
                if speed >= self.right_wave['tail_angle']:
                    out = self.pR, self.dR, self.uR
                else:
                    if speed <= self.right_wave['head_speed']:
                        out = self.pstar, self.dstarR, self.ustar
                    else:
                        p, d, u = self.right_fan_state(
                            speed, self.pR, self.dR, self.uR)
                        out = p, d, u
            else: # Right shock
                if speed >= self.right_wave['shock_speed']:
                    out = self.pR, self.dR, self.uR
                else:
                    out = self.pstar, self.dstarR, self.ustar
        return out

    def star_state(self):
        pstar = self.pressure_solve()
        ustar = self.uR + self.delta_u(pstar, self.pR, self.dR)
        dstarL = self.dL * self.d_ratio(pstar / self.pL)
        dstarR = self.dR * self.d_ratio(pstar / self.pR)
        return pstar, ustar, dstarL, dstarR

    def fan_state(self, speed, p0, u0, d0, leftright_in):
        leftright = {'left': -1, 'right': 1}[leftright_in]
        a0 = self.sound_speed(self, p0, u0)
        p = self.fan_pressure(speed, p0, u0, a0, leftright)
        u = u0 + leftright * self.delta_u(p, p0, d0)
        d = d0 * self.d_ratio(p / p0)
        out = p, d, u
        return out

    def fan_pressure(self, speed, p0, u0, a0, leftright):
        exponent = 2 * self.gamma / (self.gamma - 1)
        coeff1 = 2 / (self.gamma + 1)
        coeff2 = (self.gamma - 1) / (a0 * (self.gamma + 1))
        out = p0 * (coeff1 - leftright * coeff2 * (u0 - speed)) ** exponent
        return out

    def pressure_func(self, p):
        ustarL = self.uL - self.delta_u(p, self.pL, self.dL)
        ustarR = self.uR + self.delta_u(p, self.pR, self.dR)
        return ustarR - ustarL

    def pressure_guess(self):
        return 0.5 * (self.pL + self.pR)

    def pressure_solve(self):
        guess = self.pressure_guess()
        return opt.newton(self.pressure_func, guess)

    def d_ratio(self, alpha):
        if alpha > 1:
            temp = ((self.gamma - 1) / (self.gamma + 1))
            out = (temp + alpha) / (temp * alpha + 1)
        else:
            out = alpha ** (1 / self.gamma)
        return out

    def left_shock_speed(self):
        return self.uL - self.shock_speed_term(self.pstar, self.pL, self.dL)

    def right_shock_speed(self):
        return self.uR + self.shock_speed_term(self.pstar, self.pR, self.dR)

    def shock_mass_flow(self, pstar, p0, d0):
        coeff1 = 2 / ((self.gamma + 1) * d0)
        coeff2 = p0 * (self.gamma - 1) / (self.gamma + 1)
        return ((pstar + coeff2) / coeff1) ** 0.5

    def shock_speed_term(self, pstar, p0, d0):
        return self.shock_mass_flow(pstar, p0, d0) / d0

    def sound_speed(self, p, d):
        return (self.gamma * p / d) ** 0.5

    def delta_u(self, pstar, p0, d0):
        if pstar / p0 > 1:
            out = (pstar - p0) / self.shock_mass_flow(pstar, p0, d0)
        else:
            coeff1 = 2 / (self.gamma - 1)
            exponent = ((self.gamma - 1) / (2 * self.gamma))
            out = (coeff1 * self.sound_speed(p0, d0) *
                   ((pstar / p0) ** exponent - 1))
        return out

    def left_fan_speed(self, p, d, u):
        return u - self.sound_speed(p, d)

    def right_fan_speed(self, p, d, u):
        return u + self.sound_speed(p, d)


def test_URS(inputs_dict, star_state, speeds, precision):
    URS = UnsteadyRiemannSolution(**inputs_dict)
    star_state = URS.pstar, URS.ustar, URS.dstarL, URS.dstarR
    for ind in range(4):
        assert_almost_equal(star_state[ind], star_state_test[ind],
                            decimal=precision)
    assert(speeds[0]['type'] == URS.left_wave['type'])
    # Check left wave speeds
    if speeds[0]['type'] == 'fan':
        assert_almost_equal(speeds[0]['head_speed'],
                            URS.left_wave['head_speed'], decimal=precision)
        assert_almost_equal(speeds[0]['tail_speed'],
                            URS.left_wave['tail_speed'], decimal=precision)
    else:
        assert_almost_equal(speeds[0]['shock_speed'],
                            URS.left_wave['shock_speed'], decimal=precision)
    # Check right wave speeds
    assert(speeds[0]['type'] == URS.left_wave['type'])
    if speeds[1]['type'] == 'fan':
        assert_almost_equal(speeds[1]['head_speed'],
                            URS.right_wave['head_speed'], decimal=precision)
        assert_almost_equal(speeds[1]['tail_speed'],
                            URS.right_wave['tail_speed'], decimal=precision)
    else:
        assert_almost_equal(speeds[1]['shock_speed'],
                            URS.right_wave['shock_speed'], decimal=precision)




if __name__ == "__main__":
    # Sod's problem with solution
    test_inputs_dict = {'pL': 1, 'dL': 1, 'uL': 0,
                        'pR': 0.1, 'dR': 0.125, 'uR': 0, 'gamma': 1.4}
    star_state_test = 0.30313, 0.92745, 0.42632, 0.26557
    speeds_test = ({'type': 'fan',
                   'head_speed': -1.18322,
                    'tail_speed': -0.0702745},
                   {'type': 'shock', 'shock_speed': 1.75216})
    test_URS(test_inputs_dict, star_state_test, speeds_test, 5)

    # Left shock, right rarefaction test
    test_inputs_dict = {'pL': 0.01, 'dL': 1, 'uL': 0,
                        'pR': 100, 'dR': 1, 'uR': 0, 'gamma': 1.4}
    star_state_test = 46.095, -6.19633, 5.99242, .57511
    speeds_test = ({'type': 'shock',
                    'shock_speed': -7.43747},
                   {'type': 'fan',
                    'tail_speed': 4.39658,
                    'head_speed': 11.8322})
    test_URS(test_inputs_dict, star_state_test, speeds_test, 3)

    print "All tests passed!"
