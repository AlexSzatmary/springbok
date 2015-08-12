import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

set_name = 'n_mixed'

def DFRO(c, dcdx, ell):
    return ell * dcdx / (c + 1) ** 2


class Neutrophil(springbok.Cell):
    def __init__(
            self, length=10., sensitivity=200., persistence=10., speed=1.,
            n_t=None, xy_0=None, K_d_F=1., K_d_L=1., index=None):
        self.length = length
        self.sensitivity = sensitivity
        self.persistence = persistence
        self.speed = speed
        self.a_theta = np.zeros(n_t)
        self.a_xy = np.zeros((n_t, 2))
        self.a_xy[0] = xy_0
        self.K_d_F = K_d_F
        self.K_d_L = K_d_L
        self.n_t = n_t
        self.index = index
        self.random_state = np.random.RandomState(index)

    def orient(self, L_condition):
        F, dFdx = L_condition[0]
        L, dLdx = L_condition[1]
        kappa = self.sensitivity * (
            DFRO(F / self.K_d_F, dFdx / self.K_d_F, self.length) +
            DFRO(L / self.K_d_L, dLdx / self.K_d_L, self.length))
        if kappa >= 0.:
            self.theta = self.random_state.vonmises(0., kappa)
        else:
            self.theta = self.random_state.vonmises(np.pi, -kappa)

    def velocity(self, conditions):
        return self.speed * np.array([np.cos(self.theta), np.sin(self.theta)])

    def secrete(self, L_condition):
        return np.array([0., 6.])


class FPRpos(Neutrophil):
    def orient(self, L_condition):
        super().orient([L_condition[0], (0., 0.)])


class FPR0(Neutrophil):
    def orient(self, L_condition):
        super().orient([(0., 0.), L_condition[1]])

    def secrete(self, L_condition):
        return np.array([0., 0.])


def make_pde_stepper(ell):
    u_0 = lambda x: np.exp((x - 2000.) / ell)
    x_L = 0.
    x_r = 3e3
    n = 301
    F_pde = tiger.PDE(f=None, u_0=u_0,
                      x_L=x_L, x_r=x_r, n=n,
                      u_L=tiger.Dirichlet(u_0(x_L)),
                      u_r=tiger.Dirichlet(u_0(x_r)))
    LTB_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n,
                        u_L=tiger.Dirichlet(0.), u_r=tiger.Dirichlet(0.))

    F_pde.f_functional = (lambda *foo: lambda t, x, L_u, dudx, d2udx2: 
                          np.zeros(np.shape(F_pde.x[1:-1])))
    LTB_pde.DL = 6e4
    LTB_pde.gamma = 0.6
    LTB_pde.f_functional = (
        lambda a_F, a_L:
        lambda t, x, L_u, dudx, d2udx2: 
        LTB_pde.DL * d2udx2 + a_L[1:-1] - LTB_pde.gamma * L_u[1][1:-1])

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, LTB_pde],
        dt=10., Nt=61)
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), name='n_mixed'):
    n_neutrophil = 100
    L_FPRpos = [FPRpos(xy_0=xy_0, n_t=61)
                for xy_0 in
                1000. * np.random.random((n_neutrophil // 2, 2)) +
                np.array([1000., 0.])]
    L_FPR0 = [FPR0(xy_0=xy_0, n_t=61)
              for xy_0 in
              1000. * np.random.random((n_neutrophil // 2, 2)) +
              np.array([1000., 0.])
              ]
    L_cell_group = [springbok.CellGroup(L_FPRpos),
                    springbok.CellGroup(L_FPR0)]
    
    model = springbok.Springbok(
        L_cell_group=L_cell_group, pde_stepper=pde_stepper,
        clock_start=1, clock_end=61)
    model.name = name
    model.set_name = set_name
    return model

def make_vary_ell(L_ell=[1e2, 3e2, 1e3, 3e3, 1e4]):
    return [setup(pde_stepper=make_pde_stepper(ell), name='n_mixed-{:.2e}'.format(ell)) for ell in L_ell]
