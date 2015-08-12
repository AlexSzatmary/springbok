import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

set_name = 'neutrophil_static_gradient'

def DFRO(c, dcdx, ell):
    return ell * dcdx / (c + 1) ** 2


class Neutrophil(springbok.Cell):
    def __init__(
            self, length=10., sensitivity=200., persistence=10., speed=1.,
            n_t=None, xy_0=None, K_d=1., index=None):
        self.length = length
        self.sensitivity = sensitivity
        self.persistence = persistence
        self.speed = speed
        self.a_theta = np.zeros(n_t)
        self.a_xy = np.zeros((n_t, 2))
        self.a_xy[0] = xy_0
        self.K_d = K_d
        self.n_t = n_t
        self.random_state = np.random.RandomState(index)

    def orient(self, L_conditions):
        c, dcdx = L_conditions[0]
        c0 = c / self.K_d
        dcdx0 = (self.length / self.K_d) * dcdx
        kappa = self.sensitivity * DFRO(c, dcdx, self.length)
        if kappa >= 0.:
            self.theta = self.random_state.vonmises(0., kappa)
        else:
            self.theta = self.random_state.vonmises(np.pi, -kappa)

    def velocity(self, conditions):
        return self.speed * np.array([np.cos(self.theta), np.sin(self.theta)])

    def secrete(self, L_condition):
        return np.array([0.])

def make_pde_stepper(ell):
    u_0 = lambda x: np.exp((x - 2000.) / ell)
    x_L = 0.
    x_r = 3e3
    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[tiger.PDE(f=None, u_0=u_0,
                         x_L=x_L, x_r=x_r, n=2001,
                         u_L=tiger.Dirichlet(u_0(x_L)),
                         u_r=tiger.Dirichlet(u_0(x_r)))],
    dt=10., Nt=61)
    pde_stepper.L_pde[0].f_functional = (
        lambda foo: lambda t, x, L_u, dudx, d2udx2: 
        np.zeros(np.shape(pde_stepper.L_pde[0].x[1:-1])))
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), name='nsg'):
    n_neutrophil = 100
    L_cell = [Neutrophil(xy_0=xy_0, n_t=61, index=j)
              for (j, xy_0) in 
              enumerate(1000. * np.random.random((n_neutrophil, 2)) +
              np.array([1000., 0.]))]
    
    model = springbok.Springbok(
        L_cell_group=[springbok.CellGroup(L_cell)], pde_stepper=pde_stepper,
        clock_start=1, clock_end=61)
    model.name = name
    model.set_name = set_name
    return model

def make_vary_ell(L_ell=[1e2, 3e2, 1e3, 3e3, 1e4]):
    return [setup(pde_stepper=make_pde_stepper(ell), name='nsg-{:.2e}'.format(ell)) for ell in L_ell]
