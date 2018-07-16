# Probably obsolete
import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

set_name = 'validate_0'

def DFRO(c, dcdx, ell):
    return ell * dcdx / (c + 1) ** 2


class Neutrophil(springbok.Cell):
    def __init__(
            self, length=10., sensitivity=200., persistence=1., speed=0.,
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
        return np.array([0., 1.])

def make_pde_stepper(n_mesh):
    u_0 = lambda x: np.exp((x - 2000.) / 1e3)
    x_L = 0.
    x_r = 2e3

    F_pde = tiger.PDE(f=None, u_0=u_0,
                      x_L=x_L, x_r=x_r, n=n_mesh,
                      u_L=tiger.Dirichlet(u_0(x_L)),
                      u_r=tiger.Dirichlet(u_0(x_r)))
    LTB_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n_mesh,
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
        dt=1., Nt=61)
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), name='validate_0'):
    n_neutrophil = 1
    L_cell = [Neutrophil(xy_0=np.array([1000., 0.]), n_t=61, index=j)
              for j in range(n_neutrophil)]
    
    model = springbok.Springbok(
        L_cell_group=[springbok.CellGroup(L_cell)], pde_stepper=pde_stepper,
        clock_start=1, clock_end=60)
    model.name = name
    model.set_name = set_name
    return model

def make_vary_n_mesh(L_n_mesh=[126, 251, 501, 1001, 2001]):
    return [setup(pde_stepper=make_pde_stepper(n_mesh), name='validate_0-{:05}'.format(n_mesh)) for n_mesh in L_n_mesh]

def calculate_error(r):
    DL = r.pde_stepper.L_pde[1].DL
    gamma = r.pde_stepper.L_pde[1].gamma
    x = r.pde_stepper.L_pde[1].x
    ell = np.sqrt(DL / gamma)
    c_0 = 1. / (2. * np.sqrt(DL * gamma))
    c_analytic = lambda x: c_0 * np.exp(-x/ell) if x >= 0. else c_0 * np.exp(x/ell)
    u_analytic = np.array(list(map(c_analytic, x - 1e3)))
    u_end = r.pde_stepper.L_pde[1].u[-1]
    n = r.pde_stepper.L_pde[1].n
    return (np.max(np.abs(u_end - u_analytic)),
            u_end[(n - 1)/2] - u_analytic[(n - 1)/2])
