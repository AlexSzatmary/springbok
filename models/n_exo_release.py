import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

set_name = 'n_exo_release'

def FRO(c):
    return c / (1. + c)

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
        exo, dexodx = L_condition[1]
        L, dLdx = L_condition[2]
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
        return np.array([0., 30., 0.]) * (0.05 + 0.95 * FRO(L_condition[0][0]))


class BLTpos(Neutrophil):
    pass


class BLTposcrosstalk(BLTpos):
    def orient(self, L_condition):
        F, dFdx = L_condition[0]
        exo, dexodx = L_condition[1]
        L, dLdx = L_condition[2]
        kappa = self.sensitivity * (
            DFRO(F / self.K_d_F, dFdx / self.K_d_F, self.length) +
            DFRO(L / self.K_d_L, dLdx / self.K_d_L, self.length) * np.exp(-10. * F))
        if kappa >= 0.:
            self.theta = self.random_state.vonmises(0., kappa)
        else:
            self.theta = self.random_state.vonmises(np.pi, -kappa)


class NRelease(BLTposcrosstalk):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        n_t = kwargs['n_t']
        self.a_reservoir = np.zeros(n_t)
        self.a_reservoir[0] = 1.

    def secrete(self, L_condition, clock=None):
        q_max = 30.
        q = q_max * (0.05 + 0.95 * FRO(L_condition[0][0])) * self.a_reservoir[clock - 1]
        self.a_reservoir[clock] = (self.a_reservoir[clock - 1] - q / q_max +
                                   0.1 * (1. - self.a_reservoir[clock - 1]))
        return np.array([0., q, 0.])


class BLT0(Neutrophil):
    def orient(self, L_condition):
        super().orient([L_condition[0], (0., 0.), (0., 0.)])

    # def secrete(self, L_condition):
    #     return np.array([0., 0.])


def make_pde_stepper(ell):
    u_0 = lambda x: np.exp((x - 2000.) / ell)
    x_L = 0.
    x_r = 3e3
    n = 301
    F_pde = tiger.PDE(f=None, u_0=u_0,
                      x_L=x_L, x_r=x_r, n=n,
                      u_L=tiger.Dirichlet(u_0(x_L)),
                      u_r=tiger.Dirichlet(u_0(x_r)))
    exo_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n,
                        u_L=tiger.Dirichlet(0.), u_r=tiger.Dirichlet(0.),
                        )
    LTB_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n,
                        u_L=tiger.Dirichlet(0.), u_r=tiger.Dirichlet(0.),
                        )

    F_pde.f_functional = (lambda *foo: lambda t, x, L_u, dudx, d2udx2: 
                          np.zeros(np.shape(F_pde.x[1:-1])))
    exo_pde.gamma = 1.
    exo_pde.s = 10.
    exo_pde.f_functional = (
            lambda a_F, a_exo, a_L:
                lambda t, x, L_u, dudx, d2udx2: 
            a_exo[1:-1] - exo_pde.gamma * L_u[1][1:-1])

    LTB_pde.DL = 6e4
    LTB_pde.gamma = 0.6
    LTB_pde.f_functional = (
            lambda a_F, a_exo, a_L:
                lambda t, x, L_u, dudx, d2udx2: 
            LTB_pde.DL * d2udx2 + a_L[1:-1] - LTB_pde.gamma * L_u[2][1:-1] +
            exo_pde.s * L_u[1][1:-1])

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, exo_pde, LTB_pde],
        dt=10., Nt=61)
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), has_BLT=True,
          has_crosstalk=True, name=set_name):
    n_neutrophil = 100
    if has_BLT:
        if has_crosstalk:
            CellType = NRelease
        else:
            CellType = BLTpos
    else:
        CellType = BLT0
    cg = springbok.RectCellGroup(
        CellType,
        np.array([1000., 0.]), np.array([1800., 1000.]),
        n_neutrophil, n_t=61)
    model = springbok.Springbok(
        L_cell_group=[cg], pde_stepper=pde_stepper,
        clock_start=1, clock_end=60)
    model.name = name
    model.set_name = set_name
    return model

def make_vary_ell(L_ell=[1e2, 1e3, 1e4], L_has_BLT=[True]):
    L_runs = []
    for ell in L_ell:
        for has_BLT in L_has_BLT:
            run = setup(pde_stepper=make_pde_stepper(ell), has_BLT=has_BLT,
                        name=('n_BLT-' + ('BLT' if has_BLT else 'BLT0') +
                              '-{:.2e}').format(ell))
            run.ell = ell
            run.has_BLT = has_BLT
            L_runs.append(run)
    return L_runs

def make_3_runs():
    ell = 1e3
    L_has_BLT = [False, True, True]
    L_has_crosstalk = [False, False, True]
    L_name = ['No exo', 'No crosstalk', 'Crosstalk']
    L_runs = []
    for (has_BLT, has_crosstalk, name) in zip(
            L_has_BLT, L_has_crosstalk, L_name):
        run = setup(pde_stepper=make_pde_stepper(ell), has_BLT=has_BLT,
                    has_crosstalk=has_crosstalk,
                    name=name)
        run.ell = ell
        run.has_BLT = has_BLT
        run.has_crosstalk = has_crosstalk
        L_runs.append(run)
    return L_runs

