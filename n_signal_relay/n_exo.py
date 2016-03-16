import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

set_name = 'n_exo'

def FRO(c):
    return c / (1. + c)

def DFRO(c, dcdx, ell):
    return ell * dcdx / (c + 1) ** 2

def get_d_N_props():
    return dict(
        length=10., persistence=1., speed=10.,
        sensitivity_F=200., sensitivity_L=200., F_xt=1e10,
        sigma_CL0=30., b_L=0., sigma_CE0=0., b_E=0.,
        K_d_F=1., K_d_L=1.)

def get_d_E_props():
    return dict(sigma_EL0=1., gamma_E=0.01)

def get_d_PDE_props():
    d_N_props = get_d_N_props()
    return dict(ell=1e3, Nt=61, dt=d_N_props['persistence'])

class Neutrophil(springbok.Cell):
    def __init__(
            self, length=None,
            persistence=None, speed=None,
            sensitivity_F=None, sensitivity_L=None, F_xt=None,
            sigma_CL0=None, b_L=None, sigma_CE0=None, b_E=None,
            n_t=None, xy_0=None, K_d_F=None, K_d_L=None, index=None):
        self.length = length
        self.sensitivity_F = sensitivity_F
        self.sensitivity_L = sensitivity_L
        self.persistence = persistence
        self.speed = speed
        self.a_theta = np.zeros(n_t)
        self.a_xy = np.zeros((n_t, 2))
        self.a_xy[0] = xy_0
        self.a_secrete = np.zeros((n_t, 3))
        self.K_d_F = K_d_F
        self.K_d_L = K_d_L
        self.n_t = n_t
        self.index = index
        self.random_state = np.random.RandomState(index)
        self.F_xt = F_xt
        self.sigma_CL0 = sigma_CL0
        self.b_L = b_L
        self.sigma_CE0 = sigma_CE0
        self.b_E = b_E

    def orient(self, L_condition, clock):
        kappa = self.kappa(L_condition)
        if kappa >= 0.:
            self.theta = self.random_state.vonmises(0., kappa)
        else:
            self.theta = self.random_state.vonmises(np.pi, -kappa)
        self.a_theta[clock] = self.theta

    def kappa(self, L_condition):
        F, dFdx = L_condition[0]
        exo, dexodx = L_condition[1]
        L, dLdx = L_condition[2]
        kappa = (
            self.sensitivity_F *
            DFRO(F / self.K_d_F, dFdx / self.K_d_F, self.length) +
            self.sensitivity_L *
            DFRO(L / self.K_d_L, dLdx / self.K_d_L, self.length) *
            np.exp(-F / self.F_xt))
        return kappa

    def velocity(self, conditions):
        return self.speed * np.array([np.cos(self.theta), np.sin(self.theta)])

    def secrete(self, L_condition, clock):
        F, dFdx = L_condition[0]
        exo, dexodx = L_condition[1]
        L, dLdx = L_condition[2]
        sigma_CE = self.sigma_CE0 * (self.b_E +
                                     (1. - self.b_E) * (F / (1. + F)))
        sigma_CL = self.sigma_CL0 * (self.b_L +
                                     (1. - self.b_L) * (F / (1. + F)))
        self.a_secrete[clock] = np.array([0., sigma_CE, sigma_CL])
        return self.a_secrete[clock]


class BLT0(Neutrophil):
    def orient(self, L_condition):
        super().orient([L_condition[0], (0., 0.), (0., 0.)])


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
    exo_pde.gamma = 10.
    exo_pde.sigma_EL = 10.
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
            exo_pde.sigma_EL * L_u[1][1:-1])

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, exo_pde, LTB_pde],
        dt=10., Nt=61)
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), has_BLT=True,
          has_crosstalk=False, name='n_BLT'):
    n_neutrophil = 200
    if has_BLT:
        if has_crosstalk:
            CellType = BLTposcrosstalk
        else:
            CellType = BLTpos
    else:
        CellType = BLT0
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., 0.]), np.array([2000., 1000.]),
        n_neutrophil, n_t=61)
    model = springbok.Springbok(
        L_cell_group=[cg], pde_stepper=pde_stepper,
        clock_start=1, clock_end=60)
    model.name = name
    model.set_name = set_name
    return model

def make_vary_ell(L_ell=[1e2, 1e3, 1e4], L_has_BLT=[False, True]):
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


def new_make_pde_stepper(ell, d_E_props, d_PDE_props):
    u_0 = lambda x: np.exp((x - 3600.) / ell)
    x_L = 0.
    x_r = 7e3
    n = 701
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
    exo_pde.gamma = d_E_props['gamma_E']
    exo_pde.sigma_EL = d_E_props['sigma_EL0']
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
            exo_pde.sigma_EL * L_u[1][1:-1])

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, exo_pde, LTB_pde],
        dt=d_PDE_props['dt'], Nt=61)
    return pde_stepper

def new_setup(
        pde_stepper=None,
        name='new_setup', d_N_props=None, d_PDE_props=None, d_E_props=None):
    if pde_stepper is None:
        pde_stepper=new_make_pde_stepper(1e3, d_E_props, d_PDE_props)
    n_neutrophil = 700
    CellType = Neutrophil
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., 0.]), np.array([7000., 1000.]),
        n_neutrophil, n_t=61, **d_N_props)
    model = springbok.Springbok(
        L_cell_group=[cg], pde_stepper=pde_stepper,
        clock_start=1, clock_end=60)
    model.name = name
    model.set_name = set_name
    return model


def make_vary_r_L(L_r_L=None, phi_E=0.):
    if L_r_L is None:
        L_r_L = [0., 1e2, 1e4, 1e6, 1e8]
    L_run = []
    for r_L in L_r_L:
        d_N_props = get_d_N_props()
        d_E_props = get_d_E_props()
        d_PDE_props = get_d_PDE_props()
        d_N_props['sigma_CL0'] = r_L / 100 * (1. - phi_E)
        d_N_props['sigma_CE0'] = r_L / 100 * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
        run = new_setup(name='r_L' + str(r_L) + 'phi_E' + str(phi_E),
                        d_N_props=d_N_props, d_E_props=d_E_props, d_PDE_props=d_PDE_props)
        L_run.append(run)
    return L_run
