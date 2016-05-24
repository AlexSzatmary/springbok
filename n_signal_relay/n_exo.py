import numpy as np
from scipy.stats import vonmises
import sys
#sys.path.insert(-1, '../springbok')
import springbok
from springbok import tiger

SET_NAME = 'n_exo'

def FRO(c):
    return c / (1. + c)

def DFRO(c, dcdx, ell):
    return ell * dcdx / (c + 1) ** 2

def get_d_gen_props():
    return dict(Nt=61)

def get_d_N_props():
    return dict(
        length=10., persistence=1., speed=10.,
        sensitivity_F=200., sensitivity_L=200., F_xt=1e0,
        sigma_CL0=30., b_L=0., sigma_CE0=0., b_E=0.,
        K_d_F=1., K_d_L=1., n=500, x_max=1e4, n_t=get_d_gen_props()['Nt'])

def get_d_E_props():
    return dict(sigma_EL0=1., gamma_E=0.01)

def get_d_PDE_props(d_N_props, d_gen_props):
    return dict(ell=4e2, x_0=3600.,
                Nt=d_gen_props['Nt'], dt=d_N_props['persistence'],
                x_r=1e4, n=1001, gamma_F=0., DL=6e4, gamma_L=1/1.5,
                u_0=exp_u_0)


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


def exp_u_0(d_PDE_props):
    return lambda x: np.exp((x - d_PDE_props['x_0']) / d_PDE_props['ell']) 


def exp_bell_u_0(d_PDE_props):
    ell = d_PDE_props['ell']
    x_0 = d_PDE_props['x_0']
    a = d_PDE_props['a']
    ca = d_PDE_props['ca']
    def f(x):
        if type(x) is np.ndarray:
            a_x = x
        else:
            a_x = [x]
        L_f = []
        for xi in a_x:
            if xi - x_0 < -a:
                L_f.append(ca * np.exp((xi - (x_0 - a)) / ell))
            elif -a <= xi - x_0 < a:
                L_f.append(ca * (2. - (np.exp((xi - (x_0 - a)) / ell) + np.exp(-(xi - (x_0 + a)) / ell)) /
                             (1. + np.exp(2. * a / ell))))
            else:
                L_f.append(ca * np.exp(-(xi - (x_0 + a)) / ell))
        if type(x) is np.ndarray:
            return np.array(L_f)
        else:
            return L_f[0]
    return f

def new_make_pde_stepper(d_E_props, d_PDE_props):
    u_0 = d_PDE_props['u_0'](d_PDE_props)
    x_L = 0.
    x_r = d_PDE_props['x_r']
    n = d_PDE_props['n']
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
                        u_L=tiger.VonNeumann(0.), u_r=tiger.VonNeumann(0.),
                        )

    F_pde.gamma_F = d_PDE_props['gamma_F']
    F_pde.f_functional = (lambda *foo: lambda t, x, L_u, dudx, d2udx2: 
                          -F_pde.gamma_F * L_u[0][1:-1])
    exo_pde.gamma = d_E_props['gamma_E']
    exo_pde.sigma_EL = d_E_props['sigma_EL0']
    exo_pde.f_functional = (
            lambda a_F, a_exo, a_L:
                lambda t, x, L_u, dudx, d2udx2: 
            a_exo[1:-1] - exo_pde.gamma * L_u[1][1:-1])

    LTB_pde.DL = d_PDE_props['DL']
    LTB_pde.gamma = d_PDE_props['gamma_L']
    def outer(a_F, a_exo, a_L):
        def inner(t, x, L_u, dudx, d2udx2):
            E = L_u[1][1:-1]
            L = L_u[2][1:-1]
            return (LTB_pde.DL * d2udx2 + a_L[1:-1] - LTB_pde.gamma * L +
                    exo_pde.sigma_EL * E)
        return inner

    LTB_pde.f_functional = outer
    def bc_f(t, x, L_u, dudx, d2udx2):
        return (LTB_pde.DL * d2udx2 - LTB_pde.gamma * L_u[2] +
                exo_pde.sigma_EL * L_u[1])
    LTB_pde.u_L.f = bc_f
    LTB_pde.u_r.f = bc_f

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, exo_pde, LTB_pde],
        dt=d_PDE_props['dt'], Nt=d_PDE_props['Nt'])
    return pde_stepper

def new_setup(
        pde_stepper=None,
        name='new_setup', d_gen_props=None,
        L_cell_group=None, d_PDE_props=None, d_E_props=None, set_name=None):
    if pde_stepper is None:
        pde_stepper = new_make_pde_stepper(d_E_props, d_PDE_props)
    model = springbok.Springbok(
        L_cell_group=L_cell_group, pde_stepper=pde_stepper,
        clock_start=1, clock_end=d_gen_props['Nt'] - 1)
    model.name = name
    if set_name is None:
        model.set_name = SET_NAME
    else:
        model.set_name = set_name
    return model


def new_setup_random_N(d_N_props=None, **kwargs):
    CellType = Neutrophil
    n_neutrophil = d_N_props.pop('n')
    x_max = d_N_props.pop('x_max')
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., -x_max / 2.]), np.array([x_max, x_max / 2.]),
        n_neutrophil, **d_N_props)
    return new_setup(L_cell_group=[cg], **kwargs)


def make_r_L_phi_E(r_L=None, phi_E=0.):
    d_gen_props = get_d_gen_props()
    d_N_props = get_d_N_props()
    d_E_props = get_d_E_props()
    d_PDE_props = get_d_PDE_props(d_N_props, d_gen_props)
    d_N_props['sigma_CL0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = new_setup_random_N(name='r_L' + str(r_L) + 'phi_E' + str(phi_E),
                             d_gen_props=d_gen_props,
                             d_N_props=d_N_props, d_E_props=d_E_props,
                             d_PDE_props=d_PDE_props)
    run.d_gen_props = d_gen_props
    run.d_N_props = d_N_props
    run.d_E_props = d_E_props
    run.d_PDE_props = d_PDE_props
    return run


def make_vary_r_L(L_r_L=None, phi_E=0.):
    if L_r_L is None:
        L_r_L = [0., 1e2, 1e4, 1e6, 1e8]
    L_run = []
    for r_L in L_r_L:
        run = make_r_L_phi_E(r_L, phi_E)
        L_run.append(run)
    return L_run
