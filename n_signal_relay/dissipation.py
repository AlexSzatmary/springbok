# #meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner
import tiger


L_r_L = [1e0]
# L_r_L = [0., 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 1.5e0, 2e0, 3e0, 5e0, 7.5e0, 1e1, 1.5e1, 2e1, 3e1, 5e1, 7.5e1, 1e2, 1.5e2, 2e2, 3e2, 5e2, 7.5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4]
L_r_L = [r_L * 10. for r_L in [0., 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 1.5e0, 2e0, 3e0, 5e0, 7.5e0, 1e1, 1.5e1, 2e1, 3e1, 5e1, 7.5e1, 1e2, 1.5e2, 2e2, 3e2, 5e2, 7.5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4]]
# , 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10]
L_phi_E = [0., 0.25, 0.5, 0.75, 1.]
L_p = [0.5, 1., 2., 4.]

# Required by Meerkat
L_variables = [L_r_L, L_phi_E, L_p]

one_job = True

if one_job:
    L_variables = [[L[0]] for L in L_variables]


def job_name(r_L, phi_E, p):
    return 'n_exo-dissipation-r_L{:.2f}phi_E{:.2f}p{:.2f}'.format(
        r_L, phi_E, p)


def prototype():
    return d_runs[(1e6, 0.)]


default_suffix = '.run.pkl'


def run():
    runner.setup(L_job_name)
# setup is also required
# End required by Meerkat


def make_pde_stepper_dissipation(d_E_props, d_PDE_props):
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
    LTB_pde.p = d_PDE_props['p']

    def outer(a_F, a_exo, a_L):
        def inner(t, x, L_u, dudx, d2udx2):
            E = L_u[1][1:-1]
            L = L_u[2][1:-1]
            return (LTB_pde.DL * d2udx2 + a_L[1:-1]
                    - LTB_pde.gamma * L ** LTB_pde.p +
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


def make_r_L_phi_E(r_L=None, phi_E=0., p=1, seed=0):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_PDE_props['p'] = p
    d_N_props['sigma_CL0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (
        1. - phi_E)
    d_N_props['sigma_CE0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (
        d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E)
    pde_stepper = make_pde_stepper_dissipation(d_E_props, d_PDE_props)
    run = n_exo.new_setup_random_N(
        name='r_L' + str(r_L) + 'phi_E' + str(phi_E),
        pde_stepper=pde_stepper,
        d_gen_props=d_gen_props,
        d_N_props=d_N_props, d_E_props=d_E_props,
        d_PDE_props=d_PDE_props, seed=seed)
    run.d_gen_props = d_gen_props
    run.d_N_props = d_N_props
    run.d_E_props = d_E_props
    run.d_PDE_props = d_PDE_props
    return run


def setup(r_L, phi_E, p):
    jn = job_name(r_L, phi_E, p)
    run = make_r_L_phi_E(r_L=r_L, phi_E=phi_E, p=p)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
