#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


L_r_L = [1e1, 1e2]
L_phi_E = [0., 1.]
L_has_BLT = [False, True]

#Required by Meerkat
L_variables = [L_r_L, L_phi_E, L_has_BLT]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E, has_BLT):
    return 'n_exo_linear-r_L{:.2f}phi_E{:.2f}'.format(r_L, phi_E) + 'BLT' + ('+' if has_BLT else '-')
def prototype():
    return d_runs[(1e6, 0., True)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

class BLT0(n_exo.Neutrophil):
    def kappa(self, L_condition):
        return super().kappa([L_condition[0], L_condition[1], (0., 0.)])

    def secrete(self, L_condition, clock):
        return super().secrete([L_condition[0], L_condition[1], (0., 0.)], clock)


def new_setup_random_N_BLT0(d_N_props=None, has_BLT=True, **kwargs):
    n_neutrophil = d_N_props.pop('n')
    x_max = d_N_props.pop('x_max')
    if has_BLT:
        CellType = n_exo.Neutrophil
    else:
        CellType = BLT0
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., -x_max / 2.]), np.array([x_max, x_max / 2.]),
        n_neutrophil, **d_N_props)
    return n_exo.new_setup(L_cell_group=[cg], **kwargs)


def make_r_L_phi_E(r_L=None, phi_E=0., has_BLT=True):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_PDE_props['ell'] = 1e3
    d_PDE_props['u_0'] = lambda d_PDE_props: lambda x: (x / d_PDE_props['ell'])
    d_N_props['sigma_CL0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = new_setup_random_N_BLT0(has_BLT=has_BLT,
                                  name='r_L' + str(r_L) + 'phi_E{:.2f}'.format(phi_E),
                                  d_gen_props=d_gen_props,
                                  d_N_props=d_N_props, d_E_props=d_E_props,
                                  d_PDE_props=d_PDE_props,
                                  set_name='n_exo_linear')
    return run


def setup(r_L, phi_E, has_BLT):
    jn = job_name(r_L, phi_E, has_BLT)
    run = make_r_L_phi_E(r_L=r_L, phi_E=phi_E, has_BLT=has_BLT)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
