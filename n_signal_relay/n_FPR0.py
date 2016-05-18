#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


L_r_L = [1e4]
#L_r_L = [0., 1e4, 1e6, 1e8, 1e10, 1e12]
#L_phi_E = [0., 0.25, 0.5, 0.75, 1.]
L_phi_E = [0.]
L_case = ['WT', 'FPR0', 'Mixed']

#Required by Meerkat
L_variables = [L_r_L, L_phi_E, L_case]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E, case):
    return 'n_FPR0-r_L' + str(r_L) + 'phi_E' + str(phi_E) + '-' + case
def prototype():
    return d_runs[(1e4, 0., 'Mixed')]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

class FPR0(n_exo.Neutrophil):
    def orient(self, L_condition, clock):
        super().orient([(0., 0.), L_condition[1], L_condition[2]], clock)

    def secrete(self, L_condition, clock):
        return super().secrete([(0., 0.), L_condition[1], L_condition[2]], clock)


def new_setup_random_N_FPR0(d_N_props=None, case=None, **kwargs):
    n_neutrophil = d_N_props.pop('n')
    x_max = d_N_props.pop('x_max')
    L_cell_group = []
    if case in ['WT', 'Mixed']:
        cg = springbok.RectCellGroup(
            n_exo.Neutrophil,
            np.array([0., -x_max / 2.]), np.array([x_max, x_max / 2.]),
            n_neutrophil, **d_N_props)
        L_cell_group.append(cg)
    if case in ['FPR0', 'Mixed']:
        cg = springbok.RectCellGroup(
            FPR0,
            np.array([0., -x_max / 2.]), np.array([x_max, x_max / 2.]),
            n_neutrophil, seed=n_neutrophil, **d_N_props)
        L_cell_group.append(cg)
    return n_exo.new_setup(L_cell_group=L_cell_group, **kwargs)


def make_r_L_phi_E(r_L=None, case=None, phi_E=0.):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_N_props['n'] = 100
    d_N_props['x_max'] = 4000
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
#    d_PDE_props['x_0'] = d_N_props['x_max']
    d_N_props['sigma_CL0'] = r_L / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = new_setup_random_N_FPR0(name='r_L' + str(r_L) + 'phi_E' + str(phi_E),
                                  case=case,
                                  d_gen_props=d_gen_props,
                                  d_N_props=d_N_props, d_E_props=d_E_props,
                                  d_PDE_props=d_PDE_props, set_name='n_FPR0')
    return run


def setup(r_L, phi_E, case):
    jn = job_name(r_L, phi_E, case)
    run = make_r_L_phi_E(r_L=r_L, phi_E=phi_E, case=case)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
