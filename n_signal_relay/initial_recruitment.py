#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


L_r_L = [1e6, 1e7, 1e8]
L_phi_E = [0., 1.]

#Required by Meerkat
L_variables = [L_r_L, L_phi_E]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E):
    return 'init_rec-r_L' + str(r_L) + 'phi_E' + str(phi_E)
def prototype():
    return d_runs[(1e6, 0.)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

def make_r_L_phi_E(r_L=None, phi_E=0., seed=0):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_N_props['n'] = 1
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_N_props['sigma_CL0'] = r_L / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    n_neutrophil = d_N_props.pop('n')
    d_N_props['x_max'] = d_PDE_props['x_0']
    x_max = d_N_props.pop('x_max')
    cell_group = springbok.CellGroup([n_exo.Neutrophil(
                xy_0=np.array([x_max, 500.]), index=(j + seed), **d_N_props)
                                      for j in range(n_neutrophil)])
    run = n_exo.new_setup(name=job_name(r_L, phi_E),
                    d_gen_props=d_gen_props,
                          L_cell_group=[cell_group], d_E_props=d_E_props, d_PDE_props=d_PDE_props, set_name='init_rec')
    return run


def setup(r_L, phi_E):
    jn = job_name(r_L, phi_E)
    run = make_r_L_phi_E(r_L=r_L, phi_E=phi_E)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
