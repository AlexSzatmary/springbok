#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


#L_r_L = [1e6]
L_r_L = [1e0]
#, 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10]
L_phi_E = [0., 0.5, 1.]

#Required by Meerkat
L_variables = [L_r_L, L_phi_E]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E):
    return 'n_exo_validate-r_L' + str(r_L) + 'phi_E{:.2f}'.format(phi_E)
def prototype():
    return d_runs[(1e0, 0.)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat


def make_r_L_phi_E(r_L=None, phi_E=0.):
    if r_L is None:
        r_L = 1e6
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_PDE_props['ell'] = 1e20
    # d_PDE_props['x_r'] = 1e3
    # d_N_props['x_max'] = d_PDE_props['x_r']
    d_N_props['sigma_CL0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = n_exo.new_setup_random_N(name='validate_exo-r_L' + str(r_L) + 'phi_E' + str(phi_E),
                             d_gen_props=d_gen_props,
                             d_N_props=d_N_props, d_E_props=d_E_props,
                             d_PDE_props=d_PDE_props)
    run.d_gen_props = d_gen_props
    run.d_N_props = d_N_props
    run.d_E_props = d_E_props
    run.d_PDE_props = d_PDE_props
    return run


def setup(r_L, phi_E):
    jn = job_name(r_L, phi_E)
    run = make_r_L_phi_E(r_L=r_L, phi_E=phi_E)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
