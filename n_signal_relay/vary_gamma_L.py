#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


#L_r_L = [1e6]
#L_gamma_L = [400., 40., 4., 1/1.5, 0.4, 0.04, 0.004]
L_gamma_L = [5e2, 2e2, 1e2, 5e1, 2e1, 1e1, 5e0, 2e0, 1e0, 5e-1, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3]
#L_gamma_L = [5e-1, 2e-1, 1e-1]
L_r_L = [1e4, 1e5, 1e6, 1e7]
L_phi_E = [0., 0.5, 1.]

#Required by Meerkat
L_variables = [L_gamma_L, L_r_L, L_phi_E]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(gamma_L, r_L, phi_E):
    return 'vary_gamma_L' + str(gamma_L) + '_F-r_L' + str(r_L) + 'phi_E' + str(phi_E)
def prototype():
    return d_runs[(6e4, 1e6, 0.)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

def make_vary_gamma_L(gamma_L=None, r_L=None, phi_E=0., name=None):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_PDE_props['gamma_L'] = gamma_L
    d_N_props['sigma_CL0'] = r_L / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = n_exo.new_setup_random_N(
        name=name,
        d_gen_props=d_gen_props,
        d_N_props=d_N_props, d_E_props=d_E_props,
        d_PDE_props=d_PDE_props, set_name='vary_gamma_L')
    return run

def setup(gamma_L, r_L, phi_E):
    jn = job_name(gamma_L, r_L, phi_E)
    run = make_vary_gamma_L(gamma_L=gamma_L, r_L=r_L, phi_E=phi_E, name=jn)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
