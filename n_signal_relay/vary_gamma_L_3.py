#meerkat-copy springbok runner_vary_gamma_L.py n_signal_relay/n_exo.py
#meerkat-copy flux.py post_processing
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner_vary_gamma_L


#L_r_L = [1e6]
L_gamma_L = [5e2, 2e2, 1e2, 5e1, 2e1, 1e1, 5e0, 2e0, 1e0, 5e-1, 4 / 15, 2e-1, 1e-1, 5e-2, 2e-2, 1e-2, 5e-3, 2e-3, 1e-3]
gamma_L_0 = 4 / 15
#L_gamma_L = [5e-1, 2e-1, 1e-1]
#L_i = [0, 1, 2]
L_i = [0, 1]
kstart = 3
L_k = [kstart + i for i in L_i]
L_phi_E = [0., 0.25, 0.5, 0.75, 1.]

# rough values taken from n_exo-2016-05-27
# d_r_L = {0.: 50., 0.25: 75., 0.5: 100., 0.75: 100., 1.: 200.}
# better values taken from n_exo-2016-08-23
d_r_L = {0.: 20., 0.25: 30., 0.5: 30., 0.75: 50.0, 1.0: 75.0}
L_q = [0.01, 0.1, 1., 10., 100]
d_L_r_L = {}
for gamma_L, phi_E in zip(L_gamma_L, L_phi_E):
    d_L_r_L[gamma_L, phi_E] = [q * d_r_L[phi_E] * gamma_L / gamma_L_0
                               for q in L_q]
#Required by Meerkat
L_variables = [L_gamma_L, L_phi_E, L_k]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(gamma_L, phi_E, k):
    return 'vary_gamma_L_3-{:.4f}-phi_E{:.2f}-k{:03d}'.format(gamma_L, phi_E, k)
def prototype():
    return d_runs[(4 / 15, 1e0, 0)]
default_suffix = '.run.pkl'
def run():
    runner_vary_gamma_L.setup(L_job_name)
# setup is also required
#End required by Meerkat

def make_vary_gamma_L(gamma_L=None, r_L=None, phi_E=0., name=None):
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_PDE_props['gamma_L'] = gamma_L
    d_N_props['sigma_CL0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L * d_PDE_props['x_r'] / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = n_exo.new_setup_random_N(
        name=name,
        d_gen_props=d_gen_props,
        d_N_props=d_N_props, d_E_props=d_E_props,
        d_PDE_props=d_PDE_props, set_name='vary_gamma_L')
    run.d_gen_props = d_gen_props
    run.d_N_props = d_N_props
    run.d_E_props = d_E_props
    run.d_PDE_props = d_PDE_props
    return run

def setup(gamma_L, phi_E, k):
    jn = job_name(gamma_L, phi_E, k)
    r_L = d_L_r_L[(gamma_L, phi_E)][k - kstart]
    run = make_vary_gamma_L(gamma_L=gamma_L, r_L=r_L, phi_E=phi_E, name=jn)
    run.job_name = jn
    run.gamma_L = gamma_L
    run.phi_E = phi_E
    run.k = k
    run.r_L = r_L

    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
