#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


#L_r_L = [1e6]
L_r_L = [0., 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1, 1e0, 1.5e0, 2e0, 3e0, 5e0, 7.5e0, 1e1, 1.5e1, 2e1, 3e1, 5e1, 7.5e1, 1e2, 1.5e2, 2e2, 3e2, 5e2, 7.5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4]
#, 1e7, 2e7, 5e7, 1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10, 2e10, 5e10]
L_phi_E = [0., 0.25, 0.5, 0.75, 1.]

#Required by Meerkat
L_variables = [L_r_L, L_phi_E]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E):
    return 'n_exo_F-r_L' + str(r_L) + 'phi_E' + str(phi_E)
def prototype():
    return d_runs[(1e6, 0.)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

def setup(r_L, phi_E):
    jn = job_name(r_L, phi_E)
    run = n_exo.make_r_L_phi_E(r_L=r_L, phi_E=phi_E)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
