#meerkat-copy springbok runner.py n_signal_relay/n_exo.py
import sys
import springbok
import n_exo
import cloudpickle
import os
import numpy as np
import runner


L_r_L = [1e6]
L_phi_E = [0., 1.]

#Required by Meerkat
L_variables = [L_r_L, L_phi_E]

one_job = False

if one_job:
    L_variables = [[L[0]] for L in L_variables]

def job_name(r_L, phi_E):
    return 'decay_F-r_L' + str(r_L) + 'phi_E' + str(phi_E)
def prototype():
    return d_runs[(1e6, 0.)]
default_suffix = '.run.pkl'
def run():
    runner.setup(L_job_name)
# setup is also required
#End required by Meerkat

def setup(r_L, phi_E):
    jn = job_name(r_L, phi_E)
    run = n_exo.make_decay_F(r_L=r_L, phi_E=phi_E)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
