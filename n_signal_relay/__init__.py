from . import validate_0
import numpy as np

study_name = 'n_signal_relay'

def run_validate_0():
    L_run = validate_0.make_vary_n_mesh()
    a = np.zeros((len(L_run), 2))
    for r in L_run:
        r.run()
    a[:, 0] = np.array([r.pde_stepper.L_pde[0].dx for r in L_run])
    a[:, 1] = np.array([validate_0.calculate_error(r) for r in L_run])[:, 0]
    if not os.path.exists(os.path.join(study_name, validate_0.set_name)):
        os.mkdir(os.path.join(study_name, validate_0.set_name))
    np.savetxt(os.path.join(study_name, validate_0.set_name, validation.txt),
               a, header='dx\terror')                             
    return L_run
