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

class SpringbokPreRun(springbok.Springbok):
    def __init__(self, L_cell_group, pde_stepper, clock_start, clock_end,
                 prerun):
        self.prerun = prerun
        for cg in prerun.L_cell_group: # Keep cells in prerun immobile
            for cell in cg.L_cell:
                cell.speed = 0.
        super().__init__(L_cell_group, pde_stepper, clock_start, clock_end)

    def run(self):
        self.prerun.run()

        for pde in self.pde_stepper.L_pde:
            pde.u = np.zeros((self.pde_stepper.Nt, np.size(pde.u_0)))
        self.pde_stepper.L_pde[0].u_0 = np.zeros(np.shape(
                self.prerun.pde_stepper.L_pde[0].u_0))
        self.pde_stepper.L_pde[1].u_0 = self.prerun.pde_stepper.L_pde[1].u[-1]
        self.pde_stepper.L_pde[2].u_0 = self.prerun.pde_stepper.L_pde[2].u[-1]
        self.pde_stepper.u = np.zeros((self.pde_stepper.Nt, np.size(self.pde_stepper.u_0)))
        self.pde_stepper.u_0 = np.hstack([pde.u_0 for pde in self.pde_stepper.L_pde])
        self.pde_stepper.u[0] = self.pde_stepper.u_0

        super().run()


def new_setup(
        pde_stepper=None, pde_stepper_prerun=None,
        name='new_setup', d_gen_props=None, d_N_props=None,
        L_cell_group=None, d_PDE_props=None, d_E_props=None, set_name=None):
    if pde_stepper is None:
        pde_stepper = n_exo.new_make_pde_stepper(d_E_props, d_PDE_props)
    if pde_stepper_prerun is None:
        pde_stepper_prerun = n_exo.new_make_pde_stepper(d_E_props, d_PDE_props)
    cg_prerun = make_random_N(d_N_props=d_N_props)
    prerun = springbok.Springbok(
        L_cell_group=[cg_prerun], pde_stepper=pde_stepper_prerun,
        clock_start=1, clock_end=d_gen_props['Nt'] - 1)
    model = SpringbokPreRun(
        L_cell_group=L_cell_group, pde_stepper=pde_stepper,
        clock_start=1, clock_end=d_gen_props['Nt'] - 1,
        prerun=prerun)
    model.name = name
    if set_name is None:
        model.set_name = SET_NAME
    else:
        model.set_name = set_name
    return model


def make_random_N(d_N_props=None, **kwargs):
    d_N_props = d_N_props.copy()
    CellType = n_exo.Neutrophil
    n_neutrophil = d_N_props.pop('n')
    x_max = d_N_props.pop('x_max')
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., -x_max / 2.]), np.array([x_max, x_max / 2.]),
        n_neutrophil, **d_N_props)
    return cg


def new_setup_random_N(d_N_props=None, **kwargs):
    cg = make_random_N(d_N_props=d_N_props)
    return new_setup(L_cell_group=[cg], d_N_props=d_N_props, **kwargs)


def make_decay_F(r_L=None, phi_E=0., gamma_F=0.):
    if r_L is None:
        r_L = 1e6
    d_gen_props = n_exo.get_d_gen_props()
    d_N_props = n_exo.get_d_N_props()
    d_E_props = n_exo.get_d_E_props()
    d_PDE_props = n_exo.get_d_PDE_props(d_N_props, d_gen_props)
    d_N_props['x_max'] = d_PDE_props['x_r']
    d_PDE_props['gamma_F'] = gamma_F
    d_N_props['sigma_CL0'] = r_L / d_N_props['n'] * (1. - phi_E)
    d_N_props['sigma_CE0'] = r_L / d_N_props['n'] * d_E_props['gamma_E'] / d_E_props['sigma_EL0'] * phi_E
    run = new_setup_random_N(name='decay_F-r_L' + str(r_L) + 'phi_E' + str(phi_E),
                    d_gen_props=d_gen_props,
                    d_N_props=d_N_props, d_E_props=d_E_props, d_PDE_props=d_PDE_props, set_name='n_exo_decay')
    return run


def setup(r_L, phi_E):
    jn = job_name(r_L, phi_E)
    run = make_decay_F(r_L=r_L, phi_E=phi_E, gamma_F=0.)
    run.job_name = jn
    with open(jn + '.pkl', 'wb') as hout:
        cloudpickle.dump(run, hout)
