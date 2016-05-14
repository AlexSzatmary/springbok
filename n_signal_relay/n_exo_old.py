def make_pde_stepper(ell):
    u_0 = lambda x: np.exp((x - 2000.) / ell)
    x_L = 0.
    x_r = 3e3
    n = 301
    F_pde = tiger.PDE(f=None, u_0=u_0,
                      x_L=x_L, x_r=x_r, n=n,
                      u_L=tiger.Dirichlet(u_0(x_L)),
                      u_r=tiger.Dirichlet(u_0(x_r)))
    exo_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n,
                        u_L=tiger.Dirichlet(0.), u_r=tiger.Dirichlet(0.),
                        )
    LTB_pde = tiger.PDE(f=None,
                        u_0=0.,
                        x_L=x_L, x_r=x_r, n=n,
                        u_L=tiger.Dirichlet(0.), u_r=tiger.Dirichlet(0.),
                        )

    F_pde.f_functional = (lambda *foo: lambda t, x, L_u, dudx, d2udx2: 
                          np.zeros(np.shape(F_pde.x[1:-1])))
    exo_pde.gamma = 10.
    exo_pde.sigma_EL = 10.
    exo_pde.f_functional = (
            lambda a_F, a_exo, a_L:
                lambda t, x, L_u, dudx, d2udx2: 
            a_exo[1:-1] - exo_pde.gamma * L_u[1][1:-1])

    LTB_pde.DL = 6e4
    LTB_pde.gamma = 0.6
    LTB_pde.f_functional = (
            lambda a_F, a_exo, a_L:
                lambda t, x, L_u, dudx, d2udx2: 
            LTB_pde.DL * d2udx2 + a_L[1:-1] - LTB_pde.gamma * L_u[2][1:-1] +
            exo_pde.sigma_EL * L_u[1][1:-1])

    pde_stepper = tiger.CoupledPDEStepper2(
        L_pde=[F_pde, exo_pde, LTB_pde],
        dt=10., Nt=61)
    return pde_stepper

def setup(pde_stepper=make_pde_stepper(1000.), has_BLT=True,
          has_crosstalk=False, name='n_BLT'):
    n_neutrophil = 200
    if has_BLT:
        if has_crosstalk:
            CellType = BLTposcrosstalk
        else:
            CellType = BLTpos
    else:
        CellType = BLT0
    cg = springbok.RectCellGroup(
        CellType,
        np.array([0., 0.]), np.array([2000., 1000.]),
        n_neutrophil, n_t=61)
    model = springbok.Springbok(
        L_cell_group=[cg], pde_stepper=pde_stepper,
        clock_start=1, clock_end=60)
    model.name = name
    model.set_name = SET_NAME
    return model

def make_vary_ell(L_ell=[1e2, 1e3, 1e4], L_has_BLT=[False, True]):
    L_runs = []
    for ell in L_ell:
        for has_BLT in L_has_BLT:
            run = setup(pde_stepper=make_pde_stepper(ell), has_BLT=has_BLT,
                        name=('n_BLT-' + ('BLT' if has_BLT else 'BLT0') +
                              '-{:.2e}').format(ell))
            run.ell = ell
            run.has_BLT = has_BLT
            L_runs.append(run)
    return L_runs

def make_3_runs():
    ell = 1e3
    L_has_BLT = [False, True, True]
    L_has_crosstalk = [False, False, True]
    L_name = ['No exo', 'No crosstalk', 'Crosstalk']
    L_runs = []
    for (has_BLT, has_crosstalk, name) in zip(
            L_has_BLT, L_has_crosstalk, L_name):
        run = setup(pde_stepper=make_pde_stepper(ell), has_BLT=has_BLT,
                    has_crosstalk=has_crosstalk,
                    name=name)
        run.ell = ell
        run.has_BLT = has_BLT
        run.has_crosstalk = has_crosstalk
        L_runs.append(run)
    return L_runs
