import numpy as np
import scipy
from scipy.interpolate import interp1d
import sys
sys.path.insert(-1, 'springbok/tiger')
import tiger
import pdb


class Cell:
    def __init__(self, properties):
        pass

    def orient(self, L_condition):
        pass

    def move(self, clock, dt, L_condition):
        self.a_xy[clock] = (self.a_xy[clock - 1] +
                            self.velocity(L_condition) * dt)

    def behave(self, L_condition):
        pass


class CellGroup:
    def __init__(self, L_cell):
        self.L_cell = L_cell


class RectCellGroup(CellGroup):
    def __init__(self, CellType, xya, xyb, n_cell, seed=0, **kwargs):
        self.xya = xya
        self.xyb = xyb
        L_cell = [CellType(xy_0=xy_0, index=(j + seed), **kwargs)
                  for (j, xy_0) in
                  enumerate((xyb - xya) * np.random.random((n_cell, 2)) + xya)]
        super().__init__(L_cell)


class Exosome:
    def __init__(self, properties):
        pass

    def behave(self, L_condition):
        pass


class Molecules:
    def __init__(self):
        pass

    def step(self, sources):
        pass


class Springbok:
    def __init__(self, L_cell_group, pde_stepper, clock_start, clock_end):
        self.L_cell_group = L_cell_group
        self.pde_stepper = pde_stepper
        self.clock_start = clock_start
        self.clock_end = clock_end

    def run(self):
        for self.clock in range(self.clock_start, self.clock_end + 1):
            L_condition_interp = [
                (interp1d(pde.x, pde.u[self.clock - 1]),
                 interp1d(pde.x[1:-1], tiger.opdudx(pde.u[self.clock - 1], pde.dx)))
                for pde in self.pde_stepper.L_pde]
            L_a_secrete = [np.zeros(len(pde.x)) for pde in self.pde_stepper.L_pde]
            for g in self.L_cell_group:
                for c in g.L_cell:
                    xy = c.a_xy[self.clock - 1]
                    L_condition = [(f_c(xy[0]), f_dcdx(xy[0])) for (pde, (f_c, f_dcdx)) in zip(self.pde_stepper.L_pde, L_condition_interp)]
                    c.orient(L_condition)
                    c.move(self.clock, self.pde_stepper.dt, L_condition)
                    s = c.secrete(L_condition)
                    for (j, (a_secrete, pde)) in enumerate(
                            zip(L_a_secrete, self.pde_stepper.L_pde)):
                        i0 = int(xy[0] / pde.dx)
                        a_secrete[i0] += s[j] / pde.dx

            for pde in self.pde_stepper.L_pde:
                pde.f = pde.f_functional(*L_a_secrete)
            self.pde_stepper.step(self.clock)
