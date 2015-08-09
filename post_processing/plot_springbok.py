import numpy as np
from . import platypus
import os


def plot_tracks(model, file_name=True, path=None, style=None):
    if file_name is True:
        file_name = 'track-' + model.name
    xmin, xmax = (np.min(model.pde_stepper.L_pde[0].x), np.max(model.pde_stepper.L_pde[0].x))
    fig = platypus.multi_plot(
        [n.a_xy[:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[:, 1] for n in model.L_cell_group[0].L_cell],
        xlim=(xmin, xmax), ylim=(xmin, xmax),
        file_name=None)
    figx = platypus.multi_plot(
        [n.a_xy[-1:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[-1:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        L_marker=['o'] * len(model.L_cell_group[0].L_cell),
        markerfacecolor='none', markersize=(1e0),
        xlim=(1e3, 2e3), ylim=(0e3, 1e3),
        file_name=None)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def do_all_post_processing(L_model):
    list(map(get_a_CI, L_model))

def get_a_CI(model):
    model.a_CI = [(c.a_xy[-1, 0] - c.a_xy[0, 0]) /
                  (c.speed * c.persistence * (c.n_t - 1))
                  for c in model.L_cell_group[0].L_cell]

def plot_CI(L_model, file_name='CI', path=None, style=None):
    fig = style()
    platypus.boxplot([model.a_CI for model in L_model],
                labels=[model.name for model in L_model], fig=fig)
    fig.savefig(file_name, path=path)
    return fig

def make_all_plots(L_model, path=None, style=platypus.Print):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, style.style)
    L_fig = [plot_tracks(model, path=path, style=style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, style=style)
    L_fig.append(fig)
    return L_fig
