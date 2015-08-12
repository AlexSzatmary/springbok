import numpy as np
from . import platypus
import os


def plot_tracks(model, file_name=True, path=None, Style=None, fig=None):
    if fig is None:
        fig = Style()
    if file_name is True:
        file_name = 'track-' + model.name
    xmin, xmax = (np.min(model.pde_stepper.L_pde[0].x), np.max(model.pde_stepper.L_pde[0].x))
    fig = platypus.multi_plot(
        [n.a_xy[:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        xlim=(xmin, xmax), ylim=(xmin, xmax),
        file_name=None)
    figx = platypus.multi_plot(
        [n.a_xy[-1:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[-1:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        L_marker=['o'] * len(model.L_cell_group[0].L_cell),
        markerfacecolor='none', markersize=(4e0),
        xlim=(1e3, 2e3), ylim=(0e3, 1e3),
        file_name=None)
    if len(model.L_cell_group) == 2:
        figx = platypus.multi_plot(
            [n.a_xy[:, 0] for n in model.L_cell_group[1].L_cell],
            [n.a_xy[:, 1] for n in model.L_cell_group[1].L_cell],
            fig=fig,
            xlim=(xmin, xmax), ylim=(xmin, xmax),
            color_f=lambda x: platypus.BLACK,
            file_name=None)
        figx = platypus.multi_plot(
            [n.a_xy[-1:, 0] for n in model.L_cell_group[1].L_cell],
            [n.a_xy[-1:, 1] for n in model.L_cell_group[1].L_cell],
            fig=fig,
            L_marker=['o'] * len(model.L_cell_group[1].L_cell),
            markerfacecolor='w', markersize=(4e0),
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

def plot_CI(L_model, file_name='CI', path=None, Style=None, fig=None):
    if fig is None:
        fig = Style()
    platypus.boxplot([model.a_CI for model in L_model],
                labels=[model.name for model in L_model], fig=fig)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def make_all_plots(L_model, path=None, Style=platypus.Print):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, Style.style)
    L_fig = [plot_tracks(model, path=path, Style=Style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, Style=Style)
    L_fig.append(fig)
    return L_fig

def make_n_static_fig(L_model):
    n_row, n_col = (3, 2)
    fig = platypus.Print(subplot=(n_row, n_col, 1))
    path = os.path.join('plots', L_model[0].set_name, fig.style)
    for (i, model) in enumerate(L_model):
        if i > 0:
            fig.add_subplot(n_row, n_col, i + 1)
        figx = plot_tracks(model, file_name=None, fig=fig)
        fig.title(model.name)
    fig.add_subplot(n_row, n_col, 6)
    fig = plot_CI(L_model, file_name=None, fig=fig)
    fig.savefig('fig0', path=path)
    return fig
