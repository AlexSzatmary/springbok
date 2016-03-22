import numpy as np
from . import platypus
import springbok
import os


def plot_tracks(model, file_name=True, path=None, Style=None, fig=None):
    if fig is None:
        fig = Style()
    if file_name is True:
        file_name = 'track-' + model.name
    # xmin, xmax = (model.L_cell_group[0].xya[0], model.L_cell_group[0].xyb[0])
    # ymin, ymax = (model.L_cell_group[0].xya[1], model.L_cell_group[0].xyb[1])
    xmin, xmax = 500., 2500.
    ymin, ymax = -500., 1500.
    fig = platypus.multi_plot(
        [n.a_xy[:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        xlabel=r'x, $\mathrm{\mu m}$', ylabel=r'y, $\mathrm{\mu m}$',
        xlim=(xmin, xmax), ylim=(ymin, ymax),
        color_f=platypus.color_f_color,
        file_name=None)
    figx = platypus.multi_plot(
        [n.a_xy[-1:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[-1:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        L_marker=['o'] * len(model.L_cell_group[0].L_cell),
        markerfacecolor='none', markersize=(4e0),
        file_name=None)
    if len(model.L_cell_group) == 2:
        figx = platypus.multi_plot(
            [n.a_xy[:, 0] for n in model.L_cell_group[1].L_cell],
            [n.a_xy[:, 1] for n in model.L_cell_group[1].L_cell],
            fig=fig,
            xlim=(xmin, xmax), ylim=(ymin, ymax),
            color_f=platypus.color_f_black,
            file_name=None)
        figx = platypus.multi_plot(
            [n.a_xy[-1:, 0] for n in model.L_cell_group[1].L_cell],
            [n.a_xy[-1:, 1] for n in model.L_cell_group[1].L_cell],
            fig=fig,
            L_marker=['o'] * len(model.L_cell_group[1].L_cell),
            markerfacecolor='none', markersize=(4e0),
            file_name=None)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def do_all_post_processing(L_model):
    list(map(get_a_CI, L_model))

def CI(cell):
    return ((cell.a_xy[-1, 0] - cell.a_xy[0, 0]) /
            (cell.speed * cell.persistence * (cell.n_t - 1)))

def get_a_CI(model):
    model.a_CI = [CI(c)
                  for c in model.L_cell_group[0].L_cell]

def plot_CI(L_model, file_name='CI', path=None, Style=None, fig=None):
    if fig is None:
        fig = Style()
    platypus.boxplot([model.a_CI for model in L_model],
                labels=[model.name for model in L_model], fig=fig)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def plot_CI_nsg(L_model_nsg, file_name='CI', path='plots/neutrophil_static_gradient/print'):
#    if fig is None:
    fig = platypus.Print()
    platypus.boxplot([model.a_CI for model in L_model_nsg],
                labels=['{:.0f}'.format(float(model.name[4:])) for model in L_model_nsg],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index')
    if file_name:
        fig.savefig(file_name, path=path)
    return fig
    

def plot_CI_n_mixed(L_model_n_mixed, file_name='CI', path='plots/n_mixed/print'):
#    if fig is None:
    fig = platypus.Print()
    platypus.boxplot([model.a_CI for model in L_model_n_mixed],
                labels=['{:.0f}'.format(float(model.name[8:])) for model in L_model_n_mixed],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index', title='FPR+')
    if file_name:
        fig.savefig(file_name, path=path)
    fig = platypus.Print()
    platypus.boxplot([[CI(c)
                  for c in model.L_cell_group[1].L_cell] for model in L_model_n_mixed],
                labels=['{:.0f}'.format(float(model.name[8:])) for model in L_model_n_mixed],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index', title='FPR-')
    if file_name:
        fig.savefig(file_name + '-FPR0', path=path)
    return fig
    

def plot_CI_n_BLT(L_model_n_BLT, file_name='CI', path='plots/n_BLT/print'):
#    if fig is None:
    fig = platypus.Print(axes=[0.4, 0.25, 0.55, 0.65])
    platypus.boxplot([model.a_CI for model in L_model_n_BLT],
                labels=['{:.0f}'.format(model.ell) + r' $\mathrm{\mu m}$,' +
                  ' BLT{}'.format('+' if model.has_BLT else '-')
                        for model in L_model_n_BLT],
                     fig=fig, #xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', 
                     xlabel='Chemotactic index', vert=False)
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

def make_all_plots_n_BLT(L_model, path=None, Style=platypus.Print):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, Style.style)
    L_fig = [plot_tracks(model, path=path, Style=Style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, Style=Style)
    L_fig.append(fig)
    return L_fig

def make_all_plots_n_BLT_play(L_model, path=None, Style=platypus.Print):
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

def plot_count(L_model, L_legend=['0', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$'], fig=None):
    L_count = [np.sum(np.array([cell.a_xy[:, 0] > 2000. for cell in r.L_cell_group[0].L_cell]), axis=0) for r in L_model]
    fig = platypus.multi_plot([np.arange(0., 61.)] * len(L_count), L_count, xlabel='Time, min', ylabel='Number of cells migrated\nbeyond starting line', fig=fig); fig.legend(L_legend, loc='upper left', bbox_to_anchor=(0., 1.2)); #fig.savefig('count-t', path='plots/n_BLT_play')


def make_n_BLT_play_fig(L_model):
    L_model = [L_model[0], L_model[1], L_model[3], L_model[5], L_model[7]]
    n_row, n_col = (3, 2)
    fig = platypus.Print(subplot=(n_row, n_col, 1))
    path = os.path.join('plots', L_model[0].set_name, fig.style)
    figx = plot_count(L_model, fig=fig)
    for (i, model) in enumerate(L_model):
#        if i > 0:
        fig.add_subplot(n_row, n_col, i + 2)
        figx = plot_tracks(model, file_name=None, fig=fig)
        fig.title(model.name)
#    fig.add_subplot(n_row, n_col, 6)
#    fig = plot_CI(L_model, file_name=None, fig=fig)
    fig.savefig('fig0', path=path)
    return fig

def make_n_BLT_play_fig(L_model):
    L_model = [L_model[0], L_model[1], L_model[3], L_model[5], L_model[7]]
    n_row, n_col = (3, 2)
    fig = platypus.Print(subplot=(n_row, n_col, 1))
    path = os.path.join('plots', L_model[0].set_name, fig.style)
    figx = plot_count(L_model, fig=fig)
    for (i, model) in enumerate(L_model):
#        if i > 0:
        fig.add_subplot(n_row, n_col, i + 2)
        figx = plot_tracks(model, file_name=None, fig=fig)
        fig.title(model.name)
#    fig.add_subplot(n_row, n_col, 6)
#    fig = plot_CI(L_model, file_name=None, fig=fig)
    fig.savefig('fig0', path=path)
    return fig

def confection0(model):
    fig = platypus.Print(subplot=(2, 1, 1), panesize=(3., 3.))
    path = os.path.join('plots', model.set_name, fig.style)
    figx = plot_tracks(model, file_name=None, fig=fig)    
    if model.set_name == 'n_BLT':
        fig.title('{:.0f}'.format(model.ell) + r' $\mathrm{\mu m}$,' +
                  ' BLT{}'.format('+' if model.has_BLT else '-'))
    elif model.set_name == 'neutrophil_static_gradient':
        fig.title('{:.0f}'.format(float(model.name[4:])))
    elif model.set_name == 'n_mixed':
        fig.title('{:.0f}'.format(float(model.name[8:])))
    if hasattr(model, 'title'):
        fig.title(model.title)
                
    fig.add_subplot(4, 1, 3)
    fMLP_pde = model.pde_stepper.L_pde[0]
    fig.multi_plot([fMLP_pde.x], [fMLP_pde.u[0]],
                   color_f=platypus.color_f_black,
                   xlabel=r'x, $\mathrm{\mu m}$', ylabel='Concentration',
                   xlim=(500., 2500.), ylim=(0., 1.))
    if len(model.pde_stepper.L_pde) >= 2:
        LTB_pde = model.pde_stepper.L_pde[1]
        fig.multi_plot(
            [LTB_pde.x] * 11, LTB_pde.u[::6],
            color_f=platypus.color_f_color)
    fig.add_subplot(4, 1, 4)
    # platypus.boxplot([model.a_CI for model in L_model],
    #                  labels=[model.name for model in L_model], fig=fig)
    fig.multi_plot([[cell.a_xy[0, 0]
                    for cell in model.L_cell_group[0].L_cell]],
                   # [cell.a_xy[0, 0]
                   #  for cell in model.L_cell_group[0].L_cell],)
                   [[CI(cell)
                    for cell in model.L_cell_group[0].L_cell]], L_marker=['.'],
                   L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                   xlim=(500., 2500.),
                   ylabel='Chemotactic index', ylim=(-0.2, 1.))
    if len(model.L_cell_group) >= 2:
        fig.multi_plot([[cell.a_xy[0, 0]
                        for cell in model.L_cell_group[1].L_cell]],
                       # [cell.a_xy[0, 0]
                       #  for cell in model.L_cell_group[0].L_cell],)
                       [[CI(cell)
                        for cell in model.L_cell_group[1].L_cell]], L_marker=['.'],
                       L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                       xlim=(500., 2500.), color_f=platypus.color_f_color,
                       ylabel='Chemotactic index', ylim=(-0.2, 1.))
    fig.savefig(model.name + 'confection0', path=path)
    return fig

def confection1(model):
    fig = platypus.Print(subplot=(3, 1, 1), panesize=(3., 3.))
    path = os.path.join('plots', model.set_name, fig.style)
    figx = plot_tracks(model, file_name=None, fig=fig)    
    if model.set_name == 'n_BLT':
        fig.title('{:.0f}'.format(model.ell) + r' $\mathrm{\mu m}$,' +
                  ' BLT{}'.format('+' if model.has_BLT else '-'))
    elif model.set_name == 'neutrophil_static_gradient':
        fig.title('{:.0f}'.format(float(model.name[4:])))
    elif model.set_name == 'n_mixed':
        fig.title('{:.0f}'.format(float(model.name[8:])))
    if hasattr(model, 'title'):
        fig.title(model.title)
                
    fig.add_subplot(6, 1, 3)
    # platypus.boxplot([model.a_CI for model in L_model],
    #                  labels=[model.name for model in L_model], fig=fig)
    fig.multi_plot([[cell.a_xy[0, 0]
                    for cell in model.L_cell_group[0].L_cell]],
                   # [cell.a_xy[0, 0]
                   #  for cell in model.L_cell_group[0].L_cell],)
                   [[CI(cell)
                    for cell in model.L_cell_group[0].L_cell]], L_marker=['.'],
                   L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                   xlim=(500., 2500.),
                   ylabel='Chemotactic index', ylim=(-0.2, 1.))
    if len(model.L_cell_group) >= 2:
        fig.multi_plot([[cell.a_xy[0, 0]
                        for cell in model.L_cell_group[1].L_cell]],
                       # [cell.a_xy[0, 0]
                       #  for cell in model.L_cell_group[0].L_cell],)
                       [[CI(cell)
                        for cell in model.L_cell_group[1].L_cell]], L_marker=['.'],
                       L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                       xlim=(500., 2500.), color_f=platypus.color_f_color,
                       ylabel='Chemotactic index', ylim=(-0.2, 1.))
    fig.add_subplot(6, 1, 4)
    Nt = model.pde_stepper.Nt
    n = 11
    fMLP_pde = model.pde_stepper.L_pde[0]
    if hasattr(fMLP_pde, 'gamma_F') and fMLP_pde.gamma_F != 0.:
        fig.multi_plot([fMLP_pde.x] * n,
#                       fMLP_pde.u[0::(Nt-1) / (n - 1)],
                       fMLP_pde.u[0:n],
                       xlabel=r'x, $\mathrm{\mu m}$', ylabel='fMLP',
                       xlim=(500., 2500.), ylim=(0., 1.),
                       color_f=lambda i: platypus.greys_color_f(1. - float(i) / n))
    else:
        fig.multi_plot([fMLP_pde.x],
                       fMLP_pde.u[0:1],
                       xlabel=r'x, $\mathrm{\mu m}$', ylabel='fMLP',
                       xlim=(500., 2500.), ylim=(0., 1.),
                       color_f=platypus.color_f_black)
    fig.add_subplot(6, 1, 5)
    if len(model.pde_stepper.L_pde) >= 3:
        exo_pde = model.pde_stepper.L_pde[1]
        fig.multi_plot(
            [exo_pde.x] * n, exo_pde.u[::(Nt-1)/(n-1)],
            xlabel=r'x, $\mathrm{\mu m}$', ylabel='Exosomes',
            color_f=lambda i: platypus.blues_color_f(float(i) / n))
    fig.add_subplot(6, 1, 6)
    fMLP_pde = model.pde_stepper.L_pde[0]
    if len(model.pde_stepper.L_pde) >= 2:
        LTB_pde = model.pde_stepper.L_pde[2]
        fig.multi_plot(
            [LTB_pde.x] * n, LTB_pde.u[::(Nt-1)/(n-1)],
            xlabel=r'x, $\mathrm{\mu m}$', ylabel=r'$\mathrm{LTB_{4}}$',
            color_f=lambda i: platypus.blues_color_f(float(i) / n))

    fig.savefig(model.name + 'confection1', path=path)
    return fig

def get_range_time_averaged(model, CI_threshold):
    return np.mean([get_range(model, CI_threshold, j) for j in range(model.pde_stepper.Nt)])

def get_range(model, CI_threshold, j):
    a_kappa = get_a_kappa(model, j)
#    return np.sum(a_kappa > 1.16) * model.pde_stepper.L_pde[0].dx # 1.16 is kappa threshold for CI_threshold = 0.5
    return np.sum(a_kappa > 0.4) * model.pde_stepper.L_pde[0].dx # 0.4 is kappa threshold for CI_threshold = 0.2

def get_a_kappa(model, j):
    a_L_condition = np.array([(pde.u[j, 1:-1], springbok.tiger.opdudx(pde.u[j], pde.dx))
                              for pde in model.pde_stepper.L_pde])
    cell = model.L_cell_group[0].L_cell[0]
    return np.array([cell.kappa(a_L_condition[:, :, i]) for i in range(np.size(a_L_condition, 2))])
