import flux
import numpy as np
from . import platypus
import springbok
import os


def plot_tracks(model, file_name=True, path=None, Style=None, fig=None,
                xlim=None):
    if fig is None:
        fig = Style()
    if file_name is True:
        file_name = 'track-' + model.name
    # xmin, xmax = (model.L_cell_group[0].xya[0], model.L_cell_group[0].xyb[0])
    # ymin, ymax = (model.L_cell_group[0].xya[1], model.L_cell_group[0].xyb[1])
    if xlim:
        delta_x = xlim[1] - xlim[0]
        ymin, ymax = -delta_x / 2, delta_x / 2
    fig = platypus.multi_plot(
        [n.a_xy[:, 0] for n in model.L_cell_group[0].L_cell],
        [n.a_xy[:, 1] for n in model.L_cell_group[0].L_cell],
        fig=fig,
        xlabel=r'x, $\mathrm{\mu m}$', ylabel=r'y, $\mathrm{\mu m}$',
        xlim=xlim, ylim=(ymin, ymax),
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
            xlim=xlim, ylim=(ymin, ymax),
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
    path = os.path.join('plots', model.path_name, fig.style)
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
    path = os.path.join('plots', model.path_name, fig.style)
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


def plot_CI_x(model, fig=None):
    fig.multi_plot([[cell.a_xy[0, 0]
                    for cell in model.L_cell_group[0].L_cell]],
                   # [cell.a_xy[0, 0]
                   #  for cell in model.L_cell_group[0].L_cell],)
                   [[CI(cell)
                    for cell in model.L_cell_group[0].L_cell]], L_marker=['.'],
                   L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                   ylabel='Chemotactic index', ylim=(-0.2, 1.))


def plot_concentration(model, pde, fig=None, n=11, **kwargs):
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([pde.x] * n,
                       pde.u[0::(Nt-1) / (n - 1)],
#                       pde.u[0:n],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       color_f=lambda i: platypus.blues_color_f(float(i) / n), **kwargs)
    else:
        fig.multi_plot([pde.x],
                       pde.u[0:1],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       ylim=(0., 1.),
                       color_f=platypus.color_f_black)
                       

def plot_fMLP(model, fig=None, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = 'fMLP'
    if 'ylim' not in kwargs:
        kwargs['ylim'] = (0., 1.)
    fMLP_pde = model.pde_stepper.L_pde[0]
    plot_concentration(model, fMLP_pde, fig=fig, **kwargs)


def plot_exo(model, fig=None, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = 'Exosome\nactivity'
    exo_pde = model.pde_stepper.L_pde[1]
    plot_concentration(model, exo_pde, fig=fig, **kwargs)


def plot_LTB(model, fig=None, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$\mathrm{LTB_{4}}$'
    LTB_pde = model.pde_stepper.L_pde[2]
    plot_concentration(model, LTB_pde, fig=fig, **kwargs)


def plot_kappa(model, fig=None, n=11, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$\kappa$'
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1]] * n,
#                       pde.u[0::(Nt-1) / (n - 1)],
                       [get_a_kappa(model, j) for j in
                        range(0, Nt, int((Nt-1) / (n - 1)))],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       color_f=lambda i: platypus.blues_color_f(float(i) / n), **kwargs)
    else:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1]] * n,
                       [get_a_kappa(model, 0)],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       ylim=(0., 1.),
                       color_f=platypus.color_f_black)
                       

def plot_cos(model, fig=None, n=11, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$<\cos \theta>(x)$'
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1]] * n,
#                       pde.u[0::(Nt-1) / (n - 1)],
                       [flux.mean_velocity_fast(get_a_kappa(model, j)) for j in
                        range(0, Nt, int((Nt-1) / (n - 1)))],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       color_f=lambda i: platypus.blues_color_f(float(i) / n), **kwargs)
    else:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1]] * n,
                       [flux.mean_velocity_fast(get_a_kappa(model, 0))],
                       xlabel=r'x, $\mathrm{\mu m}$',
                       ylim=(0., 1.),
                       color_f=platypus.color_f_black)
                       

def confection_n_FPR0(n_FPR0, file_name='confection_n_FPR0',
                      Style=platypus.Print,
                      **kwargs):
    kwargs['xlim'] = (0., 3000.)
    model = n_FPR0.L_runs[0]
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3., 3.))
    fig = Style(subplot=(2, 1, 1), **kw0)
    plot_tracks(model, file_name=False, fig=fig, **kwargs)
    fig.add_subplot(4, 1, 3)
    plot_fMLP(model, fig=fig, ylog=False, **kwargs)
    fig.add_subplot(4, 1, 4)
    plot_LTB(model, fig=fig, ylog=False, **kwargs)
    print(fig.AB_font_properties.get_size())
    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_BLT0(n_BLT0, file_name='confection_n_BLT0',
                      Style=platypus.Print,
                      **kwargs):
    kwargs['xlim'] = (0., 3000.)
    model = n_BLT0.L_runs[0]
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3., 3.))
    fig = Style(subplot=(2, 2, 1), **kw0)
    fig.title('BLT-')
    plot_tracks(n_BLT0.L_runs[0], file_name=False, fig=fig, **kwargs)
    fig.add_subplot(2, 2, 2)
    plot_tracks(n_BLT0.L_runs[1], file_name=False, fig=fig, **kwargs)
    fig.title('BLT+')
    fig.add_subplot(4, 2, 5)
    plot_fMLP(n_BLT0.L_runs[0], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(4, 2, 6)
    plot_fMLP(n_BLT0.L_runs[1], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(4, 2, 7)
    plot_LTB(n_BLT0.L_runs[0], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(4, 2, 8)
    plot_LTB(n_BLT0.L_runs[1], fig=fig, ylog=False, **kwargs)
    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_BLT0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_init_rec(init_rec, r_L=1e5, Style=platypus.Print):
    if Style is platypus.Poster:
        kw0 = dict(panesize=(7., 3.5))
    else:
        kw0 = dict(panesize=(3., 1.5))
    fig = Style(subplot=(3, 2, 1), **kw0)
    path = os.path.join('plots', init_rec.set_name, fig.style)
    kwargs = {}
    kwargs['xlim'] = (600., 3600.)
    plot_fMLP(init_rec.d_runs[(r_L, 0.)], fig=fig, **kwargs)
    fig.add_subplot(3, 2, 2)
    plot_fMLP(init_rec.d_runs[(r_L, 1.)], fig=fig, **kwargs)
    fig.add_subplot(3, 2, 3)
    plot_LTB(init_rec.d_runs[(r_L, 0.)], fig=fig, ylog=False, ylim=(0, 8), yint=True, **kwargs)
    fig.add_subplot(3, 2, 4)
    plot_LTB(init_rec.d_runs[(r_L, 1.)], fig=fig, ylog=False, ylim=(0, 8), yint=True, **kwargs)
    fig.add_subplot(3, 2, 5)
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), n=1, **kwargs)
    fig.add_subplot(3, 2, 6)
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), n=1, **kwargs)
    fig.set_AB_labels()
    fig.savefig('confection_init_rec', path=path)
    return fig

def confection_decay_F(decay_F, Style=platypus.Print, **kwargs):
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3., 3.))
    fig = Style(subplot=(3, 2, 1), **kw0)
    path = os.path.join('plots', decay_F.set_name, fig.style)
    kwargs['xlim'] = (0000., 3000.)
    figx = plot_tracks(decay_F.L_runs[0], file_name=None, fig=fig, **kwargs)
    fig.add_subplot(3, 2, 2)
    figx = plot_tracks(decay_F.L_runs[1], file_name=None, fig=fig, **kwargs)
    fig.add_subplot(6, 2, 5)
    ax = fig.fig.gca()
#    ax.yaxis.get_major_formatter().set_scientific(True)
#    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plot_fMLP(decay_F.L_runs[0], fig=fig, ylim=(0., 0.011), **kwargs)
    fig.fig.canvas.draw()
    fig.add_subplot(6, 2, 6)
    plot_fMLP(decay_F.L_runs[1], fig=fig, ylim=(0., 0.011), **kwargs)
    # fig.add_subplot(6, 2, 7)
    # plot_exo(decay_F.L_runs[0], fig=fig, **kwargs)
    fig.add_subplot(6, 2, 8)
    plot_exo(decay_F.L_runs[1], fig=fig, ylim=(0., 1.), **kwargs)
    fig.add_subplot(6, 2, 9)
    plot_LTB(decay_F.L_runs[0], fig=fig, ylog=False, ylim=(0., 1e0), **kwargs)
    fig.add_subplot(6, 2, 10)
    plot_LTB(decay_F.L_runs[1], fig=fig, ylog=False, ylim=(0., 1e0), **kwargs)
    fig.add_subplot(6, 2, 11)
    plot_cos(decay_F.L_runs[0], fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.add_subplot(6, 2, 12)
    plot_cos(decay_F.L_runs[1], fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.set_AB_labels()
    fig.savefig('confection_decay_F', path=path)
    return fig

def get_range_time_averaged(model, CI_threshold):
    return np.mean([get_range(model, CI_threshold, j) for j in range(model.pde_stepper.Nt)])

def get_range(model, CI_threshold, j):
    a_kappa = get_a_kappa(model, j)
    a_cos = flux.mean_velocity_fast(a_kappa)
    return np.sum(a_cos > CI_threshold) * model.pde_stepper.L_pde[0].dx # 1.16 is kappa threshold for CI_threshold = 0.5
#    return np.sum(a_kappa > 0.4) * model.pde_stepper.L_pde[0].dx # 0.4 is kappa threshold for CI_threshold = 0.2

def get_a_kappa(model, j):
    a_L_condition = np.array([(pde.u[j, 1:-1], springbok.tiger.opdudx(pde.u[j], pde.dx))
                              for pde in model.pde_stepper.L_pde])
    cell = model.L_cell_group[0].L_cell[0]
    return np.array([cell.kappa(a_L_condition[:, :, i]) for i in range(np.size(a_L_condition, 2))])

def plot_range_vary_phi_E(n_exo, fig=None, Style=platypus.Print,
                          file_name='vary-phi_E'):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [n_exo.L_phi_E] * len(n_exo.L_r_L),
        [[get_range(n_exo.d_runs[r_L, phi_E], None, 60)
          for phi_E in n_exo.L_phi_E] for r_L in n_exo.L_r_L],
        L_legend=[r'$r_L=10^{{{}}}$'.format(int(np.log10(r_L))) if r_L != 0. else r'$r_L=0$' for r_L in n_exo.L_r_L],
        xlabel=r'Fraction of LTB$_4$ secreted via exosomes, $\phi_E$',
        ylabel='Range for directed migration, $\mu m$',
        fig=fig)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_r_L(n_exo, fig=None, Style=platypus.Print,
                          file_name='vary-r_L'):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [n_exo.L_r_L] * len(n_exo.L_r_L),
        [[get_range(n_exo.d_runs[r_L, phi_E], None, 60)
          for r_L in n_exo.L_r_L] for phi_E in n_exo.L_phi_E],
        xlog=True,
        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in n_exo.L_phi_E],
        xlabel='Characteristic LTB$_4$ secretion rate, $r_L$',
        ylabel='Range for directed migration, $\mu m$',
        fig=fig)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range(n_exo, Style=platypus.Print, file_name='range'):
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3., 3.))
    fig = Style(subplot=(2, 2, 1), **kw0)
    figx = plot_range_vary_r_L(n_exo, fig=fig, file_name=None)
    fig.add_subplot(2, 2, 3)
    figx = plot_range_vary_phi_E(n_exo, fig=fig, file_name=None)
    fig.set_AB_labels(loc='upper right')
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig
