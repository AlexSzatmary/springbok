import copy
import flux
import numpy as np
from . import platypus
import matplotlib
import springbok
import os


def plot_tracks(model, file_name=True, path=None, Style=None, fig=None,
                color_f=platypus.color_f_color,
                cell_group_i=0,
                xint=True, yint=True,
                xlim=None, use_mm=True, **kwargs):
    if use_mm:
        ratio = 1e3
    else:
        ratio = 1
    if use_mm:
        xlabel = r'x, mm'
        ylabel = r'y, mm'
    else:
        xlabel = r'x, $\mathrm{\mu m}$'
        ylabel = r'y, $\mathrm{\mu m}$'

    if fig is None:
        fig = Style()
    if file_name is True:
        file_name = 'track-' + model.name
    # xmin, xmax = (model.L_cell_group[0].xya[0], model.L_cell_group[0].xyb[0])
    # ymin, ymax = (model.L_cell_group[0].xya[1], model.L_cell_group[0].xyb[1])
    if xlim:
        delta_x = xlim[1] - xlim[0]
        ymin, ymax = -delta_x / 2, delta_x / 2
    def y_offset(a_xy):
        return np.floor((a_xy[-1, 1] - ymin) / delta_x) * delta_x
    fig = platypus.multi_plot(
        [n.a_xy[:, 0] / ratio for n in model.L_cell_group[cell_group_i].L_cell],
        [n.a_xy[:, 1] / ratio - y_offset(n.a_xy) / ratio for n in model.L_cell_group[cell_group_i].L_cell],
        fig=fig,
        xlabel=xlabel, ylabel=ylabel,
        xlim=[x / ratio for x in xlim], ylim=(ymin / ratio, ymax / ratio),
        color_f=color_f,
        xint=xint, yint=yint,
        file_name=None, **kwargs)
    figx = platypus.multi_plot(
        [np.array([n.a_xy[-1, 0] / ratio for n in model.L_cell_group[cell_group_i].L_cell])],
        [np.array([n.a_xy[-1, 1] / ratio - y_offset(n.a_xy) / ratio for n in model.L_cell_group[cell_group_i].L_cell])],
        fig=fig,
        L_linestyle='none',
        L_marker=['o'],# * len(model.L_cell_group[cell_group_i].L_cell),
        markerfacecolor='none', markersize=(4e0),
        file_name=None, **kwargs)
    # if len(model.L_cell_group) == 2:
    #     figx = platypus.multi_plot(
    #         [n.a_xy[:, 0] for n in model.L_cell_group[1].L_cell],
    #         [n.a_xy[:, 1] for n in model.L_cell_group[1].L_cell],
    #         fig=fig,
    #         xlim=xlim, ylim=(ymin, ymax),
    #         color_f=platypus.color_f_black,
    #         file_name=None)
    #     figx = platypus.multi_plot(
    #         [n.a_xy[-1:, 0] for n in model.L_cell_group[1].L_cell],
    #         [n.a_xy[-1:, 1] for n in model.L_cell_group[1].L_cell],
    #         fig=fig,
    #         L_marker=['o'] * len(model.L_cell_group[1].L_cell),
    #         markerfacecolor='none', markersize=(4e0),
    #         file_name=None)
    ax = fig.fig.gca()
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['top']._linewidth = 0.5
    ax.spines['right']._linewidth = 0.5
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
    fig = platypus.MBOC()
    platypus.boxplot([model.a_CI for model in L_model_nsg],
                labels=['{:.0f}'.format(float(model.name[4:])) for model in L_model_nsg],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index')
    if file_name:
        fig.savefig(file_name, path=path)
    return fig
    

def plot_CI_n_mixed(L_model_n_mixed, file_name='CI', path='plots/n_mixed/print'):
#    if fig is None:
    fig = platypus.MBOC()
    platypus.boxplot([model.a_CI for model in L_model_n_mixed],
                labels=['{:.0f}'.format(float(model.name[8:])) for model in L_model_n_mixed],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index', title='FPR+')
    if file_name:
        fig.savefig(file_name, path=path)
    fig = platypus.MBOC()
    platypus.boxplot([[CI(c)
                  for c in model.L_cell_group[1].L_cell] for model in L_model_n_mixed],
                labels=['{:.0f}'.format(float(model.name[8:])) for model in L_model_n_mixed],
                     fig=fig, xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', ylabel='Chemotactic index', title='FPR-')
    if file_name:
        fig.savefig(file_name + '-FPR0', path=path)
    return fig
    

def plot_CI_n_BLT(L_model_n_BLT, file_name='CI', path='plots/n_BLT/print'):
#    if fig is None:
    fig = platypus.MBOC(axes=[0.4, 0.25, 0.55, 0.65])
    platypus.boxplot([model.a_CI for model in L_model_n_BLT],
                labels=['{:.0f}'.format(model.ell) + r' $\mathrm{\mu m}$,' +
                  ' BLT{}'.format('+' if model.has_BLT else '-')
                        for model in L_model_n_BLT],
                     fig=fig, #xlabel=r'Gradient characteristic length, $\mathrm{\mu m}$', 
                     xlabel='Chemotactic index', vert=False)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig
    

def make_all_plots(L_model, path=None, Style=platypus.MBOC):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, Style.style)
    L_fig = [plot_tracks(model, path=path, Style=Style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, Style=Style)
    L_fig.append(fig)
    return L_fig

def make_all_plots_n_BLT(L_model, path=None, Style=platypus.MBOC):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, Style.style)
    L_fig = [plot_tracks(model, path=path, Style=Style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, Style=Style)
    L_fig.append(fig)
    return L_fig

def make_all_plots_n_BLT_play(L_model, path=None, Style=platypus.MBOC):
    if path is None:
        path = os.path.join('plots', L_model[0].set_name, Style.style)
    L_fig = [plot_tracks(model, path=path, Style=Style)
             for model in L_model]
    fig = plot_CI(L_model, path=path, Style=Style)
    L_fig.append(fig)
    return L_fig

def make_n_static_fig(L_model):
    n_row, n_col = (3, 2)
    fig = platypus.MBOC(subplot=(n_row, n_col, 1))
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
    fig = platypus.MBOC(subplot=(n_row, n_col, 1))
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
    fig = platypus.MBOC(subplot=(n_row, n_col, 1))
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
    fig = platypus.MBOC(subplot=(2, 1, 1), panesize=(3.3, 3.3))
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
    fig = platypus.MBOC(subplot=(3, 1, 1), panesize=(3.3, 3.3))
    path = os.path.join('plots', model.path_name, fig.style)
    kwargs = {}
    kwargs['xlim'] = (0., 3000.)
    figx = plot_tracks(model, file_name=None, fig=fig, **kwargs)
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
            color_f=lambda i: platypus.blues_color_f2(float(i) / n))
    fig.add_subplot(6, 1, 6)
    fMLP_pde = model.pde_stepper.L_pde[0]
    if len(model.pde_stepper.L_pde) >= 2:
        LTB_pde = model.pde_stepper.L_pde[2]
        fig.multi_plot(
            [LTB_pde.x] * n, LTB_pde.u[::(Nt-1)/(n-1)],
            xlabel=r'x, $\mathrm{\mu m}$', ylabel=r'$\mathrm{LTB_{4}}$',
            color_f=lambda i: platypus.blues_color_f2(float(i) / n))

    fig.savefig(model.name + 'confection1', path=path)
    return fig


def confection2(module, run, Style=platypus.MBOC, title=True, **kwargs):
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(subplot=(3, 1, 1), **kw0)
    if title is True:
        title = run.job_name
    if title:
        fig.title(title)
            
    path = os.path.join('plots', module.set_name, fig.style, 'confection')
    kwargs['xlim'] = (0000., 3000.)
    figx = plot_tracks(run, file_name=None, fig=fig, **kwargs)
    fig.add_subplot(6, 1, 3)
    ax = fig.fig.gca()
#    ax.yaxis.get_major_formatter().set_scientific(True)
#    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plot_fMLP(run, fig=fig, ylim=(0., 0.011), **kwargs)
    fig.fig.canvas.draw()
    fig.add_subplot(6, 1, 4)
    plot_exo(run, fig=fig, **kwargs)
    fig.add_subplot(6, 1, 5)
    plot_LTB(run, fig=fig, ylog=False, ylim=(0., 1e0), **kwargs)
    fig.add_subplot(6, 1, 6)
    plot_cos(run, fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.set_AB_labels()
    fig.savefig('confection2-' + run.name, path=path)
    return fig



def plot_CI_x(model, fig=None, **kwargs):
    fig.multi_plot([[cell.a_xy[0, 0]
                    for cell in model.L_cell_group[0].L_cell]],
                   # [cell.a_xy[0, 0]
                   #  for cell in model.L_cell_group[0].L_cell],)
                   [[CI(cell)
                    for cell in model.L_cell_group[0].L_cell]], L_marker=['.'],
                   L_linestyle=['None'], xlabel=r'x, $\mathrm{\mu m}$',
                   ylabel='Chemotactic index', ylim=(-0.2, 1.), **kwargs)


def plot_concentration(model, pde, fig=None, n=11, use_mm=True,
                       color_f=platypus.blues_color_f2,
                       **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
        if 'xlim' in kwargs:
            kwargs['xlim'] = (kwargs['xlim'][0] / ratio,
                              kwargs['xlim'][1] / ratio)
    else:
        ratio = 1
        units = r'$\mu m$'
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([pde.x / ratio] * n,
                       pde.u[0::(Nt-1) / (n - 1)],
#                       pde.u[0:n],
                       xlabel=r'x, ' + units,
                       color_f=lambda i: color_f(float(i) / n), **kwargs)
    else:
        fig.multi_plot([pde.x / ratio],
                       pde.u[0:1],
                       xlabel=r'x, ' + units,
                       ylim=(0., 1.),
                       color_f=platypus.color_f_black)
                       

def plot_fMLP(model, fig=None, use_mm=True, xint=True, yint=True, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = 'fMLP'
    if 'ylim' not in kwargs:
        kwargs['ylim'] = (0., 2.)
    fMLP_pde = model.pde_stepper.L_pde[0]
    plot_concentration(model, fMLP_pde, fig=fig, use_mm=use_mm,
                       xint=xint, yint=yint,
                       color_f=platypus.greys_color_f2,
                       **kwargs)


def plot_exo(model, fig=None, use_mm=True, xint=True, yint=True, **kwargs):
    if 'ylabel' not in kwargs:
        #        kwargs['ylabel'] = 'Exosome activity,\n$\mathrm{K_d/min}$'
        kwargs['ylabel'] = 'Exosome activity'
    exo_pde = model.pde_stepper.L_pde[1]
    plot_concentration(model, exo_pde, fig=fig, use_mm=use_mm,
                       xint=xint, yint=yint, **kwargs)


def plot_LTB(model, fig=None, use_mm=True, xint=True, yint=True, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$\mathrm{LTB_{4}}$'
    if 'ylim' not in kwargs:
        kwargs['ylim'] = (0., 2.)
    LTB_pde = model.pde_stepper.L_pde[2]
    plot_concentration(model, LTB_pde, fig=fig, use_mm=use_mm,
                       xint=xint, yint=yint,
                       color_f=platypus.x_color_f2('Reds'),
                       **kwargs)


def plot_kappa(model, fig=None, use_mm=True, xint=True, yint=True, n=11, f_kappa=None, **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
        if 'xlim' in kwargs:
            kwargs['xlim'] = (kwargs['xlim'][0] / ratio,
                              kwargs['xlim'][1] / ratio)
    else:
        ratio = 1
        units = r'$\mu m$'
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$\kappa$'
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1] / ratio] * n,
#                       pde.u[0::(Nt-1) / (n - 1)],
                       [get_a_kappa(model, j, f_kappa=f_kappa) for j in
                        range(0, Nt, int((Nt-1) / (n - 1)))],
                       xlabel=r'x, ' + units,
                       color_f=lambda i: platypus.blues_color_f2(float(i) / n),
                       xint=xint, yint=yint,
                       **kwargs)
    else:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1]] * n,
                       [get_a_kappa(model, 0, f_kappa=f_kappa)],
                       xlabel=r'x, ' + units,
                       ylim=(0., 1.),
                       color_f=platypus.color_f_black,
                       xint=xint, yint=yint)
                       

def plot_cos(model, fig=None, use_mm=True, xint=True, yint=True, n=11, color_f=platypus.greys_color_f2, f_kappa=None, ylim=(0., 1.), **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
        if 'xlim' in kwargs:
            kwargs['xlim'] = (kwargs['xlim'][0] / ratio,
                              kwargs['xlim'][1] / ratio)
    else:
        ratio = 1
        units = r'$\mu m$'
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$<\cos \theta>(x)$'
    Nt = model.pde_stepper.Nt
    if n != 1:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1] / ratio] * n,
#                       pde.u[0::(Nt-1) / (n - 1)],
                       [flux.mean_velocity_fast(get_a_kappa(model, j, f_kappa=f_kappa)) for j in
                        range(0, Nt, int((Nt-1) / (n - 1)))],
                       xlabel=r'x, ' + units,
                       xint=xint, yint=yint, ylim=ylim,
                       color_f=lambda i: color_f(float(i) / n), **kwargs)
    else:
        fig.multi_plot([model.pde_stepper.L_pde[0].x[1:-1] / ratio] * n,
                       [flux.mean_velocity_fast(get_a_kappa(model, 0, f_kappa=f_kappa))],
                       xlabel=r'x, ' + units,
                       ylim=ylim,
                       xint=xint, yint=yint,
                       color_f=platypus.color_f_black)
                       

def plot_cos_t(model, ix, fig=None, color_f=platypus.color_f_black, f_kappa=None, **kwargs):
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$<\cos \theta>(x)$'
    Nt = model.pde_stepper.Nt
    fig.multi_plot([np.arange(Nt)],
                   [[np.mean(flux.mean_velocity_fast(get_a_kappa(model, j, f_kappa=f_kappa)[ix]))
                    for j in np.arange(Nt)]],
                   xlabel=r't, min',
                   ylim=(0., 1.),
                   color_f=platypus.color_f_black, **kwargs)

def plot_cos_t_decay_F(decay_F, Style=platypus.MBOC, ix=350,
                       file_name='cos_t_decay_F', r_L=1e1):
    fig = Style();
    path = os.path.join('plots', decay_F.set_name, fig.style)
    figx = plot_cos_t(decay_F.d_runs[1e1, 0.], ix, fig=fig, use_mm=use_mm, color_f=lambda x: platypus.color_f_color(1))
    figx = plot_cos_t(decay_F.d_runs[1e1, 1.], ix, fig=fig, use_mm=use_mm, color_f=lambda x: platypus.color_f_color(1))
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def plot_cos_t_decay_F_vary_gamma_E(decay_F_vary_gamma_E, Style=platypus.MBOC, ix=None,
                                    file_name='cos_t_decay_F_vary_gamma_E', r_L=1e1, phi_E=1., **kwargs):
    if ix is None:
        x = decay_F_vary_gamma_E.prototype().pde_stepper.L_pde[0].x
        xmin = 3000.
        xmax = 4000.
        ix = np.where(np.logical_and(x >= xmin, x <= xmax))[0]
    fig = Style();
    path = os.path.join('plots', decay_F_vary_gamma_E.set_name, fig.style)
    if 'ylabel' not in kwargs:
        kwargs['ylabel'] = r'$<\cos \theta>(x)$'
    Nt = decay_F_vary_gamma_E.prototype().pde_stepper.Nt
    L_x = [np.arange(Nt) for gamma_E in decay_F_vary_gamma_E.L_gamma_E]
    L_y = []
    for gamma_E in decay_F_vary_gamma_E.L_gamma_E:
        run = decay_F_vary_gamma_E.d_runs[r_L, phi_E, gamma_E]
        y = [np.mean(flux.mean_velocity_fast(
                    get_a_kappa(run, j))[ix]) for j in np.arange(Nt)]
        L_y.append(y)
    print([np.shape(L_x), np.shape(L_y)])
    fig.multi_plot(L_x, L_y,
                   xlabel=r't, min',
                   ylim=(0., 1.),
                   color_f=lambda x: platypus.blues_color_f2(x / len(L_x)), **kwargs)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig


def confection_n_FPR0(n_FPR0, file_name='confection_n_FPR0',
                      Style=platypus.MBOC, xlim=(0., 7000.),
                      **kwargs):
    kwargs['xlim'] = xlim
    model = n_FPR0.L_runs[-1]
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))

    fig = Style(subplot=(2, 1, 1), **kw0)
    plot_tracks(model, file_name=False, fig=fig, **kwargs)
    fig.add_subplot(4, 1, 3)
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 4.), **kwargs)
    fig.add_subplot(4, 1, 4)
    plot_LTB(model, fig=fig, ylog=False, **kwargs)
    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_FPR0_2(
        n_FPR0, file_name='confection_n_FPR0_2',
        Style=platypus.MBOC, xlim=(0., 3000.),
        **kwargs):
    kwargs['xlim'] = xlim
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    model = n_FPR0.L_runs[0]
    fig = Style(subplot=(3, 1, 1), **kw0)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    fig.title('WT')

    fig.add_subplot(3, 1, 2)
    model = n_FPR0.L_runs[1]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(2), **kwargs)
    fig.title('FPR-')

    fig.add_subplot(3, 1, 3)
    model = n_FPR0.L_runs[2]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(2),
                cell_group_i=1,
                **kwargs)
    fig.title('Mixed')
    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_FPR0_3(
        n_FPR0, file_name='confection_n_FPR0_3',
        Style=platypus.MBOC, xlim=(2e3, 4e3),
        **kwargs):
    kwargs['xlim'] = xlim
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    model = n_FPR0.L_runs[0]
    gs = platypus.make_grid_spec(3, 3, height_ratios=[1, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    fig = Style(gs=gs, figsize=figsize, **kw0)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    fig.title('WT', y=0.98)

    fig.add_subplot(gs[1, 0])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.add_subplot(gs[2, 0])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.add_subplot(gs[0, 1])
    model = n_FPR0.L_runs[1]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(2), **kwargs)
    fig.title('FPR-', y=0.98)

    fig.add_subplot(gs[1, 1])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.add_subplot(gs[2, 1])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.add_subplot(gs[0, 2])
    model = n_FPR0.L_runs[2]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(2),
                cell_group_i=1,
                **kwargs)
    fig.title('WT & FPR-', y=0.98)

    fig.add_subplot(gs[1, 2])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.add_subplot(gs[2, 2])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, **kwargs)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig

def confection_n_FPR0_projector(
        n_FPR0, file_name='confection_n_FPR0_3',
        Style=platypus.Projector, xlim=(2e3, 4e3),
        **kwargs):
    kwargs['xlim'] = xlim
    if Style is platypus.Poster:
        kw0 = {}
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
        if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
        if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    model = n_FPR0.L_runs[0]
    gs = platypus.make_grid_spec(3, 3, height_ratios=[1, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    fig = Style(gs=gs, figsize=figsize, **kw0)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    fig.title('WT', y=0.98)

    fig.add_subplot(gs[1, 0])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Primary', **kwargs)

    fig.add_subplot(gs[2, 0])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Secondary', **kwargs)

    fig.add_subplot(gs[0, 1])
    model = n_FPR0.L_runs[1]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(4), **kwargs)
    fig.title('PR-', y=0.98)

    fig.add_subplot(gs[1, 1])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Primary', **kwargs)

    fig.add_subplot(gs[2, 1])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Secondary', **kwargs)

    fig.add_subplot(gs[0, 2])
    model = n_FPR0.L_runs[2]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(4),
                cell_group_i=1,
                **kwargs)
    fig.title('WT & PR-', y=0.98)

    fig.add_subplot(gs[1, 2])
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Primary', **kwargs)

    fig.add_subplot(gs[2, 2])
    plot_LTB(model, fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Secondary', **kwargs)

    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig

def confection_n_FPR0_4(
        n_FPR0, file_name='confection_n_FPR0_4',
        Style=platypus.MBOC, xlim=(2e3, 4e3),
        **kwargs):
    kwargs['xlim'] = xlim
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    model = n_FPR0.L_runs[0]
    gs = matplotlib.gridspec.GridSpec(6, 1)
    tup = (gs.get_geometry()[0], gs.get_geometry()[1], 1)
    fig = Style(subplot=tup, **kw0)
    plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 4.), **kwargs)

    fig.add_subplot(gs.new_subplotspec((1, 0), rowspan=2))
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    fig.title('WT')

    # fig.add_subplot(4, 3, 10)
    # plot_LTB(model, fig=fig, ylog=False, **kwargs)

    # fig.add_subplot(2, 3, 2)
    # model = n_FPR0.L_runs[1]
    # plot_tracks(model, file_name=False, fig=fig,
    #             color_f=lambda x: platypus.color_f_color(2), **kwargs)
    # fig.title('FPR-')

    # fig.add_subplot(4, 3, 8)
    # plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 4.), **kwargs)

    # fig.add_subplot(4, 3, 11)
    # plot_LTB(model, fig=fig, ylog=False, **kwargs)

    fig.add_subplot(gs.new_subplotspec((3, 0), rowspan=2))
    model = n_FPR0.L_runs[2]
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(0), **kwargs)
    plot_tracks(model, file_name=False, fig=fig,
                color_f=lambda x: platypus.color_f_color(2),
                cell_group_i=1,
                **kwargs)
    fig.title('Mixed')

    # fig.add_subplot(4, 3, 9)
    # plot_fMLP(model, fig=fig, ylog=False, ylim=(0., 4.), **kwargs)

    fig.add_subplot(gs[5])
    plot_LTB(model, fig=fig, ylog=False, **kwargs)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_FPR0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


#[post_processing.confection_n_BLT0(n_BLT0, r_L=r_L, file_name='confection_n_BLT0-r_L{:.1e}'.format(r_L)) for r_L in n_BLT0.L_r_L]
def confection_n_BLT0(n_BLT0, file_name='confection_n_BLT0',
                      Style=platypus.MBOC,
                      xlim=(1.5e3, 3.5e3),
                      r_L=1e2,
                      **kwargs):
    kwargs['xlim'] = xlim
    model = n_BLT0.L_runs[0]
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    
    gs = platypus.make_grid_spec(3, 2, height_ratios=[1, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    fig = Style(gs=gs, figsize=figsize, **kw0)
    fig.title('BLT-', y=0.98)
    plot_tracks(n_BLT0.d_runs[r_L, 0., False], file_name=False, fig=fig, **kwargs)
    fig.add_subplot(gs[0, 1])
    plot_tracks(n_BLT0.d_runs[r_L, 0., True], file_name=False, fig=fig, **kwargs)
    fig.title('BLT+', y=0.98)
    fig.add_subplot(gs[1, 0])
    plot_fMLP(n_BLT0.d_runs[r_L, 0., False], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(gs[1, 1])
    plot_fMLP(n_BLT0.d_runs[r_L, 0., True], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(gs[2, 0])
    plot_LTB(n_BLT0.d_runs[r_L, 0., False], fig=fig, ylog=False, ylim=(0., 2.), **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_LTB(n_BLT0.d_runs[r_L, 0., True], fig=fig, ylog=False, ylim=(0., 2.), **kwargs)
    if file_name:
        path = os.path.join('plots', n_BLT0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_BLT0_TmT(
        n_BLT0, file_name='confection_n_BLT0-TmT',
        Style=platypus.MBOC,
        xlim=(1.5e3, 3.5e3),
        r_L=1e2,
        **kwargs):
    kwargs['xlim'] = xlim
    model = n_BLT0.L_runs[0]
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(subplot=(1, 2, 1), axes=[0., 0., 0.8, 0.8], **kw0)
    plot_tracks(n_BLT0.d_runs[r_L, 0., False], file_name=False, fig=fig, **kwargs)
    ax = fig.fig.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for k, v in ax.spines.items():
        v.set_visible(False)
#    fig.fig.add_subplot(1, 2, 2, axes=[0.1, 0., 0.9, 0.9])
    fig.fig.add_axes([0.6, 0., 0.4, 0.8])
    plot_tracks(n_BLT0.d_runs[r_L, 0., True], file_name=False, fig=fig, **kwargs)
    ax = fig.fig.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for k, v in ax.spines.items():
        v.set_visible(False)
    if file_name:
        path = os.path.join('plots', n_BLT0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_BLT0_projector(n_BLT0, file_name='confection_n_BLT0',
                      Style=platypus.Projector,
                      xlim=(1.5e3, 3.5e3),
                      r_L=1e2,
                      **kwargs):
    kwargs['xlim'] = xlim
    model = n_BLT0.L_runs[0]
    print(Style)
    print(issubclass(Style, platypus.Projector))
    if Style is platypus.Poster:
        kw0 = {}
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
        print(kw0, kwargs)
        if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
        if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    
    gs = platypus.make_grid_spec(3, 2, height_ratios=[1, 0.5, 0.5], axes=Style.def_axes)
    print(kw0)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    fig = Style(gs=gs, figsize=figsize, **kw0)
    fig.title('SR-', y=0.98, color='red')
    plot_tracks(n_BLT0.d_runs[r_L, 0., False], file_name=False, fig=fig, **kwargs)
    fig.add_subplot(gs[0, 1])
    plot_tracks(n_BLT0.d_runs[r_L, 0., True], file_name=False, fig=fig, **kwargs)
    fig.title('SR+', y=0.98, color='red')
    fig.add_subplot(gs[1, 0])
    plot_fMLP(n_BLT0.d_runs[r_L, 0., False], fig=fig, ylog=False, yint=True, ylabel='Primary', **kwargs)
    fig.add_subplot(gs[1, 1])
    plot_fMLP(n_BLT0.d_runs[r_L, 0., True], fig=fig, ylog=False, yint=True, ylabel='Primary', **kwargs)
    fig.add_subplot(gs[2, 0])
    plot_LTB(n_BLT0.d_runs[r_L, 0., False], fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Secondary', **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_LTB(n_BLT0.d_runs[r_L, 0., True], fig=fig, ylog=False, ylim=(0., 2.), yint=True, ylabel='Secondary', **kwargs)
    if file_name:
        path = os.path.join('plots', n_BLT0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_BLT0_2(
        n_BLT0, file_name='confection_n_BLT0_2',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        r_L=1e1,
        **kwargs):
    kwargs['xlim'] = xlim
    model = n_BLT0.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(4, 1, height_ratios=[0.5, 1, 1, 0.5], axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
#    gs = matplotlib.gridspec.GridSpec(4, 1, height_ratios=[1, 2, 2, 1])
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(gs=gs, figsize=figsize, **kw0)
#    fig.add_subplot(gs.new_subplotspec((0, 0)))
    # if file_name:
    #     path = os.path.join('plots', n_BLT0.set_name, fig.style)    
    #     fig.savefig(file_name, path=path)
    # return fig
    plot_fMLP(n_BLT0.d_runs[r_L, 0., False], fig=fig, ylog=False, **kwargs)

#    fig.add_subplot(gs.new_subplotspec((1, 0), rowspan=2))
    fig.add_subplot(gs.new_subplotspec((1, 0)))
    fig.title('BLT-')
    plot_tracks(n_BLT0.d_runs[r_L, 0., False], file_name=False, fig=fig, **kwargs)

#    fig.add_subplot(gs.new_subplotspec((3, 0), rowspan=2))
    fig.add_subplot(gs.new_subplotspec((2, 0)))
    plot_tracks(n_BLT0.d_runs[r_L, 0., True], file_name=False, fig=fig, **kwargs)
    fig.title('BLT+')

#    fig.add_subplot(gs[5])
    fig.add_subplot(gs[3])
    plot_LTB(n_BLT0.d_runs[r_L, 0., True], fig=fig, ylog=False, ylim=(0., 2.), **kwargs)
    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_BLT0.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_exo(
        n_exo, file_name='confection_n_exo',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        r_L=1e1,
        **kwargs):
    kwargs['xlim'] = xlim
    model = n_exo.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(5, 2, height_ratios=[0.5, 1, 0.5, 0.5, 0.5], axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    tup = (gs.get_geometry()[0], gs.get_geometry()[1], 1)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(gs=gs, figsize=figsize, **kw0)
    plot_fMLP(n_exo.d_runs[r_L, 0.], fig=fig, ylog=False, **kwargs)
    fig.title('Direct', y=0.9)

    fig.add_subplot(gs[0, 1])
    plot_fMLP(n_exo.d_runs[r_L, 1.], fig=fig, ylog=False, **kwargs)
    fig.title('Exosome', y=0.9)

    fig.add_subplot(gs[1, 0])
    plot_tracks(n_exo.d_runs[r_L, 0.], file_name=False, fig=fig, **kwargs)

    fig.add_subplot(gs[1, 1])
    plot_tracks(n_exo.d_runs[r_L, 1.], file_name=False, fig=fig, **kwargs)

    # fig.add_subplot(gs[2, 0])
    # plot_exo(n_exo.d_runs[r_L, 0.], fig=fig, ylim=(0., 1.), **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_exo(n_exo.d_runs[r_L, 1.], fig=fig, ylim=(0., 1.), **kwargs)

    fig.add_subplot(gs[3, 0])
    plot_LTB(n_exo.d_runs[r_L, 0.], fig=fig, **kwargs)
    fig.add_subplot(gs[3, 1])
    plot_LTB(n_exo.d_runs[r_L, 1.], fig=fig, **kwargs)
    fig.add_subplot(gs[4, 0])
    plot_cos(n_exo.d_runs[r_L, 0.], fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.add_subplot(gs[4, 1])
    plot_cos(n_exo.d_runs[r_L, 1.], fig=fig, ylim=(-0.1, 1.), **kwargs)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_exo_projector(
        n_exo, file_name='confection_n_exo',
        Style=platypus.Projector,
        xlim=(1e3, 4e3),
        r_L=1e1,
        **kwargs):
    if Style is platypus.Poster:
        kw0 = {}
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
        if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
        if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    kwargs['xlim'] = xlim
    model = n_exo.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(5, 2, height_ratios=[1., 0.5, 0.5, 0.5, 0.5], axes=axes)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    tup = (gs.get_geometry()[0], gs.get_geometry()[1], 1)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(gs=gs, figsize=figsize, **kw0)
    plot_tracks(n_exo.d_runs[r_L, 0.], file_name=False, fig=fig, **kwargs)
    fig.title('Direct', y=0.98)

    fig.add_subplot(gs[0, 1])
    plot_tracks(n_exo.d_runs[r_L, 1.], file_name=False, fig=fig, **kwargs)
    fig.title('Exosome', y=0.98)

    fig.add_subplot(gs[1, 0])
    plot_fMLP(n_exo.d_runs[r_L, 0.], fig=fig, ylog=False, ylabel='Primary', **kwargs)

    fig.add_subplot(gs[1, 1])
    plot_fMLP(n_exo.d_runs[r_L, 1.], fig=fig, ylog=False, ylabel='Primary', **kwargs)

    # fig.add_subplot(gs[2, 0])
    # plot_exo(n_exo.d_runs[r_L, 0.], fig=fig, ylim=(0., 1.), **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_exo(n_exo.d_runs[r_L, 1.], fig=fig, ylim=(0., 1.), ylabel='Exosome', **kwargs)

    fig.add_subplot(gs[3, 0])
    plot_LTB(n_exo.d_runs[r_L, 0.], fig=fig, ylabel='Secondary', **kwargs)
    fig.add_subplot(gs[3, 1])
    plot_LTB(n_exo.d_runs[r_L, 1.], fig=fig, ylabel='Secondary', **kwargs)
    fig.add_subplot(gs[4, 0])
    plot_cos(n_exo.d_runs[r_L, 0.], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', **kwargs)
    fig.add_subplot(gs[4, 1])
    plot_cos(n_exo.d_runs[r_L, 1.], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', **kwargs)

    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_n_exo_vary_r_L(
        n_exo, file_name='confection_n_exo_vary_r_L',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        L_r_L=[1e0, 1e1, 1e2], phi_E=0.,
        **kwargs):
    kwargs['xlim'] = xlim
    model = n_exo.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(3, len(L_r_L), height_ratios=[0.5, 0.5, 0.5],
                                 axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))

    fig = Style(gs=gs, figsize=figsize, **kw0)
    for (i, r_L) in enumerate(L_r_L):
        if i:
            fig.add_subplot(gs[0, i])
        plot_fMLP(n_exo.d_runs[r_L, phi_E], fig=fig, ylog=False, **kwargs)
        fig.title(r'$r_L={:.2e}$'.format(r_L), y=0.9)

    
    for (i, r_L) in enumerate(L_r_L):
        fig.add_subplot(gs[1, i])
        plot_LTB(n_exo.d_runs[r_L, phi_E], fig=fig, **kwargs)

    for (i, r_L) in enumerate(L_r_L):
        fig.add_subplot(gs[2, i])
        cell = n_exo.d_runs[r_L, phi_E].L_cell_group[0].L_cell[0]
        f_kappa = lambda L_condition: cell.kappa(
            [(0., 0.), (0., 0.), L_condition[2]])
        plot_cos(n_exo.d_runs[r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.x_color_f2('Reds'),
                 f_kappa=f_kappa,
                 **kwargs)
        plot_cos(n_exo.d_runs[r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.blues_color_f2,
                 **kwargs)
        f_kappa = lambda L_condition: cell.kappa(
            [L_condition[0], (0., 0.), (0., 0.)])
        plot_cos(n_exo.d_runs[r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.greys_color_f,
                 f_kappa=f_kappa,
                 **kwargs)
        fig.multi_plot([[0., 7.]], [[0.5, 0.5]], L_linestyle=['--'],
                       color_f=platypus.color_f_black)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_vary_D_L(
        vary_D_L, file_name='confection_vary_D_L',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        L_D_L=[1e4, 6e4, 5e5], r_L=1e0, phi_E=0.,
        **kwargs):
    kwargs['xlim'] = xlim
    model = vary_D_L.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(3, len(L_D_L), height_ratios=[0.5, 0.5, 0.5],
                                 axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))

    fig = Style(gs=gs, figsize=figsize, **kw0)
    for (i, D_L) in enumerate(L_D_L):
        if i:
            fig.add_subplot(gs[0, i])
        plot_fMLP(vary_D_L.d_runs[D_L, r_L, phi_E], fig=fig, ylog=False, **kwargs)
        fig.title(r'$D_L={:.2e}$'.format(D_L), y=0.9)

    
    for (i, D_L) in enumerate(L_D_L):
        fig.add_subplot(gs[1, i])
        plot_LTB(vary_D_L.d_runs[D_L, r_L, phi_E], fig=fig, **kwargs)

    for (i, D_L) in enumerate(L_D_L):
        fig.add_subplot(gs[2, i])
        cell = vary_D_L.d_runs[D_L, r_L, phi_E].L_cell_group[0].L_cell[0]
        f_kappa = lambda L_condition: cell.kappa(
            [(0., 0.), (0., 0.), L_condition[2]])
        plot_cos(vary_D_L.d_runs[D_L, r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.x_color_f2('Reds'),
                 f_kappa=f_kappa,
                 **kwargs)
        plot_cos(vary_D_L.d_runs[D_L, r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.blues_color_f2,
                 **kwargs)
        f_kappa = lambda L_condition: cell.kappa(
            [L_condition[0], (0., 0.), (0., 0.)])
        plot_cos(vary_D_L.d_runs[D_L, r_L, phi_E], fig=fig,
                 yint=False, color_f=platypus.greys_color_f,
                 f_kappa=f_kappa,
                 **kwargs)
        fig.multi_plot([[0., 7.]], [[0.5, 0.5]], L_linestyle=['--'],
                       color_f=platypus.color_f_black)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', vary_D_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_vary_gamma_L(
        vary_gamma_L, file_name='confection_vary_gamma_L',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        L_gamma_L=[1e-1, 1., 1e1], L_0=1e0, phi_E=0.,
        **kwargs):
    kwargs['xlim'] = xlim
    model = vary_gamma_L.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(3, len(L_gamma_L), height_ratios=[0.5, 0.5, 0.5],
                                 axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))

    fig = Style(gs=gs, figsize=figsize, **kw0)
    for (i, gamma_L) in enumerate(L_gamma_L):
        if i:
            fig.add_subplot(gs[0, i])
        plot_fMLP(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], fig=fig, ylog=False, **kwargs)
        fig.title(r'$\gamma_L={:.2e}$'.format(gamma_L), y=0.9)

    
    for (i, gamma_L) in enumerate(L_gamma_L):
        fig.add_subplot(gs[1, i])
        plot_LTB(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], fig=fig, **kwargs)

    for (i, gamma_L) in enumerate(L_gamma_L):
        fig.add_subplot(gs[2, i])
        cell = vary_gamma_L.d_runs[gamma_L, L_0, phi_E].L_cell_group[0].L_cell[0]
        f_kappa = lambda L_condition: cell.kappa(
            [(0., 0.), (0., 0.), L_condition[2]])
        plot_cos(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], fig=fig,
                 yint=False, color_f=platypus.x_color_f2('Reds'),
                 f_kappa=f_kappa,
                 **kwargs)
        plot_cos(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], fig=fig,
                 yint=False, color_f=platypus.blues_color_f2,
                 **kwargs)
        f_kappa = lambda L_condition: cell.kappa(
            [L_condition[0], (0., 0.), (0., 0.)])
        plot_cos(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], fig=fig,
                 yint=False, color_f=platypus.greys_color_f,
                 f_kappa=f_kappa,
                 **kwargs)
        fig.multi_plot([[0., 7.]], [[0.5, 0.5]], L_linestyle=['--'],
                       color_f=platypus.color_f_black)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', vary_gamma_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def confection_vary_gamma_L_2(
        vary_gamma_L_2, file_name='confection_vary_gamma_L_2',
        Style=platypus.MBOC,
        xlim=(1e3, 4e3),
        L_gamma_L=[1e-1, 2/3, 1e1], L_00=1e0, phi_E=0.,
        **kwargs):
    kwargs['xlim'] = xlim
    model = vary_gamma_L_2.L_runs[0]
    axes = Style.def_axes
    gs = platypus.make_grid_spec(3, len(L_gamma_L), height_ratios=[0.5, 0.5, 0.5],
                                 axes=axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    if Style in [platypus.Poster, platypus.MBOC]:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))

    fig = Style(gs=gs, figsize=figsize, **kw0)
    for (i, gamma_L) in enumerate(L_gamma_L):
        if i:
            fig.add_subplot(gs[0, i])
        plot_fMLP(vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E], fig=fig, ylog=False, **kwargs)
        fig.title(r'$\gamma_L={:.2e}$'.format(gamma_L), y=0.9)

    
    for (i, gamma_L) in enumerate(L_gamma_L):
        fig.add_subplot(gs[1, i])
        plot_LTB(vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E], fig=fig, **kwargs)

    for (i, gamma_L) in enumerate(L_gamma_L):
        fig.add_subplot(gs[2, i])
        cell = vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E].L_cell_group[0].L_cell[0]
        f_kappa = lambda L_condition: cell.kappa(
            [(0., 0.), (0., 0.), L_condition[2]])
        plot_cos(vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E], fig=fig,
                 yint=False, color_f=platypus.x_color_f2('Reds'),
                 f_kappa=f_kappa,
                 **kwargs)
        plot_cos(vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E], fig=fig,
                 yint=False, color_f=platypus.blues_color_f2,
                 **kwargs)
        f_kappa = lambda L_condition: cell.kappa(
            [L_condition[0], (0., 0.), (0., 0.)])
        plot_cos(vary_gamma_L_2.d_runs[gamma_L, L_00, phi_E], fig=fig,
                 yint=False, color_f=platypus.greys_color_f,
                 f_kappa=f_kappa,
                 **kwargs)
        fig.multi_plot([[0., 7.]], [[0.5, 0.5]], L_linestyle=['--'],
                       color_f=platypus.color_f_black)

    fig.set_AB_labels()
    if file_name:
        path = os.path.join('plots', vary_gamma_L_2.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig


def set_init_rec(init_rec, x_0=3600.):
    for run in init_rec.L_runs:
        for pde in run.pde_stepper.L_pde:
            pde.x -= x_0

def confection_init_rec(init_rec, r_L=1e2, Style=platypus.MBOC,
                        file_name='confection_init_rec'):
    if Style is platypus.Poster:
        kw0 = dict(panesize=(7., 3.5))
    else:
        kw0 = dict(panesize=(3., 1.5))
    gs = platypus.make_grid_spec(3, 2, height_ratios=[0.5, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, Style.def_panesize)
    fig = Style(gs=gs, figsize=figsize, **kw0)

    path = os.path.join('plots', init_rec.set_name, fig.style)
    kwargs = {}
    kwargs['xlim'] = (600., 3600.)
    plot_fMLP(init_rec.d_runs[(r_L, 0.)], fig=fig, yint=True, **kwargs)
    fig.title('Direct', y=0.96)
    fig.add_subplot(gs[0, 1])
    fig.title('Exosome', y=0.96)
    plot_fMLP(init_rec.d_runs[(r_L, 1.)], fig=fig, yint=True, **kwargs)
    fig.add_subplot(gs[1, 0])
    plot_LTB(init_rec.d_runs[(r_L, 0.)], fig=fig, ylog=False, yint=True, **kwargs)
    fig.add_subplot(gs[1, 1])
    plot_LTB(init_rec.d_runs[(r_L, 1.)], fig=fig, ylog=False, yint=True, **kwargs)
    fig.add_subplot(gs[2, 0])
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), n=1, **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), n=1, **kwargs)
    fig.set_AB_labels()
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def confection_init_rec_projector(init_rec, r_L=1e2, Style=platypus.Projector,
                        file_name='confection_init_rec', **kwargs):
    print(Style)
    print(Style is platypus.Projector)
    if Style is platypus.Poster:
        kw0 = dict(panesize=(7., 3.5))
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
        if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
        if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    else:
        kw0 = dict(panesize=(3., 1.5))
    gs = platypus.make_grid_spec(3, 2, height_ratios=[0.5, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    fig = Style(gs=gs, figsize=figsize, **kw0)

    path = os.path.join('plots', init_rec.set_name, fig.style)

    kwargs['xlim'] = (-3e3, 0.)
    plot_fMLP(init_rec.d_runs[(r_L, 0.)], fig=fig, yint=True, ylabel='Primary', **kwargs)
    fig.title('Direct', y=0.96)
    fig.add_subplot(gs[0, 1])
    fig.title('Exosome', y=0.96)
    plot_fMLP(init_rec.d_runs[(r_L, 1.)], fig=fig, yint=True, ylabel='Primary', **kwargs)
    fig.add_subplot(gs[1, 0])
    plot_LTB(init_rec.d_runs[(r_L, 0.)], fig=fig, ylog=False, yint=True, ylabel='Secondary', **kwargs)
    fig.add_subplot(gs[1, 1])
    plot_LTB(init_rec.d_runs[(r_L, 1.)], fig=fig, ylog=False, yint=True, ylabel='Secondary', **kwargs)
    fig.add_subplot(gs[2, 0])
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 0.)], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', n=1, **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', **kwargs)
    plot_cos(init_rec.d_runs[(r_L, 1.)], fig=fig, ylim=(-0.1, 1.), ylabel='Directionality', n=1, **kwargs)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def confection_decay_F(decay_F, Style=platypus.MBOC,
                       file_name='confection_decay_F', r_L=1e1,
                       xlim=(2e3, 4e3),
                       **kwargs):
    if Style is platypus.Poster:
        kw0 = {}
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    gs = platypus.make_grid_spec(5, 2, height_ratios=[1, 0.5, 0.5, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    fig = Style(gs=gs, figsize=figsize, **kw0)
    path = os.path.join('plots', decay_F.set_name, fig.style)
    kwargs['xlim'] = xlim
    figx = plot_tracks(decay_F.d_runs[r_L, 0.], file_name=None, fig=fig, **kwargs)
    fig.title('Direct', y=0.98)

    fig.add_subplot(gs[0, 1])
    figx = plot_tracks(decay_F.d_runs[r_L, 1.], file_name=None, fig=fig, **kwargs)
    fig.title('Exosome', y=0.98)
    fig.add_subplot(gs[1, 0])
    ax = fig.fig.gca()
#    ax.yaxis.get_major_formatter().set_scientific(True)
#    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plot_fMLP(decay_F.d_runs[r_L, 0.], fig=fig, **kwargs)
#    fig.fig.canvas.draw()
    fig.add_subplot(gs[1, 1])
    plot_fMLP(decay_F.d_runs[r_L, 1.], fig=fig, **kwargs)
    # fig.add_subplot(6, 2, 7)
    # plot_CI_x(decay_F.d_runs[r_L, 0.], fig=fig, xlim=xlim)
    # plot_CI_x(decay_F.d_runs[r_L, 1.], fig=fig, color_f=lambda x: platypus.colo
    #              r_f_color(x + 1), xlim=xlim)
    # plot_exo(decay_F.d_runs[r_L, 0.], fig=fig, **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_exo(decay_F.d_runs[r_L, 1.], fig=fig, ylim=(0., 1.), **kwargs)
    fig.add_subplot(gs[3, 0])
    plot_LTB(decay_F.d_runs[r_L, 0.], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(gs[3, 1])
    plot_LTB(decay_F.d_runs[r_L, 1.], fig=fig, ylog=False, **kwargs)
    fig.add_subplot(gs[4, 0])
    plot_cos(decay_F.d_runs[r_L, 0.], fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.add_subplot(gs[4, 1])
    plot_cos(decay_F.d_runs[r_L, 1.], fig=fig, ylim=(-0.1, 1.), **kwargs)
    fig.set_AB_labels()
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def confection_decay_F_projector(decay_F, Style=platypus.Projector,
                                 file_name='confection_decay_F', r_L=1e1,
                                 xlim=(2e3, 4e3),
                                 **kwargs):
    if Style is platypus.Poster:
        kw0 = {}
    elif Style is platypus.Projector:
        kw0 = dict(panesize=(3., 3.))
        if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
        if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    gs = platypus.make_grid_spec(5, 2, height_ratios=[1, 0.5, 0.5, 0.5, 0.5, 0.5], axes=Style.def_axes)
    figsize = platypus.make_figsize(gs, kw0['panesize'])
    fig = Style(gs=gs, figsize=figsize, **kw0)
    path = os.path.join('plots', decay_F.set_name, fig.style)
    kwargs['xlim'] = xlim
    figx = plot_tracks(decay_F.d_runs[r_L, 0.], file_name=None, fig=fig, **kwargs)
#    fig.title('Direct', y=0.98)

    fig.add_subplot(gs[0, 1])
    figx = plot_tracks(decay_F.d_runs[r_L, 1.], file_name=None, fig=fig, **kwargs)
#    fig.title('Exosome', y=0.98)
    fig.add_subplot(gs[1, 0])
    ax = fig.fig.gca()
#    ax.yaxis.get_major_formatter().set_scientific(True)
#    ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
    plot_fMLP(decay_F.d_runs[r_L, 0.], fig=fig, ylabel='Primary', **kwargs)
#    fig.fig.canvas.draw()
    fig.add_subplot(gs[1, 1])
    plot_fMLP(decay_F.d_runs[r_L, 1.], fig=fig, ylabel='Primary', **kwargs)
    # fig.add_subplot(6, 2, 7)
    # plot_CI_x(decay_F.d_runs[r_L, 0.], fig=fig, xlim=xlim)
    # plot_CI_x(decay_F.d_runs[r_L, 1.], fig=fig, color_f=lambda x: platypus.colo
    #              r_f_color(x + 1), xlim=xlim)
    # plot_exo(decay_F.d_runs[r_L, 0.], fig=fig, **kwargs)
    fig.add_subplot(gs[2, 1])
    plot_exo(decay_F.d_runs[r_L, 1.], fig=fig, ylim=(0., 1.),
             ylabel='Exosome',
             **kwargs)
    fig.add_subplot(gs[3, 0])
    plot_LTB(decay_F.d_runs[r_L, 0.], fig=fig, ylabel='Secondary', ylog=False, **kwargs)
    fig.add_subplot(gs[3, 1])
    plot_LTB(decay_F.d_runs[r_L, 1.], fig=fig, ylabel='Secondary', ylog=False, **kwargs)
    fig.add_subplot(gs[4, 0])
    plot_cos(decay_F.d_runs[r_L, 0.], fig=fig, ylabel='Directionality', ylim=(-0.1, 1.), **kwargs)
    fig.add_subplot(gs[4, 1])
    plot_cos(decay_F.d_runs[r_L, 1.], fig=fig, ylabel='Directionality', ylim=(-0.1, 1.), **kwargs)
    if file_name:
        fig.savefig(file_name, path=path)
    return fig

def get_range_time_averaged(model, CI_threshold):
    return np.mean([get_range(model, CI_threshold, j) for j in range(model.pde_stepper.Nt)])

def get_range(model, CI_threshold, j, f_kappa=None):
    a_kappa = get_a_kappa(model, j, f_kappa=f_kappa)
    a_cos = flux.mean_velocity_fast(a_kappa)
    a_continuous = get_continuous(a_cos > CI_threshold)
    return np.sum(a_cos > CI_threshold) * model.pde_stepper.L_pde[0].dx
#    return np.sum(a_kappa > 0.4) * model.pde_stepper.L_pde[0].dx

def get_range_continuous(model, CI_threshold, j, f_kappa=None):
    a_kappa = get_a_kappa(model, j, f_kappa=f_kappa)
    a_cos = flux.mean_velocity_fast(a_kappa)
    a_continuous = get_continuous(a_cos > CI_threshold)
    print(a_continuous)
    return (a_continuous[-1, 1] - a_continuous[-1, 0]) * model.pde_stepper.L_pde[0].dx
#    return np.sum(a_kappa > 0.4) * model.pde_stepper.L_pde[0].dx

def get_range_max_continuous(model, CI_threshold, j, f_kappa=None):
    a_kappa = get_a_kappa(model, j, f_kappa=f_kappa)
    a_cos = flux.mean_velocity_fast(a_kappa)
    a_continuous = get_continuous(a_cos > CI_threshold)
#    print(a_continuous)
    return np.max(a_continuous[:, 1] - a_continuous[:, 0]) * model.pde_stepper.L_pde[0].dx
#    return np.sum(a_kappa > 0.4) * model.pde_stepper.L_pde[0].dx

def get_continuous(a):
    L_start = []
    L_end = []
    if a[0]:
        L_start.append(0)
    for i in range(1, np.size(a)):
        if a[i] and not a[i - 1]:
            L_start.append(i)
        if not a[i] and a[i - 1]:
            L_end.append(i)
    if a[-1]:
        L_end.append(np.size(a))
    return np.array((L_start, L_end)).T

def get_a_kappa(model, j, f_kappa=None):
    a_L_condition = np.array([(pde.u[j, 1:-1], springbok.tiger.opdudx(pde.u[j], pde.dx))
                              for pde in model.pde_stepper.L_pde])
    cell = model.L_cell_group[0].L_cell[0]
    if f_kappa is None:
        f_kappa = lambda L_condition: cell.kappa(L_condition)
    return np.array([f_kappa(a_L_condition[:, :, i]) for i in range(np.size(a_L_condition, 2))])

def table_about_range(n_exo, f_kappa=None):
    for F_xt in [1e-2, 3e-2, 1e-1, 3e-1, 1e0, 3e0, 1e1, 1e2, 1e3, 1e4, 1e5]:
        L_backwards = []
        L_range = []
        L_range_continuous = []
        for run in n_exo.L_runs:
            run.L_cell_group[0].L_cell[0].F_xt = F_xt
            a_kappa = get_a_kappa(run, 60, f_kappa=f_kappa)
            a_cos = flux.mean_velocity_fast(a_kappa)
            L_backwards.append(np.sum(a_cos < 0.))
            L_range.append(get_range(run, 0.5, 60))
            L_range_continuous.append(get_range_continuous(run, 0.5, 60))
        print(F_xt, max(L_backwards), max(L_range), max(L_range_continuous))


def plot_range_vary_phi_E(n_exo, fig=None, Style=platypus.MBOC,
                          file_name='vary-phi_E', use_mm=True):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    L_r_L = [n_exo.L_r_L[0]]
    L_r_L.extend(n_exo.L_r_L[13:27:3])
    figx = platypus.multi_plot(
        [n_exo.L_phi_E] * len(L_r_L),
        [[get_range_continuous(n_exo.d_runs[r_L, phi_E], 0.5, 60) / ratio
          for phi_E in n_exo.L_phi_E] for r_L in L_r_L],
        #        L_legend=[r'$r_L=10^{{{}}}$'.format(int(np.log10(r_L))) if r_L != 0. else r'$r_L=0$' for r_L in L_r_L],
        L_legend=[r'$r_L={:.2e}'.format(r_L) + 'K_d/\mathrm{min}$' if r_L != 0. else r'$r_L=0$' for r_L in L_r_L],
        xlabel=r'Fraction of LTB$_4$ secreted via exosomes, $\phi_E$',
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx


def make_n_exo_semi_analytic(n_exo):
    for phi_E in n_exo.L_phi_E:
        for r_L in n_exo.L_r_L:
            if r_L != 1.:
                n_exo.d_runs[r_L, phi_E] = copy.deepcopy(
                    n_exo.d_runs[1., phi_E])
                n_exo.d_runs[r_L, phi_E].pde_stepper.L_pde[2].u *= r_L

def plot_range_vary_r_L(n_exo, fig=None, Style=platypus.MBOC,
                        file_name='vary-r_L', use_mm=True, **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [n_exo.L_r_L] * len(n_exo.L_phi_E),
        [[get_range_continuous(n_exo.d_runs[r_L, phi_E], 0.5, 60) / ratio
          for r_L in n_exo.L_r_L] for phi_E in n_exo.L_phi_E],
        xlog=True,
        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in n_exo.L_phi_E],
        xlabel='Characteristic LTB$_4$ secretion rate, $r_L$, $K_d/\mathrm{min}$',
        xlim=(1e-2, 1e5),
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig, **kwargs)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range(n_exo, Style=platypus.MBOC, file_name='range'):
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3))
    fig = Style(subplot=(2, 2, 1), **kw0)
    figx = plot_range_vary_r_L(n_exo, fig=fig, file_name=None)
    fig.add_subplot(2, 2, 3)
    figx = plot_range_vary_phi_E(n_exo, fig=fig, file_name=None)
    fig.set_AB_labels(loc='upper right')
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig

def plot_range_2(n_exo, vary_gamma_L, vary_gamma_L_2, Style=platypus.MBOC, file_name='range_2'):
    if Style is platypus.Poster:
        kw0 = {}
    else:
        kw0 = dict(panesize=(3.3, 3.3), xlabelpad=2)
    fig = Style(subplot=(2, 2, 1), **kw0)
    figx = plot_range_vary_r_L(n_exo, fig=fig, file_name=None)
    fig.add_subplot(2, 2, 3)
    figx = plot_range_vary_gamma_L(vary_gamma_L, fig=fig, L_0=1e0, file_name=None)
    fig.add_subplot(2, 2, 4)
    figx = plot_range_vary_gamma_L_2(vary_gamma_L_2, fig=fig, file_name=None)
    fig.set_AB_labels(loc='upper right')
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return fig

def rr_vs_t(n_exo):
    fig = platypus.multi_plot([np.arange(1, 61)]*5, [np.array([get_range_continuous(run, 0.5, i) for i in range(1, 61)]) for run in [n_exo.d_runs[1.5, phi_E] for phi_E in n_exo.L_phi_E]], file_name='rr-vs-t-r_L1.5', path='plots/n_exo-2016-05-21/print/', xlabel='Time, min', ylabel=r'Recruitment range, $\mu m$')
    return fig

def plot_range_vary_D_L(vary_D_L, fig=None, Style=platypus.MBOC,
                          file_name='vary-D_L', r_L=1e6, use_mm=True):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [vary_D_L.L_D_L] * len(vary_D_L.L_r_L),
        [[get_range_continuous(vary_D_L.d_runs[D_L, r_L, phi_E], 0.5, 60) / ratio
          for D_L in vary_D_L.L_D_L] for phi_E in vary_D_L.L_phi_E],
        xlog=True,
        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in vary_D_L.L_phi_E],
        xlabel=r'LTB$_4$ diffusion coefficient, $D_L$, $\mathrm{\mu m}^2/s$',
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig)
    if file_name:
        path = os.path.join('plots', vary_D_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_gamma_L(vary_gamma_L, fig=None, Style=platypus.MBOC,
                            file_name='vary-ell_L-gamma_L', L_0=1e0, use_mm=True):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [vary_gamma_L.L_gamma_L] * len(vary_gamma_L.L_phi_E),
        [[get_range_continuous(vary_gamma_L.d_runs[gamma_L, L_0, phi_E], 0.5, 60) / ratio
          for gamma_L in vary_gamma_L.L_gamma_L] for phi_E in vary_gamma_L.L_phi_E],
        xlog=True,
#        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in vary_gamma_L.L_phi_E],
        xlabel=r'LTB$_4$ dissipation rate, $\gamma_L$, 1/min',
        xlim=(1e-2, 1e3),
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig)
    if file_name:
        path = os.path.join('plots', vary_gamma_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_gamma_L_2(vary_gamma_L_2, fig=None, Style=platypus.MBOC,
                              file_name='vary-ell_L-gamma_L_2', L_0=1e0, use_mm=True):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [vary_gamma_L_2.L_gamma_L] * len(vary_gamma_L_2.L_phi_E),
        [[get_range_continuous(vary_gamma_L_2.d_runs[gamma_L, L_0, phi_E], 0.5, 60) / ratio
          for gamma_L in vary_gamma_L_2.L_gamma_L] for phi_E in vary_gamma_L_2.L_phi_E],
        xlog=True,
#        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in vary_gamma_L_2.L_phi_E],
        xlabel=r'LTB$_4$ dissipation rate, $\gamma_L$, 1/min',
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig)
    if file_name:
        path = os.path.join('plots', vary_gamma_L_2.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_ell_L_D_L(vary_D_L, fig=None, Style=platypus.MBOC,
                          file_name='vary-D_L', r_L=1e6):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [vary_D_L.L_ell_L] * len(vary_D_L.L_r_L),
        [[get_range_continuous(vary_D_L.d_runs[D_L, r_L, phi_E], 0.5, 60)
          for D_L in vary_D_L.L_D_L] for phi_E in vary_D_L.L_phi_E],
        xlog=True,
        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in vary_D_L.L_phi_E],
        xlabel=r'LTB$_4$ characteristic length, $\ell_L$, $\mathrm{\mu m}$',
        ylabel='Recruitment range, $\mu m$',
        ylim=(0., 2000.),
        fig=fig)
    if file_name:
        path = os.path.join('plots', vary_D_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_ell_L_gamma_L(vary_gamma_L, fig=None, Style=platypus.MBOC,
                          file_name='vary-ell_L-gamma_L', r_L=1e6):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [vary_gamma_L.L_ell_L] * len(vary_gamma_L.L_r_L),
        [[get_range_continuous(vary_gamma_L.d_runs[gamma_L, r_L, phi_E], 0.5, 60)
          for gamma_L in vary_gamma_L.L_gamma_L] for phi_E in vary_gamma_L.L_phi_E],
        xlog=True,
        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in vary_gamma_L.L_phi_E],
        xlabel=r'LTB$_4$ characteristic length, $\ell_L$, $\mathrm{\mu m}$',
        ylabel='Recruitment range, $\mu m$',
        ylim=(0., 2000.),
        fig=fig)
    if file_name:
        path = os.path.join('plots', vary_gamma_L.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx


def plot_range_vary_gamma_E(decay_F_vary_gamma_E, fig=None, Style=platypus.MBOC,
                          file_name='decay_F_vary_gamma_E', phi_E=1.):
    if fig is None:
        fig = Style(subplot=(1, 2, 1))
    figx = platypus.multi_plot(
        [decay_F_vary_gamma_E.L_gamma_E] * len(decay_F_vary_gamma_E.L_r_L),
        [[get_range_max_continuous(decay_F_vary_gamma_E.d_runs[r_L, phi_E, gamma_E], 0.5, 60)
          for gamma_E in decay_F_vary_gamma_E.L_gamma_E] for r_L in decay_F_vary_gamma_E.L_r_L],
        xlog=True,
        L_legend=[r'$r_L={:.2e}$'.format(r_L) for r_L in decay_F_vary_gamma_E.L_r_L],
        xlabel=r'Exosome decay rate $\gamma_E$, 1/min',
        ylabel='Recruitment range, $\mu m$',
        ylim=(0., 2000.),
        fig=fig)
    if file_name:
        path = os.path.join('plots', decay_F_vary_gamma_E.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx


def plot_range_direct(n_exo, fig=None, Style=platypus.MBOC,
                      file_name='range-direct', use_mm=True, yint=True,
                      panesize=None, axes=None, **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(panesize=panesize, axes=None)
    figx = platypus.multi_plot(
        [n_exo.L_r_L] * len(n_exo.L_phi_E),
        [[get_range_continuous(n_exo.d_runs[r_L, 0.], 0.5, 60) / ratio
          for r_L in n_exo.L_r_L]],
        xlog=True,
#        L_legend=[r'$\phi_E={}$'.format(phi_E) for phi_E in n_exo.L_phi_E],
        xlabel='Secondary secretion rate, $K_d/\mathrm{min}$',
        xlim=(1e-2, 1e5), yint=yint,
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig, **kwargs)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx

def plot_range_vary_r_L_projector(
        n_exo, fig=None, Style=platypus.MBOC,
        file_name='vary-r_L', use_mm=True, yint=True,
        panesize=None, axes=None, **kwargs):
    if use_mm:
        ratio = 1e3
        units = 'mm'
    else:
        ratio = 1
        units = r'$\mu m$'
    if fig is None:
        fig = Style(panesize=panesize, axes=None)
    figx = platypus.multi_plot(
        [n_exo.L_r_L] * len(n_exo.L_phi_E),
        [[get_range_continuous(n_exo.d_runs[r_L, phi_E], 0.5, 60) / ratio
          for r_L in n_exo.L_r_L] for phi_E in n_exo.L_phi_E],
        xlog=True, yint=yint,
        xlabel='Secondary secretion rate, $K_d/\mathrm{min}$',
        xlim=(1e-2, 1e5),
        ylabel='Recruitment range, ' + units,
        ylim=(0., 2000. / ratio),
        fig=fig, **kwargs)
    if file_name:
        path = os.path.join('plots', n_exo.set_name, fig.style)    
        fig.savefig(file_name, path=path)
    return figx


def plot_all_projector(n_FPR0, n_BLT0, n_exo, init_rec, decay_F, **kwargs):
    if 'linewidth' not in kwargs: kwargs['linewidth'] = 2
    if 'markeredgewidth' not in kwargs: kwargs['markeredgewidth'] = 1.5
    if 'Style' not in kwargs: kwargs['Style'] = platypus.Projector
    fig = confection_n_FPR0_projector(n_FPR0, **kwargs)
    fig = confection_n_BLT0_projector(n_BLT0, **kwargs)
    fig = confection_init_rec_projector(init_rec, **kwargs)
    fig = confection_decay_F_projector(decay_F, **kwargs)
    # fig = plot_range_direct(n_exo, Style=Style, panesize=(6., 6.),
    #                         axes=[0.15, 0.15, 0.7, 0.7], **kwargs)
    # fig = plot_range_vary_r_L_projector(
    #     n_exo, Style=Style, panesize=(6., 6.),
    #     axes=[0.15, 0.15, 0.7, 0.7], **kwargs)
    return kwargs
