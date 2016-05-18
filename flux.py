import functools
import matplotlib
import numpy as np
import os
import post_processing.platypus
import scipy
from scipy.interpolate import interp1d
from scipy.stats import vonmises
import string
import sys

    
def vonmisespdf(theta, kappa):
    if kappa >= 0.:
        return vonmises.pdf(theta, kappa)
    else:
        return vonmises.pdf(np.pi - theta, -kappa)


# def spherical_vonmisespdf(theta, kappa):
#     return np.exp(kappa * np.cos(theta)) * kappa / (2. * np.pi * (np.exp(kappa) - np.exp(-kappa)))



def load_Zigmond_1977_data():
    a_po3x = np.loadtxt('data/Zigmond-1977-Fig-5-3x.txt')  # percent oriented
    a_c3x = np.array([1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5])
    a_po10x = np.loadtxt('data/Zigmond-1977-Fig-5-10x.txt')
    a_po10x = a_po10x[a_po10x[:, 0].argsort()]
    a_c10x = np.array([1e-8, 1e-7, 1e-6, 3e-6, 1e-5, 1e-4])
    return [(a_c3x, a_po3x[:, 1], 3.), (a_c10x, a_po10x[:, 1], 10.)]


def calculate_Zigmond_c_dcdx(a_c, a_po, fold):
    L_g = 1000.  # 1000 um = 1 mm
    c = lambda fold, c_0, L_g: c_0 * (1. / fold + 1.) / 2.
    dcdx = lambda fold, c_0, L_g: (1. - 1. / fold) * c_0 / L_g
    return (c(fold, a_c, L_g), dcdx(fold, a_c, L_g))


def calculate_DFRO(c, dcdx, ell, K_d):
    return ell / K_d * dcdx / (c / K_d + 1.) ** 2


def runner_Zigmond():
    L_a_c_a_po = load_Zigmond_1977_data()
    L_c_dcdx = [calculate_Zigmond_c_dcdx(a_c, a_po, fold)
                for (a_c, a_po, fold) in L_a_c_a_po]
    L_DFRO = [calculate_DFRO(c, dcdx, 10., 1e-5) for (c, dcdx) in L_c_dcdx]
    L_po = [po for (x, po, x) in L_a_c_a_po]
    L_a_kappa = [np.array(list(map(orientation_to_kappa, a_po / 100.)))[:,0]
                 for a_po in L_po]
    # kappa_by_DFRO = scipy.optimize.curve_fit(
    #     lambda x, c, b: c * x + b, np.hstack(L_DFRO),
    #     np.hstack(L_a_kappa), p0=[100., 0.])[0][0]
    #todo use the following line
    #
    kappa_by_DFRO = scipy.optimize.curve_fit(
        lambda x, c, b: 100. * v_fraction_oriented(c * x + b),
        np.hstack(L_DFRO), np.hstack(L_po), p0=[100., 1.])[0][0]

    # Fraction oriented (minus 0.5) / DFRO
    # fo_by_DFRO = scipy.optimize.curve_fit(
    #     lambda x, c: c * x, np.hstack(L_DFRO), np.hstack(L_po) / 100. - 0.5,
    #     p0=[100.])[0][0]
    # this uses L_po converted from percent to fraction, minus 0.5 (random
    # orientation)
    return (L_a_c_a_po, L_c_dcdx, L_DFRO, L_a_kappa, kappa_by_DFRO)


def plot_match_Zigmond(K_d2bK_d1=1., fK_d1=1., file_name='compare-Zigmond',
                       fig=None, L_legend=True, **kwargs):
    def b(c, dcdx, ell):
        return ell * dcdx * 1 / (c + 1) ** 2
    c_max = np.logspace(-8, -4, 201)
    L_fold = [3., 10.]
    # c_min = c_max / fold
    # c = c_min + (c_max - c_min) * x / L_g
    c = lambda fold, c_0, L_g: c_0 * (1. / fold + 1.) / 2.
    dcdx = lambda fold, c_0, L_g: (1. - 1. / fold) * c_0 / L_g
    K_d1 = 1e-5 # A guess, Zigmond puts it at 1e-6 to 1e-5
    K_d2 = K_d1 * K_d2bK_d1
    ell = 10.
    L_g = 1000.
    L_c1 = [c(fold, c_max / K_d1, L_g) for fold in L_fold]
    L_c2 = [c(fold, c_max / K_d2, L_g) for fold in L_fold]
#    L_f = [f(c) for c in L_c]
    L_dcdx1 = [dcdx(fold, c_max / K_d1, L_g)
               for (fold, c) in zip(L_fold, L_c1)]
    L_dcdx2 = [dcdx(fold, c_max / K_d2, L_g)
               for (fold, c) in zip(L_fold, L_c1)]
    L_b1 = [b(c, dcdx, ell) for (c, dcdx) in zip(L_c1, L_dcdx1)]
    L_b2 = [b(c, dcdx, ell) for (c, dcdx) in zip(L_c2, L_dcdx2)]
    L_b = [fK_d1 * b1 + (1. - fK_d1) * b2 for (b1, b2) in zip(L_b1, L_b2)]
    if L_legend is True:
        L_legend=['3x', '10x']
    fig = platypus.multi_plot(
        [c_max] * len(L_fold),
        [ L_b],
        xlabel=r'Concentration, M',
        ylabel=r'DFRO',
        xlog=True,
#        ylog=True,
        L_legend=L_legend,
        file_name=file_name, path='plot', fig=fig, **kwargs)
    return (L_fold, L_b, fig)


def plot_Zigmond_1977_fig_5(L_a_c_a_po, path=None, **kwargs):
    L_a_c = [a_c for (a_c, a_po, fold) in L_a_c_a_po]
    L_po = [po for (x, po, x) in L_a_c_a_po]
    L_c_dcdx = [calculate_Zigmond_c_dcdx(a_c, a_po, fold)
                for (a_c, a_po, fold) in L_a_c_a_po]
    fig = platypus.multi_plot(
        L_a_c, L_po, xlog=True,
        L_marker=L_marker, L_markeredgewidth=L_markeredgewidth,
        L_linestyle=['--'] * 2,
        xlabel='Concentration of fMMM, M', ylabel='% Orientation',
        xlim=(1e-8, 1e-4), ylim=(50., 100.),
        **kwargs)

    K_d = 1e-5
    c_max = np.logspace(-8, -3, 201)
    L_fold = [3., 10.]
    c = lambda fold, c_0, L_g: c_0 * (1. / fold + 1.) / 2.
    dcdx = lambda fold, c_0, L_g: (1. - 1. / fold) * c_0 / L_g
    L_c = [c(fold, c_max, 1e3) for fold in L_fold]
    L_dcdx = [dcdx(fold, c_max, 1e3)
               for (fold, c) in zip(L_fold, L_c)]
    L_DFRO = [calculate_DFRO(c, dcdx, 10., 1e-5)
              for (c, dcdx) in zip(L_c, L_dcdx)]
    figx = platypus.multi_plot(
        L_c,
        [100. * v_fraction_oriented(474. * a_DFRO) for a_DFRO in L_DFRO],
        fig=fig, **kwargs)
    ax = fig.fig.gca()
    ax.text(4e-6, 96., '10x gradient', font_properties=fig.font_properties)
    ax.text(1.3e-6, 61., '3x gradient', font_properties=fig.font_properties)
#    plot_match_Zigmond(file_name=None, fig=fig)
    fig.savefig('Zigmond_1977_fig_5', path=path)


def table_neutrophil_sensitivity(L_a_c_a_po):
    L_a_c = [a_c for (a_c, a_po, fold) in L_a_c_a_po]
    L_po = [po for (x, po, x) in L_a_c_a_po]
    L_a_kappa = [np.array(list(map(orientation_to_kappa, a_po / 100.)))[:,0]
                 for a_po in L_po]    
    L_c_dcdx = [calculate_Zigmond_c_dcdx(a_c, a_po, fold)
                for (a_c, a_po, fold) in L_a_c_a_po]
    L_DFRO = [calculate_DFRO(c, dcdx, 10., 1e-5) for (c, dcdx) in L_c_dcdx]
    print(L_a_c)
    print(L_po)
    L_S = [kappa / DFRO for (kappa, DFRO) in zip(L_a_kappa, L_DFRO)]
    print(L_S)
#po, kappa, DFRO, and S = kappa / DFRO


def orientation_to_kappa(f_oriented, my_vonmisespdf=vonmisespdf):
    return scipy.optimize.fsolve(
        lambda kappa: fraction_oriented(kappa, my_vonmisespdf=my_vonmisespdf) - f_oriented, f_oriented)


def fraction_oriented(kappa, my_vonmisespdf=None):
    fo = scipy.integrate.quad(
        lambda theta: my_vonmisespdf(theta, kappa),
        -np.pi / 2, np.pi / 2)[0]
    # if kappa < 0.:
    #     fo = 1. - fo
    return fo


def v_fraction_oriented(a_kappa, my_vonmisespdf=vonmisespdf):
    return np.array([fraction_oriented(kappa, my_vonmisespdf=my_vonmisespdf)
                     for kappa in a_kappa])


def mean_velocity(kappa, my_vonmisespdf=None):
    return scipy.integrate.quad(
        lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa),
        -np.pi, np.pi)[0]

_a_kappa = np.logspace(-4, 2)
_interp_mean_v = np.array(list(map(lambda kappa: mean_velocity(kappa, my_vonmisespdf=vonmisespdf), _a_kappa)))
mean_velocity_fast = interp1d(np.hstack(([-1e8, -10000.], -_a_kappa[::-1], [0.], _a_kappa, [10000., 1e8])), np.hstack(([-1., -1.], -_interp_mean_v[::-1], [0.], _interp_mean_v, [1., 1.])))


def flux(kappa, my_vonmisespdf=None):
    if my_vonmisespdf.__name__ is not spherical_vonmisespdf.__name__:
        return scipy.integrate.quad(
            lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa),
            -np.pi / 2, np.pi / 2)[0]
    else:
        return scipy.integrate.quad(
            lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa) *
            2. * np.pi * np.sin(theta),
            0., np.pi / 2)[0]


v_flux = np.vectorize(flux)


def normalized_flux(kappa, my_vonmisespdf=None):
    if my_vonmisespdf.__name__ is not spherical_vonmisespdf.__name__:
        return flux(kappa, my_vonmisespdf=my_vonmisespdf) * np.pi
    else:
        return flux(kappa, my_vonmisespdf=my_vonmisespdf) * 4.


v_normalized_flux = np.vectorize(normalized_flux)


def flux_up(kappa, my_vonmisespdf=None):
    if my_vonmisespdf is not spherical_vonmisespdf:
        return flux(kappa, my_vonmisespdf=my_vonmisespdf) * np.pi - 1.
    else:
        return (flux(kappa, my_vonmisespdf=my_vonmisespdf) - 0.25) / 0.25


v_flux_up = np.vectorize(flux_up)


def random_flux(kappa, my_vonmisespdf=None):
    return my_vonmisespdf(np.pi, kappa) * 2.


v_random_flux = np.vectorize(random_flux)
