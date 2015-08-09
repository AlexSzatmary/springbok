import functools
import numpy as np
import os
import scipy
from scipy.stats import vonmises
import string
import tiger_models


def vonmisespdf(theta, kappa):
    if kappa >= 0.:
        return vonmises.pdf(theta, kappa)
    else:
        return vonmises.pdf(np.pi - theta, -kappa)


def spherical_vonmisespdf(theta, kappa):
    return np.exp(kappa * np.cos(theta)) * kappa / (2. * np.pi * (np.exp(kappa) - np.exp(-kappa)))


def orientation_to_kappa(f_oriented, my_vonmisespdf=vonmisespdf):
    return scipy.optimize.fsolve(
        lambda kappa: fraction_oriented(kappa, my_vonmisespdf=my_vonmisespdf) - f_oriented, f_oriented)


def fraction_oriented(kappa, my_vonmisespdf=None):
    fo = scipy.integrate.quad(
        lambda theta: my_vonmisespdf(theta, kappa),
        -np.pi / 2, np.pi / 2)[0]
    #     fo = 1. - fo
    return fo


def v_fraction_oriented(a_kappa, my_vonmisespdf=vonmisespdf):
    return np.array([fraction_oriented(kappa, my_vonmisespdf=my_vonmisespdf)
                     for kappa in a_kappa])


def mean_velocity(kappa, my_vonmisespdf=None):
    return scipy.integrate.quad(
        lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa),
        -np.pi, np.pi)[0]
                                                   

def flux(kappa, my_vonmisespdf=None):
    if my_vonmisespdf is not spherical_vonmisespdf:
        return scipy.integrate.quad(
            lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa),
            -np.pi / 2, np.pi / 2)[0]
    else:
        return scipy.integrate.quad(
            lambda theta: np.cos(theta) * my_vonmisespdf(theta, kappa) *
            2. * np.pi * np.sin(theta),
            0., np.pi / 2)[0]


v_flux = np.vectorize(flux)


def flux_up(kappa, my_vonmisespdf=None):
    if my_vonmisespdf is not spherical_vonmisespdf:
        return flux(kappa, my_vonmisespdf=my_vonmisespdf) * np.pi - 1.
    else:
        return (flux(kappa, my_vonmisespdf=my_vonmisespdf) - 0.25) / 0.25


v_flux_up = np.vectorize(flux_up)


def random_flux(kappa, my_vonmisespdf=None):
    return my_vonmisespdf(np.pi, kappa) * 2.


v_random_flux = np.vectorize(random_flux)
