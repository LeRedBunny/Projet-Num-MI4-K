'''Approximates the wavefunction using Schrödinger's equation.'''

from PaquetOndeGauss1d4k import GaussWP
from numpy import empty, linspace


if __name__ == '__main__' :

    # Parameters
    k0 = 0
    a = 1

    # Space
    x_min = 0
    x_max = 10
    x_splits = 50
    x_int = linspace(x_min, x_max, x_splits)

    # Time
    t_min = 0
    t_max = 10
    t_splits = 50
    t_int = linspace(t_min, t_max, t_splits)

    # Wavefunction
    wavefunction = empty((x_splits, t_splits))
    wavefunction[0] = [GaussWP(k0, a, t, 0) for t in t_int]

