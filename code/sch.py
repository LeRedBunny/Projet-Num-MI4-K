'''Approximates the wavefunction using Schrödinger's equation.'''

from derivative import derivative, derivative2
from types import WaveFunction, Interval
from PaquetOndeGauss1d4k import GaussWP
from numpy import empty, linspace
from constants import H_BAR




def ApproximateWF (wavefunction: WaveFunction, x_int: Interval, t_int: Interval) -> None :
    '''Completes the array representing the wavefunction using Schrödinger's equation.'''
    
    initial_wf = wavefunction[:, 0]
    laplacian = derivative2(initial_wf, x_int)

    for i in range(len(x_int)) :
        for j in range(1, len(t_int)) :

            wavefunction[i, j] = 




if __name__ == '__main__' :

    # Parameters
    k0 = 0
    a = 1

    # Space
    x_min = 0
    x_max = 10
    nx = 50
    x_int = linspace(x_min, x_max, nx)

    # Time
    t_min = 0
    t_max = 10
    nt = 50
    t_int = linspace(t_min, t_max, nt)

    # Wavefunction
    wavefunction = empty((nx, nt))
    wavefunction[:, 0] = [GaussWP(k0, a, x, 0) for x in x_int]

