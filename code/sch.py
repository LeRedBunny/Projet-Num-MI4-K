'''Approximates the wavefunction using Schrödinger's equation.'''

from derivative import derivative, derivative2
from types import WaveFunction, Interval
from PaquetOndeGauss1d4k import GaussWP
from numpy import empty, linspace
from constants import H_BAR, ME


def ApproximateInitialPosition (wavefunction: WaveFunction, t_int: Interval, potential: float) -> None :
    '''Completes the first line of the array, corresponding to initial position, using Schrödinger's equation.'''

    laplacian = derivative2(wavefunction[:, 0], x_int)

    for i in range(len(t_int) - 1) :
        dt = t_int[i + 1] - t_int[i]
        wavefunction[0, i + 1] = (
            (1 - 1j * potential * dt / H_BAR) *  wavefunction[0, i] 
            + 0.5j * H_BAR * dt * laplacian[i] / ME
        )


def ApproximateWF (wavefunction: WaveFunction, x_int: Interval, t_int: Interval) -> None :
    '''Completes the array representing the wavefunction using Schrödinger's equation.'''
    
    initial_wf = wavefunction[:, 0]
    laplacian = derivative2(initial_wf, x_int)
    
    ApproximateInitialPosition(wavefunction, t_int)

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

