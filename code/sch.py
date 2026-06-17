'''Approximates the wavefunction using Schrödinger's equation.'''

from derivative import derivative, derivative2
from customTypes import WaveFunction, Interval
from PaquetOndeGauss1d4k import GaussWP
from numpy import empty, linspace
from constants import H_BAR, ME


def ApproximateNextPosition (position: int, wavefunction: WaveFunction, t_int: Interval, potential: float) -> None :
    '''Completes the next line of the array, corresponding to the next position, using Schrödinger's equation.'''

    laplacian = derivative2(wavefunction[:, position], x_int)

    # Proof n°1
    for i in range(len(t_int) - 1) :
        dt = t_int[i + 1] - t_int[i]
        wavefunction[position + 1, i + 1] = (
            (1 - 1j * potential * dt / H_BAR) *  wavefunction[position, i] 
            + 0.5j * H_BAR * dt * laplacian[i] / ME
        )


def ApproximateWF (wavefunction: WaveFunction, x_int: Interval, t_int: Interval, potential: Interval) -> None :
    '''Completes the array representing the wavefunction using Schrödinger's equation.'''
    
    for i in range(len(x_int) - 1) :
        ApproximateNextPosition(i, wavefunction, t_int, potential[i])




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
    for i in nx :
        wavefunction[0, i] = [GaussWP(k0, a, x, 0) for x in x_int]

