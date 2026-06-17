'''Approximates the wavefunction using Schrödinger's equation.'''

from customTypes import WaveFunction, Interval, Function
from derivative import derivative, derivative2
from PaquetOndeGauss1d4k import GaussWP
from numpy import empty, linspace
from constants import H_BAR, ME
import matplotlib.pyplot as plt


def ApproximateNextColumn (column: int, wavefunction: WaveFunction, x_int: Interval, t_int: Interval, potential: float) -> None :
    '''Completes the next column of the array, corresponding to the next instant in time, using Schrödinger's equation.'''

    laplacian = derivative2(wavefunction[:, column], x_int)

    # Proof n°1
    for i in range(len(t_int) - 1) :
        dt = t_int[i + 1] - t_int[i]
        wavefunction[i + 1, column + 1] = (
            (1 - 1j * potential * dt / H_BAR) *  wavefunction[column, i] 
            + 0.5j * H_BAR * dt * laplacian[i] / ME
        )


def ApproximateWaveFunction (wavefunction: WaveFunction, x_int: Interval, t_int: Interval, potential: Function) -> None :
    '''Completes the array representing the wavefunction using Schrödinger's equation.'''
    
    for i in range(len(t_int) - 1) :
        ApproximateNextColumn(i, wavefunction, x_int, t_int, potential[i])




if __name__ == '__main__' :

    # Parameters
    k0 = 0
    a = 1

    # Space
    x_min = -2
    x_max = 2
    nx = 100
    x_int = linspace(x_min, x_max, nx)

    # Time
    t_min = 0
    t_max = 10
    nt = 50
    t_int = linspace(t_min, t_max, nt)

    # Potential
    potential = [0] * nx

    # Wavefunction
    wavefunction = empty((nx, nt), dtype=complex)
    for i in range(nx) :
        wavefunction[i, 0] = GaussWP(k0, a, x_int[i], t_min)
    
    ApproximateWaveFunction(wavefunction, x_int, t_int, potential)


    # Display

    t_index = 1
    t = t_int[t_index]

    fig, (ax_real, ax_imag) = plt.subplots(2)
    fig.suptitle(f'Wavefunction Ψ at time t={t}s')

    ax_real.plot(x_int, [wavefunction[i, t_index] * wavefunction[i, t_index].conjugate() for i in range(nx)])

    ax_imag.plot(x_int, [GaussWP(k0, a, x_int[i], t) * GaussWP(k0, a, x_int[i], t).conjugate()  for i in range(nx)])


    plt.show()



