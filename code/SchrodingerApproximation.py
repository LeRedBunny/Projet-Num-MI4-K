'''Approximates the wavefunction using Schrödinger's equation.'''

from customTypes import WaveFunction, Interval, Function
from PaquetOndeGauss1d4k import GaussWP
from derivative import derivative2
from numpy import empty, linspace
from constants import H_BAR, ME
import matplotlib.pyplot as plt



def completeNextColumn (wavefunction: WaveFunction, t_i: int, x_int: Interval, t_int: Interval, potential: Function) -> None :
    '''Completes the column of the array, corresponding to the instant in time t_i + 1.'''

    laplacian = derivative2(wavefunction[:, t_i], x_int)
    delta_t = t_int[t_i + 1] - t_int[t_i]

    for x_i in range(len(x_int) - 1) :
        # Pour la preuve, voir annexe code 1
        wavefunction[x_i + 1, t_i + 1] = (
            (1 - 1j * potential[x_i] * delta_t / H_BAR) *  wavefunction[x_i, t_i] 
            + 0.5j * H_BAR * delta_t * laplacian[x_i] / ME
        )



def ApproximateWaveFunction (wavefunction: WaveFunction, x_int: Interval, t_int: Interval, potential: Function) -> None :
    '''Completes the array representing the wavefunction using Schrödinger's equation.'''
    
    for t_index in range(len(t_int) - 1) :
        completeNextColumn(wavefunction, t_index, x_int, t_int, potential)



def getMaximum (wavefunction: WaveFunction, t_i: int, x_int: Interval) -> int :
    '''Returns the x index of the maximum probability density at time index t_i'''
    return max([i for i in range(len(x_int))], key = lambda x_i : abs(wavefunction[x_i, t_i]) ** 2)



if __name__ == '__main__' :

    # Parameters

    k0 = 0
    a = 1


    # Space
    
    x_min = -10
    x_max = 10
    nx = 100
    x_int = linspace(x_min, x_max, nx)


    # Time

    t_min = 0
    t_max = 10
    nt = 50
    t_int = linspace(t_min, t_max, nt)


    # Potential

    potential = [1 if 5 < x < 6 else 0 for x in x_int]


    # Wavefunction

    wavefunction = empty((nx, nt), dtype=complex)
    for i in range(nx) :
        wavefunction[i, 0] = GaussWP(k0, a, x_int[i], t_min)
    
    ApproximateWaveFunction(wavefunction, x_int, t_int, potential)


    # Display

    t_index = 25
    t = t_int[t_index]

    fig, (ax_real, ax_imag) = plt.subplots(2)
    fig.suptitle(f'Wavefunction Ψ at time t={t}s')

    ax_real.plot(x_int, [abs(wavefunction[i, t_index]) ** 2 for i in range(nx)])

    ax_imag.plot(x_int, [abs(GaussWP(k0, a, x_int[i], t)) ** 2  for i in range(nx)])


    plt.show()




