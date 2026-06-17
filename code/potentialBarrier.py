'''Simulate the interaction of a Gaussian Wave Packet with a Rectangular Potential Barrier in one dimension.'''
from SchrodingerApproximation import ApproximateWaveFunction
from PaquetOndeGauss1d4k import GaussWP
from numpy import linspace, empty






if __name__ == '__main__' :

    # Parameters

    k0 = 0
    a = 1


    # Space
    
    x_min = -5
    x_max = 5
    nx = 100
    x_int = linspace(x_min, x_max, nx)


    # Time

    t_min = 0
    t_max = 10
    nt = 50
    t_int = linspace(t_min, t_max, nt)


    # Potential Barrier

    V_0 = 1 # J
    L = 1 # m ??
    potential = [V_0 if 0 <= x <= L else 0 for x in x_int]


    # Wavefunction

    wavefunction = empty((nx, nt), dtype=complex)
    for i in range(nx) :
        wavefunction[i, 0] = GaussWP(k0, a, x_int[i], t_min)
    
    ApproximateWaveFunction(wavefunction, x_int, t_int, potential)


    
