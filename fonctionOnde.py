import numpy as np
from derivee import derivee, derivee2

h_bar = 1.054571818e-34




if __name__ == '__main__' :

    # Wavefunction
    nx = 10
    nt = 10
    fonction_onde = np.empty((nx, nt))
    
    # Space interval
    x_min = 0
    x_max = 10
    x_int = np.linspace(x_min, x_max)

    # Time interval
    t_min = 0
    t_max = 10
    t_int = np.linspace(t_min, t_max)

    # Potential
    v_0 = 0

    derivee(fonction_onde[0, :])
