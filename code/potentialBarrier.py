'''Simulate the interaction of a Gaussian Wave Packet with a Rectangular Potential Barrier in one dimension.'''
from SchrodingerApproximation import ApproximateWaveFunction, getMaximum
from matplotlib import pyplot as plt, animation as anim
from PaquetOndeGauss1d4k import GaussWP
from numpy import linspace, empty, abs





if __name__ == '__main__' :

    # Parameters

    k0 = 0
    a = 1


    # Space
    
    x_min = -50
    x_max = 5
    nx = 500
    x_int = linspace(x_min, x_max, nx)


    # Time

    t_min = 0
    t_max = 10
    nt = 500
    t_int = linspace(t_min, t_max, nt)


    # Potential Barrier

    value = 0 # J
    start = 5
    length = 1 # m ??
    potential = [value if start <= x <= start + length else 0 for x in x_int]


    # Wavefunction

    wavefunction = empty((nx, nt), dtype=complex)
    for i in range(nx) :
        wavefunction[i, 0] = GaussWP(k0, a, x_int[i], t_min)
    
    ApproximateWaveFunction(wavefunction, x_int, t_int, potential)
    


    # Animation Attempt

    t_i = 0

    fig, wf = plt.subplots()
    fig.suptitle(f't={t_min}s')

    x = x_int
    y = [abs(wavefunction[i, t_i]) ** 2 for i in range(nx)]
    wf.set(xlim=(-5, 10), xlabel='Position (m)')
    line = wf.plot(x, y)[0]

    def update (frame) -> None :
        '''Moves forward in time'''
        global t_i
        if t_i < nt :
            t_i += 1
            line.set_ydata([abs(wavefunction[i, t_i]) ** 2 for i in range(nx)])
            fig.suptitle(f'Probability Density at t={t_int[t_i]}s')
            print(f'maximum probability density in x = {x_int[getMaximum(wavefunction, t_i, x_int)]}')

    animation = anim.FuncAnimation(fig, update, nt)
    #animation.save(filename="animation.gif", writer="pillow")
    plt.show()



