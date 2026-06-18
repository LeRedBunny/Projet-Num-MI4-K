'''Simulate the interaction of a Gaussian Wave Packet with a Rectangular Potential Barrier in one dimension.'''
from SchrodingerApproximation import ApproximateWaveFunction, getMaximum
from matplotlib import pyplot as plt, animation as anim
from PaquetOndeGauss1d4k import GaussWP
from numpy import linspace, empty, abs
from constants import H_BAR, ME





if __name__ == '__main__' :

    # Parameters

    k0 = 1 # m^-1
    a = 1 # m
    average_energy = H_BAR ** 2 * (1 + (a * k0) ** 2) / (2 * ME * a) # J


    # Space
    
    x_min = -60
    x_max = 15
    nx = 1000
    x_int = linspace(x_min, x_max, nx)


    # Time

    t_min = 0
    t_max = 10
    nt = 1000
    t_int = linspace(t_min, t_max, nt)


    # Potential Barrier

    energy_ratio = 0.5 # Chosen value of E/V
    potential_value = average_energy / energy_ratio # J
    barrier_start = 5
    barrier_length = 1
    potential = [potential_value if barrier_start <= x <= barrier_start + barrier_length else 0 for x in x_int]
    print(f'Average energy E = {average_energy}J\nBarrier of potential V = {potential_value}J\nRatio E/V = {energy_ratio}')


    # Wavefunction

    wavefunction = empty((nx, nt), dtype=complex)
    for i in range(nx) :
        wavefunction[i, 0] = GaussWP(k0, a, x_int[i], t_min)
    
    ApproximateWaveFunction(wavefunction, x_int, t_int, potential)
    

    # Other

    

    # Animation Attempt

    t_i = 0

    fig, wf = plt.subplots()
    fig.suptitle(f't={t_min}s')

    x = x_int
    y = [abs(wavefunction[i, t_i]) ** 2 for i in range(nx)]
    wf.set(xlim=(0, x_max - 3), xlabel='Position (m)')
    line = wf.plot(x, y)[0]

    def update (frame) -> None :
        '''Moves forward in time'''
        global t_i
        if t_i < nt :
            t_i += 1
            line.set_ydata([abs(wavefunction[i, t_i]) ** 2 for i in range(nx)])
            fig.suptitle(f'Probability Density at t={t_int[t_i]}s')

            max_i = getMaximum(wavefunction, t_i, x_int)
            #print(f'maximum probability density in x = {x_int[max_i]}')

    animation = anim.FuncAnimation(fig, update, nt, interval=50)
    # animation.save(filename="animation.gif", writer="pillow")
    plt.show()


    # for t_i in range(nt) :
    #     max_i = getMaximum(wavefunction, t_i, x_int)
    #     print(wavefunction[max_i + 3, t_i])



