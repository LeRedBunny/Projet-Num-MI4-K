

from crankNicolson import *
from collections import namedtuple


SimResults = namedtuple('SimResults', (
    'k_0',
    'a',
    'x_0',
    'E',    # Energy
    'V',    # Potential
    'L',    # Barrier length
    'time',
    'R',    # Reflection probability
    'T'     # Transmission probability
))


def simulation (k_0: float, a: float, x_0: float, potential: Array, x_int: Array, t_int: Array) -> SimResults :
    '''Creates a full simulation and displays the results.'''

    wavefunction = initWaveFunction(x_int, t_int, a, k_0, x_0)
    approximateWaveFunction(wavefunction, x_int, t_int, potential)

    barrier_start_index, barrier_end_index = getBarrierIndices(x_int, potential)

    encounter_index = 0
    traversal_index = 0 # time index at which the particle crosses the barrier

    reflection_probabilities = []
    transmission_probabilities = []

    for t_i in range(len(t_int)) :

        # Before transmission
        if encounter_index == 0 or traversal_index == 0 :
            maximums = getLocalMaximums(wavefunction, t_i, x_int)
            current_pos = maximums[-1]
            if current_pos >= barrier_start_index - 5 and encounter_index == 0 :
                encounter_index = t_i
            elif current_pos >= barrier_end_index and traversal_index == 0 :
                traversal_index = t_i
                #print(f'time {traversal_index, encounter_index}')
        
        # After transmission
        else :

            reflection_probabilities.append(integrateDensity(wavefunction, t_i, x_int, end=barrier_end_index))
            transmission_probabilities.append(integrateDensity(wavefunction, t_i, x_int, start=barrier_end_index))
    

    # Display

    energy = H_BAR ** 2 * (1 + (a * k_0) ** 2) / (2 * ME * a ** 2)
    # print(f'Energy E = {energy} J\nPotential V = {potential[barrier_start_index + 1]} J')
    # print(f'Time to traverse barrier : {t_int[traversal_index] - t_int[encounter_index]} s')
    reflection = sum(reflection_probabilities) / len(reflection_probabilities)
    transmission = sum(transmission_probabilities) / len(transmission_probabilities)
    reflection /= reflection + transmission
    transmission /= reflection + transmission
    # print(f'Probability of Reflection R = {100 * reflection} %\nProbability of transmission T = {100 * transmission} %')


    results = SimResults(
        k_0,
        a,
        x_0,
        energy,
        potential[barrier_start_index + 1],
        x_int[barrier_end_index] - x_int[barrier_start_index],
        t_int[traversal_index] - t_int[encounter_index],
        reflection, 
        transmission
    )

    return results



if __name__ == '__main__' :


    # Parameters

    a = 1  # m
    k_0 = 5  # m^-1
    x_0 = 0
    energy = H_BAR ** 2 * (1 + (a * k_0) ** 2) / (2 * ME * a ** 2)   # J


    x_min = -30
    x_max = 25
    nx = 500
    x_int = np.linspace(x_min, x_max, nx)

    t_min = 0   # Don't change that
    t_max = 3
    nt = 1000
    t_int = np.linspace(t_min, t_max, nt)


    energy_ratio = 0.8  # This is the ratio E/V
    v0 = energy / energy_ratio  # J
    barrier_start_index = 350
    


    # Simulation to plot traversal time by barrier width

    # times = []

    # values = 5, 30

    # for i in range(*values) :
    #     barrier_end = x_int[barrier_start_index + i]
    #     potential = np.array([v0 if x_int[barrier_start_index] <= x <= barrier_end else 0 for x in x_int])
    #     times.append(simulation(k_0, a, x_0, potential, x_int, t_int).time)
    
    # fig, ax = plt.subplots()
    dx = x_int[1] - x_int[0]
    # fig.suptitle(f'Time to cross the barrier proportionally to barrier length')
    # ax.set(xlabel='Length (m)', ylabel='Time (s)')
    # ax.plot([i * dx for i in range(*values)], times)
    # plt.show()


    #Simulation to plot transmission probability to barrier height

    trans = []

    for i in range(1, 30) :
        potential = np.array([i * energy / 10 if x_int[barrier_start_index] <= x <= x_int[barrier_start_index + 20] else 0 for x in x_int])
        trans.append(simulation(k_0, a, x_0, potential, x_int, t_int).T)
    
    fig, ax = plt.subplots()

    dx = x_int[1] - x_int[0]
    fig.suptitle(f'Time to cross the barrier proportionally to potential')

    ax.set(xlabel='Potential / Energy', ylabel='Probability of transmission')
    ax.plot([i / 10 for i in range(1, 30)], trans)

    plt.show()