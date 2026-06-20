

from crankNicolson import *



def simulation (k_0: float, a: float, x_0: float, potential: Array, x_int: Array, t_int: Array) -> None :
    '''Creates a full simulation and displays the results'''

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
            current_pos = x_int[maximums[-1]]
            if current_pos >= barrier_start - 5 and encounter_index == 0 :
                encounter_index = t_i
            if current_pos >= barrier_start + barrier_length and traversal_index == 0 :
                traversal_index = t_i
        
        # After transmission
        else :

            reflection_probabilities.append(integrateDensity(wavefunction, t_i, x_int, end=barrier_end_index))
            transmission_probabilities.append(integrateDensity(wavefunction, t_i, x_int, start=barrier_end_index))
    

    # Display

    energy = H_BAR ** 2 * (1 + (a * k_0) ** 2) / (2 * ME * a ** 2)
    print(f'Energy E = {energy} J\nPotential V = {potential[barrier_start_index + 1]} J')
    print(f'Time to traverse barrier : {t_int[traversal_index] - t_int[encounter_index]} s')
    reflection = sum(reflection_probabilities) / len(reflection_probabilities)
    transmission = sum(transmission_probabilities) / len(transmission_probabilities)
    reflection /= reflection + transmission
    transmission /= reflection + transmission
    print(f'Probability of Reflection R = {100 * reflection} %\nProbability of transmission T = {100 * transmission} %')



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
    barrier_start = 5
    barrier_length = 2
    potential = np.array([v0 if barrier_start <= x <= barrier_start + barrier_length else 0 for x in x_int])


    # Simulation

    simulation(k_0, a, x_0, potential, x_int, t_int)