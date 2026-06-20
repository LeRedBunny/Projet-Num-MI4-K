''''''

from crankNicolson import *





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
    t_max = 10
    nt = 1000
    t_int = np.linspace(t_min, t_max, nt)


    energy_ratio = 0.5  # This is the ratio E/V
    v0 = energy / energy_ratio  # J
    barrier_start = 5
    barrier_length = 2
    potential = np.array([v0 if barrier_start <= x <= barrier_start + barrier_length else 0 for x in x_int], dtype=complex)
    print(f'Energy of the particle E = {energy} J\nEnergy of the barrier V = {v0} J\n')

    if v0 != 0 :
        barrier_start_index = 0
        barrier_end_index = 0
        i = 0
        while barrier_start_index == 0 or barrier_end_index == 0 :
            if potential[i] == 0 and potential[i + 1] != 0 :
                barrier_start_index = i 
            if potential[i] != 0 and potential[i + 1] == 0 :
                barrier_end_index = i
            i += 1
        print(barrier_start_index, barrier_end_index)


    wavefunction = initWaveFunction(x_int, t_int, a, k_0, x_0)
    time = approximateWaveFunction(wavefunction, x_int, t_int, potential)



    choice = -1
    while choice not in (0, 1) :
        choice = input('0 : Animation\n1 : Calculs')

    # Animation

    fig, wf = plt.subplots()
    fig.suptitle(f't={t_min}s')

    x = x_int
    y1 = [abs(wavefunction[i, t_min]) ** 2 for i in range(nx)]
    y2 = [bool(potential[i]) for i in range(nx)]

    wf.set(xlim=(-5, 15), xlabel='Position (m)', ylabel='Probability Density')

    line = wf.plot(x, y1)[0]
    wf.plot(x, y2)
    points = []


    barrier_enter_time = 0
    barrier_exit_time = 0

    def update (frame) -> None :
        '''Moves forward in time'''
        global points, barrier_enter_time, barrier_exit_time

        frame = min(frame, nt - 1)
        density = probabilityDensity(wavefunction, frame, x_int)

        line.set_ydata(density)
        fig.suptitle(f'Probability Density at t={t_int[frame]}s')

        for point in points :
            point.remove()
        points = []

        maximums = getLocalMaximums(wavefunction, frame, x_int)
        for x_i in maximums :
            points.append(plt.plot(x_int[x_i], density[x_i], 'ro')[0])
        
        # current_pos = x_int[maximums[-1]]
        # if current_pos >= barrier_start and barrier_enter_time == 0 :
        #     barrier_enter_time = t_int[frame]
        # if current_pos >= barrier_start + barrier_length and barrier_exit_time == 0 :
        #     barrier_exit_time = t_int[frame]
        # print(barrier_enter_time, barrier_exit_time)
        # if barrier_enter_time and barrier_exit_time :
        #     print(f'Time to traverse barrier : {barrier_exit_time - barrier_enter_time} s')

        reflection = integrateDensity(wavefunction, frame, x_int, end=barrier_end_index)
        transmission = integrateDensity(wavefunction, frame, x_int, start=barrier_end_index)
        print(f'R = {100 * reflection} %, T = {100 * transmission} %')
        

    animation = anim.FuncAnimation(fig, update, nt, interval=30)
    #animation.save(filename="animation.gif", writer="pillow")
    plt.show()