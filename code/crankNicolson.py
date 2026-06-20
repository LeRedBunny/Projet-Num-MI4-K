'''Simulates the behavior of a gaussian wave packet interacting with a potential barrier using the Crank-Nicolson method'''


        # Imports

import numpy as np
from scipy.linalg import solve_banded
from matplotlib import pyplot as plt, animation as anim
from itertools import chain
from typing import Optional


        # Constants

TWO_PI = 6.283185307179586  # rad
H_BAR = 1 # 1.054571817e-34 / TWO_PI    # Reduced Planck's constant [Js]
ME = 1 # 9.1093837015e-31       # Mass of electron [kg]


type Array = np.array 



        # Functions

def initWaveFunction (x_int: Array, t_int: Array, a: float, k_0: float, x_0: float) -> Array :
    '''Returns the array representing the wavefunction, with the column corresponding to inital position initialized'''

    wavefunction = np.zeros((nx, nt), dtype=complex)

    for x_i in range(nx) :
        x = x_int[x_i]
        # Preuve : Annexe code 2
        wavefunction[x_i, 0] = (4 / (TWO_PI * a ** 2)) ** 0.25 * np.exp(-((x - x_0) / a) ** 2 + 1j * k_0 * x)

    # normalizeWaveFunction(wavefunction, 0, x_int)
    return wavefunction


def integrateDensity (wavefunction: Array, t_i: int, x_int: Array, *, start: Optional[int] = None, end: Optional[int] = None) -> float :
    '''Integrates the probability density from start to end. If start or end are None, they will be at infinity.'''

    dx = x_int[1] - x_int[0]
    
    mask = smoothingArray(x_int)

    if start is None :
        start = 0
    if end is None :
        end = len(x_int)
    
    integral = 0
    density = probabilityDensity(wavefunction, t_i, x_int)
    for x_i in range(start, end - 1) :
        integral += (density[x_i] + density[x_i + 1]) * dx / 2 * mask[x_i]
    
    return integral


def normalizeWaveFunction (wavefunction: Array, t_i: int, x_int: Array) -> None :
    '''Normalizes the wavefunction at a given time by approximating the integral of the probability density'''
    factor = integrateDensity(wavefunction, t_i, x_int) ** -0.5
    for x_i in range(len(x_int)) :
        wavefunction[x_i, t_i] *= factor



def smoothingArray (x_int: Array, *, width: float = 3.0) -> Array :
    '''Returns an array of values to soften the impact of approximation errors on the edges of the simulation.'''

    mask = np.ones_like(x_int)
    
    left = x_int < (x_int[0] + width)
    mask[left] = np.exp(-((x_int[left] - (x_int[0] + width)) / width) ** 2)
    # Bord droit
    zone_d = x_int > (x_int[-1] - width)
    mask[zone_d] = np.exp(-((x_int[zone_d] - (x_int[-1] - width)) / width) ** 2)

    return mask



def approximateWaveFunction (wavefunction: Array, x_int: Array, t_int: Array, potential: Array) -> float :
    '''Approximates the wave function using the Crank-Nicolson method. Returns the time it took to do so.'''

    dx = x_int[1] - x_int[0]
    dt = t_int[1] - t_int[0]

    r = 0.25j * H_BAR * dt / (ME * dx ** 2)
    
    diag_A  =  1 + 2 * r + 1j * dt / (2 * H_BAR) * potential
    off_A   = -r * np.ones(nx - 1)

    diag_B  =  1 - 2 * r - 1j * dt / (2 * H_BAR) * potential
    off_B   =  r * np.ones(nx - 1)

    A_banded = np.zeros((3, nx), dtype=complex)
    A_banded[0, 1:]  = off_A   # sur-diagonale
    A_banded[1, :]   = diag_A  # diagonale
    A_banded[2, :-1] = off_A   # sous-diagonale

    mask = smoothingArray(x_int)

    for j in range(nt - 1):
        # Calcul du membre de droite (explicite) : B @ wavefunction[j]
        rhs = diag_B * wavefunction[:, j]
        rhs[1:]  += off_B * wavefunction[:-1, j]
        rhs[:-1] += off_B * wavefunction[1:, j]

        # Résolution implicite du système tridiagonal
        wavefunction[:, j+1] = solve_banded((1, 1), A_banded, rhs) * mask


def getLocalMaximums (wavefunction: Array, t_i: int, x_int: Array, *, search_width: int = 5) -> list[int] :
    '''Returns a list of local maximums of the probability density at a given time.'''

    density = probabilityDensity(wavefunction, t_i, x_int)
    local_maximums = []
    
    for x_i in range(search_width, len(x_int) - search_width) :

        if all(density[x_i] > density[x_i + i] for i in chain(range(-search_width, 0), range(1, 1 + search_width))) :
            local_maximums.append(x_i)

    return local_maximums
        

def probabilityDensity (wavefunction: Array, t_i: int, x_int: Array) -> Array :
    '''Returns the probability density function at a given time.'''
    return np.array([abs(wavefunction[x_i, t_i]) ** 2 for x_i in range(len(x_int))])




 


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


    energy_ratio = 0.8  # This is the ratio E/V
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
    approximateWaveFunction(wavefunction, x_int, t_int, potential)



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
        
        current_pos = x_int[maximums[-1]]
        if current_pos >= barrier_start - 5 and barrier_enter_time == 0 :
            barrier_enter_time = t_int[frame]
        if current_pos >= barrier_start + barrier_length and barrier_exit_time == 0 :
            barrier_exit_time = t_int[frame]
            print(f'Time to traverse barrier : {barrier_exit_time - barrier_enter_time} s')
            

        reflection = integrateDensity(wavefunction, frame, x_int, end=barrier_end_index)
        transmission = integrateDensity(wavefunction, frame, x_int, start=barrier_end_index)
        # print(f'R = {100 * reflection} %, T = {100 * transmission} %')
        

    animation = anim.FuncAnimation(fig, update, nt, interval=30)
    #animation.save(filename="animation.gif", writer="pillow")
    plt.show()


    # t_i = nt // 3

    # fig, wf = plt.subplots()
    # wf.set(xlim=(-5, x_max - 3), xlabel='Position (m)', ylabel='Probability Density')
    # x = x_int
    # y = probabilityDensity(wavefunction, t_i, x_int)
    # wf.plot(x, y)
    # print([x_int[i] for i in getLocalMaximums(wavefunction, t_i, x_int)])
    # plt.show()
