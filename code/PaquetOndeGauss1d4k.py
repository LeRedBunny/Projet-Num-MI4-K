'''Gaussian Wave Packet.'''

from constants import H_BAR, ME, TWO_PI
from numpy import exp, linspace
import matplotlib.pyplot as plt


def GaussWP (k0: float, a: float, x: float, t: float) -> complex :
    '''Returns the value of the wavefunction of a Gaussian wave packet at position x and time t.'''
    return TWO_PI ** (-3/4) * (2 * ME * TWO_PI * a / (a ** 2 * ME + 2j * H_BAR * t)) ** 2 * exp((ME * (a ** 2 * k0 + 2j * x) ** 2 / (ME * a ** 2 + 2j * H_BAR * t) - (a * k0) ** 2) / 4)


if __name__ == '__main__' :
	
	t_0 = 1
	t = t_0
	
	k0 = 1
	a=1
	
	x_range = (0, 2)
	splits = 100
	
	x = linspace(x_range[0], x_range[1], splits)
	y = [GaussWP(k0, a, x, t).imag for x in x]
	
	fig, ax = plt.subplots()
	ax.plot(x, y)
	plt.show()
