'''Plane waves.'''

from numpy import pi, exp, sqrt, real, imag, zeros, linspace
import matplotlib.pyplot as plt
from helper import splitInterval




def PlaneWave (amp: complex, k: float, omega: float, x: float, t: float) -> complex :
	'''Salut!!!!!'''
	return amp * exp(1j * (omega * t - k * x))



if __name__ == '__main__' :

	amp = 1
	
	t_0 = 1
	t = t_0
	
	k = 1
	omega = 1
	
	x_range = (0, 100)
	splits = 100
	
	x = splitInterval(x_range[0], x_range[1], splits)
	y = [PlaneWave(amp, k, omega, x, t) for x in x]
	
	fig, ax = plt.subplots()
	ax.plot(x, y)
	plt.show()