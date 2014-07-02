try:
	
	from lenstools.topology import ConvergenceMap

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.topology import ConvergenceMap

import numpy as np
import matplotlib.pyplot as plt

test_map = ConvergenceMap("map.fit")

#Set bin edges
l_edges = np.arange(200.0,50000.0,200.0)

def test_compute():

	#Compute
	l,Pl = test_map.powerSpectrum(l_edges)
	assert type(l)==np.ndarray
	assert type(Pl)==np.ndarray

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*Pl/(2.0*np.pi))
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	plt.savefig("power_spectrum.png")