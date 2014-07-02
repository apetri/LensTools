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
thresholds = np.arange(-1.0,5.0,0.2)

def test_compute():

	#Compute
	nu,pk = test_map.peakCount(thresholds,norm=True)

	#Check if computation went OK
	assert type(nu)==np.ndarray
	assert type(pk)==np.ndarray

	#Visualize
	fig,ax = plt.subplots()
	ax.plot(nu,pk)
	ax.set_xlabel(r"$\nu=\kappa/\sigma$")
	ax.set_ylabel(r"$dN/d\nu$")

	plt.savefig("peaks.png")
