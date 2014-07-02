try:
	
	from lenstools.topology import ConvergenceMap

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.topology import ConvergenceMap

import numpy as np
import matplotlib.pyplot as plt


test_map = ConvergenceMap("map.fit")

def test_types():
	
	assert test_map.kappa.dtype == np.float
	assert type(test_map.side_angle) == np.float

def test_plot():
	
	side_x,side_y = test_map.kappa.shape

	fig,ax = plt.subplots()
	ax.imshow(test_map.kappa,origin="lower",interpolation="nearest",extent=[0,side_x*test_map.side_angle,0,side_y*test_map.side_angle])
	ax.set_xlabel(r"$x$(deg)")
	ax.set_ylabel(r"$y$(deg)")
	plt.savefig("test_map.png")