try:
	
	from lenstools import ConvergenceMap,ShearMap
	from lenstools.defaults import load_fits_default_convergence,load_fits_default_shear

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools import ConvergenceMap,ShearMap
	from lenstools.defaults import load_fits_default_convergence,load_fits_default_shear

import numpy as np
import matplotlib.pyplot as plt

test_map_conv = ConvergenceMap.fromfilename("conv.fit",loader=load_fits_default_convergence)

def test_smooth():

	test_map_conv_smoothed = test_map_conv.smooth(1.0)

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].imshow(test_map_conv.kappa,origin="lower",interpolation="nearest",extent=[0,test_map_conv.side_angle,0,test_map_conv.side_angle])
	ax[1].imshow(test_map_conv_smoothed.kappa,origin="lower",interpolation="nearest",extent=[0,test_map_conv.side_angle,0,test_map_conv.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[0].set_title("Unsmoothed")
	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")
	ax[1].set_title(r"$1^\prime$")

	fig.tight_layout()
	plt.savefig("smooth.png")

	plt.clf()