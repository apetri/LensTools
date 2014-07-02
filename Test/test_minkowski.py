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
thresholds = np.arange(-2.0,2.0,0.2)

def test_compute():

	#Compute
	nu,V0,V1,V2 = test_map.minkowskiFunctionals(thresholds,norm=True)

	#Assert computation went OK
	assert hasattr(test_map,"gradient_x")
	assert hasattr(test_map,"gradient_y")
	assert hasattr(test_map,"hessian_xx")
	assert hasattr(test_map,"hessian_yy")
	assert hasattr(test_map,"hessian_xy")

	#Visualize
	fig,ax = plt.subplots(1,3,figsize=(24,8))
	ax[0].plot(nu,V0)
	ax[1].plot(nu,V1)
	ax[2].plot(nu,V2)

	ax[0].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[0].set_ylabel(r"$V_0(\nu)$")

	ax[1].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[1].set_ylabel(r"$V_1(\nu)$")

	ax[2].set_xlabel(r"$\nu=\kappa/\sigma$")
	ax[2].set_ylabel(r"$V_2(\nu)$")

	fig.tight_layout()

	plt.savefig("minkowski.png")