try:
	
	from lenstools.shear import ShearMap

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.shear import ShearMap

import numpy as np
import matplotlib.pyplot as plt


test_map = ShearMap("shear.fit")

l_edges = np.arange(300.0,5000.0,200.0)

def test_visualize():

	assert hasattr(test_map,"gamma")
	assert hasattr(test_map,"side_angle")
	assert test_map.gamma.shape[0] == 2

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].imshow(test_map.gamma[0],origin="lower",interpolation="nearest",extent=[0,test_map.side_angle,0,test_map.side_angle])
	ax[1].imshow(test_map.gamma[1],origin="lower",interpolation="nearest",extent=[0,test_map.side_angle,0,test_map.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[0].set_title(r"$\gamma_1$")

	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")
	ax[1].set_title(r"$\gamma_2$")

	plt.savefig("shear.png")
	plt.clf()

def test_EB_decompose():

	l,EE,BB,EB = test_map.decompose(l_edges)

	assert l.shape == EE.shape == BB.shape == EB.shape

	fig,ax = plt.subplots()
	ax.plot(l,l*(l+1)*EE/(2.0*np.pi),label=r"$P_{EE}$")
	ax.plot(l,l*(l+1)*BB/(2.0*np.pi),label=r"$P_{BB}$")
	ax.plot(l,l*(l+1)*np.abs(EB)/(2.0*np.pi),label=r"$\vert P_{EB}\vert$")

	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	
	ax.legend(loc="Upper left")

	plt.savefig("EB.png")
