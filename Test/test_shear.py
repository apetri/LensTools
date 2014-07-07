try:
	
	from lenstools.shear import ShearMap
	from lenstools.topology import ConvergenceMap

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.shear import ShearMap
	from lenstools.topology import ConvergenceMap

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits

def two_file_loader(*args):

	shear_file_1 = fits.open(args[0])
	angle = shear_file_1[0].header["ANGLE"]
	gamma = shear_file_1[0].data.astype(np.float)
	shear_file_1.close()

	shear_file_2 = fits.open(args[1])
	assert shear_file_2[0].header["ANGLE"] == angle
	gamma = np.array((gamma,shear_file_2[0].data.astype(np.float)))
	shear_file_2.close()

	return angle,gamma




test_map = ShearMap.fromfilename("shear1.fit","shear2.fit",loader=two_file_loader)
test_map_conv = ConvergenceMap.fromfilename("conv.fit")

l_edges = np.arange(200.0,50000.0,200.0)

def test_visualize1():

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
	plt.clf()

def test_visualize2():

	fig,ax = plt.subplots()
	
	ax.imshow(test_map_conv.kappa,origin="lower",interpolation="nearest",cmap=plt.cm.Reds,extent=[0.0,test_map_conv.side_angle,0.0,test_map_conv.side_angle])
	test_map.sticks(ax,pixel_step=40)

	ax.set_xlabel(r"$x$(deg)")
	ax.set_ylabel(r"$y$(deg)")
	plt.savefig("sticks.png")


