try:
	
	from lenstools import ConvergenceMap
	from lenstools.defaults import load_fits_default_convergence

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools import ConvergenceMap
	from lenstools.defaults import load_fits_default_convergence

import numpy as np
import matplotlib.pyplot as plt

#Visualize the masked map
def test_visualize():

	conv_map = ConvergenceMap.fromfilename("Data/unmasked.fit",loader=load_fits_default_convergence)
	mask_profile = ConvergenceMap.fromfilename("Data/mask.fit",loader=load_fits_default_convergence)

	conv_map.mask(mask_profile)

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].imshow(mask_profile.kappa,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0.0,mask_profile.side_angle,0.0,mask_profile.side_angle])
	ax[1].imshow(conv_map.kappa,origin="lower",interpolation="nearest",extent=[0.0,mask_profile.side_angle,0.0,mask_profile.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")

	ax[0].set_title("Mask")
	ax[1].set_title("Masked map")

	fig.tight_layout()
	plt.savefig("mask.png")

	plt.clf()

#Pad the mask with zeros and see the effect on the power spectrum
def test_power():

	conv_map = ConvergenceMap.fromfilename("Data/unmasked.fit",loader=load_fits_default_convergence)
	mask_profile = ConvergenceMap.fromfilename("Data/mask.fit",loader=load_fits_default_convergence)

	l_edges = np.arange(200.0,50000.0,200.0)
	
	fig,ax = plt.subplots()

	l,P_original = conv_map.powerSpectrum(l_edges)
	l,P_masked = (conv_map*mask_profile).powerSpectrum(l_edges)
	l,P_mask = mask_profile.powerSpectrum(l_edges)

	ax.plot(l,l*(l+1)*P_original/(2*np.pi),label="Unmasked")
	ax.plot(l,l*(l+1)*P_masked/(2*np.pi),label="Masked")
	ax.plot(l,l*(l+1)*P_mask/(2*np.pi),label="Mask")
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	ax.legend(loc="lower left")

	plt.savefig("power_mask.png")
	plt.clf()