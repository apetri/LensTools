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
from scipy import integrate

#Visualize the masked map
def test_visualize():

	conv_map = ConvergenceMap.fromfilename("Data/unmasked.fit",loader=load_fits_default_convergence)
	mask_profile = ConvergenceMap.fromfilename("Data/mask.fit",loader=load_fits_default_convergence)

	masked_fraction = conv_map.mask(mask_profile,inplace=True)

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	ax[0].imshow(mask_profile.kappa,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0.0,mask_profile.side_angle,0.0,mask_profile.side_angle])
	ax[1].imshow(conv_map.kappa,origin="lower",interpolation="nearest",extent=[0.0,mask_profile.side_angle,0.0,mask_profile.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")

	ax[0].set_title("Mask")
	ax[1].set_title("Masked map: masking fraction {0:.2f}".format(masked_fraction))

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

#Check the mask boundaries, defined by gradients and hessians
def test_boundaries():

	conv_map = ConvergenceMap.fromfilename("Data/unmasked.fit",loader=load_fits_default_convergence)
	mask_profile = ConvergenceMap.fromfilename("Data/mask.fit",loader=load_fits_default_convergence)

	masked_map = conv_map.mask(mask_profile)
	assert hasattr(masked_map,"_mask")
	assert masked_map._masked
	assert masked_map.side_angle == conv_map.side_angle
	assert masked_map.kappa.shape == conv_map.kappa.shape

	#Compute boundaries
	perimeter_area = masked_map.maskBoundaries()

	fig,ax = plt.subplots(1,3,figsize=(24,8))

	#Plot gradient boundary
	ax[0].imshow(masked_map._gradient_boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle,0,conv_map.side_angle])

	#Plot hessian (but not gradient) boundary
	ax[1].imshow(masked_map._gradient_boundary ^ masked_map._hessian_boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle,0,conv_map.side_angle])

	#Plot gradient and hessian boundary
	ax[2].imshow(masked_map.boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle,0,conv_map.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[0].set_title("Gradient boundary")

	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")
	ax[1].set_title("Hessian overhead")

	ax[2].set_xlabel(r"$x$(deg)")
	ax[2].set_ylabel(r"$y$(deg)")
	ax[2].set_title("Full boundary: perimeter/area={0:.3f}".format(perimeter_area))

	fig.tight_layout()

	plt.savefig("boundaries.png")
	plt.clf()

#Check the differences in peak counts with and without masking
def test_peaks():

	th_peaks = np.ogrid[-0.04:0.12:50j]

	conv_map = ConvergenceMap.fromfilename("Data/unmasked.fit",loader=load_fits_default_convergence)
	mask_profile = ConvergenceMap.fromfilename("Data/mask.fit",loader=load_fits_default_convergence)

	masked_map = conv_map.mask(mask_profile)

	v,pk_orig = conv_map.peakCount(th_peaks)
	v,pk_masked = masked_map.peakCount(th_peaks)

	#Plot the difference
	plt.plot(v,pk_orig,label=r"Unmasked: $N_p=${0}".format(int(integrate.simps(pk_orig,x=v))))
	plt.plot(v,pk_masked,label=r"With {0:.1f}% area masking: $N_p=${1}".format(masked_map.maskedFraction * 100,int(integrate.simps(pk_masked,x=v))))

	#Labels
	plt.xlabel(r"$\kappa$")
	plt.ylabel(r"$dN/d\kappa$")
	plt.legend(loc="upper left")

	plt.savefig("masked_peaks.png")
