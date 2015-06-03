import os

from .. import ConvergenceMap,Mask
from ..utils.defaults import load_fits_default_convergence

from .. import dataExtern

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

#Visualize the masked map
# do not edit! added by PythonBreakpoints
from pdb import set_trace as _breakpoint

def test_visualize():

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_fraction = conv_map.mask(mask_profile,inplace=True)

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	mask_profile.visualize(fig,ax[0],cmap="binary")
	conv_map.visualize(fig,ax[1])

	ax[0].set_title("Mask")
	ax[1].set_title("Masked map: masking fraction {0:.2f}".format(masked_fraction))

	fig.tight_layout()
	fig.savefig("mask.png")

#Pad the mask with zeros and see the effect on the power spectrum
def test_power():

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	l_edges = np.arange(200.0,50000.0,200.0)
	
	fig,ax = plt.subplots()

	l,P_original = conv_map.powerSpectrum(l_edges)
	l,P_masked = (conv_map*mask_profile).powerSpectrum(l_edges)
	l,P_mask = mask_profile.powerSpectrum(l_edges)

	ax.plot(l,l*(l+1)*P_original/(2*np.pi),label="Unmasked")
	ax.plot(l,l*(l+1)*P_masked/(2*np.pi),label="Masked")
	ax.plot(l,l*(l+1)*P_masked/(2*np.pi*(1.0 - mask_profile.maskedFraction)**2),label="Re-scaled")
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

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_map = conv_map.mask(mask_profile)
	assert hasattr(masked_map,"_mask")
	assert masked_map._masked
	assert masked_map.side_angle == conv_map.side_angle
	assert masked_map.data.shape == conv_map.data.shape

	#Compute boundaries
	perimeter_area = masked_map.maskBoundaries()

	fig,ax = plt.subplots(1,3,figsize=(24,8))

	#Plot gradient boundary
	ax[0].imshow(masked_map._gradient_boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle.value,0,conv_map.side_angle.value])

	#Plot hessian (but not gradient) boundary
	ax[1].imshow(masked_map._gradient_boundary ^ masked_map._hessian_boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle.value,0,conv_map.side_angle.value])

	#Plot gradient and hessian boundary
	ax[2].imshow(masked_map.boundary,origin="lower",cmap=plt.cm.binary,interpolation="nearest",extent=[0,conv_map.side_angle.value,0,conv_map.side_angle.value])

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

#Check the differences in PDF with and without masking
def test_pdf():

	th_pdf = np.ogrid[-0.15:0.15:50j]

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_map = conv_map.mask(mask_profile)

	v,p_original = conv_map.pdf(th_pdf)
	v,p_masked = masked_map.pdf(th_pdf)

	#Plot the two histograms
	plt.plot(v,p_original,label="Unmasked")
	plt.plot(v,p_masked,label="Masked {0:.1f}%".format(mask_profile.maskedFraction * 100))

	#Labels
	plt.xlabel(r"$\kappa$")
	plt.ylabel(r"$P(\kappa)$")
	plt.legend(loc="upper right")

	plt.savefig("masked_pdf.png")
	plt.clf()

#Check the differences in peak counts with and without masking
def test_peaks():

	th_peaks = np.ogrid[-0.04:0.12:50j]

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_map = conv_map.mask(mask_profile)

	v,pk_orig = conv_map.peakCount(th_peaks)
	v,pk_masked = masked_map.peakCount(th_peaks)

	#Plot the difference
	plt.plot(v,pk_orig,label=r"Unmasked: $N_p=${0}".format(int(integrate.simps(pk_orig,x=v))))
	plt.plot(v,pk_masked,label=r"With {0:.1f}% area masking: $N_p=${1}".format(mask_profile.maskedFraction * 100,int(integrate.simps(pk_masked,x=v))))
	plt.plot(v,pk_masked/(1.0 - mask_profile.maskedFraction),label="Re-scaled")

	#Labels
	plt.xlabel(r"$\kappa$")
	plt.ylabel(r"$dN/d\kappa$")
	plt.legend(loc="upper left")

	plt.savefig("masked_peaks.png")
	plt.clf()

#Locate the peaks on the map
def test_peak_locations():

	th_peaks = np.arange(0.24,0.5,0.01)

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_map = conv_map.mask(mask_profile)

	#Locate the peaks on the map
	values,location = masked_map.locatePeaks(th_peaks)

	#Visualize the map and the peak locations
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	masked_map.visualize(fig=fig,ax=ax[0],colorbar=True)
	masked_map.visualize(fig=fig,ax=ax[1])

	ax[1].scatter(location[:,0].value,location[:,1].value,color="black")
	ax[1].set_xlim(0.0,masked_map.side_angle.value)
	ax[1].set_ylim(0.0,masked_map.side_angle.value)

	#Save the figure
	fig.tight_layout()
	fig.savefig("masked_peak_locations.png")



#Check the differences in minkowski functionals with and without masking
def test_minkowski():

	th_minkowski = np.ogrid[-0.15:0.15:50j]

	#Set up plots
	fig,ax = plt.subplots(1,3,figsize=(24,8))

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	#Compute and plot the MFs for the unmasked map
	v,V0,V1,V2 = conv_map.minkowskiFunctionals(th_minkowski)
	ax[0].plot(v,V0,label="Unmasked")
	ax[1].plot(v,V1,label="Unmasked")
	ax[2].plot(v,V2,label="Unmasked")

	#Compute and plot the MFs for the masked, zero padded mask
	v,V0,V1,V2 = (conv_map*mask_profile).minkowskiFunctionals(th_minkowski)
	ax[0].plot(v,V0,linestyle="--",label="Zero padded")
	ax[1].plot(v,V1,linestyle="--",label="Zero padded")
	ax[2].plot(v,V2,linestyle="--",label="Zero padded")

	#Compute and plot the MFs for the masked map
	masked_fraction = conv_map.mask(mask_profile,inplace=True)
	v,V0,V1,V2 = conv_map.minkowskiFunctionals(th_minkowski)
	ax[0].plot(v,V0,label="Masked {0:.1f}%".format(masked_fraction*100))
	ax[1].plot(v,V1,label="Masked {0:.1f}%".format(masked_fraction*100))
	ax[2].plot(v,V2,label="Masked {0:.1f}%".format(masked_fraction*100))

	#Labels
	ax[0].set_xlabel(r"$\kappa$")
	ax[0].set_ylabel(r"$V_0(\kappa)$")

	ax[1].set_xlabel(r"$\kappa$")
	ax[1].set_ylabel(r"$V_1(\kappa)$")

	ax[2].set_xlabel(r"$\kappa$")
	ax[2].set_ylabel(r"$V_2(\kappa)$")

	ax[0].legend(loc="upper right")

	fig.tight_layout()

	plt.savefig("masked_minkowski.png")
	plt.clf()

#Check the differences in moments with and without masking
def test_moments():

	conv_map = ConvergenceMap.load(os.path.join(dataExtern(),"unmasked.fit"))
	mask_profile = Mask.load(os.path.join(dataExtern(),"mask.fit"))

	masked_map = conv_map.mask(mask_profile)

	#Compute the moments in both masked and unmasked cases
	mom_original = conv_map.moments(connected=True)
	mom_masked = masked_map.moments(connected=True)
	rel_difference = np.abs(mom_masked/mom_original - 1.0)

	#Save the values and relative differences to file
	np.savetxt("masked_moments.txt",np.array([mom_original,mom_masked,rel_difference]),fmt="%.2e")


