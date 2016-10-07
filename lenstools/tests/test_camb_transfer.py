import os

from ..simulations.camb import CAMBTransferFunction
from ..simulations.raytracing import PotentialPlane

from .. import dataExtern

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

def test_transfer_plot():

	#Set up plot, load transfer function
	tfr = CAMBTransferFunction.read(os.path.join(dataExtern(),"camb","camb_tfr.pkl"))
	fig,ax = plt.subplots()

	#Plot at three different redshifts
	zp = [0.5,1.,2.]
	k = np.logspace(-3,1,100) / u.Mpc

	for z in zp:
		ax.plot(k,tfr(z,k),label=r"$z={0:.2f}$".format(z))

	#Labels
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$k({\rm Mpc}^{-1})$",fontsize=22)
	ax.set_ylabel(r"$T_{\rm CDM}(k)/k^2 (h^{-2}{\rm Mpc}^2)$",fontsize=22)
	ax.legend(loc="upper right")

	#Save 
	fig.savefig("camb_tfr.png")


def test_plane_scaling():

	#Load the transfer function, load the lens plane
	plane = PotentialPlane.load(os.path.join(dataExtern(),"plane.fits"))
	tfr = CAMBTransferFunction.read(os.path.join(dataExtern(),"camb","camb_tfr.pkl"))

	#Set up the plot
	fig,ax = plt.subplots(1,3,figsize=(24,8))

	#Visualize the original plane in the first panel
	plane.visualize(fig=fig,ax=ax[0],colorbar=True,cbar_label=r"$\Phi$")
	ax[0].set_title(r"Original ($z={0:.2f}$)".format(plane.redshift))

	#Scale to redshift 1.5 (uniform)
	zs = 1.5
	plane.scaleWithTransfer(zs,tfr,scaling_method="uniform")
	plane.visualize(fig=fig,ax=ax[1],colorbar=True,cbar_label=r"$\Phi$")
	ax[1].set_title(r"Uniform scaling ($z={0:.2f}$)".format(zs))

	#Scale to redshift 1.5 (FFT)
	zs = 1.5
	plane = PotentialPlane.load(os.path.join(dataExtern(),"plane.fits"))
	plane.scaleWithTransfer(zs,tfr,scaling_method="FFT")
	plane.visualize(fig=fig,ax=ax[2],colorbar=True,cbar_label=r"$\Phi$")
	ax[2].set_title(r"FFT scaling ($z={0:.2f}$)".format(zs))

	#Save 
	fig.tight_layout()
	fig.savefig("plane_scaling_tfr.png")