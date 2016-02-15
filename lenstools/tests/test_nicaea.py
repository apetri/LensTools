import os

import sys
from ..simulations import Nicaea,NicaeaSettings

from .. import dataExtern

import numpy as np
import matplotlib.pyplot as plt

#astropy
from astropy.cosmology import WMAP9
import astropy.units as u

#Settings
settings = NicaeaSettings.default()

#Cosmology
cosmo = Nicaea.fromCosmology(WMAP9)

def test_power_spectrum():

	#Compute power spectrum
	ell = np.arange(300.0,1.0e5,500.0)

	try:
		power = cosmo.convergencePowerSpectrum(ell,z=2.0,settings=settings)
	except ImportError:
		return

	#Plot
	fig,ax = plt.subplots(figsize=(8,8))
	ax.plot(ell,ell*(ell+1)*power/(2.0*np.pi))
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$l$",fontsize=18)
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$",fontsize=18) 

	#Save figure
	fig.savefig("power_nicaea.png")


def test_2pt_correlation():

	#Compute the correlation function
	theta = np.arange(1.0,100.0,1.0)*u.arcmin

	try:
		xi_plus = cosmo.shearTwoPoint(theta,z=2.0,settings=settings,kind="+")
		xi_minus = cosmo.shearTwoPoint(theta,z=2.0,settings=settings,kind="-")
	except ImportError:
		return

	#Plot
	fig,ax = plt.subplots(figsize=(8,8))
	ax.plot(theta,xi_plus,label=r"$\xi_+$")
	ax.plot(theta,xi_minus,label=r"$\xi_-$")
	ax.set_xscale("log")
	ax.set_xlabel(r"$\theta(\mathrm{arcmin})$",fontsize=18)
	ax.set_ylabel(r"$\xi(\theta)$",fontsize=18) 
	ax.legend()

	#Save figure
	fig.savefig("2pt_nicaea.png")