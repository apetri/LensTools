import os

from .. import ConvergenceMap,ShearMap,GaussianNoiseGenerator
from ..utils.defaults import sample_power_shape

from .. import dataExtern

import numpy as np
import matplotlib.pyplot as plt

from astropy.units import deg,arcmin

test_map_conv = ConvergenceMap.load(os.path.join(dataExtern(),"conv.fit"))

shape_noise_gen = GaussianNoiseGenerator.forMap(test_map_conv)
corr_noise_gen = GaussianNoiseGenerator.forMap(test_map_conv)

test_map_noisy = test_map_conv + shape_noise_gen.getShapeNoise(z=1.0,ngal=15.0*arcmin**-2,seed=1)

l = np.arange(200.0,50000.0,200.0)
scale = 5000.0


def test_smooth():

	test_map_conv_smoothed = test_map_conv.smooth(1.0*arcmin)

	fig,ax = plt.subplots(1,2,figsize=(16,8))
	test_map_conv.visualize(fig,ax[0])
	test_map_conv_smoothed.visualize(fig,ax[1])

	ax[0].set_title("Unsmoothed")
	ax[1].set_title(r"$1^\prime$")

	fig.tight_layout()
	fig.savefig("smooth.png")

def test_shape_noise():

	fig,ax = plt.subplots(1,3,figsize=(24,8))

	test_map_conv.visualize(fig,ax[0])
	test_map_noisy.visualize(fig,ax[1])

	test_map_conv_smoothed = test_map_noisy.smooth(1.0*arcmin)
	test_map_conv_smoothed.visualize(fig,ax[2])

	ax[0].set_title("Bare")
	ax[1].set_title("Noisy")
	ax[2].set_title(r"Noisy, $1^\prime$")

	fig.tight_layout()
	fig.savefig("shape_noise.png")

def test_correlated_convergence_power():

	fig,ax = plt.subplots()
	
	#Plot power spectral density
	ax.plot(l,l*(l+1)*sample_power_shape(l,scale=scale)/(2.0*np.pi),label="Original")

	#Generate three realizations of this power spectral density and plot power spectrum for cross check
	for i in range(3):
		noise_map = corr_noise_gen.fromConvPower(sample_power_shape,seed=i,scale=scale)
		ell,Pl = noise_map.powerSpectrum(l)

		ax.plot(ell,ell*(ell+1)*Pl/(2.0*np.pi),label="Realization {0}".format(i+1),linestyle="--")


	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	ax.legend(loc="lower right")

	ax.set_yscale("log")

	plt.savefig("correlated_power.png")
	plt.clf()

def test_correlated_convergence_maps():

	fig,ax = plt.subplots(1,3,figsize=(24,8))

	#Generate three realizations of this power spectral density and plot them for cross check
	for i in range(3):
		
		noise_map = corr_noise_gen.fromConvPower(sample_power_shape,seed=i,scale=scale)
		noise_map.visualize(fig,ax[i])

	fig.tight_layout()
	fig.savefig("correlated_maps.png")

def test_interpolated_convergence_power():

	fig,ax = plt.subplots()

	power_func = np.loadtxt(os.path.join(dataExtern(),"ee4e-7.txt"),unpack=True)
	l_in,Pl_in = power_func
	
	#Plot power spectral density
	ax.plot(l_in,l_in*(l_in+1)*Pl_in/(2.0*np.pi),label="Original")

	#Generate three realizations of this power spectral density and plot power spectrum for cross check
	for i in range(3):
		noise_map = corr_noise_gen.fromConvPower(power_func,seed=i,bounds_error=False,fill_value=0.0)
		ell,Pl = noise_map.powerSpectrum(l_in)

		ax.plot(ell,ell*(ell+1)*Pl/(2.0*np.pi),label="Realization {0}".format(i+1),linestyle="--")


	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")

	ax.legend(loc="upper right")

	ax.set_yscale("log")

	plt.savefig("interpolated_power.png")
	plt.clf()

def test_interpolated_convergence_maps():

	fig,ax = plt.subplots(1,3,figsize=(24,8))
	power_func = np.loadtxt(os.path.join(dataExtern(),"ee4e-7.txt"),unpack=True)

	#Generate three realizations of this power spectral density and plot them for cross check
	for i in range(3):
		
		noise_map = corr_noise_gen.fromConvPower(power_func,seed=i,bounds_error=False,fill_value=0.0)
		noise_map.visualize(fig,ax[i])

	fig.tight_layout()
	fig.savefig("interpolated_maps.png")







