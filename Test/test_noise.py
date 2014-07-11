try:
	
	from lenstools import ConvergenceMap,ShearMap,GaussianNoiseGenerator
	from lenstools.defaults import load_fits_default_convergence,load_fits_default_shear,sample_power_shape

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools import ConvergenceMap,ShearMap,GaussianNoiseGenerator
	from lenstools.defaults import load_fits_default_convergence,load_fits_default_shear,sample_power_shape

import numpy as np
import matplotlib.pyplot as plt

test_map_conv = ConvergenceMap.fromfilename("conv.fit",loader=load_fits_default_convergence)

shape_noise_gen = GaussianNoiseGenerator.forMap(test_map_conv)
corr_noise_gen = GaussianNoiseGenerator.forMap(test_map_conv)

test_map_noisy = test_map_conv + shape_noise_gen.getShapeNoise(z=1.0,ngal=15.0,seed=1)

l = np.arange(200.0,50000.0,200.0)
scale = 5000.0


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

def test_shape_noise():

	fig,ax = plt.subplots(1,3,figsize=(24,8))
	
	ax[0].imshow(test_map_conv.kappa,origin="lower",interpolation="nearest",extent=[0,test_map_conv.side_angle,0,test_map_conv.side_angle])
	ax[1].imshow(test_map_noisy.kappa,origin="lower",interpolation="nearest",extent=[0,test_map_conv.side_angle,0,test_map_conv.side_angle])

	test_map_conv_smoothed = test_map_noisy.smooth(1.0)

	ax[2].imshow(test_map_conv_smoothed.kappa,origin="lower",interpolation="nearest",extent=[0,test_map_conv.side_angle,0,test_map_conv.side_angle])

	ax[0].set_xlabel(r"$x$(deg)")
	ax[0].set_ylabel(r"$y$(deg)")
	ax[0].set_title("Bare")
	ax[1].set_xlabel(r"$x$(deg)")
	ax[1].set_ylabel(r"$y$(deg)")
	ax[1].set_title("Noisy")
	ax[2].set_xlabel(r"$x$(deg)")
	ax[2].set_ylabel(r"$y$(deg)")
	ax[2].set_title(r"Noisy, $1^\prime$")

	fig.tight_layout()
	plt.savefig("shape_noise.png")

	plt.clf()

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
		ax[i].imshow(noise_map.kappa,origin="lower",interpolation="nearest",extent=[0,noise_map.side_angle,0,noise_map.side_angle])
		ax[i].set_xlabel(r"$x$(deg)")
		ax[i].set_ylabel(r"$y$(deg)")

	fig.tight_layout()
	plt.savefig("correlated_maps.png")

	plt.clf()

def test_interpolated_convergence_power():

	fig,ax = plt.subplots()

	power_func = np.loadtxt("ee4e-7.txt",unpack=True)
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
	power_func = np.loadtxt("ee4e-7.txt",unpack=True)

	#Generate three realizations of this power spectral density and plot them for cross check
	for i in range(3):
		
		noise_map = corr_noise_gen.fromConvPower(power_func,seed=i,bounds_error=False,fill_value=0.0)
		ax[i].imshow(noise_map.kappa,origin="lower",interpolation="nearest",extent=[0,noise_map.side_angle,0,noise_map.side_angle])
		ax[i].set_xlabel(r"$x$(deg)")
		ax[i].set_ylabel(r"$y$(deg)")

	fig.tight_layout()
	plt.savefig("interpolated_maps.png")

	plt.clf()







