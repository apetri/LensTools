"""

.. module:: noise
	:platform: Unix
	:synopsis: This module implements various tools to generate simulated noise maps


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from .convergence import ConvergenceMap

import numpy as np

#FFT engine
from ..utils.fft import NUMPYFFTPack
fftengine = NUMPYFFTPack()

from scipy import interpolate

#Units 
from astropy.units import deg,arcmin

########################################################
########GaussianNoiseGenerator class####################
########################################################

class GaussianNoiseGenerator(object):

	"""
	A class that handles generation of Gaussian simulated noise maps

	"""

	def __init__(self,shape,side_angle,label):

		#Sanity check
		assert side_angle.unit.physical_type=="angle"
		
		self.shape = shape
		self.side_angle = side_angle
		self.label = label

	@classmethod
	def forMap(cls,conv_map):

		"""
		This class method generates a Gaussian noise generator intended to be used on a convergence map: i.e. the outputs of its methods can be added to the convergence map in question to simulate the presence of noise

		:param conv_map: The blueprint of the convergence map you want to generate the noise for
		:type conv_map: ConvergenceMap instance

		:raises: AssertionError if conv_map is not a ConvergenceMap instance

		"""

		assert isinstance(conv_map,ConvergenceMap)

		return cls(conv_map.data.shape,conv_map.side_angle,label="convergence")

	def getShapeNoise(self,z=1.0,ngal=15.0*arcmin**-2,seed=0):

		"""
		This method generates a white, gaussian shape noise map for the given redshift of the map

		:param z: single redshift of the backround sources on the map
		:type z: float.

		:param ngal: assumed angular number density of galaxies (must have units of angle^-2)
		:type ngal: float.

		:param seed: seed of the random generator
		:type seed: int.

		:returns: ConvergenceMap instance of the same exact shape as the one used as blueprint

		"""

		#Sanity check
		assert (ngal.unit**-0.5).physical_type=="angle"

		if self.label == "convergence":
		
			#Compute shape noise amplitude
			pixel_angular_side = self.side_angle / self.shape[0]
			sigma = ((0.15 + 0.035*z) / (pixel_angular_side * np.sqrt(ngal))).decompose().value

			#Generate shape noise
			np.random.seed(seed)
			noise_map = np.random.normal(loc=0.0,scale=sigma,size=self.shape) 

			#Build the ConvergenceMap object
			return ConvergenceMap(noise_map,self.side_angle)

		else:

			raise ValueError("Only convergence implemented so far!!!")

	def _fourierMap(self,power_func,**kwargs):

		#Assert the shape of the blueprint, to tune the right size for the fourier transform
		lpix = 360.0/self.side_angle.to(deg).value
		lx = fftengine.rfftfreq(self.shape[0]) * self.shape[0] * lpix
		ly = fftengine.fftfreq(self.shape[0]) * self.shape[0] * lpix

		#Compute the multipole moment of each FFT pixel
		l = np.sqrt(lx[np.newaxis,:]**2 + ly[:,np.newaxis]**2)

		#Compute the power spectrum at each l and check that it is positive 
		if isinstance(power_func,np.ndarray):
			
			#Check for correct shape
			assert power_func.shape[0] == 2,"If you want an interpolated power spectrum you should pass a (l,Pl) array!"

			#Perform the interpolation
			ell,Pell = power_func
			power_interp = interpolate.interp1d(ell,Pell,**kwargs)
			Pl = power_interp(l)

		else:
			
			Pl = power_func(l,**kwargs)
		

		assert Pl[Pl>=0.0].size == Pl.size

		#Generate real and imaginary parts
		real_part = np.sqrt(0.5*Pl) * np.random.normal(loc=0.0,scale=1.0,size=l.shape) * lpix/(2.0*np.pi)
		imaginary_part = np.sqrt(0.5*Pl) * np.random.normal(loc=0.0,scale=1.0,size=l.shape) * lpix/(2.0*np.pi)

		#Get map in real space and return
		ft_map = (real_part + imaginary_part*1.0j) * l.shape[0]**2
		ft_map[0,0] = 0.0

		return ft_map



	def fromConvPower(self,power_func,seed=0,**kwargs):

		"""
		This method uses a supplied power spectrum to generate correlated noise maps in real space via FFTs

		:param power_func: function that given a numpy array of l's returns a numpy array with the according Pl's (this is the input power spectrum); alternatively you can pass an array (l,Pl) and the power spectrum will be calculated with scipy's interpolation routines
		:type power_func: function with the above specifications, or numpy array (l,Pl) of shape (2,n) 

		:param seed: seed of the random generator 
		:type seed: int.

		:param kwargs: keyword arguments to be passed to power_func, or to the interpolate.interp1d routine

		:returns: ConvergenceMap instance of the same exact shape as the one used as blueprint

		"""
		assert self.label == "convergence"

		#Initialize random number generator
		np.random.seed(seed)

		#Generate a random Fourier realization and invert it
		ft_map = self._fourierMap(power_func,**kwargs)
		noise_map = fftengine.irfft2(ft_map)

		return ConvergenceMap(noise_map,self.side_angle)




