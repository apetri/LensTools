"""

.. module:: convergence
	:platform: Unix
	:synopsis: This module implements the tools to compute topological statistics on 2D convergence maps: measure the power spectrum, counting peaks, measure the minkowski functionals; handling of masks is also provided


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from operator import mul
from functools import reduce

from external import _topology

import numpy as np
from scipy.ndimage import filters


################################################
########ConvergenceMap class####################
################################################

class ConvergenceMap(object):

	"""
	A class that handles 2D convergence maps and allows to compute their topological descriptors (power spectrum, peak counts, minkowski functionals)

	>>> from lenstools import ConvergenceMap 
	>>> from lenstools.defaults import load_fits_default_convergence

	>>> test_map = ConvergenceMap.fromfilename("map.fit",loader=load_fits_default_convergence)
	>>> imshow(test_map.kappa)

	"""

	def __init__(self,kappa,angle,masked=False):

		self.kappa = kappa
		self.side_angle = angle
		self._masked = masked

	@classmethod
	def fromfilename(cls,*args,**kwargs):
		
		"""
		This class method allows to read the map from a data file; the details of the loading are performed by the loader function. The only restriction to this function is that it must return a tuple (angle,kappa)

		:param args: The positional arguments that are to be passed to the loader (typically the file name)

		:param kwargs: Only one keyword is accepted "loader" is a pointer to the previously defined loader method (a template is defaults.load_fits_default_convergence)

		"""

		assert "loader" in kwargs.keys(),"You must specify a loader function!"
		loader = kwargs["loader"]

		angle,kappa = loader(*args)
		return cls(kappa,angle)

	def mask(self,mask_profile,inplace=False):

		"""
		Applies a mask to the convergence map: all masked pixels are given a nan value because they cannot be used in the calculations

		:param mask_profile: profile of the mask, must be an array of 1 byte intergers that are either 0 (if the pixel is masked) or 1 (if the pixel is not masked). Must be of the same shape as the original map
		:type mask_profile: array. or ConvergenceMap instance

		:param inplace: if True the masking is performed in place and the original map is lost, otherwise a new instance of ConvergenceMap with the masked map is returned

		:returns: the masked convergence map if inplace is False, otherwise a float corresponding to the masked fraction of the map

		"""

		if inplace:
			new_map = self
		else:
			new_map = ConvergenceMap(self.kappa.copy(),self.side_angle) 

		if isinstance(mask_profile,np.ndarray):

			assert mask_profile.dtype == np.int8
			assert len(mask_profile[mask_profile==0]) + len(mask_profile[mask_profile==1]) == reduce(mul,mask_profile.shape),"The mask must be made of 0 and 1 only!"
			assert mask_profile.shape == self.kappa.shape

			new_map.kappa[mask_profile==0] = np.nan

		elif isinstance(mask_profile,ConvergenceMap):

			assert mask_profile.side_angle == self.side_angle
			assert len(mask_profile.kappa[mask_profile.kappa==0]) + len(mask_profile.kappa[mask_profile.kappa==1]) == reduce(mul,mask_profile.kappa.shape),"The mask must be made of 0 and 1 only!"
			assert mask_profile.kappa.shape == self.kappa.shape

			new_map.kappa[mask_profile.kappa==0] = np.nan

		else: 

			raise TypeError("Mask type not supported")

		new_map._masked = True
		new_map._mask = ~np.isnan(new_map.kappa)
		new_map._masked_fraction = 1.0 - new_map._mask.sum() / reduce(mul,new_map.kappa.shape)

		#Recompute gradients
		if (hasattr(new_map,"gradient_x") or hasattr(new_map,"gradient_y")):
			new_map.gradient()

		if (hasattr(new_map,"hessian_xx") or hasattr(new_map,"hessian_yy") or hasattr(new_map,"hessian_xy")):
			new_map.hessian()

		#Return
		if inplace:
			return new_map._masked_fraction
		else:
			return new_map 

	def maskBoundaries(self):

		"""
		Computes the mask boundaries defined in the following way: a boundary is a region where the convergence value is defined, but the gradients are not defined.

		:returns: float. : perimeter/area ratio of the mask boundaries
		
		"""

		if not self._masked:
			print("The map is not masked!!")
			return None

		#First check that the instance has the gradient and hessian attributes; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Check where gradient starts to have problems
		nan_gradient_pixels = np.isnan(self.gradient_x) + np.isnan(self.gradient_y)
		self._gradient_boundary = ~self._mask ^ nan_gradient_pixels

		#Check where the hessian has alsp problems
		nan_gradient_pixels = nan_gradient_pixels + np.isnan(self.hessian_xx) + np.isnan(self.hessian_yy) + np.isnan(self.hessian_xy)
		self._hessian_boundary = ~self._mask ^ nan_gradient_pixels

		#Create attribute that holds the full mask (including gradients)
		self._full_mask = self._mask * (~nan_gradient_pixels)
		assert self._full_mask.sum() < self._mask.sum()

		#Compute perimeter/area of the mask
		perimeter_area = self._hessian_boundary.sum() / (~self._full_mask).sum()

		#Return
		return perimeter_area

	@property
	def maskedFraction(self):

		"""
		Computes the masked fraction of the map

		:returns: float, the masked fraction of the map

		"""

		if not self._masked:
			return 0.0
		else:
			return self._masked_fraction

	@property
	def boundary(self):

		"""
		Computes the boundaries of the masked regions, defined as the regions in which the convergence is still well defined but the first and second derivatives are not

		:returns: array of bool. of the same shape as the map, with True values along the boundaries
		
		"""

		if not hasattr(self,"_hessian_boundary"):
			self.maskBoundaries()

		return self._hessian_boundary
			

	def gradient(self):
		
		"""
		Computes the gradient of the map and sets the gradient_x,gradient_y attributes accordingly

		:returns: tuple -- (gradient_x,gradient_y)

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> gx,gy = test_map.gradient()

		"""
		self.gradient_x, self.gradient_y = _topology.gradient(self.kappa)
		return self.gradient_x,self.gradient_y

	def hessian(self):
		
		"""
		Computes the hessian of the map and sets the hessian_xx,hessian_yy,hessian_xy attributes accordingly

		:returns: tuple -- (hessian_xx,hessian_yy,hessian_xy)

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> hxx,hyy,hxy = test_map.hessian()

		"""
		self.hessian_xx,self.hessian_yy,self.hessian_xy = _topology.hessian(self.kappa)
		return self.hessian_xx,self.hessian_yy,self.hessian_xy

	def pdf(self,thresholds,norm=False):

		"""
		Computes the one point probability distribution function of the convergence map

		:param thresholds: thresholds extremes that define the binning of the pdf
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (threshold midpoints -- array, pdf normalized at the midpoints -- array)

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> thresholds = np.arange(map.kappa.min(),map.kappa.max(),0.05)
		>>> nu,p = test_map.pdf(thresholds)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		if norm:
			sigma = self.kappa.std()
		else:
			sigma = 1.0

		#Compute the histogram
		if self._masked:
			hist,bin_edges = np.histogram(self.kappa[self._mask],bins=thresholds*sigma,density=True)
		else:
			hist,bin_edges = np.histogram(self.kappa,bins=thresholds*sigma,density=True)

		#Return
		return midpoints,hist*sigma


	def peakCount(self,thresholds,norm=False):
		
		"""
		Counts the peaks in the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (threshold midpoints -- array, differential peak counts at the midpoints -- array)

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> thresholds = np.arange(map.kappa.min(),map.kappa.max(),0.05)
		>>> nu,peaks = test_map.peakCount(thresholds)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		#Check if the map is masked
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			mask_profile = self._full_mask

		else:

			mask_profile = None 

		#Decide if normalizin thresholds
		if norm:
			
			if self._masked:
				sigma = self.kappa[self._full_mask].std()
			else:
				sigma = self.kappa.std()

		else:
			sigma = 1.0

		return midpoints,_topology.peakCount(self.kappa,mask_profile,thresholds,sigma)

	def minkowskiFunctionals(self,thresholds,norm=False):

		"""
		Measures the three Minkowski functionals (area,perimeter and genus characteristic) of the specified map excursion sets

		:param thresholds: thresholds that define the excursion sets to consider
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (nu -- array, V0 -- array, V1 -- array, V2 -- array) nu are the bins midpoints and V are the Minkowski functionals

		:raises: AssertionError if thresholds array is not provided

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> thresholds = np.arange(-2.0,2.0,0.2)
		>>> nu,V0,V1,V2 = test_map.minkowskiFunctionals(thresholds,norm=True)

		"""

		assert thresholds is not None
		midpoints = 0.5 * (thresholds[:-1] + thresholds[1:])

		#Check if the map is masked
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()

			mask_profile = self._full_mask

		else:

			mask_profile = None 

		#Decide if normalize thresholds or not
		if norm:

			if self._masked:
				sigma = self.kappa[self._full_mask].std()
			else:
				sigma = self.kappa.std()
		
		else:
			sigma = 1.0

		#Check if gradient and hessian attributes are available; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Compute the Minkowski functionals and return them as tuple
		v0,v1,v2 = _topology.minkowski(self.kappa,mask_profile,self.gradient_x,self.gradient_y,self.hessian_xx,self.hessian_yy,self.hessian_xy,thresholds,sigma)

		return midpoints,v0,v1,v2

	def moments(self,connected=False,dimensionless=False):

		"""
		Measures the first nine moments of the convergence map (two quadratic, three cubic and four quartic)

		:param connected: if set to True returns only the connected part of the moments
		:type connected: bool.

		:param dimensionless: if set to True returns the dimensionless moments, normalized by the appropriate powers of the variance
		:type dimensionless: bool. 

		:returns: array -- (sigma0,sigma1,S0,S1,S2,K0,K1,K2,K3)

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> var0,var1,sk0,sk1,sk2,kur0,kur1,kur2,kur3 = test_map.moments()
		>>> sk0,sk1,sk2 = test_map.moments(dimensionless=True)[2:5]
		>>> kur0,kur1,kur2,kur3 = test_map.moments(connected=True,dimensionless=True)[5:]

		"""

		#First check that the instance has the gradient and hessian attributes; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Decide if using the full map or only the unmasked region
		if self._masked:

			if not hasattr(self,"_full_mask"):
				self.maskBoundaries()
			
			kappa = self.kappa[self._full_mask]
			gradient_x = self.gradient_x[self._full_mask]
			gradient_y = self.gradient_y[self._full_mask]
			hessian_xx = self.hessian_xx[self._full_mask]
			hessian_yy = self.hessian_yy[self._full_mask]
			hessian_xy = self.hessian_xy[self._full_mask]

		else:

			kappa = self.kappa
			gradient_x = self.gradient_x
			gradient_y = self.gradient_y
			hessian_xx = self.hessian_xx
			hessian_yy = self.hessian_yy
			hessian_xy = self.hessian_xy

		
		#Quadratic moments
		sigma0 = kappa.std()
		sigma1 = np.sqrt((gradient_x**2 + gradient_y**2).mean())

		#Cubic moments
		S0 = (kappa**3).mean()
		S1 = ((kappa**2)*(hessian_xx + hessian_yy)).mean()
		S2 = ((gradient_x**2 + gradient_y**2)*(hessian_xx + hessian_yy)).mean()

		#Quartic moments
		K0 = (kappa**4).mean()
		K1 = ((kappa**3) * (hessian_xx + hessian_yy)).mean()
		K2 = ((kappa) * (gradient_x**2 + gradient_y**2) * (hessian_xx + hessian_yy)).mean()
		K3 = ((gradient_x**2 + gradient_y**2)**2).mean()

		#Compute connected moments (only quartic affected)
		if connected:
			K0 -= 3 * sigma0**4
			K1 += 3 * sigma0**2 * sigma1**2
			K2 += sigma1**4
			K3 -= 2 * sigma1**4

		
		#Normalize moments to make them dimensionless
		if dimensionless:
			S0 /= sigma0**3
			S1 /= (sigma0 * sigma1**2)
			S2 *= (sigma0 / sigma1**4)

			K0 /= sigma0**4
			K1 /= (sigma0**2 * sigma1**2)
			K2 /= sigma1**4
			K3 /= sigma1**4

			sigma0 /= sigma0
			sigma1 /= sigma1

		#Return the array
		return np.array([sigma0,sigma1,S0,S1,S2,K0,K1,K2,K3])

	def powerSpectrum(self,l_edges):

		"""
		Measures the power spectrum of the convergence map at the multipole moments specified in the input

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:returns: tuple -- (l -- array,Pl -- array) = (multipole moments, power spectrum at multipole moments)

		:raises: AssertionError if l_edges are not provided

		>>> test_map = ConvergenceMap.fromfilename("map.fit")
		>>> l_edges = np.arange(200.0,5000.0,200.0)
		>>> l,Pl = test_map.powerSpectrum(l_edges)

		"""

		assert not self._masked,"Power spectrum calculation for masked maps is not allowed yet!"
		assert l_edges is not None

		l = 0.5*(l_edges[:-1] + l_edges[1:])

		#Calculate the Fourier transform of the map with numpy FFT
		ft_map = np.fft.rfft2(self.kappa)

		#Compute the power spectrum with the C backend implementation
		power_spectrum = _topology.rfft2_azimuthal(ft_map,ft_map,self.side_angle,l_edges)

		#Output the power spectrum
		return l,power_spectrum

	def cross(self,other,l_edges):

		"""
		Measures the cross power spectrum between two convergence maps at the multipole moments specified in the input

		:param other: The other convergence map
		:type other: ConvergenceMap instance

		:param l_edges: Multipole bin edges
		:type l_edges: array

		:returns: tuple -- (l -- array,Pl -- array) = (multipole moments, cross power spectrum at multipole moments)

		:raises: AssertionError if l_edges are not provided or the other map has not the same shape as the input one

		>>> test_map = ConvergenceMap.fromfilename("map.fit",loader=load_fits_default_convergence)
		>>> other_map = ConvergenceMap.fromfilename("map2.fit",loader=load_fits_default_convergence)
		
		>>> l_edges = np.arange(200.0,5000.0,200.0)
		>>> l,Pl = test_map.cross(other_map,l_edges)

		"""

		assert not (self._masked or other._masked),"Power spectrum calculation for masked maps is not allowed yet!"

		assert l_edges is not None
		l = 0.5*(l_edges[:-1] + l_edges[1:])

		assert isinstance(other,ConvergenceMap)
		assert self.side_angle == other.side_angle
		assert self.kappa.shape == other.kappa.shape

		#Calculate the Fourier transform of the maps with numpy FFTs
		ft_map1 = np.fft.rfft2(self.kappa)
		ft_map2 = np.fft.rfft2(other.kappa)

		#Compute the cross power spectrum with the C backend implementation
		cross_power_spectrum = _topology.rfft2_azimuthal(ft_map1,ft_map2,self.side_angle,l_edges)

		#Output the cross power spectrum
		return l,cross_power_spectrum


	def smooth(self,angle_in_arcmin,kind="gaussian",inplace=False):

		"""
		Performs a smoothing operation on the convergence map

		:param angle_in_arcmin: size of the smoothing kernel in arcminutes
		:type angle_in_arcmin: float.

		:param kind: type of smoothing to be performed (only implemented gaussian so far)
		:type kind: str.

		:param inplace: if set to True performs the smoothing in place overwriting the old convergence map
		:type inplace: bool.

		:returns: ConvergenceMap instance (or None if inplace is True)

		"""

		assert not self._masked,"You cannot smooth a masked convergence map!!"
		assert kind == "gaussian","Only gaussian smoothing implemented!!"

		#Compute the smoothing scale in pixel units
		smoothing_scale_pixel = angle_in_arcmin * self.kappa.shape[0] / (self.side_angle*60.0)

		#Perform the smoothing
		smoothed_kappa = filters.gaussian_filter(self.kappa,smoothing_scale_pixel)

		#Return the result
		if inplace:
			
			self.kappa = smoothed_kappa
			
			#If gradient attributes are present, recompute them
			if (hasattr(self,"gradient_x") or hasattr(self,"gradient_y")):
				self.gradient()

			if (hasattr(self,"hessian_xx") or hasattr(self,"hessian_yy") or hasattr(self,"hessian_xy")):
				self.hessian()
			
			return None

		else:
			return ConvergenceMap(smoothed_kappa,self.side_angle)


	def __add__(self,rhs):

		"""
		Defines addition operator between ConvergenceMap instances; the convergence values are summed

		:returns: ConvergenceMap instance with the result of the sum

		:raises: AssertionError if the operation cannot be performed

		"""

		if isinstance(rhs,ConvergenceMap):

			assert self.side_angle == rhs.side_angle
			assert self.kappa.shape == rhs.kappa.shape

			new_kappa = self.kappa + rhs.kappa

		elif type(rhs) == np.float:

			new_kappa = self.kappa + rhs

		elif type(rhs) == np.ndarray:

			assert rhs.shape == self.kappa.shape
			new_kappa = self.kappa + rhs

		else:

			raise TypeError("The right hand side cannot be added!!")

		return ConvergenceMap(new_kappa,self.side_angle)


	def __mul__(self,rhs):

		"""
		Defines the multiplication operator between ConvergenceMap instances; the convergence values are multiplied (for example the rhs can be a mask...)

		:returns: ConvergenceMap instances with the result of the multiplication

		:raises: AssertionError if the operation cannot be performed

		""" 

		if isinstance(rhs,ConvergenceMap):

			assert self.side_angle == rhs.side_angle
			assert self.kappa.shape == rhs.kappa.shape

			new_kappa = self.kappa * rhs.kappa

		elif type(rhs) == np.float:

			new_kappa = self.kappa * rhs

		elif type(rhs) == np.ndarray:

			assert rhs.shape == self.kappa.shape
			new_kappa = self.kappa * rhs

		else:

			raise TypeError("Cannot multiply by the right hand side!!")

		return ConvergenceMap(new_kappa,self.side_angle)



###############################################
###################Mask class##################
###############################################

class Mask(ConvergenceMap):

	"""
	A class that handles the convergence masked regions

	"""

	def __init__(self,kappa,angle,masked=False):

		super(Mask,self).__init__(kappa,angle,masked)
		assert len(self.kappa[self.kappa==0]) + len(self.kappa[self.kappa==1]) == reduce(mul,self.kappa.shape),"The mask must be made of 0 and 1 only!"
		self._masked_fraction = len(self.kappa[self.kappa==0]) / reduce(mul,self.kappa.shape)

	@property
	def maskedFraction(self):

		"""
		Computes the fraction of masked pixels

		"""

		return self._masked_fraction

	def maskBoundaries(self):
		raise AttributeError("Can only call this method on a ConvergenceMap instance!")

	@property
	def boundary(self):
		raise AttributeError("This property is only defined on ConvergenceMap instances!")

