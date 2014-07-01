"""

.. module:: topology 
	:platform: Unix
	:synopsis: This module implements the tools to compute topological statistics on 2D maps: measure the power spectrum, counting peaks, measure the minkowski functionals


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from external import _topology

import numpy as np
from astropy.io import fits

################################################
########ConvergenceMap class####################
################################################

class ConvergenceMap(object):

	"""
	A class that handles 2D convergence maps and allows to compute their topological descriptors (power spectrum, peak counts, minkowski functionals)

	:param map_filename: File name of the 2D map to analyze (must be in FITS format); sets the attibute 'kappa' to a numpy array containing the convergence values
	:type map_filename: str.

	:raises: IOError if FITS file does not exist

	>>> map = ConvergenceMap("map.fit")
	>>> imshow(map.kappa)

	"""

	def __init__(self,map_filename):

		#Open the map file and store the dara in the class instance
		self._map_filename = map_filename
		fits_file = fits.open(map_filename)
		self.kappa = fits_file[0].data.astype(np.float)
		fits_file.close()

	def gradient(self):
		
		"""
		Computes the gradient of the map and sets the gradient_x,gradient_y attributes accordingly

		:returns: tuple -- (gradient_x,gradient_y)

		>>> map = ConvergenceMap("map.fit")
		>>> gx,gy = map.gradient()

		"""
		self.gradient_x, self.gradient_y = _topology.gradient(self.kappa)
		return self.gradient_x,self.gradient_y

	def hessian(self):
		
		"""
		Computes the hessian of the map and sets the hessian_xx,hessian_yy,hessian_xy attributes accordingly

		:returns: tuple -- (hessian_xx,hessian_yy,hessian_xy)

		>>> map = ConvergenceMap("map.fit")
		>>> hxx,hyy,hxy = map.hessian()

		"""
		self.hessian_xx,self.hessian_yy,self.hessian_xy = _topology.hessian(self.kappa)
		return self.hessian_xx,self.hessian_yy,self.hessian_xy


	def peakCount(self,thresholds,norm=False):
		
		"""
		Counts the peaks in the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: array -- differential peak counts at the midpoints of the specified thresholds

		:raises: AssertionError if thresholds array is not provided

		>>> map = ConvergenceMap("map.fits")
		>>> thresholds = np.arange(map.kappa.min(),map.kappa.max(),0.05)
		>>> peaks = map.peakCount(thresholds)

		"""

		assert thresholds is not None
		if norm:
			sigma = self.kappa.std()
		else:
			sigma = 1.0

		return _topology.peakCount(self.kappa,thresholds,sigma)

	def minkowskiFunctionals(self,thresholds,norm=False):

		"""
		Measures the three Minkowski functionals (area,perimeter and genus characteristic) of the specified map excursion sets

		:param thresholds: thresholds that define the excursion sets to consider
		:type thresholds: array

		:param norm: normalization; if set to a True, interprets the thresholds array as units of sigma (the map standard deviation)
		:type norm: bool.

		:returns: tuple -- (V0 -- array, V1 -- array, V2 -- array) each array represents one of the Minkowski functionals

		:raises: AssertionError if thresholds array is not provided

		>>> map = ConvergenceMap("map.fits")
		>>> thresholds = np.arange(-2.0,2.0,0.2)
		>>> V0,V1,V2 = map.minkowski(thresholds,norm=True)

		"""

		assert thresholds is not None
		if norm:
			sigma = self.kappa.std()
		else:
			sigma = 1.0

		#Check if gradient and hessian attributes are available; if not, compute them
		if not (hasattr(self,"gradient_x") and hasattr(self,"gradient_y")):
			self.gradient()

		if not (hasattr(self,"hessian_xx") and hasattr(self,"hessian_yy") and hasattr(self,"hessian_xy")):
			self.hessian()

		#Compute the Minkowski functionals and return them as tuple
		return _topology.minkowski(self.kappa,self.gradient_x,self.gradient_y,self.hessian_xx,self.hessian_yy,self.hessian_xy,thresholds,sigma)

