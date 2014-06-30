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

	:param map_filename: File name of the 2D map to analyze (must be in FITS format)
	:type map_filename: str.

	:raises: IOError if FITS file does not exist

	"""

	def __init__(self,map_filename):

		#Open the map file and store the dara in the class instance
		self._map_filename = map_filename
		fits_file = fits.open(map_filename)
		self.data = fits_file[0].data.astype(np.float)
		fits_file.close()

	def peakCount(self,thresholds,sigma=None):
		
		"""
		Counts the peaks in the map

		:param thresholds: thresholds extremes that define the binning of the peak histogram
		:type thresholds: array

		:param sigma: threshold unit; if set to a value, normalizes the thresholds by that value. Useful if one wants to compute the peak counts in units of the map variance
		:type sigma: float.

		:returns: array -- peak counts at the midpoints of the specified thresholds

		"""

		assert thresholds is not None
		if sigma is None:
			sigma = 1.0

		return _topology.peakCount(self.data,thresholds,sigma)