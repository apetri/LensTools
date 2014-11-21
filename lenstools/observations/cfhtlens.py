from __future__ import division,print_function,with_statement

import os

import numpy as np
from astropy.io import fits
from astropy.units import deg

from .. import ConvergenceMap

######################################
#######Loader function for CFHT#######
######################################

def cfht_load(self,filename):

	kappa_file = fits.open(filename)
	angle = 3.4641016151377544 * deg
	kappa = kappa_file[0].data.astype(np.float)
	kappa_file.close()

	return angle,kappa


################################################
###############CFHTLens class###################
################################################

class CFHTLens(object):

	"""
	Class handler of the CFHTLens reduced data set, already split in 13 3x3 deg^2 subfields 

	"""

	_data_loader=cfht_load


	def __init__(self,root_path=None):

		assert root_path is not None,"You must specify the root path of your CFHTLens subfield observations local copy!"
		self.root_path = root_path

	def getName(self,subfield=1,smoothing=0.5):

		"""
		Get the names of the FITS files where the subfield images are stored

		:param subfield: the subfield number (1-13)
		:type subfield: int.

		:param smoothing: the smoothing scale of the subfield image, in arcminutes
		:type smoothing: float.

		:returns: str. : the FITS file name

		"""

		full_name = os.path.join(self.root_path,"CFHT_KS")
		full_name += "_sigma{0:02d}_".format(int(smoothing*10))
		full_name += "subfield{0:02d}".format(subfield)
		full_name += ".fits"

		return full_name


	def load(self,**kwargs):

		"""
		Loads in a CFHT observation as a ConvergenceMap instance

		:param kwargs: the keyword arguments are passed to the getName method

		:returns: ConvergenceMap instance

		"""

		filename=self.getName(**kwargs)
		return ConvergenceMap.load(filename,format=self._data_loader)