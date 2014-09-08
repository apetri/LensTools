from __future__ import division,print_function,with_statement

import os

################################################
###############CFHTLens class###################
################################################

class CFHTLens(object):

	"""
	Class handler of the CFHTLens reduced data set, already split in 13 3x3 deg^2 subfields 

	"""

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