"""

.. module:: shear 
	:platform: Unix
	:synopsis: This module implements a set of operations which are usually performed on weak lensing shear maps


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from external import _topology

import numpy as np
from astropy.io import fits

##########################################
#####Default Fits loader##################
##########################################
def load_fits_default(*args):
	"""
	This is the default fits file loader, it assumes that the two components of the shear are stored in two different image FITS files, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- gamma; gamma[0] is the gamma1 map, gamma[1] is the gamma2 map)

	:raises: IOError if the FITS files cannot be opened or do not exist

	"""

	#Open the files
	gamma_file = fits.open(args[0])

	#Read the ANGLE keyword from the header
	angle = gamma_file[0].header["ANGLE"]

	#Create the array with the shear map
	gamma = gamma_file[0].data.astype(np.float)

	#Close files and return
	gamma_file.close()

	return angle,gamma


##########################################
########ShearMap class####################
##########################################

class ShearMap(object):

	"""
	A class that handles 2D shear maps and allows to perform a set of operations on them

	:param loader: FITS file loading utility, must match the signature and return type of load_fits_default
	:type loader: keyword argument, function

	:param args: positional arguments must be the exact same as the ones that loader takes

	>>> from shear import ShearMap
	>>> test = ShearMap("shear.fit",loader=load_fits_default)
	>>>test.side_angle
	1.95
	>>>test.gamma
	#The actual map values

	"""

	def __init__(self,*args,**kwargs):

		if not("loader" in kwargs.keys()):
			loader = load_fits_default
		else:
			loader = kwargs["loader"]

		self.side_angle,self.gamma = loader(*args)

