import logging

import numpy as np
from astropy.io import fits

from topology import ConvergenceMap

##########################################
#####Default Fits loader convergence######
##########################################
def load_fits_default_convergence(*args):
	"""
	This is the default convergence fits file loader, it assumes that the two components of the shear are stored in two different image FITS files, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- kappa; kappa is the convergence map)

	:raises: IOError if the FITS files cannot be opened or do not exist

	"""

	#Open the files
	kappa_file = fits.open(args[0])

	#Read the ANGLE keyword from the header
	angle = kappa_file[0].header["ANGLE"]

	#Create the array with the shear map
	kappa = kappa_file[0].data.astype(np.float)

	#Close files and return
	kappa_file.close()

	return angle,kappa

##########################################
#####Default Fits loader shear############
##########################################
def load_fits_default_shear(*args):
	"""
	This is the default shear fits file loader, it assumes that the two components of the shear are stored in a single image FITS file, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- gamma; gamma[0] is the gamma1 map, gamma[1] is the gamma2 map); the maps must follow matrix ordering, i.e. the first axis (0) is y and the second axis (1) is x. This matters for the E/B mode decomposition 

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

################################################################################
##########Default callback loader, loads in the measured power spectrum#########
################################################################################

def default_callback_loader(args):

	assert "file_name" in args.keys()
	assert "l_edges" in args.keys()

	logging.debug("Processing {0}".format(args["file_name"]))

	conv_map = ConvergenceMap.fromfilename(args["file_name"],loader=load_fits_default_convergence)
	l,Pl = conv_map.powerSpectrum(args["l_edges"])
	return Pl