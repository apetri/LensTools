import logging

import numpy as np
from astropy.io import fits

from topology import ConvergenceMap
from index import PowerSpectrum,Peaks

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
	"""
	
	Default ensemble loader: reads a FITS data file containing a convergence map and measures its power spectrum

	:param args: A dictionary that contains all the relevant parameters as keys. Must have a "file_name" key
	:type args: Dictionary

	:returns: ndarray of the measured statistics

	:raises: AssertionError if the input dictionary doesn't have the required keywords

	"""

	assert "file_name" in args.keys()
	assert "l_edges" in args.keys()

	logging.debug("Processing {0} power".format(args["file_name"]))

	conv_map = ConvergenceMap.fromfilename(args["file_name"],loader=load_fits_default_convergence)
	l,Pl = conv_map.powerSpectrum(args["l_edges"])
	return Pl

def peaks_loader(args):

	assert "file_name" in args.keys()
	assert "thresholds" in args.keys()

	logging.debug("Processing {0} peaks".format(args["file_name"]))
	conv_map = ConvergenceMap.fromfilename(args["file_name"],loader=load_fits_default_convergence)

	v,pk = conv_map.peakCount(args["thresholds"],norm=True)
	return v

def convergence_measure_all(args):

	assert "file_name" in args.keys()
	assert "index" in args.keys()

	logging.debug("Processing {0}".format(args["file_name"]))

	#Load the map
	conv_map = ConvergenceMap.fromfilename(args["file_name"],loader=load_fits_default_convergence)

	#Allocate memory for observables
	descriptors = args["index"]
	observables = np.zeros(descriptors.size)

	#Measure descriptors as directed by input
	for n in range(descriptors.num_descriptors):

		if isinstance(descriptors[n],PowerSpectrum):
			l,observables[descriptors[n].first:descriptors[n].last] = conv_map.powerSpectrum(descriptors[n].l_edges)
		elif isinstance(descriptors[n],Peaks):
			v,observables[descriptors[n].first:descriptors[n].last] = conv_map.peakCount(descriptors[n].thresholds,norm=True)
		else:
			raise ValueError("Measurement of this descriptor not implemented!!!")

	#Return
	return observables






####################################################################################
#############Default power spectrum template for testing############################
####################################################################################

def sample_power_shape(l,**kwargs):

	return np.exp(-0.5*(l/kwargs["scale"])**2)