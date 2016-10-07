"""

.. module:: defaults 
	:platform: Unix
	:synopsis: This module contains a bunch of default functions for data files loading, specifically in FITS format. This package is totally flexible and does not require that your data is in FITS format, so this module should be regarded just as an example


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import glob,os,re

import logging

import numpy as np
from astropy.io import fits
from astropy.units import deg

from ..image.convergence import ConvergenceMap

##########################################
#####Default Fits loader convergence######
##########################################
def load_fits_default_convergence(filename):
	"""
	This is the default convergence fits file loader, it assumes that the two components of the shear are stored in two different image FITS files, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- kappa; kappa is the convergence map)

	:raises: IOError if the FITS files cannot be opened or do not exist

	"""

	#Open the files
	kappa_file = fits.open(filename)

	#Read the ANGLE keyword from the header
	angle = kappa_file[0].header["ANGLE"]

	#Create the array with the shear map
	kappa = kappa_file[0].data.astype(np.float)

	#Close files and return
	kappa_file.close()

	return angle*deg,kappa

##########################################
#####Default Fits loader shear############
##########################################
def load_fits_default_shear(filename):
	"""
	This is the default shear fits file loader, it assumes that the two components of the shear are stored in a single image FITS file, which have an ANGLE keyword in the header

	:param gamma_file: Name of the FITS file that contains the shear map
	:type gamma1_file: str.

	:returns: tuple -- (angle,ndarray -- gamma; gamma[0] is the gamma1 map, gamma[1] is the gamma2 map); the maps must follow matrix ordering, i.e. the first axis (0) is y and the second axis (1) is x. This matters for the E/B mode decomposition 

	:raises: IOError if the FITS files cannot be opened or do not exist

	"""

	#Open the files
	gamma_file = fits.open(filename)

	#Read the ANGLE keyword from the header
	angle = gamma_file[0].header["ANGLE"]

	#Create the array with the shear map
	gamma = gamma_file[0].data.astype(np.float)

	#Close files and return
	gamma_file.close()

	return angle*deg,gamma

################################################################################
##########Default callback loader, loads in the measured power spectrum#########
################################################################################

def measure_power_spectrum(filename,l_edges,columns=None):
	"""
	
	Default ensemble loader: reads a FITS data file containing a convergence map and measures its power spectrum

	:param args: A dictionary that contains all the relevant parameters as keys. Must have a "map_id" key
	:type args: Dictionary

	:returns: ndarray of the measured statistics

	:raises: AssertionError if the input dictionary doesn't have the required keywords

	"""

	logging.debug("Processing {0} power".format(filename))

	conv_map = ConvergenceMap.load(filename,format=load_fits_default_convergence)
	l,Pl = conv_map.powerSpectrum(l_edges)
	return Pl

def peaks_loader(filename,thresholds,columns=None):

	logging.debug("Processing {0} peaks".format(filename))
	conv_map = ConvergenceMap.load(filename,format=load_fits_default_convergence)

	v,pk = conv_map.peakCount(thresholds,norm=True)
	return pk


####################################################################################
#############Default power spectrum template for testing############################
####################################################################################

def sample_power_shape(l,**kwargs):

	return np.exp(-0.5*(l/kwargs["scale"])**2)

#####################################################################################
###############Default loader for 3D matter power spectrum###########################
#####################################################################################

def load_power_default(path,root_name="fiducial_matterpower_"):
	"""
	This is the default matter power spectrum loader: it loads matter power spectra generated with CAMB

	"""
	zregex = re.compile(r"\_([0-9]+)\.dat")
	zfind = lambda f:int(zregex.search(f).groups()[0])

	#Find the number of power spectrum frames (one for each z)
	power_spectrum_files = list(filter(lambda f:zfind(f)>0,glob.glob(os.path.join(path,root_name) + "*.dat")))

	#Redshift array
	z_label = np.zeros(len(power_spectrum_files),dtype=np.int)

	#Check the first file
	k_try,P = np.loadtxt(power_spectrum_files[0],unpack=True)
	z_label[0] = zfind(power_spectrum_files[0])

	#Power spectrum array P(k,z)
	power = np.zeros((len(k_try),len(z_label)))
	power[:,0] = P

	#Load all the other frames
	for n,power_file in enumerate(power_spectrum_files[1:]):
		k,P = np.loadtxt(power_file,unpack=True)
		assert all(k==k_try)
		power[:,n+1] = P
		z_label[n+1] = zfind(power_file)

	#Compute the redshifts
	z = np.array([ label/100.0 for label in z_label ])

	#The redshifts might not be in ascending order, we need to make sure they are
	z_reindex = list(range(len(z)))
	z_reindex.sort(key=z.__getitem__)
	z_reindex = np.array(z_reindex)

	#Return the tuple
	return z[z_reindex],k,power[:,z_reindex]

