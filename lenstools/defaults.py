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

from convergence import ConvergenceMap
from index import PowerSpectrum,Peaks,PDF,MinkowskiAll,MinkowskiSingle,Moments

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

def default_callback_loader(filename,l_edges):
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

def peaks_loader(filename,thresholds):

	logging.debug("Processing {0} peaks".format(filename))
	conv_map = ConvergenceMap.load(filename,format=load_fits_default_convergence)

	v,pk = conv_map.peakCount(thresholds,norm=True)
	return v

def convergence_measure_all(filename,index,fits_loader=None):

	"""
	Measures all the statistical descriptors of a convergence map as indicated by the index instance
	
	"""

	logging.debug("Processing {0}".format(filename))

	#Load the map
	if fits_loader is not None:
		conv_map = ConvergenceMap.load(filename,format=fits_loader)
	else: 
		conv_map = ConvergenceMap.load(filename,format=load_fits_default_convergence)

	#Allocate memory for observables
	descriptors = index
	observables = np.zeros(descriptors.size)

	#Measure descriptors as directed by input
	for n in range(descriptors.num_descriptors):

		
		if type(descriptors[n]) == PowerSpectrum:
			
			l,observables[descriptors[n].first:descriptors[n].last] = conv_map.powerSpectrum(descriptors[n].l_edges)

		elif type(descriptors[n]) == Moments:

			observables[descriptors[n].first:descriptors[n].last] = conv_map.moments(connected=descriptors[n].connected)
		
		elif type(descriptors[n]) == Peaks:
			
			v,observables[descriptors[n].first:descriptors[n].last] = conv_map.peakCount(descriptors[n].thresholds,norm=descriptors[n].norm)

		elif type(descriptors[n]) == PDF:

			v,observables[descriptors[n].first:descriptors[n].last] = conv_map.pdf(descriptors[n].thresholds,norm=descriptors[n].norm)
		
		elif type(descriptors[n]) == MinkowskiAll:
			
			v,V0,V1,V2 = conv_map.minkowskiFunctionals(descriptors[n].thresholds,norm=descriptors[n].norm)
			observables[descriptors[n].first:descriptors[n].last] = np.hstack((V0,V1,V2))
		
		elif type(descriptors[n]) == MinkowskiSingle:
			
			raise ValueError("Due to computational performance you have to measure all Minkowski functionals at once!")
		
		else:
			
			raise ValueError("Measurement of this descriptor not implemented!!!")

	#Return
	return observables


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

	#Find the number of power spectrum frames (one for each z)
	power_spectrum_files = glob.glob(os.path.join(path,root_name) + "*.dat")[1:]

	#Redshift array
	z_label = np.zeros(len(power_spectrum_files),dtype=np.int)

	#Check the first file
	k_try,P = np.loadtxt(power_spectrum_files[0],unpack=True)
	z_label[0] = int(re.match(os.path.join(path,root_name)+r"([0-9]+).dat",power_spectrum_files[0]).group(1))

	#Power spectrum array P(k,z)
	power = np.zeros((len(k_try),len(z_label)))
	power[:,0] = P

	#Load all the other frames
	for n,power_file in enumerate(power_spectrum_files[1:]):
		k,P = np.loadtxt(power_file,unpack=True)
		assert all(k==k_try)
		power[:,n+1] = P
		z_label[n+1] = int(re.match(os.path.join(path,root_name)+r"([0-9]+).dat",power_file).group(1))

	#Compute the redshifts
	z = np.array([ label/100.0 for label in z_label ])

	#The redshifts might not be in ascending order, we need to make sure they are
	z_reindex = range(len(z))
	z_reindex.sort(key=z.__getitem__)
	z_reindex = np.array(z_reindex)

	#Return the tuple
	return z[z_reindex],k,power[:,z_reindex]