"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence)


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from topology import ConvergenceMap,load_fits_default
from shear import ShearMap

import numpy as np
from astropy.io import fits

################################################################################
##########Default callback loader, loads in the measured power spectrum#########
################################################################################

def default_callback_loader(file_name,**kwargs):

	print("Processing {0}".format(file_name))

	conv_map = ConvergenceMap.fromfilename(file_name,loader=load_fits_default)
	l,Pl = conv_map.powerSpectrum(kwargs["l_edges"])
	return Pl

##########################################
########Ensemble class####################
##########################################

class Ensemble(object):

	"""
	A class that handles statistical operations on weak lensing maps; an ensemble is a collection of different statistical realization of the same random variable. This class has an attribute 'data' that is a numpy array which first axis corresponds to the realization number.

	>>> from lenstools.statistics import Ensemble

	"""

	def __init__(self,data,num_realizations):
		
		self.data = data
		self.num_realizations = num_realizations

	@classmethod
	def fromfilelist(cls,file_list,callback_loader,pool=None,**kwargs):

		"""
		Builds the ensemble from a file list: each file corresponds to a different realization. The callback_loader parameter allows flexibility on how to popolate the ensemble

		:param file_list: List of files on which to define the ensemble
		:type file_list: list of str.

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble
		:type callback_loader: function, must take in a file name (str.) and return a numpy array with the loaded data

		:param pool: MPI pool for multiprocessing, not functional yet
		:type pool: MPI pool object

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: Keyword arguments

		"""

		#See how many realizations are there in the ensemble
		num_realizations = len(file_list)
		assert num_realizations>0,"There are no realizations in your ensemble!!"

		#Call the loader on the first file to assert the data shape
		data_first = callback_loader(file_list[0],**kwargs)
		assert type(data_first) == np.ndarray,"There is something wrong with your callback, it doesn't return numpy arrays!"

		#Allocate memory for the ensemble data and fill with the first element
		full_data_shape = (num_realizations,) + data_first.shape
		full_data = np.zeros(full_data_shape)
		full_data[0] = data_first

		#Cycle to the remaining files to fill in the ensemble data
		for n,file_name in enumerate(file_list[1:]):

			#Load
			data = callback_loader(file_name,**kwargs)

			#Safety check
			assert type(data) == np.ndarray,"There is something wrong with your callback, it doesn't return numpy arrays!"
			assert data.shape == data_first.shape,"All the loaded data arrays must have the same shape!!"

			#Add the loaded data to the ensemble
			full_data[n + 1] = data

		#Build the ensemble instance and return it
		return cls(full_data,num_realizations)  





