"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence)


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from topology import ConvergenceMap
from shear import ShearMap

import numpy as np

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

		>>> from lenstools import Ensemble
		>>> from lenstools.statistics import default_callback_loader

		>>> map_list = ["conv1.fit","conv2.fit","conv3.fit"]
		>>> l_edges = np.arange(200.0,50000.0,200.0)

		>>> conv_ensemble = Ensemble.fromfilelist(map_list,callback_loader=default_callback_loader,l_edges=l_edges)

		"""

		#############################################
		###Wrapper for MPI pool, to remain private###
		#############################################

		def _callback_wrapper(obj):
			file_name,kwargs = obj
			return callback_loader(file_name,**kwargs)

		##############################################
		##############################################

		assert callback_loader is not None

		#See how many realizations are there in the ensemble
		num_realizations = len(file_list)
		assert num_realizations>0,"There are no realizations in your ensemble!!"

		#Build list with tasks to execute
		tasks = [ (file_name,kwargs) for file_name in file_list ]

		#Execute the callback on each file in the list
		full_data = np.array(map(_callback_wrapper,tasks))
		assert type(full_data) == np.ndarray
		assert full_data.shape[0] == num_realizations

		#Build the ensemble instance and return it
		return cls(full_data,num_realizations)

	
	def mean(self):

		"""
		Computes the ensemble average over realizations 

		:returns: ndarray with the averages, has the same shape as self.data[0]

		"""

		return self.data.mean(axis=0) 






