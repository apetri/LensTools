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
	def fromfilelist(cls,file_list):

		"""
		Builds the ensemble from a file list: each file corresponds to a different realization

		:param file_list: List of files on which to define the ensemble
		:type file_list: list of str.

		"""

		#See how many realizations are there in the ensemble
		num_realizations = len(file_list)
		assert num_realizations>0,"There are no realizations in your ensemble!!"

		#Build the ensemble instance and return it
		new_ensemble = cls(None,num_realizations)
		setattr(new_ensemble,"file_list",file_list)

		return new_ensemble

	def load(self,callback_loader,pool=None,**kwargs):
		"""
		Loads the ensemble into memory, can spread the calculations on multiple processors using a MPI pool

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble; must take in a dictionary as its only parameter and must return a numpy array
		:type callback_loader: function, must take in a file name (str.) and return a numpy array with the loaded data

		:param pool: MPI pool for multiprocessing (imported from emcee https://github.com/dfm/emcee)
		:type pool: MPI pool object

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: Keyword arguments

		>>> from lenstools import Ensemble
		>>> from lenstools.statistics import default_callback_loader

		>>> map_list = ["conv1.fit","conv2.fit","conv3.fit"]
		>>> l_edges = np.arange(200.0,50000.0,200.0)

		>>> conv_ensemble = Ensemble.fromfilelist(map_list)
		>>> conv_ensemble.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

		"""

		assert callback_loader is not None

		self.pool = pool

		#Build list with tasks to execute (dirty, maybe optimize in the future?)
		tasks = list()
		for file_name in self.file_list:
			task = kwargs.copy()
			task.update({"file_name":file_name})
			tasks.append(task)

		#Execute the callback on each file in the list (spread calculations with MPI pool if it is not none)
		if pool is not None:
			M = pool.map
		else:
			M = map

		full_data = np.array(M(callback_loader,tasks))
		
		assert type(full_data) == np.ndarray
		assert full_data.shape[0] == self.num_realizations

		self.data = full_data

	
	def mean(self):

		"""
		Computes the ensemble average over realizations 

		:returns: ndarray with the averages, has the same shape as self.data[0]

		"""

		return self.data.mean(axis=0) 






