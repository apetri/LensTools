"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence)


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from topology import ConvergenceMap
from shear import ShearMap
from lenstools.index import Indexer

import numpy as np

##########################################
########Ensemble class####################
##########################################

class Ensemble(object):

	"""
	A class that handles statistical operations on weak lensing maps; an ensemble is a collection of different statistical realization of the same random variable. This class has an attribute 'data' that is a numpy array which first axis corresponds to the realization number.

	>>> from lenstools.statistics import Ensemble

	"""

	def __init__(self,file_list=None,data=None,num_realizations=0,metric="chi2"):
		
		self.file_list = file_list
		self.data = data
		self.num_realizations = num_realizations
		self.metric = metric

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
		return cls(file_list=file_list,num_realizations=num_realizations)

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
		for map_id in self.file_list:
			task = kwargs.copy()
			task.update({"map_id":map_id})
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
		
		if not hasattr(self,"_mean"):
			self._mean = self.data.mean(axis=0)
		
		return self._mean

	def covariance(self):

		"""
		Computes the ensemble covariance matrix

		:returns: ndarray with the covariance matrix, has shape (self.data[0],self.data[0]) 

		""" 

		assert self.num_realizations>1, "I can't compute a covariance matrix with one realization only!!"
		assert self.data.dtype == np.float, "This operation is unsafe with non float numbers!!"

		if not hasattr(self,"_mean"):
			self.mean()

		subtracted = self.data - self._mean[np.newaxis,:]
		return np.dot(subtracted.transpose(),subtracted) / (self.num_realizations - 1.0)

	def compare(self,rhs):

		"""
		A comparison operation between ensembles: computes a chi2-style difference between two different ensembles to assert how different they are

		"""

		assert isinstance(rhs,Ensemble)
		assert self.metric == rhs.metric,"The two ensemble instances must have the same metric!!"

		if self.metric=="chi2":

			mean1 = self.mean()
			covariance = self.covariance()
			mean2 = rhs.mean()

			return np.dot(mean1 - mean2,np.dot(np.linalg.inv(covariance),mean1 - mean2))

		else:

			raise ValueError("Only chi2 metric implemented so far!!")

	
	def __add__(self,rhs):

		"""
		Overload of the sum operator: a sum between ensembles is another ensemble whose data are vstacked. This operation is useful if you have two ensembles with different, independent realizations and you want to create an ensemble which is the union of the two. The data must be vstackable

		"""

		#Safety checks
		assert isinstance(rhs,Ensemble)
		assert self.metric == rhs.metric,"The two ensemble instances must have the same metric!!"

		try:
			
			new_data = np.vstack((self.data,rhs.data))
		
		except ValueError:
			
			print("The data of these two ensembles cannot be vstacked!")
			return None

		return Ensemble(file_list=self.file_list+rhs.file_list,data=new_data,num_realizations=self.num_realizations+rhs.num_realizations,metric=self.metric)


	def __mul__(self,rhs):

		"""
		Overload of the multiplication operator: a multiplication of two ensembles is another ensemble whose data are hstacked. This operation is useful if you have two ensembles with the same realizations but different statistical descriptor, and you want to create an ensemble in which both descriptors are included. The data must be hstackable

		"""

		#Safety checks
		assert isinstance(rhs,Ensemble)
		assert self.metric == rhs.metric
		assert self.num_realizations == rhs.num_realizations

		#Throw an error if the file lists do not coincide
		if not self.file_list == rhs.file_list:
			raise ValueError("You are try to multiply ensembles built from different files!!")

		try:

			new_data = np.hstack((self.data,rhs.data))

		except ValueError:

			print("The data of these two ensembles cannot be hstacked!")
			return None

		return Ensemble(file_list=self.file_list,data=new_data,num_realizations=self.num_realizations,metric=self.metric)

	def split(self,index):

		"""
		Inverse of the * operator: this method uses an Indexer instance to break down a multiple descriptor ensemble in many, smaller, single descriptor ensembles

		:param index: index of descriptors with which to perform the split
		:type index: Indexer instance

		:returns: list of Ensemble instances, one for each element in index

		:raises: AssertionError if shape of the ensemble data is not suitable 

		"""

		assert isinstance(index,Indexer)

		splitted = list()
		
		for n in range(index.num_descriptors):

			splitted.append(Ensemble(file_list=self.file_list,num_realizations=self.num_realizations,data=self.data[:,index[n].first:index[n].last],metric=self.metric))

		return splitted










