"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence), wrapped around a pandas DataFrame


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from operator import add
from functools import reduce
import cPickle as pickle

from ..utils import pcaHandler as PCA

import numpy as np
import scipy.io as sio
from scipy import sparse

from emcee.ensemble import _function_wrapper

try:
	import matplotlib
	import matplotlib.pyplot as plt
except ImportError:
	matplotlib = None

try:
	import pandas as pd
	pd = pd
except ImportError:
	pd = None


##############################################################
#################Series class#################################
##############################################################

class Series(pd.Series):

	@property
	def _constructor(self):
		return Series

	@property
	def _constructor_expanddim(self):
		return Ensemble

	@staticmethod
	def make_index(*indices):

		"""
		Merge a sequence of simple indices in a single multi-index

		"""

		values = np.concatenate([idx.values for idx in indices])
		names = reduce(add,[ [idx.name]*len(idx) for idx in indices ])
		return pd.MultiIndex.from_tuples(zip(names,values))

##########################################
########Ensemble class####################
##########################################

class Ensemble(pd.DataFrame):

	"""
	A class that handles statistical operations on weak lensing maps, inherits from pandas DataFrame. The rows in the Ensemble correspond to the different ensemble realizations of the same descriptor

	"""

	################################################################
	##############DataFrame subclassing#############################
	################################################################

	_metadata = ["file_list","metric"]

	@property
	def _constructor(self):
		return self.__class__

	@property
	def _constructor_sliced(self):
		return Series

	@property
	def _constructor_expanddim(self):
		return Panel

	##################################
	########Constructor###############
	##################################

	def __init__(self,data=None,file_list=list(),metric="chi2",**kwargs):

		#Call parent constructor
		super(Ensemble,self).__init__(data=data,**kwargs)
		
		#Additional attributes
		self.file_list = file_list
		self.metric = metric

	##################################
	########Serialization#############
	##################################

	def __getstate__(self):
		meta = dict((k,getattr(self,k,None)) for k in self._metadata)
		return dict(_data=self._data,_typ=self._typ,_metadata=self._metadata,**meta)

	##################
	####Properties####
	##################

	@property
	def nobs(self):
		return self.shape[0]

	###############################
	########Grouping###############
	###############################

	def groupby(self,*args,**kwargs):

		"""
		Same as pandas DataFrame groupby
		"""

		g = super(Ensemble,self).groupby(*args,**kwargs)
		g.__class__ = EnsembleGroupBy
		g._series_constructor = self._constructor_sliced
		g._ensemble_constructor = self.__class__

		return g

	####################################
	#############I/O####################
	####################################

	@classmethod
	def read_pickle(cls,filename):

		"""
		Reads in a pickled Ensemble from disk

		:param filename: name of the file to read
		:type filename: str.

		:returns: Ensemble instance read from the file

		"""

		with open(filename,"r") as fp:
			ens = pickle.load(fp)

		return cls(ens)

	@classmethod
	def read(cls,filename,callback_loader=None,**kwargs):

		"""
		Reads a numpy file into an Ensemble

		:param filename: name of the file to read
		:type filename: str.

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble. If None provided, it unpickles the specified file. Must return an acceptable input for the Ensemble constructor
		:type callback_loader: function

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: dict.

		:returns: Ensemble instance read from the file

		"""

		if callback_loader is None:
			callback_loader = cls.read_pickle

		return cls(callback_loader(filename,**kwargs),file_list=[filename])

	@classmethod
	def readall(cls,filelist,callback_loader=None,**kwargs):

		"""
		Reads a list of files into an Ensemble

		:param filelist: list of files to read
		:type filelist: list.

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble. If None provided, it performs a numpy.load on the specified file. Must return a numpy array with the loaded data
		:type callback_loader: function

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: dict.

		:returns: Ensemble instance read from the file

		"""

		return cls.concat([ cls.read(f,callback_loader,**kwargs) for f in filelist ],axis=0,ignore_index=True)

	@classmethod
	def compute(cls,file_list,callback_loader=None,pool=None,index=None,assemble=np.array,**kwargs):
		
		"""
		Computes an ensemble, can spread the calculations on multiple processors using a MPI pool

		:param file_list: list of files that will constitute the ensemble; the callback_loader is called on each of the files to produce the different realizations
		:type file_list: list. 

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble. If None provided, it performs a numpy.load on the specified files. Must return a numpy array with the loaded data
		:type callback_loader: function

		:param pool: MPI pool for multiprocessing (imported from emcee https://github.com/dfm/emcee)
		:type pool: MPI pool object

		:param index: index of the Ensemble
		:type index: pandas Index

		:param assemble: called on the list of features (one feature per file) to assemble them in an array (defaults to np.array)
		:type assemble: callable

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: dict.

		>>> from lenstools import Ensemble
		>>> from lenstools.statistics import default_callback_loader

		>>> map_list = ["conv1.fit","conv2.fit","conv3.fit"]
		>>> l_edges = np.arange(200.0,50000.0,200.0)

		>>> conv_ensemble = Ensemble.compute(map_list,callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

		"""

		#Safety checks
		assert callback_loader is not None, "You must specify a callback loader function that returns a numpy array!"
		if index is not None:
			assert len(index)==len(file_list),"The number of elements in the index hould be the same as the number of files!"

		#Build a function wrapper of the callback loader, so it becomes pickleable
		_callback_wrapper = _function_wrapper(callback_loader,args=tuple(),kwargs=kwargs)

		#Execute the callback on each file in the list (spread calculations with MPI pool if it is not none)
		if pool is not None:
			M = pool.map
		else:
			M = map

		full_data = assemble(M(_callback_wrapper,file_list))
		
		assert type(full_data) == np.ndarray

		#Check if user provided column labels
		if "columns" in kwargs.keys():
			columns = kwargs["columns"]
		else:
			columns = None

		#Return the created ensemble from the full_data array
		return cls(full_data,file_list,index=index,columns=columns)

	
	def save(self,filename,format=None,**kwargs):
		
		"""
		Save ensemble data in an external file (in arbitrary format)

		:param filename: file name of the external file
		:type filename: str.

		:format: format in which to save the ensemble; if None the format is auto detected from the filename
		:type format: str.or callable

		:param kwargs: the keyword arguments are passed to the saver (or to format if callable)
		:type kwargs: dict.

		"""

		#Auto detect format
		if format is None:
			if filename.endswith(".npy"):
				format = "numpy"
			elif filename.endswith(".mat"):
				format = "matlab"
			elif filename.endswith(".pkl"):
				format = "pickle"
			else:
				raise ValueError("Format not recognized!")


		#Proceed to the saving procedure
		if format=="numpy":
			np.save(filename,self.values)
		elif format=="matlab":
			sio.savemat(filename,{"values": self.values},**kwargs)
		elif format=="pickle":
			with open(filename,"w") as fp:
				pickle.dump(self,fp)
		else:
			format(self,filename,**kwargs)
	

	####################################
	#############Operations#############
	####################################

	@classmethod 
	def concat(cls,ensemble_list,**kwargs):
		return pd.concat(ensemble_list,**kwargs)


	def group(self,group_size,kind="sparse"):

		"""
		Sometimes it happens that different realizations in the ensemble need to be grouped together, for example when they belong to different subfields of the same observation field. With this function you can group different realizations together by taking the mean, and reduce the total number of realizations in the ensemble

		:param group_size: how many realizations to put in a group, must divide exactly the total number of realizations in the ensemble
		:type group_size: int

		:param kind: specifies how to do the grouping; if set to "sparse" the groups are formed by taking one realizations every nobs/group_size (for example ([1,1001,...,9001],[2,1002,...,9002]) if nobs=10000 and group_size=10). If set to "contiguous" then the realizations are grouped as ([1,2,...,10],[11,12,...,20]). Otherwise you can set kind to your own sparse matrix scheme 
		:type kind: str. 

		:returns: gropby object

		"""

		assert not(self.nobs%group_size),"The group size must divide exactly the number of realizations!"

		#Build the appropriate grouping scheme
		if kind=="sparse":
			return self.groupby(lambda i:i%group_size)
		elif kind=="contiguous":
			return self.groupby(lambda i:i%group_size)
		else:
			raise NotImplementedError("Grouping scheme '{0}' not implemented".format(kind))


	#TODO: reimplement bootstraping method
	def covariance(self,bootstrap=False,**kwargs):

		"""
		Computes the ensemble covariance matrix

		:param bootstrap: if True the covariance matrix is computed with a bootstrap estimate
		:type bootstrap: bool.

		:param kwargs: the keyword arguments are passed to the bootstrap method
		:type kwargs: dict.

		:returns: ndarray with the covariance matrix, has shape (self.data[1],self.data[1]) 

		""" 

		if bootstrap:
			raise NotImplementedError
		else:
			return self.cov()


	#TODO: reimplement
	def bootstrap(self,callback,bootstrap_size=10,resample=10,seed=None):

		"""
		Computes a custom statistic on the Ensemble using the bootstrap method

		:param callback: statistic to compute on the ensemble; takes the resampled Ensemble data as an input
		:type callback: callable

		:param bootstrap_size: size of the resampled ensembles used in the bootstraping; must be less than or equal to the number of realizations in the Ensemble
		:type bootstrap_size: int.

		:param resample: number of times the Ensemble is resampled
		:type resample: int.

		:param seed: if not None, this is the random seed of the random resamples 
		:type seed: int.

		:returns: the bootstraped Ensemble statistic

		"""

		#Safety check
		assert bootstrap_size<=self.nobs,"The size of the resampling cannot exceed the original number of realizations"

		#Set the random seed
		if seed is not None:
			np.random.seed(seed)

		#TODO: Parallelize
		M = map

		#Construct the randomization matrix
		randomizer = np.random.randint(self.nobs,size=(resample,bootstrap_size))

		#Compute the statistic with the callback
		statistic = np.array(M(callback,self.iloc[randomizer]))

		#Return the bootstraped statistic expectation value
		return statistic.mean(0)


	def principalComponents(self):

		"""
		Computes the principal components of the Ensemble

		:returns: pcaHandler instance

		"""

		pca = PCA(constructor_series=self._constructor_sliced,constructor_ensemble=self.__class__,columns=self.columns)
		pca.fit(self.values)
		return pca


	def compare(self,rhs,**kwargs):

		"""
		Computes a comparison score between two Ensembles: computes a chi2-style difference between two different ensembles to assert how different they are

		"""

		assert isinstance(rhs,self.__class__)
		assert self.metric == rhs.metric,"The two ensemble instances must have the same metric!!"

		if self.metric=="chi2":

			mean1 = self.mean(0).values
			
			if "covariance" in kwargs.keys():
				covariance = kwargs["covariance"]
			else:
				covariance = self.covariance().values
			
			mean2 = rhs.mean(0).values

			return np.dot(mean1 - mean2,np.dot(np.linalg.inv(covariance),mean1 - mean2))

		else:

			raise NotImplementedError("Only chi2 metric implemented so far!!")

	def selfChi2(self):

		"""
		Computes the Ensemble distribution of chi squared values, defined with respect to the same Ensemble mean and covariance

		:returns: array with the self chi squared for each realization

		"""

		#Compute mean and covariance
		mean = self.mean(0).values
		covariance = self.covariance().values
		difference = self.values - mean[None]

		#Compute the chi2 for each realization
		return Series((difference * np.linalg.solve(covariance,difference.T).T).sum(-1))


	def shuffle(self,seed=None):

		"""
		Changes the order of the realizations in the Ensemble

		:param seed: random seed for the random shuffling
		:type seed: int.

		:returns: shuffled Ensemble

		"""

		if seed is not None:
			np.random.seed(seed)

		#Shuffle
		return self.reindex(self.index[np.random.permutation(np.arange(len(self)))])


	#####################################################################################################################

	####################################
	#############Visualization##########
	####################################

	def imshow(self,fig=None,ax=None,**kwargs):

		"""
		Visualize a two dimensional map of the Ensemble, with the index as the vertical axis and the columns as the horizontal axis

		:param fig:
		:type fig:

		:param ax:
		:type ax:

		:param kwargs:
		:type kwargs: dict.

		"""

		#Matplotlib needs to be installed
		if matplotlib is None:
			raise ImportError("matplotlib needs to be installed to use this method!")

		#Create figure if one does not exist yet
		if fig is None or ax is None:
			self.fig,self.ax = plt.subplots()
		else:
			self.fig = fig
			self.ax = ax

		#Show the ensemble
		img = plt.imshow(self.values,**kwargs)

		#Set the ticks and ticklabels
		self.ax.set_yticks(range(self.shape[0]))
		self.ax.set_xticks(range(self.shape[1]))
		self.ax.set_yticklabels(self.index)
		self.ax.set_xticklabels(self.columns,rotation="vertical")

		#Maybe the axis labels if index and columns have names
		if self.index.name is not None:
			self.ax.set_ylabel(self.index.name)

		if self.columns.name is not None:
			self.ax.set_xlabel(self.columns.name)

		#Return the handle
		return self.ax


#######################################
########Panel class####################
#######################################

class Panel(pd.Panel):

	@property 
	def _constructor(self):
		return Panel

	@property
	def _constructor_sliced(self):
		return Ensemble


###################################################################################################################################################################################################


###########################
###SeriesGroupBy class#####
###########################

class SeriesGroupBy(pd.core.groupby.SeriesGroupBy):

	"""
	Useful for grouping Series

	"""

	def mean(self,*args,**kwargs):
		return self._series_constructor(super(SeriesGroupBy,self).mean(*args,**kwargs))

###########################
###EnsembleGroupBy class###
###########################

class EnsembleGroupBy(pd.core.groupby.DataFrameGroupBy):

	"""
	Useful for grouping Ensemble

	"""

	def mean(self,*args,**kwargs):
		return self._ensemble_constructor(super(EnsembleGroupBy,self).mean(*args,**kwargs))

	def apply(self,*args,**kwargs):
		return self._ensemble_constructor(super(EnsembleGroupBy,self).apply(*args,**kwargs))

	def __getitem__(self,name):
		item = super(EnsembleGroupBy,self).__getitem__(name)
		item.__class__ = SeriesGroupBy
		item._series_constructor = self._series_constructor
		return item





