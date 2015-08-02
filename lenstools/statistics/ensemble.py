"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence), wrapped around a pandas DataFrame


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from operator import add
from functools import reduce

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

##########################################################################
#########Useful for bootstrap estimates of the covariance matrix##########
##########################################################################

def _bsp_covariance(data):
	sub = data - data.mean(0)[None]
	return np.dot(sub.T,sub) / (data.shape[0] - 1.0)

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

	_metadata = ["num_realizations","file_list","metric"]

	@property
	def _constructor(self):
		return Ensemble

	@property
	def _constructor_sliced(self):
		return Series

	##################################
	########Constructor###############
	##################################

	def __init__(self,data=None,file_list=list(),metric="chi2",**kwargs):

		#Call parent constructor
		super(Ensemble,self).__init__(data=data,**kwargs)
		
		#Additional attributes
		if data is not None:
			self.num_realizations = data.shape[0]
		else:
			self.num_realizations = 0 

		self.file_list = file_list
		self.metric = metric

	####################################
	#############I/O####################
	####################################

	@classmethod
	def read(cls,filename,callback_loader=None,**kwargs):

		"""
		Reads a numpy file into an Ensemble

		:param filename: name of the file to read
		:type filename: str.

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble. If None provided, it performs a numpy.load on the specified file. Must return a numpy array with the loaded data
		:type callback_loader: function

		:param kwargs: Any additional keyword arguments to be passed to callback_loader
		:type kwargs: dict.

		:returns: Ensemble instance read from the file

		"""

		if callback_loader is None:
			callback_loader = lambda f: np.load(f)

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

		return cls.concat([ cls.read(f,callback_loader,**kwargs) for f in filelist ])

	@classmethod
	def compute(cls,file_list,callback_loader=None,pool=None,index=None,**kwargs):
		
		"""
		Computes an ensemble, can spread the calculations on multiple processors using a MPI pool

		:param file_list: list of files that will constitute the ensemble; the callback_loader is called on each of the files to produce the different realizations
		:type file_list: list. 

		:param callback_loader: This function gets executed on each of the files in the list and populates the ensemble. If None provided, it performs a numpy.load on the specified files. Must return a numpy array with the loaded data
		:type callback_loader: function

		:param pool: MPI pool for multiprocessing (imported from emcee https://github.com/dfm/emcee)
		:type pool: MPI pool object

		:param from_old: If True, the loaded data are interpreted as an old, already existing ensemble, which means that only one file (in which the old ensemble is saved) is loaded, the first dimension of the data is 1 and hence it is discarded 
		:type from_old: bool.

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

		full_data = np.array(M(_callback_wrapper,file_list))
		
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
			else:
				raise ValueError("Format not recognized!")


		#Proceed to the saving procedure
		if format=="numpy":
			np.save(filename,self.values)
		elif format=="matlab":
			sio.savemat(filename,{"values": self.values},**kwargs)
		else:
			format(self,filename,**kwargs)
	

	####################################
	#############Operations#############
	####################################

	@classmethod 
	def concat(cls,ensemble_list):

		data = np.vstack([ ens.values for ens in ensemble_list ])
		return cls(data)
	
	#TODO: remove
	def mean(self):

		"""
		Computes the ensemble average over realizations 

		:returns: ndarray with the averages, has the same shape as self.data[0]

		"""
		
		self._mean = self.data.mean(axis=0)
		
		return self._mean

	#TODO: reimplement
	def scale(self,weights):

		"""
		Set manually the units of each row of the ensemble by multiplying it by a row of weights

		:param weights: row of weights used to scale the ensemble, must have the same length as the number of rows
		:type weights: array

		"""

		if isinstance(weights,np.ndarray):
			
			assert len(weights)==self.data.shape[0]
			self.data *= weights.reshape(self.data.shape)

		else:

			self.data *= weights

	#TODO: reimplement
	def group(self,group_size,kind="sparse",inplace=True):

		"""
		Sometimes it happens that different realizations in the ensemble need to be grouped together, for example when they belong to different subfields of the same observation field. With this function you can group different realizations together by taking the mean, and reduce the total number of realizations in the ensemble

		:param group_size: how many realizations to put in a group, must divide exactly the total number of realizations in the ensemble
		:type group_size: int

		:param kind: specifies how to do the grouping; if set to "sparse" the groups are formed by taking one realizations every num_realizations/group_size (for example ([1,1001,...,9001],[2,1002,...,9002]) if num_realizations=10000 and group_size=10). If set to "contiguous" then the realizations are grouped as ([1,2,...,10],[11,12,...,20]). Otherwise you can set kind to your own sparse matrix scheme 
		:type kind: str. or scipy.sparse 

		"""

		assert not(self.num_realizations%group_size),"The group size must divide exactly the number of realizations!"
		num_groups = self.num_realizations//group_size

		#Build the appropriate grouping scheme
		if isinstance(kind,sparse.csr.csr_matrix):

			assert kind.dtype == np.bool,"The scheme type must be bool, only 0 and 1!"
			scheme = kind
		
		elif kind=="sparse":

			scheme = reduce(add,[ sparse.eye(num_groups,self.num_realizations,k=num_groups*i,dtype=np.int8) for i in range(group_size)])

		elif kind=="contiguous":
			
			row = np.array(reduce(add,[ (i,) * group_size  for i in range(num_groups) ]))
			col = np.arange(self.num_realizations,dtype=np.int)
			dat = np.ones(self.num_realizations,dtype=np.int8)

			scheme = sparse.csr_matrix((dat,(row,col)),shape=(num_groups,self.num_realizations),dtype=np.int8)

		else:
			raise TypeError("The scheme kind you inputed is not valid!")

		self._scheme = scheme

		if inplace:

			#Dot the data with the scheme to produce the groups, and take the mean in every group
			self.num_realizations = num_groups
			self.data = scheme.dot(self.data) / group_size

		else:

			#Dot the data with the scheme to produce the groups, and take the mean in every group
			return self.__class__.fromdata(scheme.dot(self.data) / group_size)


	#TODO: remove
	def differentiate(self,step=None,order=1):

		"""
		Compute a new ensemble, in which the feature is differentiated with respect to its label (e.g. the differentiation of the first minkowski functional is the PDF)

		:param step: measure unit on the feature label axis 
		:type step: float or array

		:returns: new Ensemble instance with the differentiated feature 

		"""

		#differentiate
		diff_data = np.diff(self.data,n=order,axis=1)

		#Optionally scale the measure units
		if step is not None:
			diff_data /= step

		#return the new ensemble
		return self.__class__(data=diff_data,num_realizations=diff_data.shape[0])


	#TODO: remove
	def cut(self,min=None,max=None,feature_label=None,inplace=True):

		"""
		Allows to manually cut the ensemble along the second axis, if you want to select a feature subset; you better know what you are doing if you are using this function

		:param min: left extreme of the cut, included; if a list of indices is passed, the cut is performed on those indices, on the second axis and the remaining parameters are ignored
		:type min: int or float

		:param max: right extreme of the cut, included
		:type max: int or float

		:param feature_label: if not None, serves as a reference for min and max, in which the ensemble is cut according to the position of the min and max elements in the feature_label array
		:type feature_label: array

		:param inplace: if True, the Ensemble cut is made in place, otherwise a new ensemble is returned
		:type inplace: bool.

		:returns: the min and max cut indices if feature_label is None, otherwise it returns the array of the cut feature; if min was a list of indices, these are returned back. Returns the new Ensemble if inplace is False

		"""

		#Sanity checks
		assert self.data.ndim==2,"Only one dimensional feature cuts implemented so far!"
		assert (min is not None) or (max is not None),"No cutting extremes selected!"

		if type(min)==list:
		
			new_data = self.data[:,min]

			if inplace:

				self.data = new_data

				#Return
				return min

			else:
				return self.__class__.fromdata(new_data)

		else:
		
			if feature_label is not None:

				#Look for the corresponding indices in the feature label
				min_idx = np.abs(feature_label-min).argmin()
				max_idx = np.abs(feature_label-max).argmin()

			else:

				#Interpret max and min as indices between which to perform the cut
				min_idx = min
				max_idx = max

			new_data = self.data[:,min_idx:max_idx+1]

			if inplace:

				self.data = new_data

				#Return
				if feature_label is not None:
					return feature_label[min_idx:max_idx+1]
				else:
					return min_idx,max_idx

			else:
				return self.__class__.fromdata(new_data)


	#TODO: remove
	def subset(self,realizations):

		"""
		Returns a new ensemble that contains only the selected realizations

		:param realizations: realizations to keep in the subset Ensemble, must be compatible with numpy array indexing syntax
		:type realizations: int. or array of int.

		:returns: new Ensemble instance that contains only the selected realizations 

		"""

		data = self.data[realizations]
		
		return self.__class__.fromdata(data)


	#TODO: remove
	def transform(self,transformation,inplace=False,**kwargs):

		"""
		Allows a general transformation on the Ensemble by calling a function on its data

		:param transformation: callback function called on the ensemble data
		:type transformation: callable 

		:param inplace: if True the transformation is performed in place, otherwise a new Ensemble is created
		:type inplace: bool.

		:param kwargs: the keyword arguments are passed to the transformation callable
		:type kwargs: dict.

		:returns: new Ensemble instance if the transformation is nor performed in place

		"""

		#Apply the transformation
		transformed_data = transformation(self.data,**kwargs)

		#Return the new Ensemble 
		if inplace:
			
			self.data = transformed_data
			self.num_realizations = transformed_data.shape[0]

		else:
			return self.__class__.fromdata(transformed_data)


	#TODO: reimplement
	def covariance(self,bootstrap=False,**kwargs):

		"""
		Computes the ensemble covariance matrix

		:param bootstrap: if True the covariance matrix is computed with a bootstrap estimate
		:type bootstrap: bool.

		:param kwargs: the keyword arguments are passed to the bootstrap method
		:type kwargs: dict.

		:returns: ndarray with the covariance matrix, has shape (self.data[1],self.data[1]) 

		""" 

		assert self.num_realizations>1, "I can't compute a covariance matrix with one realization only!!"
		assert self.data.dtype == np.float, "This operation is unsafe with non float numbers!!"

		if bootstrap:
			return self.bootstrap(_bsp_covariance,**kwargs)
		
		else:
			if not hasattr(self,"_mean"):
				self.mean()

			subtracted = self.data - self._mean[np.newaxis,:]
			return np.dot(subtracted.transpose(),subtracted) / (self.num_realizations - 1.0)

	#TODO: reimplement just calling self.corr()
	def correlation(self):

		"""
		Computes the ensemble correlation matrix

		:returns: ndarray with the correlation matrix, has shape (num_realizations,num_realizations)

		"""

		assert self.data.dtype == np.float, "This operation is unsafe with non float numbers!!"
		if self.num_realizations==1:
			return np.ones((1,1))
		else:	

			if not hasattr(self,"_mean"):
				self.mean()

			subtracted = self.data - self._mean[np.newaxis,:]
			std = np.sqrt((subtracted**2).sum(-1))

			return np.dot(subtracted,subtracted.T) / np.outer(std,std)

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
		assert bootstrap_size<=self.num_realizations,"The size of the resampling cannot exceed the original number of realizations"

		#Set the random seed
		if seed is not None:
			np.random.seed(seed)

		#TODO: Parallelize
		M = map

		#Construct the randomization matrix
		randomizer = np.random.randint(self.num_realizations,size=(resample,bootstrap_size))

		#Compute the statistic with the callback
		statistic = np.array(M(callback,self.iloc[randomizer]))

		#Return the bootstraped statistic expectation value
		return statistic.mean(0)


	#TODO: reimplement adding the proper constructor to PCA
	def principalComponents(self):

		"""
		Computes the principal components of the Ensemble

		:returns: pcaHandler instance

		"""

		pca = PCA()
		pca.fit(self.values)
		return pca

	

	def compare(self,rhs,**kwargs):

		"""
		Computes a comparison score between two Ensembles: computes a chi2-style difference between two different ensembles to assert how different they are

		"""

		assert isinstance(rhs,self.__class__)
		assert self.metric == rhs.metric,"The two ensemble instances must have the same metric!!"

		if self.metric=="chi2":

			mean1 = self.mean()
			
			if "covariance" in kwargs.keys():
				covariance = kwargs["covariance"]
			else:
				covariance = self.covariance()
			
			mean2 = rhs.mean()

			return np.dot(mean1 - mean2,np.dot(np.linalg.inv(covariance),mean1 - mean2))

		else:

			raise NotImplementedError("Only chi2 metric implemented so far!!")

	#TODO: reimplement
	def selfChi2(self):

		"""
		Computes the Ensemble distribution of chi squared values, defined with respect to the same Ensemble mean and covariance

		:returns: array with the self chi squared for each realization

		"""

		#Compute mean and covariance
		mean = self.mean()
		covariance = self.covariance()
		difference = self.data - mean[None]

		#Compute the chi2 for each realization
		return (difference * np.linalg.solve(covariance,difference.T).T).sum(-1)


	#TODO: remove: the reindex() method does the same thing
	def shuffle(self,seed=None):

		"""
		Changes the order of the realizations in the Ensemble

		:param seed: random seed for the random shuffling
		:type seed: int.

		"""

		if seed is not None:
			np.random.seed(seed)

		#Shuffle
		np.random.shuffle(self.data)


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









