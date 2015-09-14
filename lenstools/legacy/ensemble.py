"""

.. module:: statistics 
	:platform: Unix
	:synopsis: This module implements a set of statistical operations on ensembles of weak lensing maps (shear/convergence)


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from operator import add
from functools import reduce

from .index import Indexer
from ..utils import pcaHandler as PCA

import numpy as np
import scipy.io as sio
from scipy import sparse

from emcee.ensemble import _function_wrapper

try:
	import pandas as pd
	pd = pd
except ImportError:
	pd = None

##########################################################################
#########Useful for bootstrap estimates of the covariance matrix###########
##########################################################################

def _bsp_covariance(data):
	sub = data - data.mean(0)[None]
	return np.dot(sub.T,sub) / (data.shape[0] - 1.0)


##########################################################################
#########Useful for bootstrap estimates of the covariance matrix###########
##########################################################################

def _bsp_covariance(data):
	sub = data - data.mean(0)[None]
	return np.dot(sub.T,sub) / (data.shape[0] - 1.0)


##########################################
########Ensemble class####################
##########################################

class Ensemble(object):

	"""
	A class that handles statistical operations on weak lensing maps; an ensemble is a collection of different statistical realization of the same random variable. This class has an attribute 'data' that is a numpy array which first axis corresponds to the realization number.

	>>> from lenstools.statistics import Ensemble

	"""

	def __init__(self,file_list=list(),data=None,num_realizations=0,metric="chi2"):
		
		self.file_list = file_list
		self.data = data
		self.num_realizations = num_realizations
		self.metric = metric

	####################################
	#############I/O####################
	####################################

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

	@classmethod
	def fromdata(cls,npy_data):

		"""
		Builds the ensemble from data in numpy array format

		:param npy_data: numpy array with the data, the first dimension must be the number of realizations
		:type npy_data: array
		
		"""

		return cls(num_realizations=npy_data.shape[0],data=npy_data)


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

		new_ensemble = cls.fromfilelist([filename])
		new_ensemble.load(from_old=True,callback_loader=callback_loader,**kwargs)
		return new_ensemble

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

		return reduce(add,[ cls.read(f,callback_loader,**kwargs) for f in filelist ])


	def load(self,callback_loader=None,pool=None,from_old=False,**kwargs):
		"""
		Loads the ensemble into memory, can spread the calculations on multiple processors using a MPI pool

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

		>>> conv_ensemble = Ensemble.fromfilelist(map_list)
		>>> conv_ensemble.load(callback_loader=default_callback_loader,pool=pool,l_edges=l_edges)

		"""

		if callback_loader is None:
			callback_loader = lambda f: np.load(f)

		self.pool = pool

		#Build a function wrapper of the callback loader, so it becomes pickleable
		_callback_wrapper = _function_wrapper(callback_loader,args=tuple(),kwargs=kwargs)

		#Execute the callback on each file in the list (spread calculations with MPI pool if it is not none)
		if pool is not None:
			M = pool.map
		else:
			M = map

		full_data = np.array(M(_callback_wrapper,self.file_list))
		
		assert type(full_data) == np.ndarray
		assert full_data.shape[0] == self.num_realizations 

		if from_old:
			full_data = full_data[0]

		self.num_realizations = full_data.shape[0]
		self.data = full_data

	
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
			np.save(filename,self.data)
		elif format=="matlab":
			sio.savemat(filename,{"data": self.data},**kwargs)
		else:
			format(self,filename,**kwargs)

	#Convert into pandas DataFrame
	def toPandas(self,column_label=None):

		"""
		Convert the Ensemble into a pandas DataFrame, interpreting the first dimension as observation and the second dimension as field

		:param column_label: columns of the DataFrame
		:type column_label: list. or Index

		:returns: DataFrame

		"""

		#Check if pandas is installed
		if pd is None:
			raise ImportError("pandas needs to be installed to use this feature!")

		#Convert into DataFrame
		row_index = pd.Index(np.arange(self.num_realizations),name="realization")
		return pd.DataFrame(self.data,index=row_index,columns=column_label)


	####################################
	#############Operations#############
	####################################
	
	def mean(self):

		"""
		Computes the ensemble average over realizations 

		:returns: ndarray with the averages, has the same shape as self.data[0]

		"""
		
		self._mean = self.data.mean(axis=0)
		
		return self._mean

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



	def subset(self,realizations):

		"""
		Returns a new ensemble that contains only the selected realizations

		:param realizations: realizations to keep in the subset Ensemble, must be compatible with numpy array indexing syntax
		:type realizations: int. or array of int.

		:returns: new Ensemble instance that contains only the selected realizations 

		"""

		data = self.data[realizations]
		
		return self.__class__.fromdata(data)


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
		statistic = np.array(M(callback,self[randomizer]))

		#Return the bootstraped statistic expectation value
		return statistic.mean(0)



	def principalComponents(self):

		"""
		Computes the principal components of the Ensemble

		:returns: pcaHandler instance

		"""

		pca = PCA()
		pca.fit(self.data)
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

			raise ValueError("Only chi2 metric implemented so far!!")


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

			splitted.append(self.__class__(file_list=self.file_list,num_realizations=self.num_realizations,data=self.data[:,index[n].first:index[n].last],metric=self.metric))

		return splitted


	#####################################################################################################################

	
	def __add__(self,rhs):

		"""
		Overload of the sum operator: a sum between ensembles is another ensemble whose data are vstacked. This operation is useful if you have two ensembles with different, independent realizations and you want to create an ensemble which is the union of the two. The data must be vstackable

		"""

		#Safety checks
		assert isinstance(rhs,self.__class__)
		assert self.metric == rhs.metric,"The two ensemble instances must have the same metric!!"
			
		if self.data is not None:
			new_data = np.vstack((self.data,rhs.data))
		else:
			new_data = rhs.data

		return self.__class__(file_list=self.file_list+rhs.file_list,data=new_data,num_realizations=self.num_realizations+rhs.num_realizations,metric=self.metric)


	def __mul__(self,rhs):

		"""
		Overload of the multiplication operator: a multiplication of two ensembles is another ensemble whose data are hstacked. This operation is useful if you have two ensembles with the same realizations but different statistical descriptor, and you want to create an ensemble in which both descriptors are included. The data must be hstackable

		"""

		#Safety checks
		assert isinstance(rhs,self.__class__)
		assert self.metric == rhs.metric
		assert self.num_realizations == rhs.num_realizations

		new_data = np.hstack((self.data,rhs.data))

		return self.__class__(file_list=list(set(self.file_list + rhs.file_list)),data=new_data,num_realizations=self.num_realizations,metric=self.metric)

	
	def __getitem__(self,n):

		"""
		Retrieves the n-th realization of the ensemble (starting from 0)

		:raises: IndexError if out of bounds

		"""

		return self.data[n]

	#Protect Ensemble against changes in the data
	def __setattr__(self,a,x):

		"""
		This overload protects the Ensemble against inplace changes in its data

		"""

		#Call parent method
		super(Ensemble,self).__setattr__(a,x)

		#Recompute means if ensemble data is modified
		if a=="data" and hasattr(self,"_mean"):
			self.mean()








