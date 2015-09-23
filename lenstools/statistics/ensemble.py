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

from ..utils.algorithms import pcaHandler as PCA

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
		return self.__class__

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

	############################################################################################################################

	def combine_columns(self,combinations):

		"""
		Combine the hierarchical columns in the Series, according to a dictionary which keys are the name of the combined features

		:param combinations: mapping of combined features onto the old ones 
		:type combinations: dict.

		:returns: Series with columns combined
		:rtype: :py:class:`Series`

		"""

		combined_columns  = list()

		#Cycle over the combinations keys
		for n in combinations.keys():

			#Select and set new index 
			combined_column = pd.concat(self[c] for c in combinations[n])
			combined_column.index.name = n
			combined_column.index = self.__class__.make_index(combined_column.index) 

			#Append to the combination list
			combined_columns.append(combined_column)

		#Concatenate everything
		return pd.concat(combined_columns)

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
			callback_loader = pd.read_pickle

		#Read the Ensemble
		loaded_ensemble = callback_loader(filename,**kwargs)

		#Preserve the metadata
		if hasattr(loaded_ensemble,"_metadata"):
			meta = dict((k,getattr(loaded_ensemble,k)) for k in loaded_ensemble._metadata)
		else:
			meta = None

		#Instantiate the Ensemble
		loaded_ensemble = cls(loaded_ensemble,file_list=[filename])

		#Complete the metadata
		if meta is not None:
			for k in meta.keys():
				if not k in loaded_ensemble._metadata:
					loaded_ensemble._metadata.append(k)
				setattr(loaded_ensemble,k,meta[k])

		#Return
		return loaded_ensemble

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
	def read_sql_query(cls,sql,con,**kwargs):
		return cls(pd.read_sql_query(sql,con,**kwargs))

	@classmethod
	def read_sql_table(cls,table_name,con,**kwargs):
		return cls(pd.read_sql_table(table_name,con,**kwargs))

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

	@classmethod
	def merge(cls,*args,**kwargs):
		return cls(pd.merge(*args,**kwargs))

	@classmethod
	def combine_from_dict(cls,ensemble_dict):

		"""
		Builds an Ensemble combining the columns of smaller Ensembles; each key in the ensemble_dict dictionary becomes an element in the top level of the resulting Ensemble index

		:param ensemble_dict: dictionary that contains the Ensembles to combine; the values in the dictionary must be Ensembles
		:type ensemble_dict: dict.

		:rtype: :py:class:`Ensemble`

		"""

		combined_list = list()
		combined_columns = list()

		#Cycle over dictionary keys
		for key in ensemble_dict.keys():
			combined_list.append(ensemble_dict[key])
			combined_columns.append(ensemble_dict[key].columns)
			combined_columns[-1].name = key

		#Concatenate
		ensemble_combined = cls.concat(combined_list,axis=1,ignore_index=True)
		ensemble_combined.columns = Series.make_index(*tuple(combined_columns))

		#Return
		return ensemble_combined

	def combine_columns(self,combinations):

		"""
		Combine the hierarchical columns in the Ensemble, according to a dictionary which keys are the name of the combined features

		:param combinations: mapping of combined features onto the old ones 
		:type combinations: dict.

		:returns: Ensemble with columns combined
		:rtype: :py:class:`Ensemble`

		"""

		combined_columns  = list()

		#Cycle over the combinations keys
		for n in combinations.keys():

			#Select
			combined_column = self[combinations[n]].copy()
			
			#Merge the column names
			combined_column_index = pd.Index(np.hstack([ combined_column[c].columns.values for c in combinations[n] ]),name=n)
			combined_column_index = Series.make_index(combined_column_index)
			combined_column.columns = combined_column_index

			#Append to the combination list
			combined_columns.append(combined_column)

		#Concatenate everything
		return self.__class__.concat(combined_columns,axis=1)

	def suppress_indices(self,by,suppress,columns):

		"""
		Combine multiple rows in the Ensemble into a single row, according to a specific criterion

		:param by: list of columns that is used to group the Ensemble: the rows in each group are combined
		:type by: list.

		:param suppress: list of columns to suppress: these indices get suppressed when the rows are combined
		:type suppress: list.

		:param columns: list of columns to keep in the rows that are combined
		:type columns: list.

		:returns: suppressed indices lookup table,combined Ensemble
		:rtype: list.

		"""

		#by group
		by_group = self.groupby(by)
		self["by_group_id"] = by_group.grouper.group_info[0]
		by_labels = self[by+["by_group_id"]].drop_duplicates()

		#suppress group
		suppress_group = self.groupby(suppress)
		self["suppress_group_id"] = suppress_group.grouper.group_info[0]
		suppress_labels = self[suppress+["suppress_group_id"]].drop_duplicates()

		#Combine the rows
		ens_combined = self[["by_group_id","suppress_group_id"]+columns].pivot(index="by_group_id",columns="suppress_group_id")
		ens_combined.columns = ens_combined.columns.rename(None,level=1)
		ens_combined.index.name = None

		#Join the by_labels
		ens_combined = self.__class__.merge(ens_combined.reset_index(),by_labels,left_on="index",right_on="by_group_id")
		ens_combined.pop(("index",""))
		ens_combined.pop("by_group_id")

		#Return
		self.pop("by_group_id")
		self.pop("suppress_group_id")
		return suppress_labels,ens_combined



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
			return self.groupby(lambda i:i/group_size)
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


	def principalComponents(self,location=None,scale=None):

		"""
		Computes the principal components of the Ensemble

		:param location: compute the principal components with respect to this location; must have the same columns as the Ensemble
		:type location: :py:class:`Series`

		:param scale: compute the principal components applying this scaling on the Ensemble columns; must have the same columns as the Ensemble
		:type scale: :py:class:`Ensemble`

		:returns: pcaHandler instance

		"""

		for l in [location,scale]:
			if l is not None:
				assert (l.index==self.columns).all(),"The column names do not match!"

		if location is not None:
			location = location.values

		if scale is not None:
			scale = scale.values

		pca = PCA(constructor_series=self._constructor_sliced,constructor_ensemble=self.__class__,columns=self.columns,location=location,scale=scale)
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
		return self.__class__

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

	def _wrap_applied_output(self,*args,**kwargs):
		
		out = super(SeriesGroupBy,self)._wrap_applied_output(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		else:
			return out

	def _wrap_aggregated_output(self,*args,**kwargs):
		
		out = super(SeriesGroupBy,self)._wrap_aggregated_output(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		else:
			return out

	def aggregate(self,*args,**kwargs):

		out = super(SeriesGroupBy,self).aggregate(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		else:
			return out



###########################
###EnsembleGroupBy class###
###########################

class EnsembleGroupBy(pd.core.groupby.DataFrameGroupBy):

	"""
	Useful for grouping Ensemble

	"""

	def _wrap_applied_output(self,*args,**kwargs):
		
		out = super(EnsembleGroupBy,self)._wrap_applied_output(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		elif isinstance(out,pd.DataFrame):
			return self._ensemble_constructor(out)
		else:
			return out

	def _wrap_aggregated_output(self,*args,**kwargs):
		
		out = super(EnsembleGroupBy,self)._wrap_aggregated_output(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		elif isinstance(out,pd.DataFrame):
			return self._ensemble_constructor(out)
		else:
			return out

	def _wrap_agged_blocks(self,*args,**kwargs):

		out = super(EnsembleGroupBy,self)._wrap_agged_blocks(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		elif isinstance(out,pd.DataFrame):
			return self._ensemble_constructor(out)
		else:
			return out

	def aggregate(self,*args,**kwargs):

		out = super(EnsembleGroupBy,self).aggregate(*args,**kwargs)
		if isinstance(out,pd.Series):
			return self._series_constructor(out)
		elif isinstance(out,pd.DataFrame):
			return self._ensemble_constructor(out)
		else:
			return out


	def __getitem__(self,name):
		item = super(EnsembleGroupBy,self).__getitem__(name)
		item.__class__ = SeriesGroupBy
		item._series_constructor = self._series_constructor
		return item





