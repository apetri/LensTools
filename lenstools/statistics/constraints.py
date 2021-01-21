"""

.. module:: constraints
	:platform: Unix
	:synopsis: This module implements the usual statistical tools you need to calculate cosmological parameters confidence intervals


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division,print_function,with_statement

import sys
from operator import mul
from functools import reduce

if sys.version_info.major>=3:
	import _pickle as pickle
else:
	import cPickle as pickle

#########################################################

import numpy as np
import pandas as pd

from scipy import stats,interpolate

from emcee.ensemble import _function_wrapper

try:
	from matplotlib.patches import Ellipse
except ImportError:
	Ellipse = None

#########################################################

from ..utils.algorithms import precision_bias_correction
from .ensemble import Series,Ensemble 
from . import samplers

#########################################################
#############Default Gaussian data likelihood############
#########################################################

def gaussian_likelihood(chi2,norm=1.0):
	return norm*np.exp(-0.5*chi2)

######################################################################
##########Default chi2 calculation with the sandwich product##########
######################################################################

def chi2(parameters,*args,**kwargs):

	model_feature = _predict(parameters,kwargs["interpolator"])
	inverse_covariance = kwargs["inverse_covariance"]

	if model_feature.ndim==1:
		observed_feature = kwargs["observed_feature"]
	else: 
		observed_feature = kwargs["observed_feature"][None,:]

	inverse_covariance_dot = np.dot(observed_feature - model_feature,inverse_covariance)

	return ((observed_feature - model_feature) * inverse_covariance_dot).sum(-1)
	

#######################################################################
#############Feature prediction wrapper################################
#######################################################################

#Fast interpolation method
def _interpolate_fast(p,parameter_grid,method,weights,epsilon):
	return method(((parameter_grid[None]-p[:,None])**2).sum(-1),epsilon).dot(weights)

def _predict(parameters,interpolator):

	#Cast to higher dimension
	parameters = np.atleast_2d(parameters)

	if isinstance(interpolator,list):
		#For each feature bin, compute its interpolated value
		interpolated_feature = np.zeros((parameters.shape[0],len(interpolator)))

		for n,i in enumerate(interpolator):
			interpolated_feature[:,n] = i()(*parameters.T)

	else:
		#Compute fast interpolation
		interpolated_feature = interpolator(parameters)

	return np.squeeze(interpolated_feature)


##############################################
###########Analysis base class################
##############################################

class Analysis(Ensemble):
	
	"""
	The base class of this module; the idea in weak lensing analysis is that one has a set of simulated data, that serves as training model, and then uses that set to fit the observations for the best model parameters. Inherits from :py:class:`Ensemble`

	"""

	_analysis_type = None

	def _check_valid(self):
		assert "parameters" in self.columns.levels[0],"There are no parameters specified for this analysis!"

	@classmethod
	def from_features(cls,features,parameters=None,feature_index=None,parameter_index=None):

		#If features and parameters are already in DataFrame instances then just append them
		if isinstance(features,pd.DataFrame) and isinstance(parameters,pd.DataFrame):
			assert len(parameters.columns.levels[0])==1 and parameters.columns.levels[0][0]=="parameters"
			return cls.concat((parameters,features),axis=1)

		#Cast shapes correctly
		if len(features.shape)==1:
			features = features[None]

		if parameters is None:
			parameters = np.arange(len(features))[None]

		if len(parameters.shape)==1:
			parameters = parameters[None]

		#Make the indices
		if parameter_index is None:
			parameter_index = Series.make_index(pd.Index(range(parameters.shape[1]),name="parameters"))
		elif isinstance(parameter_index,list):
			parameter_index = Series.make_index(pd.Index(parameter_index,name="parameters"))

		if feature_index is None:
			feature_index = Series.make_index(pd.Index(range(features.shape[1]),name="features"))
		elif isinstance(feature_index,list):
			feature_index = Series.make_index(pd.Index(feature_index,name="features"))

		#Instantiate the parameter and feature part of the analysis
		analysis_param = cls(parameters,columns=parameter_index)
		analysis_features = cls(features,columns=feature_index)

		#Instantiate Analysis
		return cls.concat((analysis_param,analysis_features),axis=1)


	##################
	####Properties####
	##################

	@property
	def feature_names(self):
		all_names = list(self.columns.levels[0])
		all_names.remove("parameters")
		return all_names

	@property
	def parameter_names(self):
		return list(self["parameters"].columns)

	@property
	def parameter_set(self):
		return self["parameters"].values

	@property
	def feature_set(self):
		return self[self.feature_names].values

	def parameters(self,names=None):

		if names is None:
			return self

		parameter_names = self.parameter_names
		if isinstance(names,str):
			names = [names]

		subset = self.copy()
		exclude_names = filter(lambda n:not n in names,parameter_names)

		for n in exclude_names:
			subset.pop(("parameters",n))

		return subset
		

	def features(self,names=None):

		if names is None:
			return self
		elif isinstance(names,str):
			return self[["parameters"]+[names]].copy()
		elif isinstance(names,list):
			return self[["parameters"]+names].copy()
		elif isinstance(names,dict):

			pieces = [self[["parameters"]]]
			for key in names.keys():
				piece = self[key][names[key]]
				piece.add_name(key)
				pieces.append(piece)

			return self.__class__.concat(pieces,axis=1)

		else:
			raise TypeError("names type not supported!")


	##################
	####Operations####
	##################

	def add_models(self,parameters,feature):

		"""
		Add a model to the training set of the current analysis

		:param parameters: parameter set of the new model
		:type parameters: array

		:param feature: measured feature of the new model
		:type feature: array

		"""

		#Cast dimensions
		if len(parameters.shape)==1:
			parameters = parameters[None]

		if len(feature.shape)==1:
			feature = feature[None]

		#Check for input valudity
		assert len(parameters)==len(feature)
		assert parameters.shape[1] == self.parameter_set.shape[1]
		assert feature.shape[1:] == self.feature_set.shape[1:]

		#hstack
		parameters_and_features = np.hstack((parameters,feature))

		#Return the newly created Analysis
		return self.append(self._constructor(parameters_and_features,columns=self.columns),ignore_index=True)


	def reparametrize(self,transformation,**kwargs):

		"""
		Reparametrize the parameter set of the analysis by calling the formatter handle on the current parameter set (can be used to enlarge/shrink/relabel the parameter set)

		:param transformation: transformation function called on the parameters, must take in a row of parameters and return another row of parameters
		:type transformation: callable

		:param kwargs: the keyword arguments are passed to the transformation callable
		:type kwargs: dict.

		:returns: reparametrized Analysis

		"""

		#Apply the transformation
		new_parameters = self["parameters"].apply(transformation,axis=1,**kwargs)
		new_parameters.columns.name = "parameters"
		new_parameters.columns = Series.make_index(new_parameters.columns)

		#Return the reparametrized analysis
		reparametrized_analysis = self.copy()
		reparametrized_analysis.pop("parameters")
		return self.__class__.concat((new_parameters,reparametrized_analysis),axis=1)


	def refeaturize(self,transformation,method="apply_row",**kwargs):

		"""
		Allows a general transformation on the feature set of the analysis by calling an arbitrary transformation function

		:param transformation: callback function called on the feature_set; must take in a row of features and return a row of features. If a dictionary is passed, the keys must be the feature names
		:type transformation: callable or dict.

		:param kwargs: the keyword arguments are passed to the transformation callable
		:type kwargs: dict.

		:returns: transformed Analysis

		"""

		#Build transformation dictionary
		if isinstance(transformation,dict):
			transformation_dict = dict((n,lambda x:x) for n in self.feature_names)
			for n in transformation.keys():
				transformation_dict[n] = transformation[n]
		else:
			transformation_dict = dict((n,transformation) for n in self.feature_names)

		#Apply the transformations to each feature
		transformed_features = list()
		for n in self.feature_names:
			
			if method=="apply_row":
				transformed_feature = self[[n]].apply(transformation_dict[n],axis=1,**kwargs)
			elif method=="apply_whole":
				transformed_feature = transformation_dict[n](self[[n]],**kwargs)
				transformed_feature.add_name(n)
				transformed_feature.index = self.index
			elif method=="dot":
				transformed_feature = self[[n]].dot(transformation_dict[n])
				transformed_feature.add_name(n)
				transformed_feature.index = self.index 
			else:
				raise NotImplementedError("transformation method {0} not implemented!".format(method))

			transformed_features.append(transformed_feature)

		#Concatenate and return
		return self.__class__.concat([self[["parameters"]]]+transformed_features,axis=1)


	def combine_features(self,combinations):

		"""
		Combine features in the Analysis, according to a dictionary which keys are the name of the combined features

		:param combinations: mapping of combined features onto the old ones 
		:type combinations: dict.

		"""

		combined_features  = list()

		#Cycle over the combinations keys
		for n in combinations.keys():

			#Select
			combined_feature = self[combinations[n]].copy()
			
			#Merge the column names
			combined_feature_index = pd.Index(np.hstack([ combined_feature[c].columns.values for c in combinations[n] ]),name=n)
			combined_feature_index = Series.make_index(combined_feature_index)
			combined_feature.columns = combined_feature_index

			#Append to the combination list
			combined_features.append(combined_feature)

		#Concatenate everything
		return self.__class__.concat([self[["parameters"]]]+combined_features,axis=1)
	

	###############################################################################################################################


	def find(self,parameters,rtol=1.0e-05):

		"""
		Finds the location in the instance that has the specified combination of parameters

		:param parameters: the parameters of the model to find
		:type parameters: array.

		:param rtol: tolerance of the search (must be less than 1)
		:type rtol: float.

		:returns: array of int. with the indices of the corresponding models

		"""

		assert len(parameters)==self.parameter_set.shape[1]

		search_result = np.all(np.isclose(self.parameter_set,parameters,rtol=rtol),axis=1)
		return np.where(search_result==True)[0]

	###############################################################################################################################

	@staticmethod
	def ellipse(center,covariance,p_value=0.684,**kwargs):

		"""

		Draws a confidence ellipse using matplotlib Ellipse patch

		:param center: center of the ellipse
		:type center: tuple.

		:param covariance: parameters covariance matrix
		:type covariance: 2D-array.

		:param p_value: p-value to calculate
		:type p_value: float.
						
		:param kwargs: the keyword arguments are passed to the matplotlib Ellipse method
		:type kwargs: dict.

		:returns: matplotlib ellipse object
		:rtype: Ellipse

		"""

		#Check that ellipse patch is available
		if Ellipse is None:
			raise ImportError("The matplotlib Ellipse patch is necessary to use this method!")

		#Compute the directions and sizes of the ellipse axes
		w,v = np.linalg.eigh(covariance)
		width,height = 2*np.sqrt(w * stats.chi2(2).ppf(p_value))

		try:
			angle = 180.*np.arctan(v[1,0] / v[0,0]) / np.pi
		except ZeroDivisionError:
			angle = 90.

		#Draw the ellipse
		return Ellipse(center,width,height,angle=angle,**kwargs)


###################################################
#############Fisher matrix analysis################
###################################################

class FisherSeries(Series):

	@property
	def _constructor_expanddim(self):
		return FisherAnalysis

class FisherAnalysis(Analysis):

	################################################################
	##############DataFrame subclassing#############################
	################################################################

	@property 
	def _constructor_sliced(self):
		return FisherSeries

	@property
	def _constructor_expanddim(self):
		raise NotImplementedError("Expand dimension not supported")

	#################################################################

	_analysis_type = "Fisher"
	_fiducial = 0

	"""
	The class handler of a Fisher matrix analysis, inherits from the base class Analysis

	"""

	def set_fiducial(self,n):

		"""
		Sets the fiducial model (with respect to which to compute the derivatives), default is 0 (i.e. self.parameter_set[0])

		:param n: the parameter set you want to use as fiducial
		:type n: int.

		"""

		assert n < self.parameter_set.shape[0],"There are less than {0} models in your analysis".format(n+1)

		self._fiducial = n

	@property
	def fiducial(self):

		return self.feature_set[self._fiducial]

	@property
	def _variations(self):

		"""
		Checks the parameter variations with respect to the fiducial cosmology

		:returns: bool array (True if the parameter is varied, False otherwise)

		"""

		return self.parameter_set!=self.parameter_set[self._fiducial]

	@property
	def variations(self):

		"""
		Checks the parameter variations with respect to the fiducial cosmology

		:returns: iterable with the positions of the variations

		"""

		for n,b in enumerate(self._variations.sum(1)):
			if b:
				yield n


	def check(self):

		"""
		Asserts that the parameters are varied one at a time, and that a parameter is not varied more than once

		:raises: AssertionError

		"""

		assert (self._variations.sum(1)<2).all(),"You can vary only a parameter at a time!"

		#Check how many variations are there for each parameter
		num_par_variations = self._variations.sum(0)
		if (num_par_variations<2).all():
			return 0
		else:
			return 1

	def where(self,par=None):

		"""
		Finds the locations of the varied parameters in the parameter set

		:returns: dict. with the locations of the variations, for each parameter

		"""

		loc = dict()
		v = np.where(self._variations==1)

		#Decide if keys are lists or simple numbers
		if self.check():

			for n in range(self.parameter_set.shape[1]):
				loc[n] = list()

			for n in range(len(v[0])):
				loc[v[1][n]].append(v[0][n])

		else:

			for n in range(len(v[0])):
				loc[v[1][n]] = v[0][n]

		if par is None:
			return loc
		else:
			return loc[par]


	@property
	def varied(self):

		"""
		Returns the indices of the parameters that are varied 

		:returns: list with the indices of the varied parameters

		"""
		
		return list(sorted(self.where())) 


	def compute_derivatives(self):

		"""
		Computes the feature derivatives with respect to the parameter sets using one step finite differences; the derivatives are computed with respect to the fiducial parameter set

		:returns: array of shape (p,N), where N is the feature dimension and p is the number of varied parameters

		"""

		assert self.parameter_set.shape[0] > 1,"You need at least 2 models to proceed in a Fisher Analysis!"
		assert self.check()==0,"Finite differences implemented only at first order! Cannot compute derivatives"

		#Find the varied parameters and their locations
		loc_varied = self.where()
		par_varied = list(sorted(loc_varied))

		#Allocate space for the derivatives
		derivatives = np.zeros((len(par_varied),)+self.feature_set.shape[1:])

		#cycle to parameters to calculate derivatives
		for n,p in enumerate(par_varied):
			
			#Calculate the finite difference derivative with respect to this parameter
			derivatives[n]  = (self.feature_set[loc_varied[p]] - self.feature_set[self._fiducial]) / (self.parameter_set[loc_varied[p],p] - self.parameter_set[self._fiducial,p])

		#set the derivatives attribute and return the result
		self._derivatives = self.__class__(derivatives,index=[self.parameter_names[n] for n in par_varied],columns=self[self.feature_names].columns)

	@property 
	def derivatives(self):
		if not hasattr(self,"_derivatives"):
			self.compute_derivatives()	
		return self._derivatives


	def chi2(self,observed_feature,features_covariance,correct=None):

		"""
		Computes the chi2 between an observed feature and the fiducial feature, using the provided covariance

		:param observed_feature: observed feature to fit, its last dimension must have the same shape as self.feature_set[0] 
		:type observed_feature: array

		:param features_covariance: covariance matrix of the simulated features, must be provided for a correct fit!
		:type features_covariance: 2 dimensional array (or 1 dimensional if diagonal)

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:returns: chi2 of the comparison
		:rtype: float.

		"""

		#Cast pandas types
		assert features_covariance is not None,"No science without the covariance matrix, you must provide one!"
		
		if isinstance(observed_feature,pd.Series) or isinstance(observed_feature,pd.DataFrame):
			observed_feature = observed_feature.values

		if isinstance(features_covariance,pd.DataFrame):
			features_covariance = features_covariance.values

		#Cast the observed feature in suitable shape
		if len(observed_feature.shape)==1:
			observed_feature = observed_feature[None]
			single = True
		else:
			single = False

		#Check for correct shape of input
		assert observed_feature.shape[-1:]==self.feature_set.shape[-1:]
		assert features_covariance.shape in [self.feature_set.shape[-1:],self.feature_set.shape[-1:]*2]

		#Compute the difference
		difference = observed_feature - self.fiducial[None]

		#Compute the chi2
		if features_covariance.shape==self.feature_set.shape[-1:]:
			result = ((difference**2)/features_covariance[None]).sum(-1)
		else:
			if correct is not None:
				result = (difference * np.linalg.solve(features_covariance/precision_bias_correction(correct,len(features_covariance)),difference.transpose()).transpose()).sum(-1)
			else:
				result = (difference * np.linalg.solve(features_covariance,difference.transpose()).transpose()).sum(-1)

		#Return the result
		if single:
			return result[0]
		else:
			return result

	
	def fit(self,observed_feature,features_covariance):

		"""
		Maximizes the gaussian likelihood on which the Fisher matrix formalism is based, and returns the best fit for the parameters given the observed feature

		:param observed_feature: observed feature to fit, must have the same shape as self.feature_set[0]
		:type observed_feature: array

		:param features_covariance: covariance matrix of the simulated features, must be provided for a correct fit!
		:type features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:returns: array with the best fitted parameter values

		"""

		assert features_covariance is not None,"No science without the covariance matrix, you must provide one!"

		#Cast pandas types
		if isinstance(observed_feature,pd.Series) or isinstance(observed_feature,pd.DataFrame):
			observed_feature = observed_feature.values

		if isinstance(features_covariance,pd.DataFrame):
			features_covariance = features_covariance.values

		#Check for correct shape of input
		assert (observed_feature.shape==self.feature_set.shape[1:]) or (observed_feature.shape[1:]==self.feature_set.shape[1:]) 
		assert (features_covariance.shape==self.feature_set.shape[1:] * 2) or (features_covariance.shape==self.feature_set.shape[1:])

		#Linear algebra manipulations (parameters = M x features)
		if features_covariance.shape==self.feature_set.shape[1:] * 2:
			Y = np.linalg.solve(features_covariance,self.derivatives.values.transpose())
		else:
			Y = (1/features_covariance[:,np.newaxis]) * self.derivatives.values.transpose()

		XY = np.dot(self.derivatives.values,Y)
		M = np.linalg.solve(XY,Y.transpose())

		#Compute difference in parameters (with respect to the fiducial model)
		if observed_feature.ndim==1:
			dP = np.dot(M,observed_feature - self.feature_set[self._fiducial])
		else:
			dP = np.dot(observed_feature - self.feature_set[[self._fiducial]],M.T)

		#Return the actual best fit
		if dP.ndim==1:
			return self._constructor_sliced(self.parameter_set[self._fiducial,self.varied] + dP,index=self.derivatives.index)
		else:
			return self.__class__(self.parameter_set[self._fiducial,self.varied][None] + dP,columns=self.derivatives.index)


	def classify(self,observed_feature,features_covariance,correct=None,labels=range(2),confusion=False):

		"""
		Performs a Fisher classification of the observed feature, choosing the most probable label based on the value of the chi2
		
		:param observed_feature: observed feature to fit, the last dimenstion must have the same shape as self.feature_set[0]
		:type observed_feature: array

		:param features_covariance: covariance matrix of the simulated features, must be provided for a correct classification!
		:type features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:param labels: labels of the classification, must be the indices of the available classes (from 0 to feature_set.shape[0])
		:type labels: iterable

		:param confusion: if True, an array with the label percentage occurrences is returned; if False an array of labels is returned
		:type confusion: bool.

		:returns: array with the labels resulting from the classification
		:rtype: int.

		"""

		fiducial_original = self._fiducial

		#Compute all the chi squared values, for each observed feature and each label
		all_chi2 = list()
		for l in labels:
			self.set_fiducial(l)
			all_chi2.append(self.chi2(observed_feature,features_covariance,correct=correct))

		self.set_fiducial(fiducial_original)

		#Cast the list into an array
		all_chi2 = np.array(all_chi2)

		#Find the minima
		chi2_min = all_chi2.argmin(0)

		#Translate into the corresponding classes
		classes = np.zeros_like(chi2_min)
		for n,l in enumerate(labels):
			classes[chi2_min==n] = l

		if confusion:

			#Compute confusion array
			confusion_array = np.zeros(n+1)
			for n,l in enumerate(labels):
				confusion_array[n] = (classes==l).sum() / len(classes)

			#Return
			return confusion_array
		
		else:	
			#Return
			return classes

	###############################################################################################################################################################

	def parameter_covariance(self,simulated_features_covariance,correct=None,observed_features_covariance=None):

		"""
		Computes the parameter covariance matrix using the associated features, that in the end allows to compute the parameter confidence contours (around the fiducial value)

		:param simulated_features_covariance: covariance matrix of the simulated features, must be provided for a correct fit!
		:type simulated_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:param observed_features_covariance: covariance matrix of the simulated features, if different from the simulated one; if None the simulated feature covariance is used
		:type observed_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:returns: 2 dimensional array with the parameter covariance matrix of the analysis

		"""

		assert simulated_features_covariance is not None,"No science without the covariance matrix, you must provide one!"

		#Cast pandas types
		if isinstance(simulated_features_covariance,pd.DataFrame):
			simulated_features_covariance = simulated_features_covariance.values

		if (observed_features_covariance is not None) and (isinstance(observed_features_covariance,pd.DataFrame)):
			observed_features_covariance = observed_features_covariance.values

		#Check for correct shape of input
		assert simulated_features_covariance.shape == self.feature_set.shape[1:] * 2 or simulated_features_covariance.shape == self.feature_set.shape[1:]

		#Linear algebra manipulations (parameters = M x features)
		if simulated_features_covariance.shape ==  self.feature_set.shape[1:] * 2:
			Y = np.linalg.solve(simulated_features_covariance,self.derivatives.values.transpose())
		else:
			Y = (1/simulated_features_covariance[:,np.newaxis]) * self.derivatives.values.transpose()
		
		XY = np.dot(self.derivatives.values,Y)

		#If we are using the same covariance matrix for observations and simulations, then XY is the Fisher matrix; otherwise we need to compute M too
		if observed_features_covariance is None:
			if correct is not None:
				return self.__class__(np.linalg.inv(XY),index=self.derivatives.index,columns=self.derivatives.index) / precision_bias_correction(correct,len(simulated_features_covariance))
			else:
				return self.__class__(np.linalg.inv(XY),index=self.derivatives.index,columns=self.derivatives.index)
		else:

			assert observed_features_covariance.shape == self.feature_set.shape[1:] * 2 or observed_features_covariance.shape == self.feature_set.shape[1:]
			
			M = np.linalg.solve(XY,Y.transpose())
			
			if observed_features_covariance.shape == self.feature_set.shape[1:] * 2:
				parameter_covariance = np.dot(M,np.dot(observed_features_covariance,M.transpose()))
			else:
				parameter_covariance = np.dot(M * observed_features_covariance,M.transpose())

			return self.__class__(parameter_covariance,index=self.derivatives.index,columns=self.derivatives.index)

	
	def fisher_matrix(self,simulated_features_covariance,correct=None,observed_features_covariance=None):

		"""
		Computes the parameter Fisher matrix using the associated features, that in the end allows to compute the parameter confidence contours (around the fiducial value)

		:param simulated_features_covariance: covariance matrix of the simulated features, must be provided for a correct fit!
		:type simulated_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:param observed_features_covariance: covariance matrix of the simulated features, if different from the simulated one; if None the simulated feature covariance is used
		:type observed_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:returns: 2 dimensional array with the Fisher matrix of the analysis

		"""

		parcov = self.parameter_covariance(simulated_features_covariance,correct,observed_features_covariance)
		return self.__class__(np.linalg.inv(parcov.values),index=parcov.index,columns=parcov.columns)


	def confidence_ellipse(self,simulated_features_covariance,correct=None,observed_feature=None,observed_features_covariance=None,parameters=["Om","w"],p_value=0.684,**kwargs):

		"""

		Draws a confidence ellipse of a specified p-value in parameter space, corresponding to fit an observed feature for the cosmological parameters

		:param observed_feature: observed feature to fit, the last dimenstion must have the same shape as self.feature_set[0]
		:type observed_feature: array

		:param simulated_features_covariance: covariance matrix of the simulated features, must be provided for a correct fit!
		:type simulated_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:param observed_features_covariance: covariance matrix of the simulated features, if different from the simulated one; if None the simulated feature covariance is used
		:type observed_features_covariance: 2 dimensional array (or 1 dimensional if assumed diagonal)

		:param parameters: parameters to compute the condifence contour of
		:type parameters: list.

		:param p_value: p-value to calculate
		:type p_value: float.

		:param kwargs: the keyword arguments are passed to the matplotlib Ellipse method
		:type kwargs: dict.

		:returns: matplotlib ellipse object
		:rtype: Ellipse

		"""

		if len(parameters)!=2:
			raise ValueError("You must specify exactly two parameters to draw the ellipse of!")

		#If the observed feature is not provided, the center of the ellipse is (0,0)
		if observed_feature is None:
			center = (0,0)
		else:
			#Fit the observed feature and put the center here
			p_fit = self.fit(observed_feature,simulated_features_covariance)
			center = tuple(p_fit[parameters])

		#The parameter covariance sets the size and orientation of the ellipse
		p_cov = self.parameter_covariance(simulated_features_covariance,correct,observed_features_covariance)[parameters].loc[parameters]
		
		#Return Ellipse object to user
		return self.ellipse(center,p_cov.values,p_value,**kwargs)


	###########################################################################################################################################


	def reparametrize(self,formatter,*args,**kwargs):

		#Call the parent method
		super(FisherAnalysis,self).reparametrize(formatter,*args,**kwargs)

		#Check that the format of the parameter set is valid
		self.check()

#######################################################
#############Full analysis#############################
#######################################################

class EmulatorSeries(Series):

	@property
	def _constructor_expanddim(self):
		return Emulator

class Emulator(Analysis):

	_analysis_type = "Emulator"

	"""
	The class handler of a full likelihood analysis; the parameter likelihood function is calculated with an interpolation of various kind between simulation points

	"""

	################################################################
	##############DataFrame subclassing#############################
	################################################################

	@property 
	def _constructor_sliced(self):
		return EmulatorSeries

	@property
	def _constructor_expanddim(self):
		raise NotImplementedError("Expand dimension not supported")

	##################################
	########Constructor###############
	##################################

	def __init__(self,*args,**kwargs):
		
		super(Emulator,self).__init__(*args,**kwargs) 
		self._likelihood_function = gaussian_likelihood

		if "_likelihood_function" not in self._metadata:
			self._metadata.append("_likelihood_function")

	#######################################################################################################################################

	def approximate_linear(self,center,derivative_precision=0.1):

		"""
		Construct a FisherAnalysis by approximating the Emulator as a linear expansion along a chosen center

		:param center: center point in parameter space
		:type center: Series

		:param derivative_precision: percentage step for the finite difference derivatives
		:type derivative_precision: float.

		:returns: linearly approximated Fisher analysis
		:rtype: FisherAnalysis

		"""

		npar = len(self.parameter_names)

		#Construct the parameter set for the Fisher Analysis
		parameters_fisher = Ensemble(np.zeros((npar+1,npar)),columns=self.parameter_names)

		#Fiducial
		for n in range(len(parameters_fisher)):
			parameters_fisher.iloc[n] = center

		#Variations
		for n in range(1,len(parameters_fisher)):
			parameters_fisher.iloc[n,n-1] += np.abs(parameters_fisher.iloc[n,n-1])*derivative_precision

		#Predict the features with the emulator, build the Fisher Analysis
		features = self.predict(parameters_fisher)
		parameters_fisher.add_name("parameters")

		return FisherAnalysis.from_features(features,parameters=parameters_fisher)

	#######################################################################################################################################

	def set_likelihood(self,function=None):

		"""
		Sets the likelihood function to a custom function input by the user: the default is the usual exp(-0.5*chi^2)

		"""

		assert function is not None
		self._likelihood_function = function

	def train(self,use_parameters="all",method="Rbf",**kwargs):

		"""
		Builds the interpolators for each of the feature bins using a radial basis function approach

		:param use_parameters: which parameters actually vary in the supplied parameter set (it doesn't make sense to interpolate over the constant ones)
		:type use_parameters: list. or "all"

		:param method: interpolation method; can be 'Rbf' or callable. If callable, it must take two arguments, a square distance and a square length smoothing scale
		:type method: str. or callable

		:param kwargs: keyword arguments to be passed to the interpolator constructor

		"""

		#input sanity check
		if use_parameters != "all":
			assert type(use_parameters) == list
			used_parameters = self.parameter_set[:,use_parameters]
		else:
			used_parameters = self.parameter_set

		#Compute total number of feature bins and reshape the training set accordingly
		if "_num_bins" not in self._metadata:
			self._metadata.append("_num_bins")
		self._num_bins = reduce(mul,self.feature_set.shape[1:])

		flattened_feature_set = self.feature_set.reshape((self.feature_set.shape[0],self._num_bins))

		#Build the interpolator
		if "_interpolator" not in self._metadata:
			self._metadata.append("_interpolator")

		if method=="Rbf":

			#Scipy Rbf method
			self._interpolator = list()

			for n in range(self._num_bins):
				self._interpolator.append(_interpolate_wrapper(interpolate.Rbf,args=(tuple(used_parameters.T) + (flattened_feature_set[:,n],)),kwargs=kwargs))

		else:

			#Compute pairwise square distance between points
			distances = ((used_parameters[None] - used_parameters[:,None])**2).sum(-1)
			epsilon = distances[np.triu_indices(len(distances),k=1)].mean()
			kernel = method(distances,epsilon)
			weights = np.linalg.solve(kernel,self.feature_set)

			#Wrap interpolator
			self._interpolator = _function_wrapper(_interpolate_fast,args=[],kwargs={"parameter_grid":used_parameters,"method":method,"weights":weights,"epsilon":epsilon})


	###############################################################################################################################################################

	def predict(self,parameters,raw=False):

		"""
		Predicts the feature at a new point in parameter space using the bin interpolators, trained with the simulated features

		:param parameters: new points in parameter space on which to compute the chi2 statistic; it'a (N,p) array where N is the number of points and p the number of parameters, or array of size p if there is only one point
		:type parameters: array  

		:param raw: if True returns raw numpy arrays
		:type raw: bool.

		:returns: predicted features
		:rtype: array or :py:class:`Ensemble`

		"""

		#If you didn't do training before, train now with the default settings
		if not hasattr(self,"_interpolator"):
			self.train()

		#Cast DataFrames to numpy arrays
		if isinstance(parameters,pd.DataFrame):
			assert (parameters.columns==self["parameters"].columns).all(),"Parameters do not match!"
			parameters = parameters.values
		elif isinstance(parameters,pd.Series):
			assert (parameters.index==self["parameters"].columns).all(),"Parameters do not match!"
			parameters = parameters.values

		#####################################
		#Interpolate to compute the features#
		#####################################

		interpolated_feature = _predict(parameters,self._interpolator)

		############################################################################################

		#Return the result
		if raw:
			return interpolated_feature
		else:
			if parameters.ndim==1:
				return Series(interpolated_feature.reshape(self.feature_set.shape[1:]),index=self[self.feature_names].columns)
			else:
				return Ensemble(interpolated_feature.reshape((parameters.shape[0],) + self.feature_set.shape[1:]),columns=self[self.feature_names].columns)


	###############################################################################################################################################################


	def chi2(self,parameters,observed_feature,features_covariance,correct=None,split_chunks=None,pool=None):

		"""
		Computes the chi2 part of the parameter likelihood with the usual sandwich product with the covariance matrix; the model features are computed with the interpolators

		:param parameters: new points in parameter space on which to compute the chi2 statistic
		:type parameters: (N,p) array where N is the number of points and p the number of parameters

		:param observed_feature: observed feature on which to condition the parameter likelihood
		:type observed_feature: array

		:param features_covariance: covariance matrix of the features, must be supplied
		:type features_covariance: array

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:param split_chunks: if set to an integer bigger than 0, splits the calculation of the chi2 into subsequent chunks, each that takes care of an equal number of points. Each chunk could be taken care of by a different processor
		:type split_chunks: int.

		:returns: array with the chi2 values, with the same shape of the parameters input

		"""

		#Sanity checks
		assert observed_feature is not None 
		assert features_covariance is not None,"No science without the covariance matrix, you must provide one!"
		assert observed_feature.shape == self.feature_set.shape[1:]
		assert features_covariance.shape == observed_feature.shape * 2

		#If you didn't do training before, train now with the default settings
		if not hasattr(self,"_interpolator"):
			self.train()

		#Reformat the parameter input into a list of chunks
		if parameters.ndim==1:
			num_points = 1
		else:
			num_points = parameters.shape[0]

		if split_chunks is None:
			
			parameter_chunks = [parameters]
		
		elif split_chunks > 0:
			
			assert num_points%split_chunks == 0,"split_chunks must divide exactly the number of points!!"
			chunk_length = num_points//split_chunks
			parameter_chunks = [ parameters[n*chunk_length:(n+1)*chunk_length] for n in range(split_chunks) ]

		else:

			raise ValueError("split_chunks must be >0!!")

		#Compute the inverse of the covariance matrix once and for all
		covinv = np.linalg.inv(features_covariance)
		if correct is not None:
			covinv *= precision_bias_correction(correct,len(covinv))

		#Build the keyword argument dictionary to be passed to the chi2 calculator
		kwargs = {"interpolator":self._interpolator,"inverse_covariance":covinv,"observed_feature":observed_feature}

		#Hack to make the chi2 pickleable (from emcee)
		chi2_wrapper = _function_wrapper(chi2,tuple(),kwargs)

		#Finally map chi2 calculator on the list of chunks
		if pool is not None:
			M = pool.map
		else:
			M = map
		
		chi2_list = list(M(chi2_wrapper,parameter_chunks))

		return np.array(chi2_list).reshape(num_points)


	def chi2Contributions(self,parameters,observed_feature,features_covariance,correct=None): 

		"""
		Computes the individual contributions of each feature bin to the chi2; the model features are computed with the interpolators. The full chi2 is the sum of the individual contributions

		:param parameters: new points in parameter space on which to compute the chi2 statistic
		:type parameters: (N,p) array where N is the number of points and p the number of parameters

		:param observed_feature: observed feature on which to condition the parameter likelihood
		:type observed_feature: array

		:param features_covariance: covariance matrix of the features, must be supplied
		:type features_covariance: array

		:param correct: if not None, correct for the bias in the inverse covariance estimator assuming the covariance was estimated by 'correct' simulations
		:type correct: int.

		:returns: numpy 2D array with the contributions to the chi2 (off diagonal elements are the contributions of the cross correlation between bins)

		"""

		#Sanity checks
		assert observed_feature is not None 
		assert features_covariance is not None,"No science without the covariance matrix, you must provide one!"
		assert observed_feature.shape == self.feature_set.shape[1:]
		assert features_covariance.shape == observed_feature.shape * 2

		#If you didn't do training before, train now with the default settings
		if not hasattr(self,"_interpolator"):
			self.train()

		#Compute each bin contribution to the chi2
		residuals = observed_feature - self.predict(parameters)

		#Compute the inverse covariance
		covinv = np.linalg.inv(features_covariance)
		if correct is not None:
			covinv *= precision_bias_correction(correct,len(inverse_covariance_dot))

		#Compute the hits map
		return np.outer(residuals,residuals) * covinv



	def likelihood(self,chi2_value,**kwargs):

		"""
		Computes the likelihood value with the selected likelihood function, given the pre-computed chi2 value

		:param chi2_value: chi squared values 
		:type chi2_value: array

		:param kwargs: keyword arguments to be passed to your likelihood function

		"""

		return self._likelihood_function(chi2_value,**kwargs)

	############################################################################################################################################

	def score(self,parameters,observed_feature,method="chi2",**kwargs):

		"""
		Compute the score for an observed feature for each combination of the proposed parameters

		:param parameters: parameter combinations to score
		:type parameters: DataFrame or array

		:param observed_feature: observed feature to score
		:type observed_feature: Series

		:param method: scoring method to use (defaults to chi2): if callable, must take in the current instance, the parameter array and the observed feature and return a score for each parameter combination
		:type method: str. or callable

		:param kwargs: keyword arguments passed to the callable method
		:type kwargs: dict.

		:returns: ensemble with the scores, for each feature
		:rtype: :py:class:`Ensemble`

		"""

		#Get the names of the features to use
		feature_names = list(observed_feature.index.levels[0])
		try:
			feature_names.remove("parameters")
		except ValueError:
			pass

		#Check that the observed feature columns and the Emulator columns correspond
		for c in feature_names:
			assert c in self.feature_names,"Feature '{0}' is not present in the Emulator!".format(c)

		#Reorder the parameters according to the ones in the Emulator
		parameters = parameters[self.parameter_names]

		#If the method is chi2, a covariance matrix must be provided
		if method=="chi2":
			assert "features_covariance" in kwargs.keys()
			features_covariance = kwargs["features_covariance"]
			del(kwargs["features_covariance"])
			assert (features_covariance.index==observed_feature.index).all()
			assert (features_covariance.columns==observed_feature.index).all()

		#Build an Ensemble with the parameters
		score_ensemble = Ensemble(parameters)

		#For each feature, compute the score
		for c in feature_names:

			#Isolate the Emulator that concerns this feature only
			sub_emulator = self.features(c)
			sub_emulator_columns = sub_emulator[c].columns

			#Isolate the observed sub_feature
			sub_feature = observed_feature[c][sub_emulator_columns].values

			#If the method is chi2, use the already implemented version of it
			if method=="chi2":
				sub_emulator.train()
				sub_feature_covariance = features_covariance[c][sub_emulator_columns].loc[c].loc[sub_emulator_columns].values
				score_ensemble[c] = sub_emulator.chi2(parameters=parameters.values,observed_feature=sub_feature,features_covariance=sub_feature_covariance,**kwargs)
			else:
				score_ensemble[c] = method(sub_emulator,parameters.values,sub_feature,**kwargs)

		#Return the score ensemble
		return score_ensemble

	############################################################################################################################################

	def sample_posterior(self,observed_feature,sample="emcee",**kwargs):

		"""
		Sample the parameter posterior distribution

		:param observed_feature: observed feature to score
		:type observed_feature: Series

		:param sample: posterior sampling method
		:type sample: str. or callable

		:returns: samples from the posterior distribution
		:rtype: dict. 

		"""

		if sample=="emcee":
			sample = samplers.emcee_sampler

		#Get the names of the features to use
		feature_names = list(observed_feature.index.levels[0])
		try:
			feature_names.remove("parameters")
		except ValueError:
			pass

		#Check that the observed feature columns and the Emulator columns correspond
		for c in feature_names:
			assert c in self.feature_names,"Feature '{0}' is not present in the Emulator!".format(c)

		#Check if the user provides a covariance matrix
		if "features_covariance" in kwargs.keys():
			features_covariance = kwargs["features_covariance"]
			del(kwargs["features_covariance"])
			assert (features_covariance.index==observed_feature.index).all()
			assert (features_covariance.columns==observed_feature.index).all()
		else:
			features_covariance = None

		#Parameter samples
		samples = dict()

		#For each feature, compute the score
		for c in feature_names:

			#Isolate the Emulator that concerns this feature only
			sub_emulator = self.features(c)
			sub_emulator_columns = sub_emulator[c].columns

			#Isolate the observed sub_feature
			sub_feature = observed_feature[c][sub_emulator_columns].values

			#Train the emulator
			sub_emulator.train()

			#Isolate the sub feature covariance matrix, proceed with the sampling
			if features_covariance is not None:
				sub_feature_covariance = features_covariance[c][sub_emulator_columns].loc[c].loc[sub_emulator_columns].values
				samples[c] = sample(emulator=sub_emulator,observed_feature=sub_feature,features_covariance=sub_feature_covariance,**kwargs)
			else:
				samples[c] = sample(emulator=sub_emulator,observed_feature=sub_feature,**kwargs)

		#Done, return the sampled points
		return samples
				

	############################################################################################################################################

	def set_to_model(self,parameters):

		"""
		Set the current model of the emulator to the one specified by the parameter set

		:param parameters: parameters for which the feature will be emulated
		:type parameters: array.

		"""

		#assert parameters.shape[0]==self.parameter_set.shape[1]
		
		#if not hasattr(self,"_interpolator"):
		#	self.train()
		
		#self._current_model_parameters = parameters
		#self._current_predicted_feature = self.predict(parameters)
		#self._current_interpolated_feature = interp1d(self.feature_label,self._current_predicted_feature)

		raise NotImplementedError


	def emulate(self,new_feature_label):

		"""
		Compute an emulated feature at the new feature label specified (multipoles, thresholds, ...) for the current model, using a linear interpolation between bins

		:param new_feature_label: new feature label for which you want to emulate the feature
		:type new_feature_label: array.

		:returns: the emulated feature

		""" 

		#return self._current_interpolated_feature(new_feature_label)
		raise NotImplementedError

#########################################################################################################################################################################################

###########################################################################
###########Hack to make scipy interpolate objects pickleable###############
###########################################################################

class _interpolate_wrapper(object):

	def __init__(self,f,args,kwargs):
		self.f = f
		self.args = args
		self.kwargs = kwargs
	
	def __call__(self):
		try:
			return self.f(*self.args,**self.kwargs)
		except:
			import traceback
			print("lenstools: Exception while building the interpolators")
			print(" exception:")
			traceback.print_exc()
			raise 
