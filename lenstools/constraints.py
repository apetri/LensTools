"""

.. module:: constraints
	:platform: Unix
	:synopsis: This module implements the usual statistical tools you need to calculate cosmological parameters confidence intervals


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

import numpy as np

##############################################
###########Analysis base class################
##############################################

class Analysis(object):
	
	"""
	The base class of this module; the idea in weak lensing analysis is that one has a set of simulated data, that serves as training model, and then uses that set to fit the observations for the best model parameters.

	:param parameter_set: the points in parameter space that the simulated data cover; the first axis refers to the model, the second to the number of parameters 
	:type parameter_set: array

	:param training_set: the measured feature in each model; the first axis refers to the model, the others to the feature indices (or bins)
	:type training_set: array

	:param observed_set: the measured feature in the data, should be a one dimensional array
	:type observed_set: array

	"""

	_analysis_type = None

	def __init__(self,parameter_set=None,training_set=None,observed_set=None):

		assert self._analysis_type is not None,"Don't instantiate this class directly, use one of its subclasses!"
		
		if parameter_set is not None and training_set is not None:
			assert parameter_set.shape[0] == training_set.shape[0],"There should be one feature for each of the simulated models!"

		self.parameter_set = parameter_set
		self.training_set = training_set
		self.observed_set = observed_set

	def __repr__(self):

		try:
			return "{0} type analysis, based on {1} models spanning a {2}-dimensional parameter space".format(self._analysis_type,self.parameter_set.shape[0],self.parameter_set.shape[1])
		except AttributeError:
			return "{0} type analysis, no models in it yet!".format(self._analysis_type)

	def add_model(self,parameters,feature):

		"""
		Add a model to the training set of the current analysis

		:param parameters: parameter set of the new model
		:type parameters: array

		:param feature: measured feature of the new model
		:type feature: array

		"""

		#If the analysis doesn't have any models, add the first, otherwise simply vstack them
		if self.parameter_set is None:
			
			assert self.training_set is None

			self.parameter_set = parameters.copy()[np.newaxis,:]
			self.training_set = feature.copy()[np.newaxis,:]

		else:

			#Check for input valudity
			assert parameters.shape[0] == self.parameter_set.shape[1]
			assert feature.shape == self.training_set.shape[1:]

			self.parameter_set = np.vstack((self.parameter_set,parameters))
			self.training_set = np.vstack((self.training_set,feature))




###################################################
#############Fisher matrix analysis################
###################################################

class FisherAnalysis(Analysis):

	_analysis_type = "Fisher"

	"""
	The class handler of a Fisher matrix analysis, inherits from the base class Analysis

	"""



