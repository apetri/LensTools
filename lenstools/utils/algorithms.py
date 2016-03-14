from __future__ import division
from operator import add
from functools import reduce

import numpy as np
import pandas as pd

#######################################################################################################################
#############################Correction factors for inverse covariance matrix estimator################################
#######################################################################################################################

precision_bias_correction = lambda nr,nb: (nr-nb-2)/(nr-1) 

#################################################################################################
##################Convenient definition of step function (fast implementation)###################
#################################################################################################

def step(x,intervals,vals):
	if vals.dtype==np.int:
		return reduce(add,[vals[n]*(np.sign(x-i[0])+np.sign(i[1]-x)).astype(np.int)//2 for n,i in enumerate(intervals)])
	else:
		return reduce(add,[0.5*vals[n]*(np.sign(x-i[0])+np.sign(i[1]-x)) for n,i in enumerate(intervals)])

#################################################################################################
#############################Principal Component Analysis handler################################
#################################################################################################

def pca_transform(data,pca,n_components):
	assert n_components<=pca.components_.shape[0]
	return pca.transform(data).T[:n_components].T

class pcaHandler(object):

	"""
	Handles principal component analysis

	"""

	def __init__(self,constructor_series,constructor_ensemble,columns,location,scale):
		
		self._constructor_series = constructor_series
		self._constructor_ensemble = constructor_ensemble
		self._columns = columns
		self._location = location
		self._scale = scale

	def fit(self,data):

		#Subtract the mean 
		if self._location is not None:
			self._pca_mean = self._location
		else:
			self._pca_mean = data.mean(0)

		#Scale by the measure units
		if self._scale is not None:
			self._pca_std = self._scale
		else:
			self._pca_std = data.std(0)

		#Whiten the data
		self._data_scaled = data.copy()
		self._data_scaled -= self._pca_mean[None]
		self._data_scaled /= self._pca_std[None]

		#Scale by sqrt(N-1)
		self._data_scaled /= np.sqrt(self._data_scaled.shape[0] - 1)

		#Perform singular value decomposition
		left,eigenvalues,right = np.linalg.svd(self._data_scaled,full_matrices=False)

		#Assign eigenvalues and eigenvectors as attributes
		self.components_ = right
		self.explained_variance_ = eigenvalues**2 

	@property
	def eigenvalues(self):
		return self._constructor_series(self.explained_variance_)

	@property
	def eigenvectors(self):
		e = self._constructor_ensemble(self.components_*np.sqrt(self._data_scaled.shape[0] - 1)*self._pca_std[None] + self._pca_mean[None],columns=self._columns)
		e.index.name = "eigenvector"
		e.columns.name = "component"
		return e

	@property
	def directions(self):
		e = self._constructor_ensemble(self.components_,columns=self._columns)
		e.index.name = "eigenvector"
		e.columns.name = "component"
		return e

	def transform(self,X):

		if isinstance(X,pd.Series) or isinstance(X,pd.DataFrame):
			X = X[self._columns].values

		#Cast X to the right dimensions
		if len(X.shape)==1:
			X_copy = X.copy()[None]
		else:
			X_copy = X.copy()

		#Subtract mean and scale by variance
		X_copy -= self._pca_mean[None]
		X_copy /= (self._pca_std[None]*np.sqrt(self._data_scaled.shape[0] - 1))

		#Compute the projection via dot product
		components = X_copy.dot(self.components_.transpose())
		if len(X.shape)==1:
			return self._constructor_series(components[0])
		else:
			return self._constructor_ensemble(components)

	def inverse_transform(self,X,n_components=None):

		try:
			X = X.values
		except AttributeError:
			pass

		#Cast X to the right dimensions
		if len(X.shape)==1:
			X_copy = X.copy()[None]
		else:
			X_copy = X.copy()

		#Use the PCA basis vectors to project back to the original space
		if n_components is not None:
			basis_vectors = self.components_[:n_components]
			X_copy = X_copy[:,:n_components]
		else:
			basis_vectors = self.components_

		#Original space
		original_components = X_copy.dot(basis_vectors)

		#De-whitening
		original_components *= (self._pca_std[None]*np.sqrt(self._data_scaled.shape[0] - 1))
		original_components += self._pca_mean[None]  

		if original_components.shape[0]==1:
			return self._constructor_series(original_components[0],index=self._columns)
		else:
			return self._constructor_ensemble(original_components,columns=self._columns)

	def select_components(self,X,n_components):

		all_components = self.transform(X)
		return self.inverse_transform(all_components,n_components=n_components)