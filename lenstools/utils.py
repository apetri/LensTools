from __future__ import division

import numpy as np
from scipy import special as sp,integrate

from .mpi import MPIWhirlPool
MPIWhirlPool = MPIWhirlPool

####################################################################
#################Hankel transforms##################################
####################################################################

def fht(n,l_binned,func,**kwargs):

	if(kwargs.has_key('theta')):
		theta = kwargs['theta']
	else:
		theta_min = 1.0/l_binned.max()
		theta = l_binned*(theta_min/l_binned.min())

	h_kernel = sp.jn(n,np.outer(l_binned,theta))
	
	integrand = np.dot(np.diag(l_binned*func),h_kernel) * (2*np.pi)
	transform = integrate.simps(integrand,l_binned,axis=0)
	
	return theta,transform

def ifht(n,l_binned,func,**kwargs):
	
	if(kwargs.has_key('theta')):
		theta = kwargs['theta']
	else:
		theta_min = 1.0/l_binned.max()
		theta = l_binned*(theta_min/l_binned.min())
	
	h_kernel = sp.jn(n,np.outer(l_binned,theta))
	
	integrand = np.dot(np.diag(l_binned*func),h_kernel) / (2*np.pi)
	transform = integrate.simps(integrand,l_binned,axis=0)
	
	return theta,transform

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

	def fit(self,data):

		#Scale the data to zero mean and unit variance
		self._pca_mean = data.mean(0)
		self._pca_std = data.std(0)
		self._data_scaled = data.copy()
		self._data_scaled -= self._pca_mean[None]
		self._data_scaled /= self._pca_std[None]
		self._data_scaled /= np.sqrt(self._data_scaled.shape[0] - 1)

		#Perform singular value decomposition
		left,eigenvalues,right = np.linalg.svd(self._data_scaled,full_matrices=False)

		#Assign eigenvalues and eigenvectors as attributes
		self.components_ = right
		self.explained_variance_ = eigenvalues**2 

	@property
	def eigenvalues(self):
		return self.explained_variance_

	@property
	def eigenvectors(self):

		return self.components_*np.sqrt(self._data_scaled.shape[0] - 1)*self._pca_std[None] + self._pca_mean[None]

	def transform(self,X):

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
			return components[0]
		else:
			return components

	def inverse_transform(self,X,n_components=None):

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
			return original_components[0]
		else:
			return original_components

	def select_components(self,X,n_components):

		all_components = self.transform(X)
		return self.inverse_transform(all_components,n_components=n_components) 
