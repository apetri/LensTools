from __future__ import division
import numpy as np

from emcee.utils import MPIPool

try:
	
	from mpi4py import MPI
	MPI=MPI
	default_op = MPI.SUM

except ImportError:
	
	MPI=None
	default_op=None
	print("Warning! mpi4py installation not found or broken!")

#####################################################################################
#######Supplying to the lack of rfftfreq implementation in numpy<1.8#################
#####################################################################################

def rfftfreq(n, d=1.0):

    if not (isinstance(n,int) or isinstance(n, integer)):
        raise ValueError("n should be an integer")
    val = 1.0/(n*d)
    N = n//2 + 1
    results = np.arange(0, N, dtype=int)
    return results * val

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


#################################################################################################
###################MPIWhirlPool: should handle one sided communications too######################
#################################################################################################

class MPIWhirlPool(MPIPool):

	"""
	MPI class handler, inherits from MPI pool and adds one sided communications utilities (using RMA windows)

	"""

	def openWindow(self,memory):

		"""
		Create a RMA window that looks from the master process onto all the other workers

		:param memory: memory buffer on which to open the window
		:type memory: numpy nd array

		"""

		#Stats of the memory to open a window onto
		assert isinstance(memory,np.ndarray)
		self.memory = memory

		#Create the window
		self.win = MPI.Win.Create(memory=memory,comm=self.comm)
		self.win.Fence()

	
	def get(self,process):

		"""
		Read data from an RMA window open on a particular process

		"""

		read_buffer = np.zeros(self.memory.shape,dtype=self.memory.dtype)

		self.win.Fence()

		if self.is_master():
			self.win.Get(read_buffer,process)

		self.win.Fence()

		if self.is_master():
			return read_buffer
		else:
			return None

	def accumulate(self,op=default_op):

		"""
		Accumulates the all the window data on the master, performing a custom operation (default is sum)

		"""

		for n in range(1,self.size+1):

			self.win.Fence()

			if(self.rank==n):
				self.win.Accumulate(self.memory,0,op=op)

			self.win.Fence()


	def closeWindow(self):

		"""
		Closes a previously opened RMA window

		"""
		self.win.Fence()
		self.win.Free()
