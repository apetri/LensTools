from __future__ import division
import numpy as np

from emcee.utils import MPIPool

try:
	from mpi4py import MPI
except:
	pass

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
		self.mpi_data_type = MPI.FLOAT
		self.memory = memory

		#Create the window
		self.win = MPI.Win.Create(memory=memory,comm=self.comm)
		self.win.Fence()

	
	def get(self,process):

		"""
		Read data from an RMA window open on a particular process

		"""

		read_buffer = np.zeros(self.memory.shape,dtype=self.memory.dtype)

		if self.is_master():
			self.win.Get([read_buffer,self.mpi_data_type],process)

		self.win.Fence()

		if self.is_master():
			return read_buffer
		else:
			return None

	def accumulate(self,op=MPI.SUM):

		"""
		Accumulates the all the window data on the master, performing a custom operation (default is sum)

		"""

		for n in range(1,self.size+1):

			if(self.rank==n):
				self.win.Accumulate(self.memory,0,op=op)

			self.win.Fence()


	def closeWindow(self):

		"""
		Closes a previously opened RMA window

		"""

		self.win.Free()