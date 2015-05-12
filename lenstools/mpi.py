try:
	
	from mpi4py import MPI
	MPI=MPI
	default_op = MPI.SUM

except ImportError:
	
	MPI=None
	default_op=None
	print("Warning! mpi4py installation not found or broken!")

from emcee.utils import MPIPool
import numpy as np

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

		#TODO: This can run in log(N) time

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