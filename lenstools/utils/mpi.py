from __future__ import division
import sys,warnings

try:
	
	from mpi4py import MPI
	MPI=MPI
	default_op = MPI.SUM

except ImportError:
	
	MPI=None
	default_op=None
	wmsg = "Could not import mpi4py! (if you set sys.modules['mpi4py']=None please disregard this message)"
	warnings.warn(wmsg)

from emcee.utils import MPIPool
import numpy as np

#################################################################################################
###################MPIWhirlPool: should handle one sided communications too######################
#################################################################################################

class MPIWhirlPool(MPIPool):

	"""
	MPI class handler, inherits from MPI pool and adds one sided communications utilities (using RMA windows)

	"""

	#######################################################################################################################
	
	def openWindow(self,memory,window_type="sendrecv"):

		"""
		Create a RMA window that looks from the master process onto all the other workers

		:param memory: memory buffer on which to open the window
		:type memory: numpy nd array

		"""

		self._window_type = window_type

		#Stats of the memory to open a window onto
		assert isinstance(memory,np.ndarray)
		self.memory = memory

		#Create the window
		if window_type=="RMA":

			self.win = MPI.Win.Create(memory=memory,comm=self.comm)
			self.win.Fence()

		elif window_type=="sendrecv":
			self._buffer = np.zeros_like(memory)

		else:
			raise NotImplementedError("Window of type {0} not implemented!".format(window_type))

	#######################################################################################################################
	
	def get(self,process):

		"""
		Read data from an RMA window open on a particular process

		"""

		if self._window_type=="RMA":
		
			read_buffer = np.zeros(self.memory.shape,dtype=self.memory.dtype)

			self.win.Fence()

			if self.is_master():
				self.win.Get(read_buffer,process)

			self.win.Fence()

			return read_buffer

		elif self._window_type=="sendrecv":
			raise NotImplementedError

	#######################################################################################################################

	def accumulate(self,op=default_op):

		"""
		Accumulates the all the window data on the master, performing a custom operation (default is sum)

		"""

		#All the tasks that participate in the communication
		tasks = range(self.size+1)

		#Cycle until only master is left
		while len(tasks)>1:

			if self._window_type=="RMA":
				self.win.Fence()
				
			#Odd tasks communicate the info to the even ones
			try:
				n = tasks.index(self.rank)
				if n%2:

					if self._window_type=="RMA":
						self.win.Accumulate(self.memory,tasks[n-1],op=op)
					elif self._window_type=="sendrecv":
						self.comm.Send(self.memory,dest=tasks[n-1],tag=tasks[n])

				elif n!=len(tasks)-1:

					if self._window_type=="sendrecv":
						self.comm.Recv(self._buffer,source=tasks[n+1],tag=tasks[n+1])

						if op==default_op:
							self.memory += self._buffer
						else:
							raise NotImplementedError

			except ValueError:
				pass

			finally:

				if self._window_type=="RMA":
					self.win.Fence()
				elif self._window_type=="sendrecv":
					self.comm.Barrier()

				#Remove all tasks in odd positions (which already communicated)
				purge = list()
				for n in range(len(tasks)):
					if n%2:
						purge.append(tasks[n])

				for t in purge:
					tasks.remove(t)

				#Safety barrier
				self.comm.Barrier()
				
	
	#######################################################################################################################

	def closeWindow(self):

		"""
		Closes a previously opened RMA window

		"""
		if self._window_type=="RMA":
			
			self.win.Fence()
			self.win.Free()
		
		elif self._window_type=="sendrecv":
			pass