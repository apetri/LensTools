import sys
from .mpi import MPIWhirlPool

try:
	from mpi4py import MPI
except ImportError:
	MPI = None

class Parallelize(object):

	@classmethod
	def masterworker(cls,func):

		def spreaded_func(*args,**kwargs):

			#MPI Pool
			try:
				pool = MPIWhirlPool()
			except ValueError:
				pool = None

			if (pool is not None) and (not pool.is_master()):
				pool.wait()
				pool.comm.Barrier()
				MPI.Finalize()
				sys.exit(0)

			#Replace the pool in the arguments with the newly created one
			if "pool" in kwargs.keys():
				kwargs["pool"] = pool

			#Execute
			result = func(*args,**kwargs)

			#Finish
			if pool is not None:
				pool.close()
				pool.comm.Barrier()
				MPI.Finalize()

			#Return the result
			if pool is None or pool.is_master():
				return result
		
		#Restore the documentation
		spreaded_func.__doc__ = func.__doc__

		#Return
		return spreaded_func