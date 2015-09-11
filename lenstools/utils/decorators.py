import sys
from .mpi import MPIWhirlPool

from mpi4py import MPI

def spread_across_pool(func):

	def spreaded_func():

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

		#Execute
		func(pool)

		#Finish
		if pool is not None:
			pool.close()
			pool.comm.Barrier()
			MPI.Finalize()
		

	return spreaded_func