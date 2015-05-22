.. _parallel_operations::

Distribute operations across different processors
=================================================

You can take advantage of the :py:class:`~lenstools.statistics.Ensemble` class to distribute custom operations (that do not require inter--process communication) across different processors. Suppose we want to compute the square of a bunch of numbers contained in a file and save the result to another file. Suppose also we want to do this on 100 different files and spread the computations across independent processes. The following code will do it:

::

	import numpy as np
	
	#Define the function that computes the square
	def square(filename):

		data = np.loadtxt(filename)
		try:
			np.savetxt("square_"+filename,data**2)
			return 0
		except IOError:
			return 1

	from lenstools import Ensemble

	#MPI 
	try:

		from emcee.utils import MPIPool
		MPIPool = MPIPool

	except ImportError:

		MPIPool = None

	#Create the pool
	try:
		pool = MPIPool()
	except (ValueError,TypeError):
		pool = None

	#Workers only
	if (pool is not None) and not(pool.is_master()):

		pool.wait()
		sys.exit(0)

	#Operate on these files
	filelist = [ "file{0}.txt".format(n+1) for n in range(100) ]
	ens = Ensemble.fromfilelist(filelist)

	#Call square on every file in the list
	ens.load(square,pool=pool)

	#ens.data will be an array of 0 and 1 depending on which files have been written succesfully and which not

	#Close pool and exit
	if pool is not None:
		pool.close()


