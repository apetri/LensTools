.. _parallel_operations::

Distribute operations across different processors
=================================================

You can take advantage of the :py:class:`~lenstools.statistics.ensemble.Ensemble` class to distribute custom operations (that do not require inter--process communication) across different processors. Suppose we want to compute the square of a bunch of numbers contained in a file and save the result to another file. Suppose also we want to do this on 100 different files and spread the computations across independent processes. The following code will do it:

::

	import numpy as np
	from lenstools.utils.decorators import Parallelize
	from lenstools.statistics.ensemble import Ensemble
	
	#Define the function that computes the square
	def square(filename):

		data = np.loadtxt(filename)
		try:
			np.savetxt("square_"+filename,data**2)
			return 0
		except IOError:
			return 1

	#Main execution
	@Parallelize.masterworker
	def main(pool):

		#Operate on these files
		filelist = [ "file{0}.txt".format(n+1) for n in range(100) ]

		#ens will be an instance of :py:class:`~lenstools.statistics.ensemble.Ensemble` with only one column, whose elements will be 0 and 1 depending on which files have been written succesfully and which not
		ens = Ensemble.compute(filelist,callback_loader=square,pool=pool)

	if __name__=="__main__":
		main(None)
	


