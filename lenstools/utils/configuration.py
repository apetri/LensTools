from ..pipeline.remote import LocalSystem
from ..simulations import Nicaea
from ..simulations.gadget2 import Gadget2SnapshotPipe
from .fft import FFTEngine,NUMPYFFTPack

#Import all the modules that use FFT operations
from ..image import convergence,shear,noise
from ..simulations import nbody,raytracing

modules_with_fft = [convergence,shear,noise,nbody,raytracing]

###################
#Default cosmology#
###################

class LensToolsCosmology(Nicaea):

	"""
	LensTools pipeline cosmology handler

	"""


#####################
#Configuration class#
#####################

class Configuration(object):

	def __init__(self):

		###################
		#Default cosmology#
		###################
		
		self.CosmoDefault = LensToolsCosmology() 

		########################
		#Default system handler#
		########################

		self.syshandler = LocalSystem()

		##########################
		#Default snapshot handler#
		##########################

		self.snapshot_handler = Gadget2SnapshotPipe

		#########################
		#Fast Fourier Transforms#
		#########################

		self.fftengine = NUMPYFFTPack

	def __setattr__(self,a,v):
		
		super(Configuration,self).__setattr__(a,v)

		#Update FFT engine in all the modules that make use of it
		if a=="fftengine":
			fftengine = v()
			assert isinstance(fftengine,FFTEngine)
			for module in modules_with_fft:
				module.fftengine = fftengine


#######################
#Default configuration#
#######################

configuration = Configuration()