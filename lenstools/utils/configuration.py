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
		
		self.CosmoDefault = LensToolsCosmology() 

		##############################################
		#Parametere name to attribute name dictionary#
		##############################################

		self.name2attr = dict()
		self.name2attr["Om"] = "Om0"
		self.name2attr["Ol"] = "Ode0"
		self.name2attr["w"] = "w0"
		self.name2attr["wa"] = "wa"
		self.name2attr["h"] = "h"
		self.name2attr["Ob"] = "Ob0"
		self.name2attr["si"] = "sigma8"
		self.name2attr["ns"] = "ns"

		############################################
		#Number of digits of precision for cosmo_id#
		############################################

		self.cosmo_id_digits = 3

		########################
		#Default system handler#
		########################

		self.syshandler = LocalSystem()

		###################################
		#Default simulation tree json file#
		###################################

		self.json_tree_file = ".tree.json"

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