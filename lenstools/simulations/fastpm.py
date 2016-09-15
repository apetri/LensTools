from .nbody import NbodySnapshot
import bigfile

import numpy as np
import astropy.units as u
import astropy.constants as cnst

######################
#FastPMSnapshot class#
######################

class FastPMSnapshot(NbodySnapshot):

	"""
	A class that handles FastPM simulation snapshots

	"""

	_header_keys = ['masses','num_particles_file','num_particles_total','box_size','num_files','Om0','Ode0','h']

	############################
	#Open the file with bigfile#
	############################

	@classmethod
	def open(cls,filename,pool=None,header_kwargs=dict(),**kwargs):
		fp = bigfile.BigFile(cls.buildFilename(filename,pool,**kwargs))
		return cls(fp,pool,header_kwargs=header_kwargs)

	###################################################################################
	######################Abstract method implementation###############################
	###################################################################################

	@classmethod
	def buildFilename(cls,root,pool):
		return root

	@classmethod
	def int2root(cls,name,n):
		return name

	def getHeader(self):

		#Initialize header
		header = dict()
		bf_header = self.fp["."].attrs

		###############################################
		#Translate fastPM header into lenstools header#
		###############################################

		#Number of particles/files
		header["num_particles_file"] = bigfile.BigData(self.fp).size
		header["num_particles_total"] = header["num_particles_file"]
		header["num_files"] = 1

		#Cosmology
		header["Om0"] = bf_header["OmegaM"][0]
		header["Ode0"] = 1. - header["Om0"]
		header["w0"] = -1.
		header["wa"] = 0.
		header["h"] = 0.72

		#Box size in kpc/h
		header["box_size"] = bf_header["BoxSize"][0]*1.0e3

		#Masses
		header["masses"] = np.array([0.,bf_header["M0"][0],0.,0.,0.,0.])

		#################

		return header

	def setLimits(self):
		self._first = None
		self._last = None

	def getPositions(self,first=None,last=None,save=True):

		#Get data pointer
		data = bigfile.BigData(self.fp)
		
		#Read in positions in Mpc/h
		if (first is None) or (last is None):
			positions = data["Position"][:]*self.Mpc_over_h
			aemit = data["Aemit"][:]
		else:
			positions = data["Position"][first:last]*self.Mpc_over_h
			aemit = data["Aemit"][first:last]

		#Maybe save
		if save:
			self.positions = positions
			self.aemit = aemit

		#Return
		return positions 

	###########################################################################################

	def getVelocities(self,first=None,last=None,save=True):
		raise NotImplementedError

	def getID(self,first=None,last=None,save=True):
		raise NotImplementedError

	def write(self,filename,files=1):
		raise NotImplementedError