from __future__ import division

from .nbody import NbodySnapshot

import numpy as np
import astropy.constants as cnst
import astropy.units as u

###############################################
#########AmigaHalos class######################
###############################################

class AmigaHalos(NbodySnapshot):

	"""
	A class that handles the Amiga Halo Finder (AHF) halo output files. Inherits from the abstract NbodySnapshot

	"""

	###############################################################################################
	#########################Abstract methods implementation#######################################
	###############################################################################################

	@classmethod
	def buildFilename(cls,root,pool,**kwargs):
		pass

	def getHeader(self):
		pass

	def getPositions(self,first=None,last=None,save=True):

		#Omega_m and h are necessary
		Om = self._header["Om0"]
		h = self._header["h"]

		if first is None:
			first = 0

		if last is None:
			m,x,y,z,rv,c = np.loadtxt(self.fp,usecols=(3,5,6,7,11,42))[first:].T
		else:
			m,x,y,z,rv,c = np.loadtxt(self.fp,usecols=(3,5,6,7,11,42))[first:last].T
		
		positions = np.array((x,y,z)).astype(np.float).T * self.kpc_over_h
		self.virial_radius = rv * self.kpc_over_h * self._header["scale_factor"]
		self.concentration = c
		self.weights = None

		if save:
			self.positions = positions
			return self.positions

		return positions


	############################################
	###########These are not necessary##########
	############################################

	def getVelocities(self,first=None,last=None,save=True):
		raise NotImplementedError

	def getID(self,first=None,last=None,save=True):
		raise NotImplementedError

	def write(self,filename,files=1):
		raise NotImplementedError