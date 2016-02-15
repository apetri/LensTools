from __future__ import division
import re

from .nbody import NbodySnapshot

import numpy as np

from astropy.cosmology import w0waCDM
import astropy.constants as cnst
import astropy.units as u

###############################################
#########AmigaHalos class######################
###############################################

class AmigaHalos(NbodySnapshot):

	"""
	A class that handles the Amiga Halo Finder (AHF) halo output files. Inherits from the abstract NbodySnapshot

	:warning: Not tested yet!

	"""

	###############################################################################################
	#########################Abstract methods implementation#######################################
	###############################################################################################

	@classmethod
	def buildFilename(cls,root,pool,**kwargs):
		
		#Make sure this is an AHF halo catalog
		assert root.endswith(".AHF_halos")

		#Substitute xxxx with the task number
		if pool is not None:
			return re.sub(r"\.[0-9]{4}\.",".{0:04d}.".format(pool.rank),root)
		else:
			return root

	@classmethod
	def int2root(cls,name,n):
		return re.sub(r"\.[0-9]+\.",".{0:04d}.".format(n),name)

	def getHeader(self,h=0.72,w0=-1.,wa=0.):

		#Create header dictionary
		header = dict()
		header["h"] = h
		header["w0"] = w0
		header["wa"] = wa

		#Read the redshift from the filename
		header["redshift"] = float(re.search(r"\.z([0-9\.]+)\.",self.fp.name).groups()[0])
		header["scale_factor"] = 1./(1.+header["redshift"])
		
		#Look for the log file that contains all the information we need in the header
		task_number = int(re.search(r"\.([0-9]{4})\.",self.fp.name).groups()[0])
		logfile = re.sub(r"\.[0-9]{4}\..*",".{0:02d}.log".format(task_number),self.fp.name)

		#Read the log file and fill the header
		with open(logfile,"r") as logfp:
			ahflog = logfp.read()

		#Cosmological parameters and box size
		Om0,Ode0,box_size = re.search(r"simu\.omega0\s*:\s*([0-9\.]+)\nsimu\.lambda0\s*:\s*([0-9\.]+)\nsimu\.boxsize\s*:\s*([0-9\.]+)",ahflog).groups()
		header["Om0"] = float(Om0)
		header["Ode0"] = float(Ode0)
		header["box_size"] = float(box_size)*1.0e3

		#These are not important
		header["masses"] = np.zeros(6)
		header["num_particles_file"] = 1.
		header["num_particles_total"] = 1.
		header["num_files"] = 1

		#Finally compute the comoving distance
		header["comoving_distance"] = w0waCDM(H0=h*100,Om0=header["Om0"],Ode0=header["Ode0"],w0=header["w0"],wa=header["wa"]).comoving_distance(header["redshift"]).to(u.kpc).value / h

		#Return to user
		return header

	def getPositions(self,first=None,last=None,save=True):

		#Matter density today
		rhoM = self.cosmology.critical_density0 * self.cosmology.Om0

		if first is None:
			first = 0

		if last is None:
			m,x,y,z,rv,c = np.loadtxt(self.fp,usecols=(3,5,6,7,11,42))[first:].T
		else:
			m,x,y,z,rv,c = np.loadtxt(self.fp,usecols=(3,5,6,7,11,42))[first:last].T
		
		positions = np.array((x,y,z)).astype(np.float32).T * self.kpc_over_h
		self.virial_radius = rv * self.kpc_over_h 
		self.concentration = c
		self.weights = ((1./(4*np.pi)) * (c**3/(np.log(1.+c)-c/(1.+c))) * m*(u.Msun/self.header["h"]) / (rhoM*(self.virial_radius**3))).decompose().value

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