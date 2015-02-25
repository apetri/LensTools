import os

############################################################
#############EnvironmentSettings class######################
############################################################

class EnvironmentSettings(object):

	"""
	This class handles the system specific environment settings, such as directory paths, modules, etc...

	"""

	def __init__(self,home="SimTest/Home",storage="SimTest/Storage"):

		"""
		Creates the home (meant to store small files like execution scripts) and storage (meant to store large files like simulation outputs) directories

		:param home: name of the simulation home directory
		:type home: str.

		:param storage: name of the simulation mass storage directory
		:type storage: str.

		"""

		self.home = home
		self.storage = storage

		#Create directories if they do not exist yet
		if not os.path.isdir(home):
			os.mkdir(home)

		if not os.path.isdir(storage):
			os.mkdir(storage)

#################################################
###########NGenICSettings class##################
#################################################

class NGenICSettings(object):

	"""
	Class handler of NGenIC settings
	
	"""

	def __init__(self,**kwargs):

		self.GlassFile = "dummy_glass_little_endian.dat"
		self.Redshift = 100.0

		#TODO these need to be computed with Lam's F77 routines
		self.GrowthFactor = 76.078362
		self.VelocityPrefactor = 5.15001

		self.SphereMode = 1 
		self.WhichSpectrum = 2
		self.FileWithInputSpectrum = "powerSpectrum.txt"

		self.InputSpectrum_UnitLength_in_cm = 3.085678e24
		self.ReNormalizeInputSpectrum = 1
		self.ShapeGamma = 0.21
		self.PrimordialIndex = 1.0

		self.NumFilesWrittenInParallel = 1

		self.UnitLength_in_cm = 3.085678e21
		self.UnitMass_in_g = 1.989e43
		self.UnitVelocity_in_cm_per_s = 1e5

		#Typically do not touch these, needed for prefactors calculations
		self._zminact = 0.0
		self._zmaxact = 110.0
		self._iwmode = 3
		self._delz = 0.000001

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])


#################################################
###########NGenICSettings class##################
#################################################

class PlaneSettings(object):

	"""
	Class handler of plane generation settings

	"""

	def __init__(self,**kwargs):

		self.directory_name = "Planes"
		self.format = "fits"

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])



