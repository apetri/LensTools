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