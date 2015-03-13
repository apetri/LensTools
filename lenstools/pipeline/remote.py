from abc import ABCMeta,abstractproperty,abstractmethod
import os,glob

###################################################
###########SystemHandler class#####################
###################################################

class SystemHandler(object):

	__metaclass__ = ABCMeta

	##################################
	######Abstract methods############
	##################################

	@abstractmethod
	def mkdir(self,d):
		pass

	@abstractmethod
	def exists(self,d):
		pass

	@abstractmethod
	def glob(self,n):
		pass

	@abstractmethod
	def open(self,f,mode):
		pass


#############################################
#########Local filesystem ###################
#############################################

class LocalSystem(SystemHandler):

	"""
	Local system handler

	"""

	def __init__(self,readonly=False):
		self.readonly = readonly

	#############################################
	######Abstract method definitions############
	#############################################

	def mkdir(self,d):

		if self.readonly:
			raise IOError("Simulation batch is read only!")

		os.mkdir(d)

	def exists(self,d):
		return os.path.exists(d)

	def glob(self,n):
		return glob.glob(n)

	def open(self,f,mode):

		if (self.readonly) and ("w" in mode or "a" in mode):
			raise IOError("Simulation batch is read only!")

		return open(f,mode)
