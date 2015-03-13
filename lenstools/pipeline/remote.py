from abc import ABCMeta,abstractproperty,abstractmethod
import os,glob

###################################################
###########SystemHandler class#####################
###################################################

class SystemHandler(object):

	__metaclass__ = ABCMeta

	"""
	SystemHandler class

	"""

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

	def mkdir(self,d):
		os.mkdir(d)

	def exists(self,d):
		return os.path.exists(d)

	def glob(self,n):
		return glob.glob(n)

	def open(self,f,mode):
		return open(f,mode)
