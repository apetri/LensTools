import os
import ast
from abc import ABCMeta,abstractproperty,abstractmethod

from distutils import config
from ConfigParser import NoOptionError

import cPickle as pkl
import json

import numpy as np
import astropy.units as u


###################################
#Select parser from file extension#
###################################

def select_parser(filename):

	if filename.endswith(".ini"):

		options = config.ConfigParser()
		options.read([filename])
		return options
	
	elif filename.endswith(".pkl") or filename.endswith(".p"):
		
		options = PickleParser()
		with open(filename,"r") as fp:
			options._buffer = pkl.loads(fp.read())
		return options

	elif filename.endswith(".json"):

		options = JSONParser()
		with open(filename,"r") as fp:
			options._buffer = json.loads(fp.read())
		return options

	else:
		raise NotImplementedError("Config file type not supported!")

#########################
#Types of option parsers#
#########################

#Generic parser
class GenericParser(object):

	__metaclass__ = ABCMeta

	@abstractmethod
	def has_section(self,section):
		pass

	@abstractmethod
	def get(self,section,option):
		pass

	@abstractmethod
	def getint(self,section,option):
		pass

	@abstractmethod
	def getfloat(self,section,option):
		pass

	@abstractmethod
	def getboolean(self,section,option):
		pass

#Pickle parser
class PickleParser(GenericParser):
	 
	def has_section(self,section):
		return True

	def get(self,section,option):
		if hasattr(self._buffer,option):
			
			parsed = getattr(self._buffer,option) 
			
			if type(parsed) in [list,np.ndarray]:
				return ",".join([str(p) for p in parsed])
			elif type(parsed)==u.quantity.Quantity:
				try:
					return ",".join([str(p) for p in parsed.value])
				except TypeError:
					return parsed.value
			elif isinstance(parsed,u.core.UnitBase):
				return parsed.to_string()

			return parsed

		else:
			raise NoOptionError(section,option)

	def getint(self,section,option):
		return self.get(section,option)

	def getfloat(self,section,option):
		return self.get(section,option)

	def getboolean(self,section,option):
		return self.get(section,option)


#JSON parser
class JSONParser(GenericParser):
	
	def has_section(self,section):
		return section in self._buffer

	def get(self,section,option):
		if option in self._buffer[section]:
			parsed = self._buffer[section][option]
			if isinstance(parsed,dict):
				return parsed
			else:
				return str(parsed)
		else:
			raise NoOptionError(section,option)

	def getint(self,section,option):
		return int(self.get(section,option))

	def getfloat(self,section,option):
		return float(self.get(section,option))

	def getboolean(self,section,option):
		parsed = self.get(section,option).lower()
		if parsed=="true":
			return True
		if parsed=="false":
			return False

		raise ValueError("Cannot cast option {0} to boolean!".format(option))


###############################################################################################################################################