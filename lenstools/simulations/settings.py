import sys,os
import ast
from abc import ABCMeta,abstractproperty,abstractmethod

from distutils import config

if sys.version_info.major>=3:
	import _pickle as pkl
	from configparser import NoOptionError
else:
	import cPickle as pkl
	from ConfigParser import NoOptionError

import json

import numpy as np
import astropy.units as u


###################################
#Select parser from file extension#
###################################

def select_parser(filename,read=True):

	if filename.endswith(".ini"):

		options = config.ConfigParser()
		if read:
			options.read([filename])
			options.filename = filename
		
		return options
	
	elif filename.endswith(".pkl") or filename.endswith(".p"):
		
		options = PickleParser()
		if read:
			with open(filename,"r") as fp:
				options._buffer = options.load(fp)
				options.filename = filename
		
		return options

	elif filename.endswith(".json"):

		options = JSONParser()
		if read:
			with open(filename,"r") as fp:
				options._buffer = options.load(fp)
				options.filename = filename
		
		return options

	else:
		raise NotImplementedError("Config file type not supported!")

############################
#Generic LensTools settings#
############################

class LTSettings(object):

	__metaclass__ = ABCMeta

	@abstractmethod
	def __init__(self,*args,**kwargs):
		pass

	#Read
	@classmethod
	def read(cls,config_file,*args):
		options = select_parser(config_file)
		return cls.get(options,*args)

	#Read from dictionary
	@classmethod
	def from_dict(cls,d):
		options = DictParser()
		options._buffer = d
		options.filename = "buffer"
		return cls.get(options)

	#Convert into dictionary
	def to_dict(self):

		#Create object dictionary
		obj_dict = dict()

		#Fill in the values and make JSON serialization possible
		self_dict = self.__dict__
		for key in self_dict:

			if isinstance(self_dict[key],u.core.UnitBase):
				obj_dict[key] = self_dict[key].to_string()
			
			elif type(self_dict[key])==u.quantity.Quantity:
				
				numeric_value = self_dict[key].value
				if isinstance(numeric_value,np.ndarray):
					obj_dict[key] = list(numeric_value)
				else:
					obj_dict[key] = numeric_value
			
			elif type(self_dict[key])==np.ndarray:
				obj_dict[key] = list(self_dict[key])
			
			else:
				obj_dict[key] = self_dict[key]

		#Return
		return obj_dict


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

	@abstractmethod
	def dumps(obj):
		pass

	@abstractmethod
	def loads(s):
		pass

	################################################

	@classmethod
	def dump(obj,fp):
		fp.write(cls.dumps(obj))

	@classmethod
	def load(cls,fp):
		return cls.loads(fp.read())


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

	@staticmethod
	def dumps(obj):
		return pkl.dumps(obj)

	@staticmethod
	def loads(s):
		return pkl.loads(s)


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
		parsed = self.get(section,option)
		
		if type(parsed)==bool:
			return parsed	
		if parsed=="True":
			return True
		if parsed=="False":
			return False

		raise ValueError("Cannot cast option {0} to boolean!".format(option))

	@staticmethod
	def dumps(obj):
		return json.dumps(obj.to_dict())

	@staticmethod
	def loads(s):
		return json.loads(s)


#Parser from dictionary
class DictParser(JSONParser):

	def has_section(self,section):
		return True

	def get(self,section,option):
		if option in self._buffer:
			parsed = self._buffer[option]
			if isinstance(parsed,dict):
				return parsed
			elif isinstance(parsed,list):
				return ",".join([str(p) for p in parsed])
			else:
				return str(parsed)
		else:
			raise NoOptionError(section,option)


###############################################################################################################################################