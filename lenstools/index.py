"""

.. module:: index 
	:platform: Unix
	:synopsis: This module implements a set of useful indexing operations useful for analyses with multiple statistical descriptors


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from abc import ABCMeta
from abc import abstractproperty

#######################################################
######Individual descriptor classes####################
#######################################################

class Descriptor(object):
	__metaclass__ = ABCMeta

	def __init__(self):
		self.label = None
		self._start = 0

	@abstractproperty
	def first(self):
		pass

	@abstractproperty
	def last(self):
		pass

class PowerSpectrum(Descriptor):
	
	def __init__(self,l_edges):
		
		super(PowerSpectrum,self).__init__()
		self.l_edges = l_edges
		self.l = 0.5*(l_edges[:-1] + l_edges[1:])

	def __str__(self):

		return "Power Spectrum"

	def first(self):
		return self._start

	def last(self):
		return self._start + len(self.l) - 1

###################################
#####Global indexer################
###################################

class Indexer(object):

	def __init__(self,descriptor_list):

		#Initialize index with descriptor list
		self.descriptor_list = descriptor_list
		self._last = 0

		#Create an index of all the descriptor bounds
		for descriptor in descriptor_list:

			increment = descriptor.last() + 1
			descriptor._start = self._last
			self._last += increment

	def __getitem__(self,n):
		return self.descriptor_list[n]
			




