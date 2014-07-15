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
		self._first = None
		self._last = None
		self._start = 0

	@abstractproperty
	def first(self):
		pass

	@abstractproperty
	def last(self):
		pass

class PowerSpectrum(Descriptor):
	"""
	Power spectrum type descriptor class
	"""
	
	def __init__(self,l_edges):
		
		super(PowerSpectrum,self).__init__()
		self.l_edges = l_edges
		self.l = 0.5*(l_edges[:-1] + l_edges[1:])

	def __str__(self):

		return "Power Spectrum descriptor: {0} bins\nEdges: {1}\nMidpoints{2}".format(len(self.l),self.l_edges.__str__(),self.l.__str__())

	@property
	def first(self):

		self._first = self._start
		return self._first

	@property
	def last(self):

		self._last = self._start + len(self.l)
		return self._last

###################################
#####Global indexer################
###################################

class Indexer(object):
	"""
	This class is useful for indexing array of statistical descriptors that are hstacked together

	:param descriptor_list: list of Descriptor subinstances, such as PowerSpectrum; each of these sub instances must be a subclass of Descriptor with a 'first' and 'last' getter methods implemented
	:type descriptor_list: list.

	"""

	def __init__(self,descriptor_list):

		#Initialize index with descriptor list
		self.descriptor_list = descriptor_list
		self._last = 0
		self.num_descriptors = len(descriptor_list)

		#Create an index of all the descriptor bounds
		for descriptor in descriptor_list:

			increment = descriptor.last
			descriptor._start = self._last
			self._last += increment

	def __getitem__(self,n):
		"""Allows to access the n-th descriptor by indexing"""

		return self.descriptor_list[n]
			




