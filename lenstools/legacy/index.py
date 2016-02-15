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


##############################################################################
##############################################################################

class PowerSpectrum(Descriptor):
	"""
	Power spectrum indexing class

	:param l_edges: bin edges of multipole moments
	:type l_edges: array

	"""

	name = "power_spectrum"
	
	def __init__(self,l_edges):
		
		super(PowerSpectrum,self).__init__()
		self.l_edges = l_edges
		self.l = 0.5*(l_edges[:-1] + l_edges[1:])

	def __repr__(self):

		return "Power Spectrum descriptor: {0} bins\nEdges: {1}\nMidpoints{2}".format(len(self.l),self.l_edges.__str__(),self.l.__str__())

	@property
	def first(self):

		self._first = self._start
		return self._first

	@property
	def last(self):

		self._last = self._start + len(self.l)
		return self._last

##############################################################################
##############################################################################

class Moments(Descriptor):
	"""
	Moments indexing class

	"""

	name = "moments"

	def __init__(self,connected=True,dimensionless=False):

		super(Moments,self).__init__()
		self.connected = connected
		self.dimensionless = dimensionless

	def __repr__(self):

		repr_str = "Moments descriptor: 2 quadratic moments, 3 cubic, 4 quartic"
		
		if self.connected:
			repr_str += "\nconnected: yes"
		else:
			repr_str += "\nconnected: no"

		if self.dimensionless:
			repr_str += "\ndimensionless: yes"
		else:
			repr_str += "\ndimensionless: no"

		return repr_str

	@property
	def first(self):

		self._first = self._start
		return self._first

	@property
	def last(self):

		self._last = self._start + 9
		return self._last

##############################################################################
##############################################################################

class Peaks(Descriptor):
	"""
	Peaks histogram indexing class

	:param thresholds: peak histogram bin edges
	:type thresholds: array

	"""

	name = "peaks"

	def __init__(self,thresholds,norm=False):

		super(Peaks,self).__init__()
		self.thresholds = thresholds
		self.midpoints = 0.5*(thresholds[:-1] + thresholds[1:])
		self.norm = norm

	def __repr__(self):

		if self.norm:
			answer = "yes"
		else:
			answer = "no"

		return "Peak count descriptor: {0} bins\nEdges: {1}\nMidpoints{2}\nThresholds in sigma units: {3}".format(len(self.midpoints),self.thresholds.__str__(),self.midpoints.__str__(),answer)

	@property
	def first(self):

		self._first = self._start
		return self._first

	@property
	def last(self):

		self._last = self._start + len(self.midpoints)
		return self._last

##############################################################################
##############################################################################

class MinkowskiSingle(Peaks):
	"""
	Single Minkowski functional indexing class, identical to Peaks

	"""

	name = "minkowski_single"

	def __repr__(self):

		if self.norm:
			answer = "yes"
		else:
			answer = "no"

		return "Minkowski functionals descriptor (single): {0} bins\nEdges: {1}\nMidpoints{2}\nThresholds in sigma units: {3}".format(len(self.midpoints),self.thresholds.__str__(),self.midpoints.__str__(),answer)


##############################################################################
##############################################################################

class PDF(Peaks):
	"""
	Probability density function indexing class, identical to Peaks
	 
	"""

	name = "pdf"

	def __repr__(self):

		if self.norm:
			answer = "yes"
		else:
			answer = "no"

		return "Probability density function descriptor: {0} bins\nEdges: {1}\nMidpoints{2}\nThresholds in sigma units: {3}".format(len(self.midpoints),self.thresholds.__str__(),self.midpoints.__str__(),answer)


##############################################################################
##############################################################################

class MinkowskiAll(Peaks):
	"""
	Minkowski functionals indexing class, inherits from Peaks; the add-ons are some methods that deal with the fact that there are 3 Minkowski functionals

	"""

	name = "minkowski_all"

	def __repr__(self):

		if self.norm:
			answer = "yes"
		else:
			answer = "no"

		return "Minkowski functionals descriptor (all): {0} bins\nEdges: {1}\nMidpoints{2}\nThresholds in sigma units: {3}".format(len(self.midpoints),self.thresholds.__str__(),self.midpoints.__str__(),answer)

	@property
	def last(self):

		self._last = self._start + 3*len(self.midpoints)
		return self._last

	def separate(self):
		"""
		Separates an MinkowskiAll index into 3 separate MinkowsiSingle indices, one for each Minkowski functional

		:returns: list of three MinkowskiSingle instances

		"""

		mink0 = MinkowskiSingle(self.thresholds)
		mink0._start = self._start

		mink1 = MinkowskiSingle(self.thresholds)
		mink1._start = self._start + len(self.midpoints)

		mink2 = MinkowskiSingle(self.thresholds)
		mink2._start = self._start + 2*len(self.midpoints) 

		return [mink0,mink1,mink2]

##############################################################################
##############################################################################

###################################
#####Global indexer################
###################################

class Indexer(object):
	"""
	This class is useful for indexing an array of statistical descriptors that are hstacked together; the Indexer instance keeps track of the memory regions where the different descriptors are stored

	:param descriptor_list: list of Descriptor subinstances, such as PowerSpectrum; each of these sub instances must be a subclass of Descriptor with a 'first' and 'last' getter methods implemented
	:type descriptor_list: list.

	"""

	def __init__(self,descriptor_list):

		#Initialize index with descriptor list
		self.descriptor_list = descriptor_list
		self.num_descriptors = len(descriptor_list)

	@classmethod
	def stack(cls,descriptor_list):

		_last = 0

		#Create an index of all the descriptor bounds by stacking the observables vectors
		for descriptor in descriptor_list:

			increment = descriptor.last - descriptor._start
			descriptor._start = _last
			_last += increment

		new_idx = cls(descriptor_list)
		setattr(new_idx,"_last",_last)

		return new_idx

	@property
	def size(self):

		if hasattr(self,"_last"):
			return self._last
		else: 
			raise AttributeError("This index has not been created by stack() method, cannot determine size!")

	def __getitem__(self,n):
		"""Allows to access the n-th descriptor by indexing"""

		return self.descriptor_list[n]
			




