from lenstools.external import _gadget

from astropy.units import kpc,g

############################################################
#################Gadget2Snapshot class######################
############################################################

class Gadget2Snapshot(object):

	"""
	A class that handles I/O from Gadget2 snapshots: it was mainly designed to parse physical information from binary gadget snapshots

	"""

	def __init__(self,fp):

		assert type(fp)==file,"Call the open() method instead!!"

		self.fp = fp
		self._header = _gadget.getHeader(fp)

		#Scale box to kpc
		self._header["box_size"] *= kpc

		#Scale masses to correct units
		self._header["masses"] *= (1.989e43 / self._header["h"])
		self._header["masses"] *= g 

	@classmethod
	def open(cls,filename):

		"""
		Opens a gadget snapshot at filename

		:param filename: file name of the gadget snapshot
		:type filename: str. or file.
		"""

		if type(filename)==str:
			fp = open(filename,"r")
		elif type(filename)==file:
			fp = filename
		else:
			raise TypeError("filename must be string or file!")
		
		return cls(fp)

	def header(self):

		"""
		Displays the snapshot header information

		:returns: the snapshot header information in dictionary form
		:rtype: dict.

		"""

		return self._header

	def close(self):

		"""
		Closes the snapshot file

		"""

		self.fp.close()