from __future__ import division

from lenstools.external import _gadget

import numpy as np
from astropy.units import kpc,Mpc,g

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
		#Convert to Mpc
		self._header["box_size"] = self._header["box_size"].to(Mpc)

		#Scale masses to correct units
		self._header["masses"] *= (1.989e43 / self._header["h"])
		self._header["masses"] *= g 

		#Update the dictionary with the number of particles per side
		self._header["num_particles_total_side"] = int(np.round(self._header["num_particles_total"]**(1/3)))

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

	@property
	def header(self):

		"""
		Displays the snapshot header information

		:returns: the snapshot header information in dictionary form
		:rtype: dict.

		"""

		return self._header

	def getPositions(self):

		"""
		Reads in the particles positions

		"""

		#Calculate the offset from the beginning of the file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
		offset = 4 + 256 + 8

		#Read in the particles positions and return the corresponding array
		return (_gadget.getPosVel(self.fp,offset,self._header["num_particles_file"]) * kpc).to(Mpc)

	def close(self):

		"""
		Closes the snapshot file

		"""

		self.fp.close()