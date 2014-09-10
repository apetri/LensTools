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

	def getPositions(self,first=None,last=None):

		"""
		Reads in the particles positions (read in of a subset is allowed): when first and last are specified, the numpy array convention is followed (i.e. getPositions(first=a,last=b)=getPositions()[a:b])

		:param first: first particle in the file to be read, if None 0 is assumed
		:type first: int. or None

		:param last: last particle in the file to be read, if None the total number of particles is assumed
		:type last: int. or None

		:returns: numpy array with the particle positions

		"""

		numPart = self._header["num_particles_file"]

		#Calculate the offset from the beginning of the file: 4 bytes (endianness) + 256 bytes (header) + 8 bytes (void)
		offset = 4 + 256 + 8

		#If first is specified, offset the file pointer by that amount
		if first is not None:
			
			assert first>=0
			offset += 4 * 3 * first
			numPart -= first

		if last is not None:

			if first is not None:
				
				assert last>=first and last<=self._header["num_particles_file"]
				numPart = last - first

			else:

				assert last<=self._header["num_particles_file"]
				numPart = last


		#Read in the particles positions and return the corresponding array
		return (_gadget.getPosVel(self.fp,offset,numPart) * kpc).to(Mpc)

	def close(self):

		"""
		Closes the snapshot file

		"""

		self.fp.close()