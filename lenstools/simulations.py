"""

.. module:: simulations 
	:platform: Unix
	:synopsis: This module handles the book keeping of a simulation set map names, cosmological parameters, etc...


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import division

from astropy.cosmology import FlatwCDM

######################################
###########IGS1 class#################
######################################

class IGS1(FlatwCDM):

	"""
	Class handler of the IGS1 simulations set, inherits the cosmological parameters from the astropy.cosmology.FlatwCDM class; the default parameter values are the fiducial ones

	"""

	def __init__(self,H0=70.0,Om0=0.26,w0=-1.0,sigma8=0.798,ns=0.960,root_path=None,name=None):

		super(IGS1,self).__init__(H0,Om0,w0=w0,name=name)
		self.sigma8 = sigma8
		self.ns = ns

		assert root_path is not None,"You must specify the root path of your IGS1 local copy!"

		self.root_path = root_path

		#Don't touch these!
		self._series_name = "m"
		self._num_particles = 512
		self._box_size_mpc = 240
		self._lens_plane_size = 4096 
		self._full_path = self.root_path.rstrip("/") + "/"+self._series_name+"-"+str(self._num_particles)+"b"+str(self._box_size_mpc)+"_Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}".format(self.Om0,1.0-self.Om0,self.w0,self.ns,self.sigma8)

	def __repr__(self):

		astropy_string = super(IGS1,self).__repr__()
		pieces = astropy_string.split(",")
		si8_piece = u" sigma8={0}".format(self.sigma8)
		ns_piece = u" ns={0}".format(self.ns)

		return ",".join(pieces[:3] + [si8_piece,ns_piece] + pieces[3:])

	def _realization_id(self,n):

		d1 = n//1000
		d2 = (n - 1000*d1)//100
		d3 = (n - 1000*d1 - 100*d2)//10
		d4 = n - 1000*d1 - 100*d2 - 10*d3

		return "{0}{1}{2}{3}".format(d1,d2,d3,d4)

	def _redshift_id(self,z):

		return self._realization_id(int(z*100))

	def _plane_id(self,z):

		if z==1.0:
			return "0029p"
		elif z==1.5:
			return "0038p"
		elif z==2.0:
			return "0046p"
		else:
			raise ValueError("IGS1 doesn't have maps at redshift {0}".format(z))

	def getNames(self,z,realizations,kind="convergence",big_fiducial_set=False):

		"""
		Get the full name of the IGS1 maps, once a redshift, realization identificator and kind are specified

		:param z: redshift plane of the maps, must be one of [1.0,1.5,2.0]
		:type z: float.

		:param realizations: list of realizations to get the names of, the elements must be in [1,1000]
		:type realizations: list.

		:param kind: decide if retrieve convergence or shear maps, must be one of [convergence,shear1,shear2]
		:type kind: str.

		:param big_fiducial_set: set to True if you want to get the names of the bigger fiducial simulation based on 45 N-body simulations
		:type big_fiducial_set: bool.

		"""

		assert type(realizations) == list
		assert z in [1.0,1.5,2.0],"IGS1 doesn't have maps at redshift {0}".format(z)
		assert kind in ["convergence","shear1","shear2"],"You must select one of these: convergence,shear1,shear2"

		if kind=="convergence":
			prefix = "WL-conv"
			direct = "Maps"
		elif kind=="shear1":
			prefix = "Wl-shear1"
			direct = "shear"
		elif kind=="shear2":
			prefix = "Wl-shear2"
			direct = "shear"

		full_path = self._full_path

		if big_fiducial_set:
			assert self.Om0==0.26 and self.w0==-1.0 and self.sigma8==0.798 and self.ns==0.96
			full_path += "_f"

		full_path += "/{0}".format(direct)

		return [full_path + "/{0}_".format(prefix)+self._series_name+"-"+str(self._num_particles)+"b"+str(self._box_size_mpc)+"_Om{0:.3f}_Ol{1:.3f}_w{2:.3f}_ns{3:.3f}_si{4:.3f}_".format(self.Om0,1.0-self.Om0,self.w0,self.ns,self.sigma8)+str(self._lens_plane_size)+"xy_{0}r_{1}_{2}z_og.gre.fit".format(self._realization_id(n),self._plane_id(z),self._redshift_id(z)) for n in realizations]






