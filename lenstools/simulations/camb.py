import sys,os
import re

if sys.version_info.major>=3:
	import _pickle as pkl
	from io import StringIO
else:
	import cPickle as pkl
	from StringIO import StringIO

import numpy as np
from scipy import interpolate
import astropy.units as u
from astropy.cosmology import FLRW

#Option parsing method
from .settings import select_parser,LTSettings

#Parse CAMB output log
def parseLog(fname):

	"""
	Parse CAMB output log

	:param fname: file name or file descriptor
	:type fname: str. or file.

	:returns: parsed log
	:rtype: dict.

	"""

	#Get the filehandle
	if type(fname)==file:
		fp = fname
	else:
		fp = open(fname,"r")

	#Dictionary with parsed log
	parsed = dict()
	parsed["sigma8"] = dict()

	#Cycle over the lines in the log
	for line in fp.readlines():

		#w0/wa
		match = re.match(r"\(w0, wa\) = \(([-\.0-9]+),[\s]+([-\.0-9]+)\)",line) 
		if match:
			parsed["w0"],parsed["wa"] = [ float(v) for v in match.groups() ]
			continue

		#Parameters
		match = re.match(r"Reion redshift[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["z_ion"] = float(match.groups()[0])

		match = re.match(r"Om_b h\^2[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Obh2"] = float(match.groups()[0])

		match = re.match(r"Om_c h\^2[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Omch2"] = float(match.groups()[0])

		match = re.match(r"Om_nu h\^2[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Onuh2"] = float(match.groups()[0])

		match = re.match(r"Om_Lambda[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Ode"] = float(match.groups()[0])

		match = re.match(r"Om_K[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Omk"] = float(match.groups()[0])

		match = re.match(r"Om_m \(1-Om_K-Om_L\)[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Om"] = float(match.groups()[0])

		match = re.match(r"100 theta \(CosmoMC\)[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["100thetaMC"] = float(match.groups()[0])

		match = re.match(r"Reion opt depth[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["tau_ion"] = float(match.groups()[0])

		match = re.match(r"Age of universe\/GYr[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["Age"] = float(match.groups()[0]) * u.Gyr

		match = re.match(r"zstar[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["zstar"] = float(match.groups()[0])

		match = re.match(r"r_s\(zstar\)/Mpc[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["rs"] = float(match.groups()[0]) * u.Mpc

		match = re.match(r"zdrag[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["zdrag"] = float(match.groups()[0])

		match = re.match(r"r_s\(zdrag\)/Mpc[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["rs(zdrag)"] = float(match.groups()[0]) * u.Mpc

		match = re.match(r"k_D\(zstar\) Mpc[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["kD(zstar)"] = float(match.groups()[0]) / u.Mpc

		match = re.match(r"100\*theta_D[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["100thetaD"] = float(match.groups()[0])

		match = re.match(r"z_EQ \(if v_nu=1\)[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["zEQ"] = float(match.groups()[0])

		match = re.match(r"100\*theta_EQ[\s]+=[\s]+([0-9\.]+)",line)
		if match:
			parsed["100thetaEQ"] = float(match.groups()[0])

		match = re.match(r"tau_recomb/Mpc[\s]+=[\s]+([0-9\.]+)[\s]+tau_now/Mpc =[\s]+([0-9\.]+)",line)
		if match:
			parsed["tau_rec"],parsed["tau_now"] = [float(v)*u.Mpc for v in match.groups()]

		match = re.match(r"[\s]+at z =[\s]+([0-9E\-\+\.]+)[\s]+sigma8 \(all matter\)=[\s]+([0-9\.]+)",line)
		if match:
			z,sigma8 = [ float(v) for v in match.groups() ]
			parsed["sigma8"][z] = sigma8

	#Return
	if type(fname)!=file:
		fp.close()
	
	return parsed



##################################################################################################

class CAMBSettings(LTSettings):

	def __init__(self,**kwargs):

		self.get_scalar_cls = True
		self.get_vector_cls = False
		self.get_tensor_cls = False
		self.get_transfer = True
		self.do_lensing = False
		self.do_nonlinear = 0
		self.l_max_scalar = 8000
		self.k_eta_max_scalar = 16000
		self.l_max_tensor = 1500
		self.k_eta_max_tensor = 3000
		self.use_physical = True

		#####################################
		
		self.cs2_lam = 1

		#####################################

		self.helium_fraction = 0.24
		self.nu_mass_eigenstates = 0
		self.nu_mass_degeneracies = 0
		self.share_delta_neff = True

		self.scalar_amplitude = 2.41e-9
		self.pivot_scalar = 0.002 * u.Mpc**-1
		self.pivot_tensor = 0.002 * u.Mpc**-1

		#################################################
		
		self.reionization = True
		self.re_use_optical_depth = True
		self.re_optical_depth = 0.087
		self.re_redshift = 11
		self.re_delta_redshift = 0.5
		self.re_ionization_frac = -1
		self.recfast_fudge = 1.14
		self.recfast_fudge_he = 0.86
		self.recfast_heswitch = 6
		self.initial_condition = 1

		###########################################################
		
		self.initial_vector = np.array([-1,0,0,0,0])

		###########################################################

		self.vector_mode = 0
		self.cobe_normalize = False
		self.cmb_outputscale = 7.4311e12
		self.transfer_high_precision = True
		self.transfer_kmax = 1000
		self.transfer_k_per_logint = 100

		################################################################
		
		self.transfer_interp_matterpower = True

		#############################################################

		self.transfer_power_var = 7
		self.scalar_output_file = "scalCls.dat"
		self.vector_output_file = "vecCls.dat"
		self.tensor_output_file = "tensCls.dat"
		self.total_output_file = "totCls.dat"
		self.lensed_output_file = "lensedCls.dat"
		self.lensed_total_output_file = "lensedtotCls.dat"
		self.fits_filename = "scalCls.fits"

		###############################################################

		self.feedback_level = 1
		self.lensing_method = 1
		self.accurate_bb = True
		self.massive_nu_approx = 3
		self.accurate_polarization = True
		self.accurate_reionization = True
		self.do_tensor_neutrinos = False
		self.do_late_rad_truncation = False
		self.number_of_threads = 0
		self.accuracy_boost = 3
		self.l_accuracy_boost = 3
		self.l_sample_boost = 3

		###############################################################

		#Allow for kwargs override
		for key in kwargs.keys():
			setattr(self,key,kwargs[key])


	def write(self,output_root,cosmology,redshifts):

		"""
		Writes a CAMB parameter file

		:param output_root: output_root for the files that CAMB will produce in output
		:type output_root: str. 

		:param cosmology: cosmological model to generate the parameter file for
		:type cosmology: FLRW

		:param redshifts: redshifts on which to compute the matter power spectrum and transfer function
		:type redshifts: array.

		:returns: string object
		:rtype: StringIO

		"""

		#Safety type check
		assert isinstance(cosmology,FLRW)

		#Sort the redshifts in chronological order
		z = -1.0*np.sort(-1.0*redshifts)

		#Instantiate StringIO object
		s = StringIO()

		s.write("output_root = {0}\n".format(output_root))

		s.write("\n\n#####################################\n\n")

		s.write('get_scalar_cls = {0}\n'.format(self.get_scalar_cls.__str__()[0]))
		s.write('get_vector_cls = {0}\n'.format(self.get_vector_cls.__str__()[0]))
		s.write('get_tensor_cls = {0}\n'.format(self.get_tensor_cls.__str__()[0]))
		s.write('get_transfer = {0}\n'.format(self.get_transfer.__str__()[0]))
		s.write('do_lensing = {0}\n'.format(self.do_lensing.__str__()[0]))
		s.write('do_nonlinear = {0}\n'.format(self.do_nonlinear))
		s.write('l_max_scalar = {0}\n'.format(self.l_max_scalar))
		s.write('k_eta_max_scalar = {0}\n'.format(self.k_eta_max_scalar))
		s.write('l_max_tensor = {0}\n'.format(self.l_max_tensor))
		s.write('k_eta_max_tensor = {0}\n'.format(self.k_eta_max_tensor))
		
		s.write("\n\n#####################################\n\n")
		s.write('use_physical = {0}\n'.format(self.use_physical.__str__()[0]))

		#############################################
		#######Cosmological parameters###############
		#############################################

		#Baryon and dark matter densities
		s.write("ombh2 = {0:.6f}\n".format(cosmology.Ob0*(cosmology.h**2)))
		s.write("omch2 = {0:.6f}\n".format((cosmology.Om0 - cosmology.Ob0)*(cosmology.h**2)))

		#Neutrino density
		if cosmology._nmassivenu==0:
			omnuh2 = 0.0
		else:
			omnuh2 = cosmology.Onu0 * (cosmology.h**2)
		s.write("omnuh2 = {0:.6f}\n".format(omnuh2))
		
		#Curvature parameter (enforce Om+Ol=1-Ok)
		s.write("omk = {0:.6f}\n".format(1-cosmology.Om0-cosmology.Ode0))

		#Hubble constant
		s.write("hubble = {0:.6f}\n".format(cosmology.h * 100))

		#Dark energy parameters
		if hasattr(cosmology,"w0"):
			w0 = cosmology.w0
		else:
			w0 = -1.

		if hasattr(cosmology,"wa"):
			wa = cosmology.wa
		else:
			wa = 0.

		s.write("w = {0:.6f}\n".format(w0))
		s.write("wa = {0:.6f}\n".format(wa))

		s.write("\n\n#####################################\n\n")

		s.write('cs2_lam = {0}\n'.format(self.cs2_lam))
		s.write('temp_cmb = {0:.3f}\n'.format(cosmology.Tcmb0.to(u.K).value))

		s.write("\n\n#####################################\n\n")
			
		s.write('helium_fraction = {0}\n'.format(self.helium_fraction))
		s.write('massless_neutrinos = {0}\n'.format(cosmology.Neff-cosmology._nmassivenu))
		s.write('massive_neutrinos = {0}\n'.format(cosmology._nmassivenu))
		s.write('nu_mass_eigenstates = {0}\n'.format(self.nu_mass_eigenstates))
		s.write('nu_mass_degeneracies = {0}\n'.format(self.nu_mass_degeneracies))
		
		#Compute the mass fractions of the massive species
		if cosmology._nmassivenu:
			fractions = (cosmology.m_nu / cosmology.m_nu.sum()).decompose().value
		else:
			fractions = 1
		s.write('nu_mass_fractions = {0}\n'.format(fractions.__str__().strip("[").strip("]")))
		
		s.write('share_delta_neff = {0}\n'.format(self.share_delta_neff.__str__()[0]))

		s.write("\n\n#####################################\n\n")

		#############################################
		#######Spectral index tilt ##################
		#############################################
			
		s.write('pivot_scalar = {0:.3f}\n'.format(self.pivot_scalar.to(u.Mpc**-1).value))
		s.write('pivot_tensor = {0:.3f}\n'.format(self.pivot_tensor.to(u.Mpc**-1).value))

		s.write('initial_power_num = {0}\n'.format(1))
		s.write('scalar_amp(1) = {0:.6e}\n'.format(self.scalar_amplitude))
		
		if hasattr(cosmology,"ns"):
			ns = cosmology.ns
		else:
			ns = 1.0
		
		s.write('scalar_spectral_index(1) = {0:.6f}\n'.format(ns))
		
		s.write('scalar_nrun(1) = {0}\n'.format(0))
		s.write('tensor_spectral_index(1) = {0}\n'.format(0))
		s.write('initial_ratio(1) = {0}\n'.format(0))

		s.write("\n\n#####################################\n\n")

		s.write('reionization = {0}\n'.format(self.reionization.__str__()[0]))
		s.write('re_use_optical_depth = {0}\n'.format(self.re_use_optical_depth.__str__()[0]))
		s.write('re_optical_depth = {0:.3f}\n'.format(self.re_optical_depth))
		s.write('re_redshift = {0}\n'.format(self.re_redshift))
		s.write('re_delta_redshift = {0:.2f}\n'.format(self.re_delta_redshift))
		s.write('re_ionization_frac = {0}\n'.format(self.re_ionization_frac))

		s.write("\n\n#####################################\n\n")
			
		s.write('RECFAST_fudge = {0}\n'.format(self.recfast_fudge))
		s.write('RECFAST_fudge_he = {0}\n'.format(self.recfast_fudge_he))
		s.write('RECFAST_heswitch = {0}\n'.format(self.recfast_heswitch))

		s.write("\n\n#####################################\n\n")
			
		s.write('initial_condition = {0}\n'.format(self.initial_condition))
		s.write('initial_vector = {0}\n'.format(self.initial_vector.__str__().strip("[").strip("]")))

		s.write("\n\n#####################################\n\n")
			
		s.write('vector_mode = {0}\n'.format(self.vector_mode))
		s.write('cobe_normalize = {0}\n'.format(self.cobe_normalize.__str__()[0]))

		s.write("\n\n#####################################\n\n")

		s.write('cmb_outputscale = {0:4e}\n'.format(self.cmb_outputscale))

		s.write("\n\n#####################################\n\n")
			
		s.write('transfer_high_precision = {0}\n'.format(self.transfer_high_precision.__str__()[0]))
		s.write('transfer_kmax = {0}\n'.format(self.transfer_kmax))
		s.write('transfer_k_per_logint = {0}\n'.format(self.transfer_k_per_logint))
		s.write('transfer_interp_matterpower = {0}\n'.format(self.transfer_interp_matterpower.__str__()[0]))

		#############################################
		#######Transfer function ####################
		#############################################

		s.write("transfer_num_redshifts = {0}\n\n".format(len(z)))
		for n in range(len(z)):
			s.write("transfer_redshift({0}) = {1:.6f}\n".format(n+1,z[n]))
			s.write("transfer_filename({0}) = transferfunc_z{1:.6f}.dat\n".format(n+1,z[n]))
			s.write("transfer_matterpower({0}) = matterpower_z{1:.6f}.dat\n".format(n+1,z[n]))

		s.write("\n\n#####################################\n\n")

		s.write('transfer_power_var = {0}\n'.format(self.transfer_power_var))

		s.write("\n\n#####################################\n\n")
			
		s.write('scalar_output_file = {0}\n'.format(self.scalar_output_file))
		s.write('vector_output_file = {0}\n'.format(self.vector_output_file))
		s.write('tensor_output_file = {0}\n'.format(self.tensor_output_file))
		s.write('total_output_file = {0}\n'.format(self.total_output_file))
		s.write('lensed_output_file = {0}\n'.format(self.lensed_output_file))
		s.write('lensed_total_output_file = {0}\n'.format(self.lensed_total_output_file))
		s.write('fits_filename = {0}\n'.format(self.fits_filename))

		s.write("\n\n#####################################\n\n")
			
		s.write('feedback_level = {0}\n'.format(self.feedback_level))
		s.write('lensing_method = {0}\n'.format(self.lensing_method))
		s.write('accurate_bb = {0}\n'.format(self.accurate_bb.__str__()[0]))

		s.write("\n\n#####################################\n\n")

		s.write('massive_nu_approx = {0}\n'.format(self.massive_nu_approx))

		s.write("\n\n#####################################\n\n")

		s.write('accurate_polarization = {0}\n'.format(self.accurate_polarization.__str__()[0]))
		s.write('accurate_reionization = {0}\n'.format(self.accurate_reionization.__str__()[0]))
		s.write('do_tensor_neutrinos = {0}\n'.format(self.do_tensor_neutrinos.__str__()[0]))
		s.write('do_late_rad_truncation = {0}\n'.format(self.do_late_rad_truncation.__str__()[0]))

		s.write("\n\n#####################################\n\n")

		s.write('number_of_threads = {0}\n'.format(self.number_of_threads))
		s.write('accuracy_boost = {0}\n'.format(self.accuracy_boost))
		s.write('l_accuracy_boost = {0}\n'.format(self.l_accuracy_boost))
		s.write('l_sample_boost = {0}\n'.format(self.l_sample_boost))

		s.seek(0)
		return s.read()

#########################################################################################################################################

#CAMB transfer function
class TransferFunction(object):

	def __init__(self,k):

		"""
		:param k: wavenumbers at which the transfer function is computed at
		:type k: quantity

		"""

		assert k.unit.physical_type=="wavenumber"
		self._k = k.to((u.Mpc)**-1)
		self._transfer = dict()
		self._interpolated = dict()

	def add(self,z,T):

		"""
		Add transfer function information at redshift z

		:param z: redshift
		:type z: float.

		:param T: CDM transfer function from CAMB output
		:type T: array 

		"""

		if hasattr(self,"_sorted_z"):
			del(self._sorted_z)

		assert T.shape==self._k.shape,"There should be exactly one transfer function value for each wavenumber! len(T)={0} len(k)={1}".format(len(T),len(self._k))
		self._transfer[z] = T

	def __getitem__(self,z):

		"""
		Returns the tabulated transfer function at z. If z is not in the table, returns the tabulated transfer function at the closest z available

		:param z: redshift at which to output the tabulated transfer function
		:type z: float.

		:returns: (tabulated z,k,tabulated transfer function)
		:rtype: tuple.

		"""

		#If the transfer function is not tabulated with z, use the closest z in the table
		if not hasattr(self,"_sorted_z"):
			self._sorted_z = np.sort(np.array(list(self._transfer.keys())))

		if z in self._transfer:
			zt = z
		else:
			zt = self._sorted_z[np.abs(self._sorted_z - z).argmin()] 

		#Return
		return zt,self._k,self._transfer[zt]



	def __call__(self,z,k):


		"""
		Compute the transfer function at redshift z by linear interpolation

		:param z: redshift
		:type z: float.

		:param k: wavenumbers at which to compute the transfer function (linearly interpolated with scipy.interp1d)
		:type k: quantity

		:returns: transfer function at k
		:rtype: array 

		"""

		assert k.unit.physical_type=="wavenumber"

		#If the transfer function is not tabulated with z, use the closest z in the table
		if not hasattr(self,"_sorted_z"):
			self._sorted_z = np.sort(np.array(list(self._transfer.keys())))

		if z in self._transfer:
			zt = z
		else:
			zt = self._sorted_z[np.abs(self._sorted_z - z).argmin()] 

		#If interpolator has not been built yet for the current redshift, build it
		if zt not in self._interpolated:
			self._interpolated[zt] = interpolate.interp1d(self._k.value,self._transfer[zt],fill_value=1,bounds_error=False)

		#Use interpolator to compute the transfer function
		return self._interpolated[zt](k.to((u.Mpc)**-1).value)

	#I/O
	def save(self,filename):

		"""
		Pickle the TransferFunction instance

		:param filename: name of the file to save the instance to
		:type filename: str.

		"""

		with open(filename,"wb") as fp:
			pkl.dump(self,fp,protocol=2)

	@classmethod
	def read(cls,filename):

		"""
		Load a previously pickled TransferFunction instance

		:param filename: name of the file from which the instance should be read
		:type filename: str.

		:rtype: :py:class:`TransferFunction`

		"""

		with open(filename,"rb") as fp:
			tfr = pkl.load(fp)

		if isinstance(tfr,cls):
			return tfr
		else:
			raise TypeError("Pickled instance is not of type {0}".format(cls.__name__))


##############################################################################################################################

class CAMBTransferFunction(TransferFunction):
	pass

class CAMBTransferFromPower(TransferFunction):

	def add(self,z,T):

		"""
		Add transfer function information at redshift z

		:param z: redshift
		:type z: float.

		:param T: CDM transfer function from CAMB output
		:type T: array 

		"""

		if hasattr(self,"_sorted_z"):
			del(self._sorted_z)

		assert T.shape==self._k.shape,"There should be exactly one transfer function value for each wavenumber! len(T)={0} len(k)={1}".format(len(T),len(self._k))
		self._transfer[z] = np.sqrt(T) 

##############################################################################################################################

#k independent transfer function for testing D(z,k) = 1/1+z
class TestTransferFunction(TransferFunction):

	def __call__(self,z,k):
		return np.ones(k.shape)/(1+z)