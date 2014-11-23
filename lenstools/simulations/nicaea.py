from astropy.cosmology import w0waCDM

try:
	from ..extern import _nicaea
	_nicaea=_nicaea
except ImportError:
	_nicaea=None


##########################################
##########NicaeaSettings class############
##########################################

class NicaeaSettings(dict):

	"""
	Class handler of the code settings (non linear modelling, tomography, transfer function, etc...)
	"""




##########################################
##########Nicaea class####################
##########################################

class Nicaea(w0waCDM):

	"""
	Main class handler for the python bindings of the NICAEA cosmological code, written by M. Kilbinger & collaborators

	"""

	def __init__(self,H0=72.0,Om0=0.26,Ode0=0.74,w0=-1.0,wa=0.0,sigma8=0.798,ns=0.960,name=None):

		super(Nicaea,self).__init__(H0,Om0,Ode0,w0=w0,wa=wa,name=name)
		self.sigma8=sigma8
		self.ns=ns

	def __repr__(self):

		astropy_string = super(Nicaea,self).__repr__()
		pieces = astropy_string.split(",")
		si8_piece = u" sigma8={0}".format(self.sigma8)
		ns_piece = u" ns={0}".format(self.ns)

		return ",".join(pieces[:3] + [si8_piece,ns_piece] + pieces[3:])