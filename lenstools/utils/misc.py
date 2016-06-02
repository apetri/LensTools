from __future__ import division

import numpy as np
from scipy import special as sp,integrate

####################################################################
#################Hankel transforms##################################
####################################################################

def fht(n,l_binned,func,**kwargs):

	if(kwargs.has_key('theta')):
		theta = kwargs['theta']
	else:
		theta_min = 1.0/l_binned.max()
		theta = l_binned*(theta_min/l_binned.min())

	h_kernel = sp.jn(n,np.outer(l_binned,theta))
	
	integrand = np.dot(np.diag(l_binned*func),h_kernel) * (2*np.pi)
	transform = integrate.simps(integrand,l_binned,axis=0)
	
	return theta,transform

def ifht(n,l_binned,func,**kwargs):
	
	if(kwargs.has_key('theta')):
		theta = kwargs['theta']
	else:
		theta_min = 1.0/l_binned.max()
		theta = l_binned*(theta_min/l_binned.min())
	
	h_kernel = sp.jn(n,np.outer(l_binned,theta))
	
	integrand = np.dot(np.diag(l_binned*func),h_kernel) / (2*np.pi)
	transform = integrate.simps(integrand,l_binned,axis=0)
	
	return theta,transform


######################################################################################
#################Approximate key matching dictionary##################################
######################################################################################

class ApproxDict(dict):

	"""
	Dictionary that looks up the closest key in the keychain

	"""

	def __init__(self,*args,**kwargs):
		super(ApproxDict,self).__init__(*args,**kwargs)
		self._sorted_keys = np.sort(self.keys())

	def __getitem__(self,key):
		closest_key = self._sorted_keys[np.abs(self._sorted_keys - key).argmin()]
		return super(ApproxDict,self).__getitem__(closest_key)
