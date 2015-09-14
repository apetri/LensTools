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

###########################################################################
###########Hack to make scipy interpolate objects pickleable###############
###########################################################################

class _interpolate_wrapper(object):

	def __init__(self,f,args,kwargs):
		self.f = f
		self.args = args
		self.kwargs = kwargs
	
	def __call__(self):
		try:
			return self.f(*self.args,**self.kwargs)
		except:
			import traceback
			print("lenstools: Exception while building the interpolators")
			print(" exception:")
			traceback.print_exc()
			raise 
