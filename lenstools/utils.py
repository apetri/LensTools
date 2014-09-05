from __future__ import division
import numpy as np

#####################################################################################
#######Disclaimer, I do not own the rights on this implementation####################
#######it is mainly for supplying to the lack of its implementation in numpy<1.8#####
#####################################################################################

def rfftfreq(n, d=1.0):
    """
Return the Discrete Fourier Transform sample frequencies
(for usage with rfft, irfft).

The returned float array `f` contains the frequency bin centers in cycles
per unit of the sample spacing (with zero at the start). For instance, if
the sample spacing is in seconds, then the frequency unit is cycles/second.

Given a window length `n` and a sample spacing `d`::

f = [0, 1, ..., n/2-1, n/2] / (d*n) if n is even
f = [0, 1, ..., (n-1)/2-1, (n-1)/2] / (d*n) if n is odd

Unlike `fftfreq` (but like `scipy.fftpack.rfftfreq`)
the Nyquist frequency component is considered to be positive.

Parameters
----------
n : int
Window length.
d : scalar, optional
Sample spacing (inverse of the sampling rate). Defaults to 1.

Returns
-------
f : ndarray
Array of length ``n//2 + 1`` containing the sample frequencies.

Examples
--------
>>> signal = np.array([-2, 8, 6, 4, 1, 0, 3, 5, -3, 4], dtype=float)
>>> fourier = np.fft.rfft(signal)
>>> n = signal.size
>>> sample_rate = 100
>>> freq = np.fft.fftfreq(n, d=1./sample_rate)
>>> freq
array([ 0., 10., 20., 30., 40., -50., -40., -30., -20., -10.])
>>> freq = np.fft.rfftfreq(n, d=1./sample_rate)
>>> freq
array([ 0., 10., 20., 30., 40., 50.])

"""
    if not (isinstance(n,int) or isinstance(n, integer)):
        raise ValueError("n should be an integer")
    val = 1.0/(n*d)
    N = n//2 + 1
    results = np.arange(0, N, dtype=int)
    return results * val

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