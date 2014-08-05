import sys

import numpy as np
import matplotlib.pyplot as plt

try:
	
	from lenstools.constraints import FisherAnalysis

except ImportError:

	sys.path.append("..")
	from lenstools.constraints import FisherAnalysis


#Test Fisher analysis with power spectrum
def test_fisher():

	f = FisherAnalysis()

	#Load the power spectral features and parameters into the analysis instance
	
	#Fiducial
	l,P,g = np.loadtxt("Data/fiducial",unpack=True)
	f.add_model(parameters=np.array([0.26,-1.0,0.798,0.960]),feature=P)

	#High Om
	l,P,g = np.loadtxt("Data/v_Om_0.29",unpack=True)
	f.add_model(parameters=np.array([0.29,-1.0,0.798,0.960]),feature=P)

	#High w
	l,P,g = np.loadtxt("Data/v_w_-0.8",unpack=True)
	f.add_model(parameters=np.array([0.26,-0.8,0.798,0.960]),feature=P)

	#High sigma8
	l,P,g = np.loadtxt("Data/v_si8_0.850",unpack=True)
	f.add_model(parameters=np.array([0.26,-1.0,0.850,0.960]),feature=P)

	#Sanity checks
	assert f.parameter_set.shape[0] == 4
	assert f.training_set.shape[0] == 4

	#Compute the derivatives
	derivatives = f.compute_derivatives()

	#Sanity check
	assert derivatives.shape == (3,len(l))
	assert f.get_varied() == range(3)

	#Plot the results
	for n in range(3): 
		plt.plot(l,np.abs(derivatives[n])*2*np.pi/(l*(l+1)),label="Parameter {0}".format(n+1))

	plt.yscale("log")
	plt.xlim(500,5000)
	plt.ylim(1.0e-11,1.0e-8)
	plt.xlabel(r"$l$")
	plt.ylabel(r"$\partial P_l/\partial p$")
	plt.legend(loc="upper right")

	plt.savefig("derivatives.png")
	plt.clf()

	return f


