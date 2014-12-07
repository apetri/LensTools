import os,sys

import numpy as np
import matplotlib.pyplot as plt

try:
	
	from lenstools import Ensemble
	from lenstools.constraints import FisherAnalysis,LikelihoodAnalysis,Emulator
	from lenstools.simulations import CFHTemu1

except ImportError:

	sys.path.append("..")

	from lenstools import Ensemble
	from lenstools.constraints import FisherAnalysis,LikelihoodAnalysis,Emulator
	from lenstools.simulations import CFHTemu1


#Test Fisher analysis with power spectrum
def test_fisher():

	f = FisherAnalysis()

	#Load the power spectral features and parameters into the analysis instance
	
	#Fiducial
	l,P,g = np.loadtxt("Data/fiducial",unpack=True)
	covariance = P**2 / (l + 0.5)
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
	assert f.varied == range(3)

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

	#Sample fit with diagonal covariance matrix
	fitted_parameters = f.fit(P,features_covariance=covariance)
	np.savetxt("fitted_parameters.txt",fitted_parameters)

	#Sample fisher matrix
	fisher = f.fisher_matrix(simulated_features_covariance=covariance)
	np.savetxt("fisher_constraints.txt",np.sqrt(np.linalg.inv(fisher).diagonal()))

	return f

#Test interpolation of power spectrum
def test_interpolation():

	root_path = "Data/all"
	analysis = LikelihoodAnalysis()

	#Read in model names
	models = CFHTemu1.getModels()[:17]
	assert len(models) == 17

	#Shuffle the models
	np.random.seed(1)
	np.random.shuffle(models)

	#Divide into training and testing
	training_models = models[:-1]
	testing_model = models[-1]

	#Read multipoles
	ell = np.load(os.path.join(root_path,"ell.npy"))

	#Load in the means of the power spectra of the 17 models, and populate the analysis instance
	for model in training_models:

		ens = Ensemble.fromfilelist([os.path.join(root_path,model._cosmo_id_string,"subfield1","sigma05","power_spectrum.npy")])
		ens.load(from_old=True)

		analysis.add_model(parameters=model.squeeze(with_ns=True),feature=ens.mean())

	#Add the multipoles to the analysis
	analysis.add_feature_label(ell)
	l = analysis.feature_label

	ens = Ensemble.fromfilelist([os.path.join(root_path,testing_model._cosmo_id_string,"subfield1","sigma05","power_spectrum.npy")])
	ens.load(from_old=True)
	testing_Pl = ens.mean()

	#Load in also the observed power spectrum
	ens = Ensemble.fromfilelist([os.path.join(root_path,"observations","subfield1","sigma05","power_spectrum.npy")])
	ens.load(from_old=True)
	observed_Pl = ens.mean() 

	#Output the analysis stats
	np.savetxt("16_parameter_points.txt",analysis.parameter_set)

	for n in range(len(training_models)):

		plt.plot(l,l*(l+1)*analysis.training_set[n]/(2*np.pi))

	plt.plot(l,l*(l+1)*observed_Pl/(2*np.pi),linestyle="--",label="Observation")	

	plt.xlabel(r"$l$")
	plt.ylabel(r"$l(l+1)P_l/2\pi$")
	plt.yscale("log")

	plt.legend(loc="upper left")

	plt.savefig("16_power_spectra.png")
	plt.clf()

	#Train the interpolators
	analysis.train(use_parameters=range(3))
	assert hasattr(analysis,"_interpolator")
	assert hasattr(analysis,"_num_bins")

	#Emulator portability test with pickle/unpickle
	analysis.save("analysis.p")
	emulator = LikelihoodAnalysis.load("analysis.p")

	#Predict the power spectrum at the remaining point
	predicted_Pl = emulator.predict(testing_model.squeeze())

	#Plot it against the measured one
	fig,ax = plt.subplots(2,1,figsize=(16,8))

	#Measured
	ax[0].plot(l,l*(l+1)*testing_Pl/(2*np.pi),label="measured")

	#Predicted
	ax[0].plot(l,l*(l+1)*predicted_Pl/(2*np.pi),label="interpolated")
	
	#Fractional difference
	ax[1].plot(l,(predicted_Pl - testing_Pl)/testing_Pl)

	ax[1].set_xlabel(r"$l$")
	ax[0].set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax[1].set_ylabel(r"$P_l^I-P_l^M/P_l^M$")
	
	ax[0].set_yscale("log")
	ax[0].legend(loc="upper left")

	plt.savefig("power_interpolator_test.png")
	plt.clf()

	#Give it a shot with two points in parameter space to test vectorization
	two_parameter_points = np.array((training_models[0].squeeze(),testing_model.squeeze()))
	two_predicted_Pl = emulator.predict(two_parameter_points)

	fig,ax = plt.subplots(2,1,figsize=(16,8))

	#Predicted
	ax[0].plot(l,l*(l+1)*two_predicted_Pl[0]/(2*np.pi),color="red",linestyle="--")
	ax[0].plot(l,l*(l+1)*two_predicted_Pl[1]/(2*np.pi),color="green",linestyle="--")

	#Measured
	ax[0].plot(l,l*(l+1)*emulator.training_set[0]/(2*np.pi),color="red",linestyle="-")
	ax[0].plot(l,l*(l+1)*testing_Pl/(2*np.pi),color="green",linestyle="-")

	#Fractional difference
	ax[1].plot(l,(two_predicted_Pl[0] - emulator.training_set[0])/emulator.training_set[0],color="red")
	ax[1].plot(l,(two_predicted_Pl[1] - testing_Pl)/testing_Pl,color="green")

	ax[1].set_xlabel(r"$l$")
	ax[0].set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax[1].set_ylabel(r"$P_l^I-P_l^M/P_l^M$")
	
	ax[0].set_yscale("log")

	plt.savefig("power_interpolator_test_2.png")
	plt.clf()

	#Generate a fudge power spectrum covariance matrix
	covariance = np.diag(testing_Pl**2/(0.5 + l))

	#Generate a fudge observation by wiggling the testing power spectrum
	observation = testing_Pl + np.random.uniform(low=-testing_Pl*0.1,high=testing_Pl*0.1)

	#Choose a bunch of points in parameter space
	points = emulator.parameter_set[:,:-1]

	#Compute the chi2
	chi2_values_1 = emulator.chi2(points,observation,covariance)
	chi2_values_2 = emulator.chi2(points,observation,covariance,split_chunks=4)

	assert chi2_values_1.shape == chi2_values_2.shape

	#Compute the individual contributions
	chi2_contributions = emulator.chi2Contributions(points[0],observation,covariance)

	#Plot
	plt.imshow(chi2_contributions,interpolation="nearest")
	plt.colorbar()
	plt.xlabel(r"$j$")
	plt.ylabel(r"$i$")
	plt.savefig("chi2_contributions.png")
	plt.clf()

	return chi2_values_1,chi2_values_2

def test_remove():

	emulator = LikelihoodAnalysis.load("analysis.p")
	emulator.remove_model([8,10])
	emulator.train(use_parameters=range(3))

	assert emulator.parameter_set.shape[0] == 14
	assert emulator.training_set.shape[0] == 14 

def test_find():

	emulator = LikelihoodAnalysis.load("analysis.p")
	parameters_to_find = emulator.parameter_set[7]

	n = emulator.find(parameters_to_find)
	assert len(n)==1
	assert n[0] == 7

def test_emulator():

	#Unpickle the emulator
	emulator = Emulator.load("analysis.p")
	emulator.train()

	#Set the model
	emulator.set_to_model(np.array([ 0.26,-2.66,1.31,0.96]))
	ell = emulator.feature_label
	Pell = emulator._current_predicted_feature

	#Select the new multipoles
	l = np.arange(900.0,3000.0,200.0)
	#Emulate the power spectrum
	Pl = emulator.emulate(l)

	#Plot
	fig,ax = plt.subplots()
	ax.plot(ell,ell*(ell+1)*Pell/(2*np.pi),label="Fully emulated")
	ax.plot(l,l*(l+1)*Pl/(2.0*np.pi),label="New multipoles",color="yellow")

	ax.set_xlabel(r"$l$")
	ax.set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax.set_yscale("log")
	ax.legend(loc="upper left")

	fig.savefig("emulated_power.png")


def test_reparametrize():

	#Unpickle the emulator
	emulator = Emulator.load("analysis.p")

	#Get the old parametrizations
	Om = emulator.parameter_set[:,0]
	w = emulator.parameter_set[:,1]
	si8 = emulator.parameter_set[:,2]

	#Define the reparametrizer ( (Om,si8) -> (si8 x Om^0.5))
	def formatter(data,n):
		print("Omega_m exponent={0:.1f}".format(n))
		return np.array([data[:,1],data[:,2]*(data[:,0]**n)]).transpose()

	#Reparametrize the parameter space
	emulator.reparametrize(formatter,0.5)

	#Check that everything worked
	assert (emulator.parameter_set[:,0]==w).all()
	assert (emulator.parameter_set[:,1]==(si8*(Om**0.5))).all()








