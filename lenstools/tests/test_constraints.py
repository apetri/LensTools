import sys,os
sys.modules["mpi4py"] = None

from .. import dataExtern

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

	
from .. import Ensemble
from ..statistics.constraints import FisherAnalysis,FisherSeries,Emulator,EmulatorSeries
from ..statistics.contours import ContourPlot
from ..simulations import CFHTemu1


#Test Fisher analysis with power spectrum
def test_fisher():

	#Load the power spectral features and parameters
	parameters = list()
	features = list()
	
	#Fiducial
	parameters.append(np.array([0.26,-1.0,0.798,0.960]))
	l,P,g = np.loadtxt(os.path.join(dataExtern(),"fiducial"),unpack=True)
	features.append(P)
	
	covariance = P**2 / (l + 0.5)

	#High Om
	parameters.append(np.array([0.29,-1.0,0.798,0.960]))
	l,P,g = np.loadtxt(os.path.join(dataExtern(),"v_Om_0.29"),unpack=True)
	features.append(P)

	#High w
	parameters.append(np.array([0.26,-0.8,0.798,0.960]))
	l,P,g = np.loadtxt(os.path.join(dataExtern(),"v_w_-0.8"),unpack=True)
	features.append(P)

	#High sigma8
	parameters.append(np.array([0.26,-1.0,0.850,0.960]))
	l,P,g = np.loadtxt(os.path.join(dataExtern(),"v_si8_0.850"),unpack=True)
	features.append(P)

	#Instantiate the FisherAnalysis
	feature_index = FisherSeries.make_index(pd.Index(l,name="ell"))
	f = FisherAnalysis.from_features(np.array(features),parameters=np.array(parameters),parameter_index=["Om","w","si8","ns"],feature_index=feature_index)

	#Sanity checks
	assert f.parameter_set.shape[0] == 4
	assert f.feature_set.shape[0] == 4

	#Compute the derivatives
	derivatives = f.derivatives.values

	#Sanity check
	assert derivatives.shape == (3,len(l))
	assert f.varied == [0,1,2]

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
	np.savetxt("fitted_parameters.txt",fitted_parameters.values)

	#Sample fisher matrix
	fisher = f.fisher_matrix(simulated_features_covariance=covariance)
	np.savetxt("fisher_constraints.txt",np.sqrt(np.linalg.inv(fisher.values).diagonal()))

	#Confidence ellipse
	fig,ax = plt.subplots()
	e = f.confidence_ellipse(covariance,observed_feature=P,parameters=["Om","si8"],fill=True,edgecolor="black",alpha=0.3)
	ax.add_artist(e)
	ax.set_xlim(0.26-5.0e-3,0.26+5.0e-3)
	ax.set_ylim(0.85-5.0e-3,0.85+5.0e-3)
	ax.legend([e],["contour"],loc="upper left")
	ax.set_xlabel(r"$\Omega_m$")
	ax.set_ylabel(r"$\sigma_8$")
	fig.savefig("fisher_contours.png")

#Test interpolation of power spectrum
def test_interpolation():

	root_path = os.path.join(dataExtern(),"all")
	parameters = list()
	features = list()

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
		ens = Ensemble(np.load(os.path.join(root_path,model._cosmo_id_string,"subfield1","sigma05","power_spectrum.npy")))
		parameters.append(model.squeeze(with_ns=True))
		features.append(ens.mean(0).values)

	#Instantiate the Emulator
	feature_index = EmulatorSeries.make_index(pd.Index(ell,name="ell"))
	analysis = Emulator.from_features(np.array(features),parameters=np.array(parameters),parameter_index=["Om","w","si8","ns"],feature_index=feature_index)

	ens = Ensemble(np.load(os.path.join(root_path,testing_model._cosmo_id_string,"subfield1","sigma05","power_spectrum.npy")))
	testing_Pl = ens.mean(0).values

	#Load in also the observed power spectrum
	ens = Ensemble(np.load(os.path.join(root_path,"observations","subfield1","sigma05","power_spectrum.npy")))
	observed_Pl = ens.mean(0) 

	#Output the analysis stats
	np.savetxt("16_parameter_points.txt",analysis.parameter_set)

	for n in range(len(training_models)):
		plt.plot(ell,ell*(ell+1)*analysis.feature_set[n]/(2*np.pi))

	plt.plot(ell,ell*(ell+1)*observed_Pl/(2*np.pi),linestyle="--",label="Observation")	

	plt.xlabel(r"$l$")
	plt.ylabel(r"$l(l+1)P_l/2\pi$")
	plt.yscale("log")

	plt.legend(loc="upper left")

	plt.savefig("16_power_spectra.png")
	plt.clf()

	#Train the interpolators
	analysis.train(use_parameters=[0,1,2])
	assert hasattr(analysis,"_interpolator")
	assert hasattr(analysis,"_num_bins")

	#Emulator portability test with pickle/unpickle
	analysis.save("analysis.pkl")
	emulator = Emulator.read("analysis.pkl")
	emulator.train(use_parameters=[0,1,2])

	#Predict the power spectrum at the remaining point
	predicted_Pl = emulator.predict(testing_model.squeeze())

	#Plot it against the measured one
	fig,ax = plt.subplots(2,1,figsize=(16,8))

	#Measured
	ax[0].plot(ell,ell*(ell+1)*testing_Pl/(2*np.pi),label="measured")

	#Predicted
	ax[0].plot(ell,ell*(ell+1)*predicted_Pl/(2*np.pi),label="interpolated")
	
	#Fractional difference
	ax[1].plot(ell,(predicted_Pl - testing_Pl)/testing_Pl)

	ax[1].set_xlabel(r"$l$")
	ax[0].set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax[1].set_ylabel(r"$P_l^I-P_l^M/P_l^M$")
	
	ax[0].set_yscale("log")
	ax[0].legend(loc="upper left")

	plt.savefig("power_interpolator_test.png")
	plt.clf()

	#Give it a shot with two points in parameter space to test vectorization
	two_parameter_points = np.array((training_models[0].squeeze(),testing_model.squeeze()))
	two_predicted_Pl = emulator.predict(two_parameter_points).values

	fig,ax = plt.subplots(2,1,figsize=(16,8))

	#Predicted
	ax[0].plot(ell,ell*(ell+1)*two_predicted_Pl[0]/(2*np.pi),color="red",linestyle="--")
	ax[0].plot(ell,ell*(ell+1)*two_predicted_Pl[1]/(2*np.pi),color="green",linestyle="--")

	#Measured
	ax[0].plot(ell,ell*(ell+1)*emulator.feature_set[0]/(2*np.pi),color="red",linestyle="-")
	ax[0].plot(ell,ell*(ell+1)*testing_Pl/(2*np.pi),color="green",linestyle="-")

	#Fractional difference
	ax[1].plot(ell,(two_predicted_Pl[0] - emulator.feature_set[0])/emulator.feature_set[0],color="red")
	ax[1].plot(ell,(two_predicted_Pl[1] - testing_Pl)/testing_Pl,color="green")

	ax[1].set_xlabel(r"$l$")
	ax[0].set_ylabel(r"$l(l+1)P_l/2\pi$")
	ax[1].set_ylabel(r"$P_l^I-P_l^M/P_l^M$")
	
	ax[0].set_yscale("log")

	plt.savefig("power_interpolator_test_2.png")
	plt.clf()

	#Generate a fudge power spectrum covariance matrix
	covariance = np.diag(testing_Pl**2/(0.5 + ell))

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

def test_find():

	emulator = Emulator.read("analysis.pkl")
	parameters_to_find = emulator.parameter_set[7]

	n = emulator.find(parameters_to_find)
	assert len(n)==1
	assert n[0] == 7


def test_reparametrize():

	#Unpickle the emulator
	emulator = Emulator.read("analysis.pkl")

	#Get the old parametrizations
	Om = emulator[("parameters","Om")].values
	w = emulator[("parameters","w")].values
	si8 = emulator[("parameters","si8")].values

	#Reparametrize the parameter space ((Om,w,si8) -> (w,si8 x Om^0.5))
	emulator = emulator.reparametrize(lambda p:pd.Series([p["w"],p["si8"]*p["Om"]**0.5],index=["w","Si8"]))

	#Check that everything worked
	assert (emulator[("parameters","w")]==w).all()
	assert (emulator[("parameters","Si8")]==(si8*(Om**0.5))).all()


#Test various methods of parameter sampling: Fisher Matrix, grid emulator, MCMC chain
def test_sampling(p_value=0.684):

	#Plot setup
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	
	#Unpickle the emulator, data and covariance matrix
	emulator = Emulator.read(os.path.join(dataExtern(),"sample","emulator.pkl"))
	test_data = pd.read_pickle(os.path.join(dataExtern(),"sample","data.pkl"))
	covariance = Ensemble.read(os.path.join(dataExtern(),"sample","covariance.pkl"))

	#Map the likelihood in the OmegaM-sigma8 plane
	p = Ensemble.meshgrid({"Om":np.linspace(0.2,0.5,50),"sigma8":np.linspace(0.6,0.9,50)})
	p["w"] = -1.
	scores = emulator.score(p,test_data,features_covariance=covariance,correct=1000)
	scores["likelihood"] = np.exp(-0.5*scores[emulator.feature_names[0]])

	for n in range(2):
		contour = ContourPlot.from_scores(scores,parameters=["Om","sigma8"],feature_names=["likelihood"],plot_labels=[r"$\Omega_m$",r"$\sigma_8$"],fig=fig,ax=ax[n])
		contour.show()
		contour.getLikelihoodValues([p_value],precision=0.01)
		contour.plotContours(colors=["red"])
		contour.labels()

	#Approximate the emulator linearly around the maximum (Fisher matrix)
	fisher = emulator.approximate_linear(center=(0.26,-1.,0.8))

	#Consider (OmegaM,sigma8) only
	fisher.pop(("parameters","w"))
	fisher = fisher.iloc[[0,1,3]]

	#Fisher confidence ellipse
	for n in range(2):
		ellipse = fisher.confidence_ellipse(covariance,correct=1000,observed_feature=test_data,parameters=["Om","sigma8"],p_value=p_value,fill=False,edgecolor="blue")
		ax[n].add_artist(ellipse)

	#MCMC sampling of (OmegaM,sigma8)
	samples = emulator.sample_posterior(test_data,features_covariance=covariance,correct=1000,pslice={"w":-1})[emulator.feature_names[0]]
	ax[1].scatter(samples["Om"],samples["sigma8"],marker=".",color="yellow")

	#Save the figure
	fig.tight_layout()
	fig.savefig("parameter_sampling.png")








