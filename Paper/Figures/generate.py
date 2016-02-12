#!/usr/bin/env python

import sys,os,argparse
import re

sys.modules["mpi4py"] = None

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u

from matplotlib import rc
import daft

###############################################

from lenstools import dataExtern

from lenstools.simulations import PotentialPlane
from lenstools.simulations.nicaea import Nicaea

from lenstools import ConvergenceMap,ShearMap

from lenstools.statistics.ensemble import Ensemble
from lenstools.statistics.constraints import Emulator
from lenstools.statistics.contours import ContourPlot

#Options
parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",dest="type",default="png",help="format of the figure to save")
parser.add_argument("fig",nargs="*")

###########################################################################################################################################

#Show a lens plane
def density_plane(cmd_args):

	#Set up plot
	fig,ax = plt.subplots()

	#Load the plane and compute the density
	potential = PotentialPlane.load(os.path.join(dataExtern(),"plane.fits"))
	density = potential.density()

	#Show the result
	density.visualize(colorbar=True,cbar_label=r"$\sigma$",fig=fig,ax=ax)

	#Save
	fig.tight_layout()
	fig.savefig("lens_plane_density."+cmd_args.type)

def potential_plane(cmd_args):

	#Set up plot
	fig,ax = plt.subplots()

	#Load the plane and compute the density
	potential = PotentialPlane.load(os.path.join(dataExtern(),"plane.fits"))

	#Show the result
	potential.visualize(colorbar=True,cbar_label=r"$\phi$",fig=fig,ax=ax)

	#Save
	fig.tight_layout()
	fig.savefig("lens_plane_potential."+cmd_args.type)

###########################################################################################################################################

#Parse logs into a pandas DataFrame
def parse_log(fp):

	#Regex to parse lines
	linefmt = re.compile(r"([0-9\-\:\.\s]+)\:lenstools\.stderr\:(INFO|DEBUG)\:(.+)\:[a-zA-Z\s]*([0-9\.]+) Gbyte \(task\)")
	
	#Keep track of this information
	timestamp = list()
	log_level = list()
	step_type = list()
	peak_memory = list()

	#Cycle through the lines
	for line in fp.readlines():
		match = linefmt.search(line)
		if match:
			match_results = match.groups()
			timestamp.append(pd.to_datetime(_fill(match_results[0]),format="%m-%d %H:%M:%S.%f"))
			log_level.append(match_results[1])
			step_type.append(match_results[2])
			peak_memory.append(float(match_results[3]))

	#Construct the DataFrame
	df = pd.DataFrame.from_dict({"timestamp":timestamp,"level":log_level,"step":step_type,"peak_memory(GB)":peak_memory})
	df["delta_timestamp"] = df.timestamp.diff()
	dt_s = df.delta_timestamp.values.astype(np.float) / 1.0e9
	dt_s[0] = 0
	df["delta_timestamp_s"] = dt_s
	df["timestamp_s"] = dt_s.cumsum()

	#Return to user
	return df

#Fill milliseconds
def _fill(s):
	last = s.split('.')[-1]
	nzeros = 3-len(last)
	if nzeros:
		return s.replace('.'+last,'.'+'0'*nzeros+last)
	else:
		return s

#Memory usage figure
def memory_usage(cmd_args):

	#Setup plot
	fig,ax = plt.subplots()

	#Plot memory usage for planes operations
	with open("Logs/planes.err","r") as fp:
		df_planes = parse_log(fp)

	ax.plot(df_planes["timestamp_s"].values,df_planes["peak_memory(GB)"].values,label="Lens Planes",color="black")
	ax.set_ylim(0,2.)

	#Plot a black line after each plane is completed
	planes_completed = df_planes.timestamp_s[df_planes.step.str.contains("Plane")].values
	planes_completed_memory = df_planes["peak_memory(GB)"][df_planes.step.str.contains("Plane")].values
	for n,t in enumerate(planes_completed):
		_bline(ax,t,planes_completed_memory[n],color="black",linestyle="--")

	#Plot memory usage for raytracing operations
	with open("Logs/ray.err","r") as fp:
		df_ray = parse_log(fp)

	ax_top = ax.twiny()
	ax_top.plot(df_ray["timestamp_s"].values,df_ray["peak_memory(GB)"].values,label="Ray--tracing",color="red")
	ax_top.set_xscale("log")
	ax_top.spines["top"].set_color("red")
	ax_top.tick_params(axis="x",colors="red")

	#Plot a red line after each lens crossing
	lens_crossed = df_ray.timestamp_s[df_ray.step.str.contains("Lens")].values
	lens_crossed_memory = df_ray["peak_memory(GB)"][df_ray.step.str.contains("Lens")].values
	for n,t in enumerate(lens_crossed):
		_tline(ax_top,t,lens_crossed_memory[n],color="red",linestyle="--")

	#Labels
	ax.set_xlabel(r"${\rm Runtime(s)}$",fontsize=22)
	ax.set_ylabel(r"${\rm Peak\,\, memory(GB)}$",fontsize=22)

	#Save
	fig.savefig("memory_usage."+cmd_args.type)

#Plot vertical line
def _bline(ax,x,top,**kwargs):
	down,up = ax.get_ylim()
	xp = np.ones(100)*x
	yp = np.linspace(down,top,100)
	ax.plot(xp,yp,**kwargs)

def _tline(ax,x,bottom,**kwargs):
	down,up = ax.get_ylim()
	xp = np.ones(100)*x
	yp = np.linspace(bottom,up,100)
	ax.plot(xp,yp,**kwargs)

###########################################################################################################################################

def flow(cmd_args):

	rc("font", family="serif", size=12)
	rc("text", usetex=False)

	r_color = {"ec" : "red"}
	g_color = {"ec" : "green"}

	#Instantiate PGM
	pgm = daft.PGM([16,7],origin=[0,0])

	#Nodes
	pgm.add_node(daft.Node("parameters","Parameters",1,2.5,aspect=3.,plot_params=r_color))
	pgm.add_node(daft.Node("geometry","geometry",3,2.5,aspect=3.))

	#Seeds
	pgm.add_node(daft.Node("seed1","seed 1",5,3.25,aspect=2.))
	pgm.add_node(daft.Node("seed2","seed 2",5,2.5,aspect=2.))
	pgm.add_node(daft.Node("seedN",r"seed $N$",5,1.75,aspect=2.))
	
	#ICS
	pgm.add_node(daft.Node("ic1","IC 1",6.5,3.25,aspect=2.))
	pgm.add_node(daft.Node("ic2","IC 2",6.5,2.5,aspect=2.))
	pgm.add_node(daft.Node("icN","IC N",6.5,1.75,aspect=2.))

	#Evolution
	pgm.add_node(daft.Node("e1",r"$\rho^1(\mathbf{x},z)$",8.5,3.25,aspect=3.))
	pgm.add_node(daft.Node("e2",r"$\rho^2(\mathbf{x},z)$",8.5,2.5,aspect=3.))
	pgm.add_node(daft.Node("eN",r"$\rho^N(\mathbf{x},z)$",8.5,1.75,aspect=3.))

	#Lens Planes
	pgm.add_node(daft.Node("p1","Lens Planes 1",10.5,3.25,aspect=3.5))
	pgm.add_node(daft.Node("p2","Lens Planes 2",10.5,2.5,aspect=3.5))
	pgm.add_node(daft.Node("pN","Lens Planes N",10.5,1.75,aspect=3.5))

	#Mix planes
	pgm.add_plate(daft.Plate([9.4,1.0,2.1,3.0],label="Mix ICs"))

	#Lensing maps
	pgm.add_node(daft.Node("lens","Lensing maps " + r"$(\kappa,\gamma)$",13.0,2.5,aspect=4.5,plot_params=g_color))

	#Executables
	pgm.add_node(daft.Node("camb","CAMB",3,0.5,aspect=2.,observed=True))
	pgm.add_node(daft.Node("ngenic","NGen-IC",5,0.5,aspect=2.,observed=True))
	pgm.add_node(daft.Node("gadget","Gadget2",6.5,0.5,aspect=2.,observed=True))
	pgm.add_node(daft.Node("planes","lenstools.planes",8.5,0.5,aspect=4.,observed=True))
	pgm.add_node(daft.Node("ray","lenstools.raytracing",10.5,4.5,aspect=5.,observed=True))

	#Edges
	pgm.add_edge("parameters","geometry")
	pgm.add_edge("geometry","seed1")
	pgm.add_edge("geometry","seed2")
	pgm.add_edge("geometry","seedN")
	pgm.add_edge("seed1","ic1")
	pgm.add_edge("seed2","ic2")
	pgm.add_edge("seedN","icN")
	pgm.add_edge("ic1","e1")
	pgm.add_edge("ic2","e2")
	pgm.add_edge("icN","eN")
	pgm.add_edge("e1","p1")
	pgm.add_edge("e2","p2")
	pgm.add_edge("eN","pN")
	pgm.add_edge("p2",'lens')
	pgm.add_edge("camb","geometry")
	pgm.add_edge("ngenic","seedN")
	pgm.add_edge("gadget","icN")
	pgm.add_edge("planes","eN")
	pgm.add_edge("ray","p1")

	#Render and save
	pgm.render()
	pgm.figure.savefig("flow."+cmd_args.type)

def inheritance(cmd_args):

	rc("font", family="serif", size=12)
	rc("text", usetex=False)

	r_color = {"ec" : "red"}
	g_color = {"ec" : "green"}

	#Instantiate PGM
	pgm = daft.PGM([17,7],origin=[0,0])

	#Batch
	pgm.add_node(daft.Node("batch","SimulationBatch",2,2.5,aspect=5.))

	#Model
	pgm.add_node(daft.Node("model","SimulationModel",5,2.5,aspect=5.,plot_params=r_color))

	#Collection
	pgm.add_node(daft.Node("collection","SimulationCollection",8,2.5,aspect=5.))

	#IC
	pgm.add_node(daft.Node("ic","SimulationIC",11,2.5,aspect=5.))

	#Planes
	pgm.add_node(daft.Node("planes","SimulationPlanes",14,2.5,aspect=5.))

	#Maps and catalogs
	pgm.add_node(daft.Node("maps","SimulationMaps",11,3.5,aspect=5.,plot_params=g_color))
	pgm.add_node(daft.Node("catalog","SimulationCatalog",11,1.5,aspect=5.,plot_params=g_color))
	pgm.add_node(daft.Node("subcatalog","SimulationSubCatalog",14,1.5,aspect=5.,plot_params=g_color))

	#Plate
	pgm.add_plate(daft.Plate([3.4,1.0,12,3.0]))

	#Edges
	pgm.add_edge("model","collection")
	pgm.add_edge("collection","ic")
	pgm.add_edge("ic","planes")
	pgm.add_edge("collection","maps")
	pgm.add_edge("collection","catalog")
	pgm.add_edge("catalog","subcatalog")

	#Render and save
	pgm.render()
	pgm.figure.savefig("inheritance."+cmd_args.type)

###########################################################################################################################################

def convergence_visualize(cmd_args):

	#Plot setup
	fig,ax = plt.subplots()
 
	#Load the convergence map and smooth on 1 arcmin
	conv = ConvergenceMap.load(os.path.join(dataExtern(),"conv1.fit"))
	conv.smooth(1.0*u.arcmin,kind="gaussianFFT",inplace=True)

	#Find the peak locations and height
	sigma_peaks = np.linspace(-2.,11.,101)
	height,positions = conv.locatePeaks(sigma_peaks,norm=True)

	#Show the map and the peaks on it (left panel)
	conv.visualize(fig=fig,ax=ax,colorbar=True,cbar_label=r"$\kappa$")
	ax.scatter(*positions[height>2.].to(u.deg).value.T,color="black",marker="x")
	ax.set_xlim(0,conv.side_angle.to(u.deg).value)
	ax.set_ylim(0,conv.side_angle.to(u.deg).value)

	#Save the figure
	fig.tight_layout()
	fig.savefig("convergence_visualize."+cmd_args.type)

def convergence_stats(cmd_args):

	#Plot setup
	fig,ax = plt.subplots()

	#Load the convergence map and smooth on 1 arcmin
	conv = ConvergenceMap.load(os.path.join(dataExtern(),"conv1.fit"))
	conv.smooth(1.0*u.arcmin,kind="gaussianFFT",inplace=True)

	#Find the peak locations and height
	sigma = np.linspace(-2.,11.,101)

	#Show the peak histogram and the PDF
	conv.peakHistogram(sigma,norm=True,fig=fig,ax=ax)

	ax_right = ax.twinx()
	conv.plotPDF(sigma,norm=True,fig=fig,ax=ax_right,color="red")
	
	#All PDF quantities are shown in red
	ax_right.spines["right"].set_color("red")
	ax_right.tick_params(axis="y",colors="red")
	ax_right.yaxis.label.set_color("red")

	#Save
	fig.savefig("convergence_stats."+cmd_args.type)


def eb_modes(cmd_args):

	#Plot setup
	fig,ax = plt.subplots()

	#Load in the shear map, compute E and B modes power spectrum
	shear = ShearMap.load(os.path.join(dataExtern(),"WLshear_z2.00_0001r.fits"))
	l_edges = np.linspace(200.,50000.,50)
	l,ee,bb,eb = shear.decompose(l_edges)

	#Plot the power spectra and prediction from NICAEA
	ax.plot(l,l*(l+1)*ee/(2.*np.pi),label=r"$P^{EE}$",color="red")
	cosmo = Nicaea(Om0=0.26,Ode0=0.74,w0=-1,sigma8=0.8)
	ax.plot(l,l*(l+1)*cosmo.convergencePowerSpectrum(l,z=2.0)/(2.*np.pi),label=r"$P^{\kappa\kappa}{\rm (NICAEA)}$",linestyle="--",color="red")
	ax.plot(l,l*(l+1)*bb/(2.*np.pi),label=r"$P^{BB}$",color="blue")

	#Labels
	ax.set_xscale("log")
	ax.set_yscale("log")
	ax.set_xlabel(r"$\ell$",fontsize=22)
	ax.set_ylabel(r"$\ell(\ell+1)P_\ell/2\pi$",fontsize=22)
	ax.legend(loc="upper left")

	#Save
	fig.savefig("eb_modes."+cmd_args.type)


###########################################################################################################################################

def parameter_sampling(cmd_args,p_value=0.684):

	#Plot setup
	fig,ax = plt.subplots(figsize=(8,8))
	
	#Unpickle the emulator, data and covariance matrix
	emulator = Emulator.read(os.path.join(dataExtern(),"sample","emulator.pkl"))
	test_data = pd.read_pickle(os.path.join(dataExtern(),"sample","data.pkl"))
	covariance = Ensemble.read(os.path.join(dataExtern(),"sample","covariance.pkl"))

	#Map the likelihood in the OmegaM-sigma8 plane
	p = Ensemble.meshgrid({"Om":np.linspace(0.2,0.5,50),"sigma8":np.linspace(0.6,0.9,50)})
	p["w"] = -1.
	scores = emulator.score(p,test_data,features_covariance=covariance,correct=1000)
	scores["likelihood"] = np.exp(-0.5*scores[emulator.feature_names[0]])

	contour = ContourPlot.from_scores(scores,parameters=["Om","sigma8"],feature_names=["likelihood"],plot_labels=[r"$\Omega_m$",r"$\sigma_8$"],fig=fig,ax=ax)
	#contour.show(cmap=plt.get_cmap("binary"))
	contour.getLikelihoodValues([p_value],precision=0.01)
	contour.plotContours(colors=["red"])
	contour.labels()

	#Approximate the emulator linearly around the maximum (Fisher matrix)
	fisher = emulator.approximate_linear(center=(0.26,-1.,0.8))

	#Consider (OmegaM,sigma8) only
	fisher.pop(("parameters","w"))
	fisher = fisher.iloc[[0,1,3]]

	#Fisher confidence ellipse
	ellipse = fisher.confidence_ellipse(covariance,correct=1000,observed_feature=test_data,parameters=["Om","sigma8"],p_value=p_value,fill=False,edgecolor="blue")
	ax.add_artist(ellipse)

	#MCMC sampling of (OmegaM,sigma8)
	samples = emulator.sample_posterior(test_data,features_covariance=covariance,correct=1000,pslice={"w":-1})[emulator.feature_names[0]]
	ax.scatter(samples["Om"],samples["sigma8"],marker=".",color="black",s=1)
	ax.set_xlim(0.2,0.5)
	ax.set_ylim(0.6,0.9)

	#Save the figure
	fig.tight_layout()
	fig.savefig("parameter_sampling."+cmd_args.type)



###########################################################################################################################################

#Method dictionary
method = dict()
method["1"] = density_plane
method["1b"] = potential_plane
method["2"] = flow
method["3"] = inheritance
method["4"] = memory_usage
method["5"] = convergence_visualize
method["5b"] = convergence_stats
method["6"] = eb_modes
method["7"] = parameter_sampling

#Main
def main():
	cmd_args = parser.parse_args()
	for fig in cmd_args.fig:
		method[fig](cmd_args)

if __name__=="__main__":
	main()