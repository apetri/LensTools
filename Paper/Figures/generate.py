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
from lenstools.image.convergence import ConvergenceMap
from lenstools.statistics.ensemble import Ensemble
from lenstools.statistics.constraints import Emulator
from lenstools.statistics.contours import ContourPlot

#Options
parser = argparse.ArgumentParser()
parser.add_argument("-t","--type",dest="type",default="png",help="format of the figure to save")
parser.add_argument("fig",nargs="*")

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
	pgm.add_plate(daft.Plate([9.4,1.0,2.1,3.0],label="Mix realizations"))

	#Lensing maps
	pgm.add_node(daft.Node("lens","Lensing maps",13.0,2.5,aspect=3.5,plot_params=g_color))

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

###########################################################################################################################################

def convergence_stats(cmd_args):

	#Plot setup
	fig,ax = plt.subplots()
 
	#Load the convergence map and smooth on 1 arcmin
	conv = ConvergenceMap.load(os.path.join(dataExtern(),"conv1.fit"))
	conv.smooth(1.0*u.arcmin,kind="gaussianFFT",inplace=True)

	#Find the peak locations and height
	sigma_peaks = np.linspace(-2.,15.,101)
	height,positions = conv.locatePeaks(sigma_peaks,norm=True)

	#Show the map and the peaks on it (left panel)
	conv.visualize(fig=fig,ax=ax,colorbar=True,cbar_label=r"$\kappa$")
	ax.scatter(*positions[height>2.].to(u.deg).value.T,color="black",marker="x")
	ax.set_xlim(0,conv.side_angle.to(u.deg).value)
	ax.set_ylim(0,conv.side_angle.to(u.deg).value)

	#Save the figure
	fig.tight_layout()
	fig.savefig("convergence_stats."+cmd_args.type)



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
method["1"] = flow
method["2"] = convergence_stats
method["3"] = parameter_sampling

#Main
def main():
	cmd_args = parser.parse_args()
	for fig in cmd_args.fig:
		method[fig](cmd_args)

if __name__=="__main__":
	main()