import sys,os

from .. import dataExtern
from ..statistics.database import Database,ScoreDatabase

import numpy as np
import matplotlib.pyplot as plt

##############################
#Plot styles for each feature#
##############################

#Markers
markers = {
"power_logb_large" : "x",
"power_logb_small" : "s",
"power_logb_all" : "o",
"power_large" : "+",
"power_small" : "x",
"power_large+small" : "o",
"power_all" : "o",
"peaks_low" : "+",
"peaks_intermediate" : "*",
"peaks_high" : "d",
"peaks_low+intermediate" : "x",
"peaks_intermediate+high" : "s",
"peaks_all" : "s",
} 

#Colors
colors = {
"power_logb_large" : "black",
"power_logb_small" : "black",
"power_logb_all" : "red",
"power_large" : "red",
"power_small" : "red",
"power_large+small" : "green",
"power_all" : "blue",
"peaks_low" : "red",
"peaks_intermediate" : "red",
"peaks_high" : "red",
"peaks_low+intermediate" : "green",
"peaks_intermediate+high" : "green",
"peaks_all" : "magenta",
} 

#Labels
labels = {
"power_logb_large" : r"$\ell\in[100,800],N_b=8(\mathrm{log})$",
"power_logb_small" : r"$\ell\in[1000,6000],N_b=7(\mathrm{log})$",
"power_logb_all" : r"$\ell\in[100,6000],N_b=15(\mathrm{log})$",
"power_logb_lowest_ell" : r"$\ell\in[100,250],N_b=4(\mathrm{log})$",
"power_large" : r"$\ell\in[100,2000],N_b=15(\mathrm{lin})$",
"power_small" : r"$\ell\in[2500,4500],N_b=15(\mathrm{lin})$",
"power_large+small" : r"$\ell\in[100,4500],N_b=30(\mathrm{lin})$",
"power_all" : r"$\ell\in[2500,6000],N_b=39(\mathrm{lin})$",
"peaks_low" : r"$\kappa_0\in[-0.06,0.09],N_b=15$",
"peaks_intermediate" : r"$\kappa_0\in[0.1,0.27],N_b=15$",
"peaks_high" : r"$\kappa_0\in[0.28,0.45],N_b=15$",
"peaks_low+intermediate" : r"$\kappa_0\in[-0.06,0.27],N_b=30$",
"peaks_intermediate+high" : r"$\kappa_0\in[0.1,0.45],N_b=30$",
"peaks_all" : r"$\kappa_0\in[-0.06,0.45],N_b=45$",
"peaks_highest_kappa" : r"$\kappa_0\in[0.44,0.48],N_b=4$",
"peaks_highest_kappa_s1" : r"$\kappa_0(\theta_G=1^\prime)>0.15,N_b=20$",
} 

#Plot order
order = {
"power_logb_large" : 8,
"power_logb_small" : 7,
"power_logb_all" : 15,
"power_large" : 15,
"power_small" : 15,
"power_large+small" : 30,
"power_all" : 39,
"peaks_low" : 15,
"peaks_intermediate" : 15,
"peaks_high" : 15,
"peaks_low+intermediate" : 30,
"peaks_intermediate+high" : 30,
"peaks_all" : 45,
} 

#Curving effect of the variance versus 1/Nr fir different Nb
def test_curving_nb():

	#Options
	db_filename = os.path.join(dataExtern(),"database","variance_scaling_nb_expected.sqlite")
	parameter = "w"
	nsim = 200
	xlim = (0,1./65)
	ylim = (1,2.5)
	nr_top = [1000,500,300,200,150,100,90,70]

	#Plot panel
	fig,ax = plt.subplots() 
	
	#################################################################################################################

	#Load the database and fit for the effective dimensionality of each feature space
	with Database(db_filename) as db:

		features  = db.tables
		features.sort(key=order.get)
		
		for f in features:

			#Read the table corresponding to each feature
			v = db.read_table(f).query("nsim=={0}".format(nsim))
			v["1/nreal"] = v.eval("1.0/nreal")
			v = v.sort_values("1/nreal")

			#Find the variance in the limit of large Nr
			s0 = v[parameter].values[0]

			#Nb,Np
			Nb = v["bins"].mean()
			Np = 3

			#Plot the variance versus 1/nreal
			ax.scatter(v["1/nreal"],v[parameter]/s0,color=colors[f],marker=markers[f],label=labels[f],s=10+(100-10)*(Nb-1)/(200-1))

			#Plot the theory predictions
			x = 1./np.linspace(1000,65,100)
			ax.plot(x,1+x*(Nb-Np),linestyle="--",color=colors[f])
			ax.plot(x,1+(Nb-Np)*x+(Nb-Np)*(Nb-Np+2)*(x**2),linestyle="-",color=colors[f])
			ax.plot(x,(1-2*x)/(1-(Nb-Np+2)*x),linestyle="-",linewidth=3,color=colors[f])


	#Axis bounds
	ax.set_xlim(*xlim)
	ax.set_ylim(*ylim)

	#Axis labels and legends
	ax.set_xlabel(r"$1/N_r$",fontsize=22)
	ax.set_ylabel(r"$\langle\hat{\sigma}^2_w\rangle/\sigma^2_{w,\infty}$",fontsize=22)
	ax.legend(loc="upper left",prop={"size":13})

	#Mirror x axis to show Nr on top
	ax1 = ax.twiny()
	ax1.set_xlim(*xlim)
	ax1.set_xticks([1./n for n in nr_top])
	ax1.set_xticklabels([str(n) for n in nr_top])
	ax1.set_xlabel(r"$N_r$",fontsize=22)

	#Save the figure
	fig.savefig("database_nb.png")

#Test the sigma8 likelihood (at fixed Om,w)
def test_sigma8():

	#Set up plot
	fig,ax = plt.subplots()

	db_name = os.path.join(dataExtern(),"database","scores_power.sqlite")
	with ScoreDatabase(db_name) as db:
		likelihood = db.pull_features(["coarse_150sim"],table_name="coarse",score_type="likelihood")

	#Plot
	likelihood_slice = likelihood.query("Om==0.2 and w==-1.0")
	ax.plot(likelihood_slice["sigma8"],likelihood_slice["coarse_150sim"])
	ax.set_xlabel(r"$\sigma_8$",fontsize=22)
	ax.set_ylabel(r"$\mathcal{L}(\sigma_8,\Omega_m=0.2,w=-1)$",fontsize=22)

	#Save
	fig.savefig("si8_likelihood.png")



