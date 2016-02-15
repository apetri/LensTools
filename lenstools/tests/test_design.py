import os

from ..simulations import Design

from .. import dataExtern

import numpy as np
import matplotlib.pyplot as plt


#Test the visualization of a design
def test_visualize():

	#This fails if GSL is not installed
	try:
		design = Design.from_specs(npoints=50,parameters=[("Om",r"$\Omega_m$",0.1,0.9),("w",r"$w$",-2.0,-1.0),("si8",r"$\sigma_8$",0.01,1.6)])
	except ImportError:
		return

	print(design)

	#The cost function should be the diagonal one
	np.testing.assert_approx_equal(design.diagonalCost(Lambda=1.0),design.cost(p=2.0,Lambda=1.0))

	#Visualize the 3d diagonal configuration
	design.visualize(color="blue")
	design.set_title("Cost={0:.2f}".format(design.diagonalCost(Lambda=1.0)))
	design.savefig("design_diagonal_3d.png")

	#Visualize the 2d (Om,si8) projection
	fig,ax = plt.subplots()
	design.visualize(fig=fig,ax=ax,parameters=["Om","si8"],color="red",marker=".")
	design.savefig("design_diagonal_2d.png")

	#Now perform the sampling of the parameter space
	deltaPerc = design.sample(Lambda=1.0,p=2.0,seed=1,maxIterations=100000)

	#Visualize the 3d configuration
	design.visualize(color="blue")
	design.set_title("Cost={0:.2f}".format(design.cost(Lambda=1.0,p=2.0)))
	design.savefig("design_3d.png")

	#Visualize the 2d (Om,si8) projection
	fig,ax = plt.subplots()
	design.visualize(fig=fig,ax=ax,parameters=["Om","si8"],color="red",marker=".")
	design.savefig("design_2d.png")

	fig,ax = plt.subplots()
	design.plot(ax=ax,kind="scatter",x="Om",y="si8")
	fig.savefig("design_2d_pandas.png")

	#Visualize the changes in the cost function
	fig,ax = plt.subplots()
	ax.plot(design.cost_values)
	ax.set_xlabel(r"$N$")
	ax.set_ylabel("cost")
	ax.set_title("Last change={0:.1e}%".format(deltaPerc*100))
	ax.set_xscale("log")
	fig.savefig("cost.png")
