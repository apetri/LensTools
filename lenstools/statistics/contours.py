"""

.. module:: contours 
	:platform: Unix
	:synopsis: This module implements a confidence contour plotting engine


.. moduleauthor:: Andrea Petri <apetri@phys.columbia.edu>


"""

from __future__ import print_function,division,with_statement

import os
import logging

from .constraints import FisherAnalysis

import numpy as np
from scipy import stats
from scipy import integrate

try: 
	import matplotlib
	import matplotlib.pyplot as plt
	from matplotlib import rc
except ImportError:
	matplotlib = None

import pandas as pd

#############################################################
############Find confidence levels in 1D likelihood##########
#############################################################

def _1d_level_values(p,l,level=0.684,quantity=2):

	"""
	Find the parameter extremes that correspons to the likelihood N--sigma level

	"""

	#Find the maximum of the likelihood
	maximum = np.where(l==l.max())[0][0]
	parmax = p[maximum]

	all_levels = np.zeros_like(l)
	
	for n in range(l.shape[0]):
		all_levels[n] = l[l>=l[n]].sum() / l.sum()

	#Find the closest level
	closest = np.argmin(np.abs(all_levels - level))

	#Find the n corresponding parameter values
	ranks = stats.rankdata(np.abs(l-l[closest])).astype(np.int) - 1

	par = list()
	for n in range(quantity):
		par.append(p[np.where(ranks==n)[0][0]])

	#Sort from left to right
	par.sort()

	return par


#############################################################
###########Find confidence levels in N-dim likelihood########
#############################################################

def _nd_level_value(likelihood,level,low,high,precision=0.01):

	middle = (low+high)/2
	current_integral = likelihood[likelihood>middle].sum()

	if np.abs((current_integral-level)/level)<precision:
		return middle
	
	#Proceed with bisection method
	if current_integral>level:
		return _nd_level_value(likelihood,level,middle,high,precision=precision)
	else:
		return _nd_level_value(likelihood,level,low,middle,precision=precision)

#############################################################
##################ContourPlot class##########################
#############################################################

class ContourPlot(object):

	"""
	A class handler for contour plots

	"""


	##############################################################
	#######Build ContourPlot instance out of score ensemble#######
	##############################################################

	@classmethod 
	def from_scores(cls,score_ensemble,parameters,feature_names=None,plot_labels=None,fig=None,ax=None,figsize=(8,8)):

		"""
		Build a ContourPlot instance out of an ensemble of scores

		:param score_ensemble: ensemble of the scores
		:type score_ensemble: :py:class:`Ensemble`

		:param parameters: columns names that contain the parameters
		:type parameters: list.

		:param feature_names: name of the features to generate the contour plot for
		:type feature_names: list.

		:param plot_labels: plot labels for the parameters
		:type plot_labels: list.

		:param figsize: size of the plot
		:type figsize: tuple.

		:returns: list of ContourPlot

		"""

		contour_plots = list()

		npar = len(parameters)
		nside = int(np.ceil(len(score_ensemble)**(1/npar)))

		#Type check
		assert isinstance(score_ensemble,pd.DataFrame)

		#Check that the ensemble is the right size
		if nside**npar!=len(score_ensemble):
			raise ValueError("Parameters cannot be cast in a meshgrid!")

		#Check that there is regular spacing between parameters
		parameter_values_dict = dict()
		for p in parameters:
			parameter_values = np.sort(score_ensemble[p].drop_duplicates().values)
			
			if len(parameter_values)!=nside:
				raise ValueError("Parameters cannot be cast in a meshgrid!")

			if not(np.isclose(parameter_values[1]-parameter_values[0],(parameter_values[1:]-parameter_values[:-1]).mean())):
				raise ValueError("Parameters cannot be cast in a meshgrid!")

			parameter_values_dict[p] = parameter_values

		#Sort the score_ensemble so that the parameters ordered like a meshgrid
		score_ensemble = score_ensemble.sort_values(by=parameters)

		#Check the feature names used to plot the contours (one contour plot for each feature name)
		if type(feature_names)==str:
			feature_names = [feature_names]

		if feature_names is not None:
			feature_labels = feature_names
		else: 
			feature_labels = filter(lambda l:l not in parameters,score_ensemble.columns)

		#One contour plot for each label
		if (matplotlib is not None) and ((fig is None) or (ax is None)):
			fig,ax = plt.subplots(figsize=figsize)
		
		for feature in feature_labels:
			
			contour = cls(fig,ax)
			contour.parameter_axes = dict()
			contour.parameter_labels = dict()
			
			#Units and labels
			for n,p in enumerate(parameters):
				contour.parameter_axes[p] = n
				contour.min[p] = parameter_values_dict[p].min()
				contour.max[p] = parameter_values_dict[p].max()
				contour.npoints[p] = len(parameter_values_dict[p])
				contour.unit[p] = parameter_values_dict[p][1] - parameter_values_dict[p][0]

				if plot_labels is not None:
					contour.parameter_labels[p] = plot_labels[n]
				else:
					contour.parameter_labels[p] = p

			#Finally load in the probability meshgrid
			contour.likelihood = score_ensemble[feature].values.reshape(*(nside,)*len(parameters))
			contour.likelihood /= contour.likelihood.sum()

			#If the likelihood is two dimensional then we can update the plot extent information
			if len(parameters)==2:
				contour.extent = (contour.min[parameters[0]],contour.max[parameters[0]],contour.min[parameters[1]],contour.max[parameters[1]])

			#Append to the list
			contour_plots.append(contour)

		#Return
		if len(contour_plots)>1:
			return contour_plots
		else:
			return contour_plots[0]


	##############################################################

	def __init__(self,fig=None,ax=None):

		try:
			
			if (fig is None) or (ax is None):
				self.fig,self.ax = plt.subplots()
				self.ax.proxy = list()
			else:
				self.fig = fig
				self.ax = ax

				if not hasattr(self.ax,"proxy"):
					self.ax.proxy = list()

		except:
			pass
		
		self.min = dict()
		self.max = dict()
		self.npoints = dict()
		self.unit = dict()

	def savefig(self,figname):

		"""
		Save the plot to file

		"""

		self.fig.savefig(figname)


	def close(self):

		"""
		Closes the figure

		"""

		plt.close(self.fig)

	def window(self):

		plt.ion()
		plt.show()

	def getUnitsFromOptions(self,options):
		
		"""
		Parse options file to get physical units of axes

		"""

		assert hasattr(self,"parameter_axes"),"You have to load in the likelihood first!"
		parameters = self.parameter_axes.keys()

		for parameter in parameters:
			
			self.min[parameter],self.max[parameter],self.npoints[parameter] = options.getfloat(parameter,"min"),options.getfloat(parameter,"max"),options.getint(parameter,"num_points")
			assert self.npoints[parameter] == self.likelihood.shape[self.parameter_axes[parameter]]
			self.unit[parameter] = (self.max[parameter] - self.min[parameter]) / (self.npoints[parameter] - 1)

	def setUnits(self,parameter,parameter_min,parameter_max,parameter_unit):

		"""
		Set manually the physical units for each of the likelihood axes

		"""
		assert hasattr(self,"parameter_axes"),"You have to load in the likelihood first!"
		assert parameter in self.parameter_axes.keys(),"You are trying to set units for a parameter that doesn't exist!"

		self.min[parameter] = parameter_min
		self.max[parameter] = parameter_max
		self.unit[parameter] = parameter_unit

		print("Units set for {0}; min={1:.3f} max={2:.3f} unit={3:.3f}".format(parameter,parameter_min,parameter_max,parameter_unit))

	def value(self,*coordinates):

		"""
		Compute the (un-normalized) likelihood value at the specified point in parameter space

		"""

		assert len(coordinates) == self.likelihood.ndim,"You must specify a coordinate (and only one) for each axis"

		#Compute the physical values of the pixels
		pix = np.zeros(len(coordinates))
		for parameter in self.parameter_axes.keys():

			assert parameter in self.unit.keys() and parameter in self.min.keys()
			axis = self.parameter_axes[parameter]
			pix[axis] = int((coordinates[axis] - self.min[parameter])/(self.unit[parameter]))

		#Return the found likelihood value
		try:
			return self.likelihood[tuple(pix)]
		except IndexError:
			print("Out of bounds!")
			return None


	def getLikelihood(self,likelihood_filename,parameter_axes={"Omega_m":0,"w":1,"sigma8":2},parameter_labels={"Omega_m":r"$\Omega_m$","w":r"$w$","sigma8":r"$\sigma_8$"}):
		
		"""
		Load the likelihood function from a numpy file

		"""

		self.parameter_axes = parameter_axes
		self.parameter_labels = parameter_labels

		if type(likelihood_filename)==str:
			
			self.likelihood = np.load(likelihood_filename)
			#Construct title label
			self.title_label = os.path.split(likelihood_filename)[1].lstrip("likelihood_").rstrip(".npy")
		
		elif type(likelihood_filename)==np.ndarray:
			
			self.likelihood = likelihood_filename
			#Construct title label
			self.title_label = "Default"

		assert len(self.parameter_axes.keys()) == self.likelihood.ndim,"The number of parameters should be the same as the number of dimensions of the likelihood!"

		#Normalize
		self.likelihood /= self.likelihood.sum()

	def getMaximum(self,which="full"):

		"""
		Find the point in parameter space on which the likelihood is maximum

		"""
		max_parameters = dict()

		if which=="full":
			
			max_loc = np.where(self.likelihood==self.likelihood.max())
			for parameter in self.parameter_axes.keys():
				max_parameters[parameter] = max_loc[self.parameter_axes[parameter]][0] * self.unit[parameter] + self.min[parameter]
		
		elif which=="reduced":
			
			max_loc = np.where(self.reduced_likelihood==self.reduced_likelihood.max())
			for n,parameter in enumerate(self.remaining_parameters):
				max_parameters[parameter] = max_loc[n][0] * self.unit[parameter] + self.min[parameter]
		
		else:
			raise ValueError("which must be either 'full' or 'reduced'")

		return max_parameters


	def expectationValue(self,function,**kwargs):

		"""
		Computes the expectation value of a function of the parameters over the current parameter likelihood

		"""

		assert hasattr(self,"likelihood"),"You have to load in the likelihood first!"

		#Parameters
		parameters = self.parameter_axes.keys()
		parameters.sort(key=self.parameter_axes.__getitem__)

		#Initialize the parameter mesh
		mesh_axes = [ np.linspace(self.min[par],self.max[par],self.npoints[par]) for par in parameters ]
		parameter_mesh = np.meshgrid(*tuple(mesh_axes),indexing="ij")

		#Compute the expectation value
		expectation = (function(parameter_mesh,**kwargs)*self.likelihood).sum() / self.likelihood.sum()  

		#Return
		return expectation


	def variance(self,function,**kwargs):

		"""
		Computes the variance of a function of the parameters over the current parameter likelihood

		"""

		expectation = self.expectationValue(function,**kwargs)

		#Parameters
		parameters = self.parameter_axes.keys()
		parameters.sort(key=self.parameter_axes.__getitem__)

		#Initialize the parameter mesh
		mesh_axes = [ np.linspace(self.min[par],self.max[par],self.npoints[par]) for par in parameters ]
		parameter_mesh = np.meshgrid(*tuple(mesh_axes),indexing="ij")

		#Compute the variance
		variance = (self.likelihood*(function(parameter_mesh,**kwargs) - expectation)**2).sum() / self.likelihood.sum()

		#Return 
		return variance


	def marginalize(self,parameter_name="w"):

		"""
		Marginalize the likelihood over the indicated parameters

		"""

		#Parse all the parameters to marginalize over
		marginalize_parameters = parameter_name.split(",")

		assert hasattr(self,"likelihood"),"You have to load in the likelihood first!"

		for par in marginalize_parameters:
			assert par in self.parameter_axes.keys(),"You are trying to marginalize over a parameter {0}, that does not exist!".format(par)

		marginalize_indices = [ self.parameter_axes[par] for par in marginalize_parameters ]
		self.reduced_likelihood = self.likelihood.sum(tuple(marginalize_indices))

		#Normalize
		self.reduced_likelihood /= self.reduced_likelihood.sum()

		#Find the remaining parameters
		self.remaining_parameters = self.parameter_axes.keys()

		for par in marginalize_parameters:
			self.remaining_parameters.pop(self.remaining_parameters.index(par))

		#Sort the remaining parameter names so that the corresponding axes are in increasing order
		self.remaining_parameters.sort(key=self.parameter_axes.get)
		
		if len(self.remaining_parameters)==2:
			
			self.extent = (self.min[self.remaining_parameters[0]],self.max[self.remaining_parameters[0]],self.min[self.remaining_parameters[1]],self.max[self.remaining_parameters[1]])
			self.ax.set_xlim(self.extent[0],self.extent[1])
			self.ax.set_ylim(self.extent[2],self.extent[3])


	def marginal(self,parameter_name="w",levels=None):

		"""
		Marginalize the likelihood over all parameters but one

		"""

		assert hasattr(self,"likelihood"),"You have to load in the likelihood first!"
		assert parameter_name in self.parameter_axes.keys(),"You are trying to compute a marginal likelihood of a parameter that does not exist!"

		remaining_parameters = self.parameter_axes.keys()
		remaining_parameters.pop(remaining_parameters.index(parameter_name))
		remaining_parameter_axes = [ self.parameter_axes[par] for par in remaining_parameters ]

		#Marginalize the likelihood
		parameter_range = np.linspace(self.min[parameter_name],self.max[parameter_name],self.npoints[parameter_name])
		marginal_likelihood = self.likelihood.sum(axis=tuple(remaining_parameter_axes))

		#Compute the normalization
		normalization = integrate.simps(marginal_likelihood,x=parameter_range)
		marginal_likelihood /= normalization

		#Compute the maximum
		par_max = parameter_range[np.where(marginal_likelihood==marginal_likelihood.max())[0][0]]

		#Compute also the contour extremes if levels 
		if levels is not None:

			par_extremes = list()
			for level in levels:

				pL = _1d_level_values(parameter_range,marginal_likelihood,level=level,quantity=3)
				par_extremes.append((pL[0],pL[-1]))

			#Return the normalized single parameter likelihood, along with the contour extremes
			return parameter_range,marginal_likelihood,par_max,par_extremes

		else:
			
			#Return the normalized single parameter likelihood
			return parameter_range,marginal_likelihood,par_max


	def slice(self,parameter_name="w",parameter_value=-1.0):

		"""
		Slice the likelihood cube by fixing one of the parameters

		"""

		assert hasattr(self,"likelihood"),"You have to load in the likelihood first!"
		assert parameter_name in self.parameter_axes.keys(),"You are trying to get a slice with a parameter that does not exist!"

			
		#Select the slice
		slice_axis = self.parameter_axes[parameter_name]
		slice_index = int((parameter_value - self.min[parameter_name]) / self.unit[parameter_name])
		assert slice_index<self.npoints[parameter_name],"Out of bounds!"

		#Get the slice
		self.reduced_likelihood = np.split(self.likelihood,self.npoints[parameter_name],axis=slice_axis)[slice_index].squeeze()
			
		#Normalize
		self.reduced_likelihood /= self.reduced_likelihood.sum()

		#Find the remaining parameters
		self.remaining_parameters = self.parameter_axes.keys()
		self.remaining_parameters.pop(self.remaining_parameters.index(parameter_name))
		#Sort the remaining parameter names so that the corresponding axes are in increasing order
		self.remaining_parameters.sort(key=self.parameter_axes.get)
		
		self.extent = (self.min[self.remaining_parameters[0]],self.max[self.remaining_parameters[0]],self.min[self.remaining_parameters[1]],self.max[self.remaining_parameters[1]])
		self.ax.set_xlim(self.extent[0],self.extent[1])
		self.ax.set_ylim(self.extent[2],self.extent[3])


	def show(self,cmap=None,interpolation="nearest"):

		"""
		Show the 2D marginalized likelihood

		"""

		if cmap is None:
			cmap = plt.get_cmap("binary_r")

		if hasattr(self,"reduced_likelihood"):
			reduced_likelihood = self.reduced_likelihood
		else:
			reduced_likelihood = self.likelihood

		assert reduced_likelihood.ndim == 2,"Can show only 2 dimensional likelihoods in the figure!!"
		
		self.likelihood_image = self.ax.imshow(reduced_likelihood.transpose(),origin="lower",cmap=cmap,extent=self.extent,interpolation=interpolation,aspect="auto")
		self.colorbar = plt.colorbar(self.likelihood_image,ax=self.ax)
		

	def labels(self,contour_label=None,fontsize=22,**kwargs):

		"""
		Put the labels on the plot

		"""

		if not hasattr(self,"remaining_parameters"):
			self.remaining_parameters = list(self.parameter_axes.keys())
			self.remaining_parameters.sort(key=self.parameter_axes.__getitem__)

		self.ax.set_xlabel(self.parameter_labels[self.remaining_parameters[0]],fontsize=fontsize)
		self.ax.set_ylabel(self.parameter_labels[self.remaining_parameters[1]],fontsize=fontsize)

		if hasattr(self,"title_label"):
			self.ax.set_title(self.title_label,fontsize=fontsize)

		if contour_label is not None:
			self.ax.legend(self.ax.proxy,contour_label,**kwargs)

	def point(self,coordinate_x,coordinate_y,color="green",marker="o"):

		"""
		Draws a point in parameter space at the specified physical coordinates

		"""

		if not hasattr(self,"remaining_parameters"):
			self.remaining_parameters = self.parameter_axes.keys()
			self.remaining_parameters.sort(key=self.parameter_axes.__getitem__)

		#First translate the physical coordinates into pixels, to obtain the likelihood value
		px = int((coordinate_x - self.min[self.remaining_parameters[0]]) / self.unit[self.remaining_parameters[0]])
		py = int((coordinate_y - self.min[self.remaining_parameters[1]]) / self.unit[self.remaining_parameters[1]])

		#Draw the point
		self.ax.plot(coordinate_x,coordinate_y,color=color,marker=marker)

		#Return the likelihood value at the specified point
		if hasattr(self,"reduced_likelihood"):
			return self.reduced_likelihood[px,py]
		else:
			return self.likelihood[px,py]


	#################################################################################################
	###############Find the likelihood values that correspond to the confidence contours#############
	#################################################################################################

	def getLikelihoodValues(self,levels,precision=0.001):

		"""
		Find the likelihood values that correspond to the selected p_values
		"""

		if hasattr(self,"reduced_likelihood"):
			likelihood = self.reduced_likelihood
		else:
			likelihood = self.likelihood

		self.original_p_values = levels

		#Check sanity of input, likelihood must be normalized
		np.testing.assert_approx_equal(likelihood.sum(),1.0)

		#Initialize list of likelihood values
		values = list()
		p_values = list()

		#Loop through levels to find corresponding likelihood values
		for level in levels:

			#Call the recursive bisection method
			value = _nd_level_value(likelihood,level,likelihood.min(),likelihood.max(),precision=precision)
			confidence_integral = likelihood[likelihood>value].sum()

			#Append the found likelihood value to the output
			values.append(value)
			p_values.append(confidence_integral)

		#Return
		self.computed_p_values = p_values
		self.likelihood_values = values
		
		return values

	################################################################################
	###############Compute the area enclosed by the confidence contours#############
	################################################################################

	def confidenceArea(self):

		"""
		Computes the area enclosed by the parameter confidence contours

		:returns: area enclosed by the confidence contours
		:rtype: dict.

		"""

		if hasattr(self,"reduced_likelihood"):
			likelihood = self.reduced_likelihood
		else:
			likelihood = self.likelihood

		#Find the maximum of the likelihood and compute the hessian matrix around it
		assert likelihood.ndim==2,"The likelihood must be 2-dimensional"

		#Find the names of the relevant parameters
		if hasattr(self,"remaining_parameters"):
			parameters = self.remaining_parameters
		else:
			parameters = self.parameter_axes.keys()

		assert len(parameters)==2,"There should be only 2 parameters!"

		parameters.sort(key=self.parameter_axes.get)

		#This is the area of a pixel in parameter space
		pixel_area = self.unit[parameters[0]] * self.unit[parameters[1]]

		#Check that the values of the likelihood levels have been computed already
		if not hasattr(self,"computed_p_values"):
			raise AttributeError("No likelihood values computed, please run getLikelihoodValues() first to compute the desired likelihood values!!")

		#Compute the enclosed area for each contour
		confidence_area = dict()
		for n,p_value in enumerate(self.computed_p_values):
			confidence_area[p_value] = (likelihood>self.likelihood_values[n]).sum() * pixel_area

		#Return to user
		return confidence_area


	######################################################################################
	##############Plot the Gaussian approximation to the confidence contours##############
	######################################################################################

	def plotEllipse(self,colors=["red"],levels=[0.683],**kwargs):

		"""
		Plot the Gaussian approximation to the confidence contours

		:param colors: colors of the confidence ellipses
		:type colors: list.

		:param levels: p_values of the confidence ellipses
		:type levels: list.

		:param kwargs: additional keyword arguments are passed to the Ellipse patch method
		:type kwargs: dict.

		"""

		if hasattr(self,"reduced_likelihood"):
			likelihood = np.log(self.reduced_likelihood)
		else:
			likelihood = np.log(self.likelihood)

		#Find the maximum of the likelihood and compute the hessian matrix around it
		assert likelihood.ndim==2,"The likelihood must be 2-dimensional"

		#Find the names of the relevant parameters
		if hasattr(self,"remaining_parameters"):
			parameters = self.remaining_parameters
		else:
			parameters = self.parameter_axes.keys()

		assert len(parameters)==2,"There should be only 2 parameters!"

		parameters.sort(key=self.parameter_axes.get)

		#Find the maximum and the center of the ellipse
		imax,jmax = [ l[0] for l in np.where(likelihood==likelihood.max()) ]
		center = tuple([ self.min[p] + self.unit[p]*(imax,jmax)[n] for n,p in enumerate(parameters) ])

		#Compute the hessian
		hessian = np.zeros((2,2))
		hessian[0,0] = (likelihood[imax+2,jmax] + likelihood[imax-2,jmax] - 2*likelihood[imax,jmax]) / 4.
		hessian[1,1] = (likelihood[imax,jmax+2] + likelihood[imax,jmax-2] - 2*likelihood[imax,jmax]) / 4.
		hessian[1,0] = (likelihood[imax+1,jmax+1] + likelihood[imax-1,jmax-1] - likelihood[imax+1,jmax-1] - likelihood[imax-1,jmax+1]) / 4.

		#Scale the hessian to parameter units
		hessian[0,0] /= (self.unit[parameters[0]]**2)
		hessian[1,1] /= (self.unit[parameters[1]]**2)
		hessian[1,0] /= (self.unit[parameters[0]]*self.unit[parameters[1]])
		hessian[0,1] = hessian[1,0]

		#Compute the parameter covariance matrix
		pcov = np.linalg.inv(-hessian)

		#Plot an ellipse for each p-value
		for n,level in enumerate(levels):
			self.ax.add_artist(FisherAnalysis.ellipse(center,pcov,level,edgecolor=colors[n],**kwargs))


	######################################################################
	##############Plot the contours on top of the likelihood##############
	######################################################################

	def plotContours(self,colors=["red","green","blue"],display_percentages=True,display_maximum=True,color_maximum="green",fill=False,**kwargs):

		"""
		Display the confidence likelihood contours

		"""

		if not hasattr(self,"likelihood_values"):
			self.getLikelihoodValues(levels=[0.683,0.95,0.997])

		if hasattr(self,"reduced_likelihood"):
			reduced_likelihood = self.reduced_likelihood
		else:
			reduced_likelihood = self.likelihood

		assert len(colors) >= len(self.likelihood_values)
		assert reduced_likelihood.ndim==2,"this routine plots 2D contours only!!"

		extent = self.extent
		likelihood = reduced_likelihood.transpose()
		values = self.likelihood_values

		unit_j = (extent[1] - extent[0])/(likelihood.shape[1] - 1)
		unit_i = (extent[3] - extent[2])/(likelihood.shape[0] - 1) 

		#Build contour levels
		fmt = dict()
		
		for n,value in enumerate(values):
			fmt[value] = "{0:.1f}%".format(self.computed_p_values[n]*100)

		if fill:
			self.contour = self.ax.contourf(likelihood,values,colors=colors,origin="lower",extent=extent,aspect="auto",**kwargs)
		else:
			self.contour = self.ax.contour(likelihood,values,colors=colors,origin="lower",extent=extent,aspect="auto",**kwargs)

		#Contour labels
		self.ax.proxy += [ plt.Rectangle((0,0),1,1,fc=color) for color in colors if color!=rc.__globals__["rcParams"]["axes.facecolor"] ]
		
		if display_percentages:
			plt.clabel(self.contour,fmt=fmt,inline=1,fontsize=9)

		if display_maximum:
			
			#Find the maximum
			likelihood_max = likelihood.max()
			imax,jmax = np.where(likelihood==likelihood_max)

			#Plot scaling to physical values
			self.ax.plot(extent[0] + np.arange(likelihood.shape[1])*unit_j,np.ones(likelihood.shape[1])*imax[0]*unit_i + extent[2],linestyle="--",color=color_maximum)
			self.ax.plot(extent[0] + np.ones(likelihood.shape[0])*jmax[0]*unit_j,extent[2] + np.arange(likelihood.shape[0])*unit_i,linestyle="--",color=color_maximum)


	##################################################################################################
	#################Plot the likelihood marginalized over all parameters except one##################
	##################################################################################################

	def plotMarginal(self,parameter,levels=[0.684],colors=["red","blue","green"],alpha=0.5,fill=False):

		"""
		Plot the likelihood function marginalized over all parameters except one

		"""

		#Compute marginalized likelihood
		p,l,par_max,par_extremes = self.marginal(parameter,levels=levels)

		#Plot the likelihood
		self.ax.plot(p,l)

		#Plot the confidence contours
		for n,level in enumerate(levels):

			relevant_indices = np.where((p>=par_extremes[n][0])*(p<=par_extremes[n][1]))[0]
			if fill:
				self.ax.fill_between(p[relevant_indices],np.ones_like(relevant_indices)*l.min(),l[relevant_indices],facecolor=colors[n],alpha=alpha)
			else:
				self.ax.plot(np.ones(100)*p[relevant_indices[0]],np.linspace(l.min(),l[relevant_indices[0]],100),color=colors[n])
				self.ax.plot(np.ones(100)*p[relevant_indices[-1]],np.linspace(l.min(),l[relevant_indices[-1]],100),color=colors[n])


		#Labels
		self.ax.set_xlabel(self.parameter_labels[parameter],fontsize=22)
		self.ax.set_ylabel(r"$\mathcal{L}$"+"$($"+self.parameter_labels[parameter]+"$)$",fontsize=22)