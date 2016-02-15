import ast
from ..statistics.contours import ContourPlot

default_colors = ["#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000"]

def main(filename,cmd_args,options):

	#Parse the parameter names
	parameter_axes = dict()
	for n,par in enumerate(cmd_args.axes.split(",")):
		parameter_axes[par] = n

	#Get the parameter labels from the options
	cosmo_labels = dict()
	for par in parameter_axes.keys():
		cosmo_labels[par] = ast.literal_eval(options.get("labels",par))

	#Build the contour plot with the ContourPlot class handler
	contour = ContourPlot()
	#Load the likelihood
	contour.getLikelihood(filename,parameter_axes=parameter_axes,parameter_labels=cosmo_labels)
	#Set the physical units
	contour.getUnitsFromOptions(options)

	#Parse levels
	levels = [ float(l) for l in cmd_args.levels.split(",") ]

	#Choose the colors
	if cmd_args.colors is not None:
		colors = cmd_args.colors.split(",")
	else:
		colors = default_colors

	#Decide if marginalize over one of the parameters or take a slice over it
	assert (cmd_args.marginalize is None)+(cmd_args.slice is None)+(cmd_args.marginal is None)==2,"You cannot both marginalize and slice over a parameter!"

	if cmd_args.marginalize is not None:
		contour.marginalize(cmd_args.marginalize)

	if cmd_args.slice is not None:
		slice_over,value = cmd_args.slice.split("=")
		value = float(value)
		contour.slice(slice_over,value)

	if cmd_args.marginal is None:

		#Show the full likelihood
		contour.show()

		#Compute the likelihood levels
		contour.getLikelihoodValues(levels=levels)

		#Display the contours
		contour.plotContours(colors=colors,fill=cmd_args.fill,display_percentages=cmd_args.display_percentages,display_maximum=cmd_args.display_maximum)

		#Optionally give a title to the figure
		if cmd_args.title is not None:
			contour.title_label = cmd_args.title
		else:
			contour.title_label = ""

		#Labels
		contour.labels()

	else:

		#Plot the marginal likelihood
		contour.plotMarginal(cmd_args.marginal,levels=levels,colors=colors,fill=cmd_args.fill)

	
	#Return the contour object
	return contour

