try:
	
	from lenstools.simulations import Design

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.simulations import Design

import matplotlib.pyplot as plt


#Test the visualization of a design
def test_visualize():

	design = Design()

	#Add the dimensions
	design.add_parameter("Omega_m",min=0.1,max=0.9,label=r"$\Omega_m$")
	design.add_parameter("w",min=-2.0,max=-0.1,label=r"$w$")
	design.add_parameter("sigma8",min=0.01,max=1.6,label=r"$\sigma_8$")

	#Lay down 50 points
	design.put_points(50)
	print(design)

	#Visualize the 3d diagonal configuration
	design.visualize(color="blue")
	design.set_title("Cost={0:.2f}".format(design.diagonalCost(1.0)))
	design.savefig("design_diagonal_3d.png")

	#Visualize the 2d (Om,si8) projection
	fig,ax = plt.subplots()
	design.visualize(fig=fig,ax=ax,parameters=["Omega_m","sigma8"],color="red",marker=".")
	design.savefig("design_diagonal_2d.png")
