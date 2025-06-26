from lenstools.simulations import PotentialPlane
import matplotlib.pyplot as plt

for i,n in enumerate(range(11,59)):

	plane_file = "snap{0}_potentialPlane0_normal0.fits".format(n)

	print("Plotting plane at {0}".format(plane_file))
	pl = PotentialPlane.load(plane_file)
	dens = pl.density()
	
	dens.visualize(colorbar=True)
	dens.ax.set_title(r"$z={0:.2f}$  $N={1:d}$".format(pl.redshift,int(pl.num_particles)))
	dens.savefig("plane{0}.png".format(i))
	plt.close(dens.fig)
