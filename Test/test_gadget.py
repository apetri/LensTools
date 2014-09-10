try:
	
	from lenstools.simulations import Gadget2Snapshot

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools.simulations import Gadget2Snapshot

import numpy as np

def test_positions():

	#Open the gadget snapshot
	snapshot = Gadget2Snapshot.open("Data/gadget/snapshot_001")

	#Get the particles positions and play with the indices
	pospart = snapshot.getPositions(first=500,last=1000)
	pos = snapshot.getPositions()
	assert np.all(pospart==pos[500:1000])

	#Visualize the snapshot
	snapshot.visualize(s=1)

	#Save the visualization
	snapshot.ax.set_title(r"${0}^3$ particles".format(snapshot.header["num_particles_total_side"]))
	snapshot.savefig("gadget_snapshot.png")

	#Close the snapshot
	snapshot.close()