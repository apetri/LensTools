import os

from ..simulations import Gadget2SnapshotDE
from ..pipeline.settings import Gadget2Settings

from .. import dataExtern

import numpy as np
from astropy.units import Mpc,m,s

def test_read():

	#Open the gadget snapshot
	snapshot = Gadget2SnapshotDE.open(os.path.join(dataExtern(),"gadget/snapshot_001"))

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

def test_write():

	#Create an empty gadget snapshot
	snap = Gadget2SnapshotDE()

	#Generate random positions and velocities
	NumPart = 32**3
	x = np.random.normal(loc=7.0,scale=5.0,size=(NumPart,3)) * Mpc
	v = np.random.uniform(-1,1,size=(NumPart,3)) * m / s

	#Put the particles in the snapshot
	snap.setPositions(x)
	snap.setVelocities(v)

	#Generate minimal header
	snap.setHeaderInfo()

	#Visualize
	snap.visualize(s=1)
	snap.savefig("gadget_initial_condition.png")

	#Write the snapshot
	snap.write("gadget_ic")


def test_paramfile():

	#Create an empty gadget snapshot
	snap = Gadget2SnapshotDE()

	#Generate random positions and velocities
	NumPart = 32**3
	x = np.random.normal(loc=7.0,scale=5.0,size=(NumPart,3)) * Mpc
	v = np.zeros((NumPart,3)) * m / s

	#Put the particles in the snapshot
	snap.setPositions(x)
	snap.setVelocities(v)

	#Generate minimal header
	snap.setHeaderInfo()

	#Split the particles between two files
	snap.write("gadget_sphere",files=2)

	#Generate the parameter file that will determine the evolution
	settings = Gadget2Settings.default()
	snap.writeParameterFile("gadget_sphere.param",settings)

	#Generate a file with the scale factor of the output snapshots
	z = np.arange(90.0,0.0,-10.0)
	a = 1.0 / (1 + z)
	np.savetxt("outputs.txt",a)

