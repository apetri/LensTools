from lenstools.simulations import Gadget2Snapshot

import numpy as np
from astropy.units import Mpc,m,s

########################Write#################################

#Create an empty gadget snapshot
snap = Gadget2Snapshot()

#Generate random positions and velocities
NumPart = 32**3
x = np.random.normal(loc=7.0,scale=5.0,size=(NumPart,3)) * Mpc
v = np.random.uniform(-1,1,size=(NumPart,3)) * m / s

#Put the particles in the snapshot
snap.setPositions(x)
snap.setVelocities(v)

#Generate minimal header
snap.setHeaderInfo()

#Write the snapshot
snap.write("gadget_ic")

######################Read and visualize#########################

#Open the snapshot
snap = Gadget2Snapshot.open("gadget_ic")

#Visualize the header
print(snap.header)

#Visualize the snapshot
snap.visualize(s=1)
snap.savefig("snapshot.png")
snap.close()

#################Actual Gadget snapshot##########################
###################Read and visualize############################

#Open the snapshot
snap = Gadget2Snapshot.open("../Test/Data/gadget/snapshot_001")

#Visualize the header
print(snap.header)

#Get positions and velocities
snap.getPositions()
snap.getVelocities()

#Visualize the snapshot
snap.visualize(s=1)
snap.savefig("snapshot_gadget.png")
snap.close()
