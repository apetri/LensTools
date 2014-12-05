from ..simulations import Gadget2Snapshot

def main(filename):

	#Open snapshot
	snap = Gadget2Snapshot.open(filename)

	#If velocities are not present, skip this snapshot
	try:
		snap.getVelocities(save=False)
	except IOError:
		print("{0} does not contain velocity information".format(filename))
		snap.close()
		return 

	#If velocities are present we are here; strip them
	snap.getPositions(save=True)
	snap.close()

	#Write the new snapshot
	print("{0} has been stripped of velocity information".format(filename))
	snap.write(filename)