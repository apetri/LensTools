import sys
sys.path.append("..")

from lenstools.simulations import Gadget2Snapshot
from lenstools.utils import MPIWhirlPool
from scipy.ndimage import filters
from mayavi import mlab

#Initialize MPIWhirlPool
try:
	pool = MPIWhirlPool()
except:
	pool = None


if pool is None:

	print("Processing ICs in series...")

	#If we run on one core only the snapshots need to be processed in series
	snap = Gadget2Snapshot.open("Test/Data/gadget/ic1.0")
	n1,r = snap.numberDensity(resolution=64,left_corner=np.array([0.0,0.0,0.0])*snap.Mpc_over_h)
	snap.close()

	snap = Gadget2Snapshot.open("Test/Data/gadget/ic1.1")
	n2,r = snap.numberDensity(resolution=64,left_corner=np.array([0.0,0.0,0.0])*snap.Mpc_over_h)
	snap.close()

	n = n1+n2
	ns = filters.gaussian_filter(n,2)
	scene = mlab.pipeline.volume(mlab.pipeline.scalar_field(ns))

else:

	if pool.is_master():
		print("Processing ICs in parallel on {0} cores".format(pool.size+1))

	#If we run on 2 cores, we can process the snapshots in parallel
	snap = Gadget2Snapshot.open("Test/Data/gadget/ic1",pool=pool)
	n,r = snap.numberDensity(resolution=64,left_corner=np.array([0.0,0.0,0.0])*snap.Mpc_over_h)
	snap.close()

	if pool.is_master():
		ns = filters.gaussian_filter(n,2)
		scene = mlab.pipeline.volume(mlab.pipeline.scalar_field(ns))
