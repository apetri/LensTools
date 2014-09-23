import numpy as np
from scipy.ndimage import filters

from mayavi import mlab

n = np.load("density_parallel.npy")
ns = filters.gaussian_filter(n,2)

scene = mlab.pipeline.volume(mlab.pipeline.scalar_field(ns))