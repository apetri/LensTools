from distutils.core import setup,Extension
import numpy.distutils.misc_util

setup(
	install_requires = ["numpy"],
	ext_modules=[Extension("_topology",["_topology.c","peaks.c","coordinates.c"])],
	include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)