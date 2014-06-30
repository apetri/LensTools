import os,sys,re

try:
	import numpy.distutils.misc_util 
except ImportError:
	print("Please install numpy!")
	sys.exit(1)


from distutils.core import setup,Extension

def rd(filename):
	
	f = file(filename,"r")
	r = f.read()
	f.close()

	return r


vre = re.compile("__version__ = \"(.*?)\"")
m = rd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "lenstools", "__init__.py"))
version = vre.findall(m)[0]

setup(
	name="lenstools",
	version=version,
	author="Andrea Petri",
	author_email="apetri@phys.columbia.edu",
	packages=["lenstools"],
	url="https://github.com/apetri/LensTools",
	license="?",
	description="Toolkit for Weak Gravitational Lensing analysis",
	long_description=rd("README.md"),
	classifiers=[
		"Development Status :: 2 - Pre-Alpha",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Programming Language :: C",
		"License :: Public Domain"
	],
	ext_package="lenstools/external",
	ext_modules=[Extension("_topology",["lenstools/external/_topology.c","lenstools/external/differentials.c","lenstools/external/peaks.c","lenstools/external/minkmom.c","lenstools/external/coordinates.c"])],
	include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)