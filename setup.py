import os,sys,re

name = "lenstools"
me = "Andrea Petri"
email = "apetri@phys.columbia.edu"
url = "https://github.com/apetri/LensTools"

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

classifiers = [
		"Development Status :: 2 - Pre-Alpha",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Programming Language :: C",
		"License :: Public Domain"
	]

external_sources = dict()
external_dir = "external"

#List external package sources here
external_sources["_topology"] = ["_topology.c","differentials.c","peaks.c","minkowski.c","coordinates.c"]

#Extra compiler options
link_flags = ['-lm']

#Extension objects
ext = list()

for ext_module in external_sources.keys():

	sources = list()
	for source in external_sources[ext_module]:
		sources.append("{0}/{1}/{2}".format(name,external_dir,source))

	ext.append(Extension(ext_module,sources,extra_link_args=link_flags))

setup(
	name=name,
	version=version,
	author=me,
	author_email=email,
	packages=[name,"{0}.{1}".format(name,external_dir)],
	url=url,
	license="?",
	description="Toolkit for Weak Gravitational Lensing analysis",
	long_description=rd("README.md"),
	classifiers=classifiers,
	ext_package="{0}/{1}".format(name,external_dir),
	ext_modules=ext,
	include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)