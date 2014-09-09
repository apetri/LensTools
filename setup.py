import os,sys,re

name = "lenstools"
me = "Andrea Petri"
email = "apetri@phys.columbia.edu"
url = "http://www.columbia.edu/~ap3020/LensTools/html"
default_gsl = "/usr/local"

try:
	import numpy.distutils.misc_util 
except ImportError:
	print("Please install numpy!")
	sys.exit(1)


try:
	from setuptools import setup,Extension
except ImportError:
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
simulations_dir = "simulations"
observations_dir = "observations"

#List external package sources here
external_sources["_topology"] = ["_topology.c","differentials.c","peaks.c","minkowski.c","coordinates.c","azimuth.c"]
external_sources["_gadget"] = ["_gadget.c","read_gadget_header.c"]

#Check for GSL installation, necessary for using the Design feature
gsl_required_includes = ["gsl_permutation.h","gsl_randist.h","gsl_rng.h","gsl_matrix.h"]
gsl_required_links = ["libgsl.a","libgslcblas.a"]
gsl_location = raw_input("Enter the location of your GSL installation (default '{0}'): ".format(default_gsl))
gsl_ok = True

if gsl_location == "":
	gsl_location = default_gsl

#Check for required GSL includes and links
for include in gsl_required_includes:
	
	include_filename = os.path.join(gsl_location,"include","gsl",include)
	sys.stderr.write("Checking if {0} exists... ".format(include_filename))

	if os.path.isfile(include_filename):
		sys.stderr.write("[OK]\n")
	else:
		sys.stderr.write("[FAIL]\n")
		gsl_ok = False
		break

for lib in gsl_required_links:

	lib_filename = os.path.join(gsl_location,"lib",lib)
	sys.stderr.write("Checking if {0} exists... ".format(lib_filename))

	if os.path.isfile(lib_filename):
		sys.stderr.write("[OK]\n")
	else:
		sys.stderr.write("[FAIL]\n")
		gsl_ok = False
		break

#Decide if we can install the Design feature, if not throw a warning
if gsl_ok:
	print("[OK] Checked GSL installation, the Design feature will be installed")
	lenstools_includes = [ os.path.join(gsl_location,"include") ]
	lenstools_link = ["-lm","-L {0}".format(os.path.join(gsl_location,"lib")),"-lgsl","-lgslcblas"]
	external_sources["_design"] = ["_design.c","design.c"] 
else:
	raw_input("[FAIL] GSL installation not found, the Design feature will not be installed, please press a key to continue: ")
	lenstools_includes = list()
	lenstools_link = ["-lm"]

#Extension objects
ext = list()

for ext_module in external_sources.keys():

	sources = list()
	for source in external_sources[ext_module]:
		sources.append(os.path.join(name,external_dir,source))

	ext.append(Extension(ext_module,sources,extra_link_args=lenstools_link))

#Data files on which the package depends on
package_data = {name:[os.path.join(simulations_dir,"book","CFHTemu1.txt")]}

#Append numpy includes
lenstools_includes += numpy.distutils.misc_util.get_numpy_include_dirs()

setup(
	name=name,
	version=version,
	author=me,
	author_email=email,
	packages=[name,"{0}.{1}".format(name,external_dir),"{0}.{1}".format(name,simulations_dir),"{0}.{1}".format(name,observations_dir)],
	package_data=package_data,
	install_requires=["numpy","scipy","astropy"],
	url=url,
	license="?",
	description="Toolkit for Weak Gravitational Lensing analysis",
	long_description=rd(os.path.join("docs","source","index.rst")),
	classifiers=classifiers,
	ext_package=os.path.join(name,external_dir),
	ext_modules=ext,
	include_dirs=lenstools_includes,
)
