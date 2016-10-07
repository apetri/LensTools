import os,sys,glob,re
import platform

if sys.version_info.major>=3:
	import configparser as cfg
else:
	import ConfigParser as cfg


#Names
name = "lenstools"
me = "Andrea Petri"
email = "apetri@phys.columbia.edu"
url = "http://lenstools.readthedocs.org"
default_cfg = "setup.cfg"

#Sub-packages
sub_package_names = ["image","catalog","statistics","utils","simulations","observations","pipeline","legacy","scripts","extern","tests"]
packages = [ name ]
for sub_name in sub_package_names:
	packages.append("{0}.{1}".format(name,sub_name))

#External sub-packages
external_dir = "extern"
external_support_dir = "cextern"

try:
	import numpy.distutils.misc_util 
except ImportError:
	print("Please install numpy!")
	sys.exit(1)

try:
	from setuptools import setup,Extension
except ImportError:
	from distutils.core import setup,Extension

#Print in color
def red(s):
	return '\033[31m' + s + '\033[39m'

def green(s):
	return '\033[32m' + s + '\033[39m'

def yellow(s):
	return '\033[33m' + s + '\033[39m'


#Shortcut for reading file
def rd(filename):
	with open(filename,"r") as fp:
		return fp.read()

#Check GSL installation, necessary for using the Design feature
def check_gsl(conf):
	
	gsl_location = conf.get("gsl","installation_path")
	gsl_required_includes = ["gsl_permutation.h","gsl_randist.h","gsl_rng.h","gsl_matrix.h"]
	gsl_required_links = ["libgsl.a","libgslcblas.a"]

	#Check for required GSL includes and links
	for include in gsl_required_includes:
	
		include_filename = os.path.join(gsl_location,"include","gsl",include)
		sys.stderr.write("Checking if {0} exists... ".format(include_filename))

		if os.path.isfile(include_filename):
			sys.stderr.write(green("[OK]\n"))
		else:
			sys.stderr.write(red("[FAIL]\n"))
			return None

	for lib in gsl_required_links:

		lib_filename = os.path.join(gsl_location,"lib",lib)
		sys.stderr.write("Checking if {0} exists... ".format(lib_filename))

		if os.path.isfile(lib_filename):
			sys.stderr.write(green("[OK]\n"))
		else:
			sys.stderr.write(red("[FAIL]\n"))
			return None

	return gsl_location


#Check fftw3 installation, required from NICAEA
def check_fftw3(conf):

	fftw3_location = conf.get("fftw3","installation_path")
	fftw3_required_includes = ["fftw3.h"]
	fftw3_required_links = ["libfftw3.a"]

	#Check for required GSL includes and links
	for include in fftw3_required_includes:
	
		include_filename = os.path.join(fftw3_location,"include",include)
		sys.stderr.write("Checking if {0} exists... ".format(include_filename))

		if os.path.isfile(include_filename):
			sys.stderr.write(green("[OK]\n"))
		else:
			sys.stderr.write(red("[FAIL]\n"))
			return None

	for lib in fftw3_required_links:

		lib_filename = os.path.join(fftw3_location,"lib",lib)
		sys.stderr.write("Checking if {0} exists... ".format(lib_filename))

		if os.path.isfile(lib_filename):
			sys.stderr.write(green("[OK]\n"))
		else:
			sys.stderr.write(red("[FAIL]\n"))
			return None

	return fftw3_location


#Check correctness of NICAEA installation
def check_nicaea(conf):

	nicaea_root = conf.get("nicaea","installation_path")

	#These are the include and lib directory
	nicaea_include = os.path.join(nicaea_root,"include","nicaea")
	nicaea_lib = os.path.join(nicaea_root,"lib")

	#Check for their existence
	sys.stderr.write("Checking for {0}...".format(nicaea_include))
	if os.path.isdir(nicaea_include):
		sys.stderr.write(green("[OK]\n"))
	else:
		sys.stderr.write(red("[FAIL]\n"))
		return None

	sys.stderr.write("Checking for {0}...".format(os.path.join(nicaea_lib,"libnicaea.a")))
	if os.path.isfile(os.path.join(nicaea_lib,"libnicaea.a")):
		sys.stderr.write(green("[OK]\n"))
	else:
		sys.stderr.write(red("[FAIL]\n"))
		return None
	
	return nicaea_include,nicaea_lib


############################################################
#####################Execution##############################
############################################################

#Read system dependent configuration file
conf = cfg.ConfigParser()
this = os.getenv("THIS")

if (this is not None) and (os.path.isfile(default_cfg+"."+this)):
	cfg_file = default_cfg+"."+this
else:
	cfg_file = default_cfg

print("Reading system dependent configuration from {0}".format(cfg_file))
conf.read([cfg_file])

#Override library locations with the ones provided from the command line
for l in ["fftw3","gsl","nicaea"]:
	for arg in sys.argv:
		
		if "--{0}".format(l) in arg:
			
			loc = arg.split("=")[-1]
			print(yellow("Command line override: looking for {0} in {1}".format(l,loc)))
			conf.set(l,"installation_path",loc)

			if l=="nicaea":
				conf.set(l,"install_python_bindings","True")

			#Purge from command line arguments
			sys.argv.remove(arg)

#########################################################################################

vre = re.compile("__version__ = \"(.*?)\"")
m = rd(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "lenstools", "__init__.py"))
version = vre.findall(m)[0]
download_url = "https://github.com/apetri/LensTools/archive/{0}.tar.gz".format(version)

classifiers = [
		"Development Status :: 4 - Beta",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Programming Language :: C",
		"License :: OSI Approved :: MIT License"
	]

external_sources = dict()
external_support = dict()

#Includes
lenstools_includes = list()

#List external package sources here
external_sources["_topology"] = ["_topology.c","differentials.c","peaks.c","minkowski.c","coordinates.c","azimuth.c"]
external_sources["_gadget2"] = ["_gadget2.c","read_gadget_header.c","read_gadget_particles.c","write_gadget_particles.c"]
external_sources["_nbody"] = ["_nbody.c","grid.c","coordinates.c"]
external_sources["_pixelize"] = ["_pixelize.c","grid.c","coordinates.c"]

######################################################################################################################################

#Decide if we can install the Design feature, if not throw a warning
gsl_location = check_gsl(conf)

if gsl_location is not None:
	print(green("[OK] Checked GSL installation, the Design feature will be installed"))
	lenstools_includes.append(os.path.join(gsl_location,"include")) 
	lenstools_link = ["-lm","-L{0}".format(os.path.join(gsl_location,"lib")),"-lgsl","-lgslcblas"]
	external_sources["_design"] = ["_design.c","design.c"] 
else:
	print(red("[FAIL] GSL installation not found, the Design feature will not be installed"))
	lenstools_link = ["-lm"]


######################################################################################################################################

if conf.getboolean("nicaea","install_python_bindings"):

	#Decide if we can install the NICAEA bindings
	fftw3_location = check_fftw3(conf)
	nicaea_location = check_nicaea(conf)
	
	if (gsl_location is not None) and (fftw3_location is not None) and (nicaea_location is not None):
	
		print(green("[OK] Checked GSL,FFTW3 and NICAEA installations, the NICAEA bindings will be installed"))
		
		#Add necessary includes and links
		lenstools_includes += [ os.path.join(fftw3_location,"include") , nicaea_location[0] ]
		lenstools_link += ["-L{0}".format(os.path.join(fftw3_location,"lib")) , "-lfftw3", "-L{0}".format(nicaea_location[1]) , "-lnicaea"]

		#Specify the NICAEA extension
		external_sources["_nicaea"] = ["_nicaea.c"]

	else:
		print(red("[FAIL] NICAEA bindings will not be installed (either enable option or check GSL/FFTW3/NICAEA installations)"))


#################################################################################################
#############################Package data########################################################
#################################################################################################

#Data files on which the package depends on
package_data = dict()
package_data[name] = [ os.path.join("data",filename) for filename in os.listdir(os.path.join(name,"data")) if os.path.isfile(os.path.join(name,"data",filename)) ]
package_data[name] += [ os.path.join(external_dir,filename) for filename in os.listdir(os.path.join(name,external_dir)) if filename.endswith(".h") ]

#################################################################################################
#############################Additional includes#################################################
#################################################################################################

#Append numpy includes
lenstools_includes += numpy.distutils.misc_util.get_numpy_include_dirs()

#Append system includes (fix OSX clang includes)
if platform.system() in ["Darwin","darwin"]:
	lenstools_includes += ["/usr/local/include","/usr/include"]

#################################################################################################
#############################Scripts#############################################################
#################################################################################################

#package scripts
scripts = [ fname for fname in glob.glob(os.path.join("scripts","*")) if os.path.basename(fname)!="README.rst" ]

#################################################################################################
#############################Extensions##########################################################
#################################################################################################

#Extension objects
ext = list()

for ext_module in external_sources.keys():

	sources = list()
	for source in external_sources[ext_module]:
		sources.append(os.path.join(name,external_dir,source))

	#Append external support sources too
	if ext_module in external_support.keys():
		sources += external_support[ext_module]

	ext.append(Extension(ext_module,
                             sources,
                             extra_link_args=lenstools_link,
                             include_dirs=lenstools_includes))

#################################################################################################
#############################Setup###############################################################
#################################################################################################

setup(
	name=name,
	version=version,
	author=me,
	author_email=email,
	packages=packages,
	package_data=package_data,
	install_requires=["numpy","scipy","pandas","matplotlib","astropy","emcee"],
	url=url,
	download_url=download_url,
	license="MIT",
	description="Toolkit for Weak Gravitational Lensing analysis",
	long_description=rd("README.rst") + "\n\n" + "Changelog\n" + "---------\n\n" + rd("HISTORY.rst"),
	scripts=scripts,
	classifiers=classifiers,
	ext_package=os.path.join(name,external_dir),
	ext_modules=ext,
	include_dirs=lenstools_includes,
	zip_safe=False,
)
