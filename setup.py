import os,sys,re

try:
	from setuptools import setup
	setup
except ImportError:
	from distutils.core import setup
	setup

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
	install_requires=["numpy","scipy","astropy"],
	classifiers=[
		"Development Status :: 2 - Pre-Alpha",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Programming Language :: Python",
		"Programming Language :: C",
		"License :: Public Domain"
	],
)
