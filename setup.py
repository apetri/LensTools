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


version = "0.1"


setup(
	name="lenstools",
	version=version,
	author="Andrea Petri",
	author_email="apetri@phys.columbia.edu",
	packages=["limber"],
	url="https://github.com/apetri/LensTools",
	license="None",
	description="Toolkit for Weak Gravitational Lensing analysis",
	long_description=rd("README.md"),
	install_requires=["numpy","scipy","astropy"],
	classifiers=[
		"Development Status :: 1 - Planning",
		"Intended Audience :: Science/Research",
		"Operating System :: OS Independent",
		"Programming Language :: Python"
	],
)
