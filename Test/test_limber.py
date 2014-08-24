try:
	
	from lenstools import LimberIntegrator
	from lenstools.defaults import load_power_default

except ImportError:
	
	import sys
	sys.path.append("..")
	from lenstools import LimberIntegrator
	from lenstools.defaults import load_power_default

import numpy as np
import matplotlib.pyplot as plt

from astropy.cosmology import WMAP9


def test_convergence_power():

	l = np.logspace(0.0,5.0,100.0)

	integrator = LimberIntegrator(cosmoModel=WMAP9)
	integrator.load3DPowerSpectrum(load_power_default,"Data/camb_output","fiducial_matterpower_")

	Cl = integrator.convergencePowerSpectrum(l)

	plt.plot(l,l*(l+1)*Cl/(2.0*np.pi))
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel(r"$l$")
	plt.ylabel(r"$l(l+1)C_l/2\pi$")

	plt.savefig("limber_power.png") 