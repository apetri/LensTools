import os

from ..simulations.limber import LimberIntegrator
from ..utils.defaults import load_power_default

from .. import dataExtern

import numpy as np
import matplotlib.pyplot as plt

from astropy.cosmology import WMAP9


def test_convergence_power():

	l = np.logspace(0.0,5.0,100)

	integrator = LimberIntegrator(cosmoModel=WMAP9)
	integrator.load3DPowerSpectrum(load_power_default,os.path.join(dataExtern(),"camb_output"),"fiducial_matterpower_")

	Cl = integrator.convergencePowerSpectrum(l)

	plt.plot(l,l*(l+1)*Cl/(2.0*np.pi))
	plt.xscale("log")
	plt.yscale("log")

	plt.xlabel("l")
	plt.ylabel("l(l+1)C_l/2*pi")

	try:
		plt.savefig("limber_power.png")
	except:
		pass