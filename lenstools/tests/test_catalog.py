import sys,os

from .. import dataExtern
from ..catalog import ShearCatalog

import matplotlib.pyplot as plt
import astropy.units as u

#Reconstruct shear map from catalogs
def test_reconstruct():

	#Options
	map_size = 3.5*u.deg
	npixel = 512
	smooth = 0.1*u.arcmin
	zbins = [(0.0,0.5),(0.5,0.7),(0.7,0.9),(0.9,1.2),(1.2,3.0)]

	#Files
	pos_files = [os.path.join(dataExtern(),"catalog","positions_bin{0}.fits".format(n)) for n in range(1,6)]
	shear_files = [os.path.join(dataExtern(),"catalog","1-2","WLshear_positions_bin{0}_0001r.fits".format(n)) for n in range(1,6)]

	#Set up plot
	fig,ax = plt.subplots(2,3,figsize=(24,16))
	ax = ax.reshape(6)

	#Read in the full catalog
	full_catalog = ShearCatalog.readall(shear_files,pos_files)

	#Plot in the last panel
	full_convergence = full_catalog.toMap(map_size=map_size,npixel=npixel,smooth=smooth).convergence()
	full_convergence.visualize(fig=fig,ax=ax[-1],colorbar=True,cbar_label=r"$\kappa$")
	ax[-1].set_title("All",fontsize=30)

	#Rebin galaxies according to bins
	catalog_rebinned = full_catalog.rebin(zbins)
	for n,ct in enumerate(catalog_rebinned):
		convergence = ct.toMap(map_size=map_size,npixel=npixel,smooth=smooth).convergence()
		convergence.visualize(fig=fig,ax=ax[n],colorbar=True,cbar_label=r"$\kappa$")
		ax[n].set_title(r"$z\in[{0:.1f},{1:.1f}]$".format(*zbins[n]),fontsize=30)

	#Save 
	fig.tight_layout()
	fig.savefig("catalog_to_convergence.png")




