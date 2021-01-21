import os
import pandas as pd
from .. import dataExtern

from .. import ShearCatalog,FlexionCatalog
from ..image.flexion import Spin1

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.table import Table


# Import an example of a simulated catalog and convert to astropy table
data = pd.read_pickle(os.path.join(dataExtern(),"flexionCatalog_SIS.pkl"))
catalog = Table.from_pandas(data)

# Map Initialization
mapSize = 3600 # in arcsec
mapOrigin = [-1800.,-1800.]
nPixel = 512
smoothFactor = 0.7 # in arcmin

#Config catalog for shear
catalog['gamma1'].name = 'shear1'
catalog['gamma2'].name = 'shear2'
shearCat = ShearCatalog(catalog)

#Use ShearCatalog class to create shear map
shearCat.setSpatialInfo('x','y',u.arcsec)
shearCat.setRedshiftInfo('z')
shearMap = shearCat.toMap(map_size=mapSize*u.arcsec, npixel=nPixel, origin=mapOrigin*u.arcsec, smooth=smoothFactor*u.arcmin)

# Use FlexionCatalog class to create flexion map
flexionCat = FlexionCatalog(catalog)
flexionCat.setSpatialInfo('x','y',u.arcsec)
flexionCat.setRedshiftInfo('z')
flexionMap = flexionCat.toMap(map_size=mapSize*u.arcsec, npixel=nPixel, origin=mapOrigin*u.arcsec, smooth=smoothFactor*u.arcmin)



def test_visualizeShear():
	
	shearMap.visualize(colorbar=True)
	shearMap.savefig('catalogToShear.png')
    
def test_visualizeFlexion():

	flexionMap.visualize(colorbar=True)
	flexionMap.savefig('catalogToFlexion.png')

def test_gradient():

	grad = Spin1(flexionMap.gradient(),angle=flexionMap.side_angle)
	fig,ax = plt.subplots(2,2,figsize=(16,16))

	#Plot the components
	grad.visualize(fig=fig,ax=ax.reshape(4),component_labels=(r"$F_{1,x}$",r"$F_{1,y}$",r"$F_{2,x}$",r"$F_{2,y}$"),colorbar=True)

	#Save
	fig.tight_layout()
	fig.savefig("flexiongradient.png")
	plt.clf()

def test_reconstruct():

	# Convergence plot setup
	fig,ax = plt.subplots(1,2,figsize=(16,8))
	
	# Create convergence map from shear map and visualize
	convMapFromShear = shearMap.convergence()
	convMapFromShear.visualize(fig=fig,ax=ax[0],colorbar=True)
	ax[0].set_title("Convergence Map from Shear")
    
	# Create convergence map from flexion map and visualize
	convMapFromFlexion = flexionMap.convergence()
	convMapFromFlexion.visualize(fig=fig,ax=ax[1],colorbar=True,cbar_label=r"$\kappa$")
	ax[1].set_title("Convergence Map from Flexion")
	ax[1].set_ylabel("")

	fig.tight_layout()
	plt.savefig("shear_flexion_reconstruction.png")
	plt.clf()
