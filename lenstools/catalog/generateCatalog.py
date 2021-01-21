"""

.. module:: generateCatalog
	:platform: Unix
	:synopsis: This module can generate example catalogs consisting of lensing fields.
                Precisely useful for testing out flexion analysis

.. moduleauthor:: Evan J. Arena <eja55@drexel.edu>
                ... and edited by Brij Patel <brp53@drexel.edu>

"""

import numpy as np
from astropy.cosmology import Planck15
import pandas as pd
from astropy import units as u

"""

Description:
.. We want to simulate a small galaxy survey that consists of some relatively small angle on the sky.
   This galaxy "catalog" is simply a dataframe, where each "galaxy" has a (random) position on the
   sky, a (random) distance from the observer, a (random) intrinsic ellipticity, a (random) intrinsic flexion, 
   and an angular size.  Lensing fields are then applied to each galaxy to create lensing signals, 
   which can be used to test convergence reconstruction.  The user can choose a few lensing flexion signals on the sky.  
   The default in the examples below is the softened Isothermal Sphere (sIS).

.. So, what one would observe when looking at this patch of sky is an observed ellipticity/shear, 
   which is the vector sum of the intrinsic ellipticity and the shear, and an observed flexion, 
   which is the vector sum of the intrinsic flexion and the lensing flexion.

.. This catalog is stored in a pandas dataframe and exported as a pickle file.

*********************************** Example 1 *********************************
.. Generate a catalog with no lensing signal, only intrinstic ellipticity and flexion.  
   Therefore, this is just noise meaning it's not useful for convergence reconstruction.

cat = MockCatalog(cat_name='toy_cf_test_mock_catalog_intrinsic_flexion')
cat.generateSourceGalaxies()
cat.exportCatalog()

*********************************** Example 2 *********************************
.. Generate a catalog with a sIS.  We, unrealistically, increase n_eff, and make sigma_v very large, 
   so that we can actually see the sIS shear flexion signals above the intrinsic noise.
   This is useful to test out convergence reconstruction from shear and flexion.

sIS_cat = MockCatalog(cat_name='toy_cf_test_mock_catalog_intrinsic_flexion_w_sIS', n_eff = 10)
sIS_cat.generateSourceGalaxies()
sIS_cat.sIScorrSig(sigma_v = 2500)
sIS_cat.exportCatalog()

"""

class MockCatalog:
   
    def __init__(self, cat_path='',cat_name='toy_cf_test_mock_catalog',
                  solid_angle=1, zmin=0.20, zmax=0.43, n_eff = 1.52, sigma_aF=0.0370):
        
        """The default parameters for redshift and number density come from redshift bin 1 of
           the DES Y1 shape catalog. The default scatter in intrinsic flexion comes from that of
           early-type objects found in Fabritius, Arena, and Goldberg arXiv:2006.03506
        """
        
        self.cosmo = Planck15                             # Make default cosmology Planck15

        self.cat_path = cat_path                          # Path to catalog. Default is working directory
        self.cat_name = cat_name                          # Name of catalog

        self.solid_angle = solid_angle*(u.deg**2.)        # Solid angle of survey in square degrees.
        self.zmin = zmin                                  # Minimum redshift of catalog.
        self.zmax = zmax                                  # Maximum redshift of catalog
        self.n_eff = n_eff                                # Effective number density of galaxies/arcmin^2
        self.sigma_aF = sigma_aF                          # Scatter in intrinsic flexion

        self.kappa_list = 0                               # ---
        self.gamma1_list = 0                              # Make all lensing fields zero by default.
        self.gamma2_list = 0                              # This can be overwritten if a function is called
        self.F1_list = 0                                  #  that generates a correlated lensing signal.
        self.F2_list = 0                                  # -- 

    def generateSourceGalaxies(self, seeds=np.array((0,1,2,3,4,5))):
        
        """Generate the background galaxies with no lensing fields. we will create N galaxies described by the following parameters:
             (i).   An (x,y) coordinate (we are just treating the galaxies as points here) in units of arcsec
             (ii).  A redshift drawn from the redshift range of the chosen bin
             (iii). An angular size, a
             (iv).  An intrinsic ellipticity eps1^(s), eps2^(s).
             (v).   An intrinsic first-flexion F1^(s), F2^(s).
        """
        
        # Based on the solid angle and number density, we can calculate the total number of galaxies
        solid_angle_arcmin2 = self.solid_angle.to(u.arcmin**2).value
        N = int(solid_angle_arcmin2*self.n_eff**2)

        # Next we will create N galaxies described by the parameters
        # .. Lay down grid coordinates in arcsec. We place the origin at the center of the square patch of sky
        xlim = (np.sqrt(self.solid_angle).to(u.arcsec).value)/2
        ylim = xlim
        self.xlim = xlim
        self.ylim = ylim
        np.random.seed(seeds[0])
        x_list = np.random.uniform(-xlim, xlim, N)
        np.random.seed(seeds[1])
        y_list = np.random.uniform(-ylim, ylim, N)
        
        # .. Assign redshifts
        z_list = np.random.uniform(self.zmin, self.zmax, N)
        
        # .. Given the redshifts, get the corresponding angular diameter distances (in Mpc) to each galaxy
        D_s_list = self.cosmo.angular_diameter_distance(z_list).value
        
        # .. Assume a typical physical galaxy size of 10 kpc, then we can get the angular size of each galaxy in arcsec
        R_gal = (10*u.kpc).to(u.Mpc).value
        a_list = ((R_gal/D_s_list)*u.rad).to(u.arcsec).value
        
        # .. Generate the intrinsic ellipticities
        # .. .. We will make use of the intrisic ellipticity distribution given in Schneider 1996:
        np.random.seed(seeds[2])
        sigma_eps = 0.2
        
        # .. .. Generate a list of instrinsic ellipticities twice as long as needed, so we can truncate anywhere where |eps_s| > 0.5
        eps_s_list = (1/(np.sqrt(np.pi)*sigma_eps*(1-np.exp(-1/sigma_eps**2.))))*np.random.normal(loc=0,scale=sigma_eps/np.sqrt(2),size=(2*N,))
        id_eps_s = np.where(abs(eps_s_list) < 0.5)
        eps_s_list = eps_s_list[id_eps_s]
        eps_s_list = eps_s_list[0:N]
        
        # .. .. Draw orientation angles from a uniform distribution
        phi_low, phi_up = 0., 2*np.pi
        np.random.seed(seeds[3])
        phi_list = np.random.uniform(low=phi_low, high=phi_up, size=(N,))
        
        # .. .. Get intrinsic ellipticity components
        eps1_s_list, eps2_s_list = eps_s_list*np.cos(2.*phi_list), eps_s_list*np.sin(2.*phi_list)
        
        # .. Draw the intrinsic flexion aF from a Gaussian distribution
        np.random.seed(seeds[4])
        aF1_list = np.random.normal(loc=0, scale=self.sigma_aF/np.sqrt(2), size=N)
        np.random.seed(seeds[5])
        aF2_list = np.random.normal(loc=0, scale=self.sigma_aF/np.sqrt(2), size=N)
        
        # .. Get intrisic flexion F1 and F2
        F1_s_list = aF1_list/a_list
        F2_s_list = aF2_list/a_list

        self.x_list = x_list
        self.y_list = y_list
        self.z_list = z_list
        self.D_s_list = D_s_list
        self.a_list = a_list
        self.eps1_s_list = eps1_s_list
        self.eps2_s_list = eps2_s_list
        self.F1_s_list = F1_s_list
        self.F2_s_list = F2_s_list

    def sinusoidalCorrSig1D(self, x=None, kappa0=10., lam=(10*u.arcmin).to(u.arcsec).value):
        
        """Lensing fields varying in x only.  The kappa field is simply a cosine wave varying in x.
           The characteristic scale is given by lambda, which has a default value of 10 arcminutes.
        """
        
        if x == None:
            x = self.x_list
        
        k = 2*np.pi/lam
        kappa = kappa0*np.cos(k*x)
        
        #  Get the shear components.  
        #  .. We can (cleverly) use the fact that psi,22 = kappa - gamma1
        #     (and noting that if kappa only varies in x then the lensing term psi,22 = 0)
        #     then we clearly have gamma1 = kappa.  We should also have gamma2 = 0 by the same token
        #     that kappa only varies in x.  Then,
        gamma1 = kappa
        gamma2 = 0*gamma1
        
        #  Now, we can create a correlated flexion signal, which is simply the derivative of the kappa field
        #   (We simply take the analytical derivative of the above kappa function)
        F1 = -k*kappa0*np.sin(k*x)
        F2 = 0*F1

        self.kappa_list = kappa
        self.gamma1_list = gamma1
        self.gamma2_list = gamma2
        self.F1_list = F1
        self.F2_list = F2

    def SIScorrSig(self, x_l='default', y_l='default', z_l=0.01, sigma_v=200):
        
        """Place a Singular Isothermal Sphere (SIS) at some redshift between the galaxies is this mock catalog
           and the observer (at z=0), the default is z=0.01.  For each source galaxy in the catalog, we first
           compute the Einstein radius, and then the lensing fields at the location of the source galaxies on the sky.
           The velocity dispersion default value is taken to be 200 km/s
        """
        
        # First, get location of the SIS
        if x_l == 'default':
          x_l =self.xlim/2
        if y_l == 'default':
          y_l =self.ylim/2
       
        
       # Get the angular diameter distance to the SIS
        cosmo = Planck15
        D_l = cosmo.angular_diameter_distance(z_l).value
        
        # Next get the distance between the source and the SIS
        D_s = self.D_s_list
        D_ls = D_s - D_l
        
        # Compute the Einstein radius
        theta_E = ((4*np.pi*(sigma_v/3.e5)**2.*D_ls/D_s)*u.rad).to(u.arcsec).value

        # Get the source galaxy's angular distance from the SIS lens center
        x_s = self.x_list
        y_s = self.y_list
        theta = np.sqrt((x_s-x_l)**2. + (y_s-y_l)**2.)
        
        # Get the position angle around the lens
        phi = np.arctan2((y_s-y_l),(x_s-x_l))

        # Get the lensing fields
        kappa = theta_E/(2*theta)
        gamma1 = -kappa*np.cos(2*phi)
        gamma2 = -kappa*np.sin(2*phi)
        F1 = -(theta_E)/(2*theta**2.)*np.cos(phi)
        F2 = -(theta_E)/(2*theta**2.)*np.sin(phi)

        self.kappa_list = kappa
        self.gamma1_list = gamma1
        self.gamma2_list = gamma2
        self.F1_list = F1
        self.F2_list = F2

    def sIScorrSig(self, x_l='default', y_l='default', z_l=0.01, sigma_v=200, theta_c=1):
        
        """Place a Softened Isothermal Sphere (sIS) at some redshift between the galaxies is this mock catalog
           and the observer (at z=0), the default is z=0.01.  For each source galaxy in the catalog, we first
           compute the Einstein radius, and then the lensing fields at the location of the source galaxies on the sky.
           The velocity dispersion default value is taken to be 200 km/s.  The default value for the core radius is taken
           to be 1 arcsec.
        """
        
        # First, get location of the SIS
        if x_l == 'default':
          x_l = self.xlim/2
        if y_l == 'default':
          y_l = self.ylim/2
        
        # First get the angular diameter distance to the SIS
        D_l = self.cosmo.angular_diameter_distance(z_l).value
        
        # Next get the distance between the source and the SIS
        D_s = self.D_s_list
        D_ls = D_s - D_l
        
        # Compute the Einstein radius
        theta_E = ((4*np.pi*(sigma_v/3.e5)**2.*D_ls/D_s)*u.rad).to(u.arcsec).value

        # Get the source galaxy's angular distance from the SIS lens center
        x_s = self.x_list
        y_s = self.y_list
        theta = np.sqrt((x_l-x_s)**2. + (y_l-y_s)**2.)
        
        # Get the position angle around the lens
        phi = np.arctan2((y_s-y_l),(x_s-x_l))

        # Get the lensing fields: convergence, shear, and flexion
        kappa = theta_E/(2*np.sqrt(theta**2.+theta_c**2.))
        gamma1 = -kappa*np.cos(2*phi)
        gamma2 = -kappa*np.sin(2*phi)
        F1 = -(theta_E)/(2*(theta**2.+theta_c**2.)**(3./2))*theta*np.cos(phi)
        F2 = -(theta_E)/(2*(theta**2.+theta_c**2.)**(3./2))*theta*np.sin(phi)

        self.kappa_list = kappa
        self.gamma1_list = gamma1
        self.gamma2_list = gamma2
        self.F1_list = F1
        self.F2_list = F2

    def exportCatalog(self, kappaFlag=False):
        
        """Export the mock catalog to a pickle file in the form of a pandas dataframe.
            The kappaFlag allows you to export the convergence along with lensing fields.
        """
        
        # First, get the observed ellipticity/shear and flexion by adding the intrinsic and lensing components
        gamma1_obs_list = self.eps1_s_list + self.gamma1_list
        gamma2_obs_list = self.eps2_s_list + self.gamma2_list
        F1_obs_list = self.F1_s_list + self.F1_list
        F2_obs_list = self.F2_s_list + self.F2_list
        
        # Let us create a pandas dataframe for this mock galaxy catalog
        if kappaFlag:
            
            col_list = ['x', 'y', 'z', 'a', 'shear1', 'shear2', 'F1', 'F2', 'kappa']
            arrs = [self.x_list, self.y_list, self.z_list, self.a_list, gamma1_obs_list, gamma2_obs_list, F1_obs_list, F2_obs_list, self.kappa_list]
        
        else:
            
            col_list = ['x', 'y', 'z', 'a', 'shear1', 'shear2', 'F1', 'F2']
            arrs = [self.x_list, self.y_list, self.z_list, self.a_list, gamma1_obs_list, gamma2_obs_list, F1_obs_list, F2_obs_list]
        
        dat = {i:arrs[j] for i,j in zip(col_list, range(len(col_list)))}
        out_frame = pd.DataFrame(data = dat, columns = col_list)
        out_frame.to_pickle(self.cat_path+self.cat_name+'.pkl')
