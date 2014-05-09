import numpy as np
from scipy import interpolate,integrate
from astropy.cosmology import wCDM 

##################################################
#############Options##############################
##################################################

#Cosmological parameters
h = 0.7
Omega_m = 0.26
w = -1.0
c = 3.0e5

#Binning options
min_redshift = 0.0
max_redshift = 2.0
num_redshifts = 20
l_min = 1
l_max = 100000

#Power spectra filename prefix
power_prefix = 'camb_output/highsi8_matterpower_'

#Convergence power spectrum file
conv_file_name = 'linear_spectra/highsi8_conv_power.txt'

##################################################
#############End of options, begin code###########
##################################################

#Power spectrum normalization
normalization = (9.0/4)*(Omega_m)**2*(100*h/c)**4

#Import cosmological model from astropy
cosmo = wCDM(H0=100*h, Om0=Omega_m, Ode0=1.0-Omega_m, w0=w)

#l bins
l = np.ogrid[l_min:l_max:1]
C = np.zeros(len(l))

#Compute comoving distances
redshift = np.ogrid[min_redshift:max_redshift:num_redshifts*1j]
chi = cosmo.comoving_distance(redshift)
chi0 = chi[num_redshifts-1]
kernel = (1.0 - chi/chi0)**2

#Load matter power spectra from camb output files#
#See how many k's are stored
kappa,try_power = (np.loadtxt(power_prefix + '0.dat')).transpose()
num_kappa = len(kappa)

#Load power spectrum
power_spectrum = np.zeros([num_kappa,num_redshifts])
for i in range(num_redshifts):
	try_power = np.loadtxt(power_prefix + ('%d.dat'%(int(redshift[i]*100))))
	power_spectrum[:,i] = try_power[:,1] / (h**3)

#Compute the integral for lensing power spectrum#
power_interpolation = interpolate.interp1d(kappa,power_spectrum,axis=0)

conv_file = open(conv_file_name,'w')

for j in range(len(l)):
	if(l[j]%1000==0):
		print l[j]/1000
	power_integrand = np.zeros(num_redshifts)
	power_integrand[1:] = power_interpolation(l[j]/chi[1:]).diagonal()
	full_integrand = kernel * (1.0 + redshift)**2 * power_integrand
	#Finally compute the integral
	C[j] = integrate.simps(full_integrand,chi) * normalization
	#Write on the file
	conv_file.write('%d %e\n'%(l[j],C[j]))

conv_file.close()













