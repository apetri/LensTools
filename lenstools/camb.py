import numpy as np

#Options
transfer_high_precision = False
transfer_kmax = 250.0
transfer_k_per_logint = 0
transfer_redshift_min = 0
transfer_redshift_max = 2.0
transfer_num_redshifts = 20
transfer_interp_matterpower = True
camb_parameters_file = 'redshifts.txt'

#Write options on file
file=open(camb_parameters_file,'w')

if(transfer_high_precision):
	file.write('transfer_high_precision = T\n')
else:
	file.write('transfer_high_precision = F\n')

file.write('transfer_kmax = %f\n'%(transfer_kmax))
file.write('transfer_k_per_logint = %d\n'%(transfer_k_per_logint))
file.write('transfer_num_redshifts = %d\n'%(transfer_num_redshifts))

if(transfer_interp_matterpower):
	file.write('transfer_interp_matterpower = T\n')
else:
	file.write('transfer_interp_matterpower = F\n')

redshift = np.ogrid[transfer_redshift_max:transfer_redshift_min:transfer_num_redshifts*1j]

for i in range(transfer_num_redshifts):
	file.write('transfer_redshift(%d) = %f\n'%(i+1,redshift[i]))
	file.write('transfer_filename(%d) = transfer_out_%d.dat\n'%(i+1,int(redshift[i]*100)))

file.write('#Matter power spectrum output against k/h in units of h^{-3} Mpc^3\n')

for i in range(transfer_num_redshifts):
	file.write('transfer_matterpower(%d) = matterpower_%d.dat\n'%(i+1,int(redshift[i]*100)))

