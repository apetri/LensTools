#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "coordinates.h"

/*Compute power spectral azimuthal averages of 2D real Fourier transforms of images*/
int azimuthal_rfft2(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *power_l){

	//Define the pixel physical size in fourier space
	double lpix = 360.0/map_angle_degrees;
	double lx,ly,l;

	//Take care of Fourier transforms normalization
	double normalization = pow((map_angle_degrees * M_PI/180.0)/(size_x*size_x),2); 

	//Counters
	long i,j,pixid;

	//Binning
	int Nbins = Nvalues - 1;
	int k, *hits;

	//Allocate memory for hits counter
	hits = (int *)malloc(sizeof(int)*Nbins);
	if(hits==NULL){
		return 1;
	}

	//Initialize hits to 0
	for(k=0;k<Nbins;k++){
		hits[k] = 0;
	}

	//Select bins
	for(i=0;i<size_x;i++){

		lx = min_long(i,size_x-i) * lpix;

		for(j=0;j<size_y;j++){

			ly = j*lpix;
			l = sqrt(lx*lx + ly*ly);

			pixid = fourier_coordinate(i,j,size_x);

			//decide in which l bin this pixel falls into
			for(k=0;k<Nbins;k++){
				
				if(l>lvalues[k] && l<=lvalues[k+1]){

					power_l[k] += creal(ft_map1[pixid])*creal(ft_map2[pixid]) + cimag(ft_map1[pixid])*cimag(ft_map2[pixid]);
					hits[k]++; 

				}

			}


		}

	}

	//Compute average
	for(k=0;k<Nbins;k++){
		if(hits[k]>0){
			power_l[k] = normalization * power_l[k]/hits[k];
		}
	}

	//Free allocated memory
	free(hits);

	return 0;

}

/*Compute power spectral azimuthal averages of 3D real Fourier transforms of scalar fields (i.e. density)*/
int azimuthal_rfft3(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,long size_z,double kpixX,double kpixY,double kpixZ,int Nvalues,double *kvalues,double *power_k,long *hits){

	//Counters
	long x,y,z,p;
	int b;
	double k,kx,ky,kz;

	//Binning
	int Nbins = Nvalues - 1;

	//Loop over all the pixels in the fourier map
	for(x=0;x<size_x;x++){
		for(y=0;y<size_y;y++){
			for(z=0;z<size_z;z++){

				kx = min_long(x,size_x-x) * kpixX;
				ky = min_long(y,size_y-y) * kpixY;
				kz = z * kpixZ;
				k = sqrt(kx*kx + ky*ky + kz*kz);
				p = x*size_y*size_z + y*size_z + z;

				//Decide in which bin this pixel falls into
				for(b=0;b<Nbins;b++){
					if(k>kvalues[b] && k<=kvalues[b+1]){
						power_k[b] += creal(ft_map1[p])*creal(ft_map2[p]) + cimag(ft_map1[p])*cimag(ft_map2[p]);
						hits[b]++;
					}
				}



			}
		}
	}


	//Return
	return 0;


}


/* Compute angular averaged bispectrum of 3D real Fourier transforms of scalar fields (i.e. density)*/
int azimuthal_bispectrum_rfft3(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,long size_z,int size_small_x, int size_small_y, int size_small_z, double kpixX,double kpixY,double kpixZ,int Nsmall,int Nlarge, double *ksmall_edges, double *klarge_edges ,double *bispectrum_k,long *hits){
	// ft_map1, ft_map2, ft_map3 : input spectra
	// size_x, size_y, size_z : the grid size of ft_map
	// size_small_x, size_small_y, size_small_z : the maximum grid points that belong to ksmall_edges array; needed to reduced calculation time
	// kpixX, kpixY, kpixZ = 2 Pi / box_size : minimum momentum
	// Nsmall, Nlarge : the length of ksmall_edges or klarge_edges correspondingly
	// ksmall_edges, klarge_edges : small or large momentum at which bispectrum are computed
	// bispectrum_k, hits : results


	//Counters
	long x1,y1,z1,x2,y2,z2,x3,y3,z3,p1,p2,p3;
	int b_large,b_small;
	double k1,k1x,k1y,k1z,k2,k2x,k2y,k2z,k3,k3x,k3y,k3z;
	int odd_or_even;


	//Binning
	int largebins = Nlarge - 1;
	int smallbins = Nsmall - 1;

	if(size_x % 2 == 0){
		odd_or_even = 1; //even
	}
	else{
		odd_or_even = 0; //odd
	}

	//Loop over all the pixels in the fourier map
	for(x1=0;x1<size_small_x;x1++){
		for(y1=0;y1<size_small_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for(y1=size_y-size_small_y+1;y1<size_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	for(x1=size_x - size_small_x + 1;x1<size_x;x1++){
		for(y1=0;y1<size_small_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		for(y1=size_y-size_small_y+1;y1<size_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													bispectrum_k[b_small*largebins + b_large] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
													hits[b_small*largebins + b_large]++;
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}	
	//Return
	return 0;
}

/* Compute bispectrum of 3D real Fourier transforms of scalar fields (i.e. density)*/
int bispectrum_rfft3(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y, long size_z, int size_small_x, int size_small_y, int size_small_z, double kpixX,double kpixY,double kpixZ,int Nsmall, int Nlarge, int Ncosine, double *ksmall_edges, double *klarge_edges, double *cosines, double *bispectrum_k,long *hits){
	// ft_map1, ft_map2, ft_map3 : input spectra
	// size_x, size_y, size_z : the grid size of ft_map
	// size_small_x, size_small_y, size_small_z : the maximum grid points that belong to ksmall_edges array; needed to reduced calculation time
	// kpixX, kpixY, kpixZ = 2 Pi / box_size : minimum momentum
	// Nsmall, Nlarge, Ncosine : the length of ksmall_edges, klarge_edges or cosines correspondingly
	// ksmall_edges, klarge_edges : small or large momentum at which bispectrum are computed
	// cosines : cosine of angle between momenta at which bispectrum are computed
	// bispectrum_k, hits : results


	//Counters
	long x1,y1,z1,x2,y2,z2,x3,y3,z3,p1,p2,p3;
	int b_large,b_small;
	long b_cosine;
	double k1,k1x,k1y,k1z,k2,k2x,k2y,k2z,k3,k3x,k3y,k3z,cos_theta;
	int odd_or_even;

	//Binning
	int largebins = Nlarge - 1;
	int smallbins = Nsmall - 1;
	int cosinebins = Ncosine - 1;

	if(size_x % 2 == 0){
		odd_or_even = 1; //even
	}
	else{
		odd_or_even = 0; //odd
	}
	//Loop over all the pixels in the fourier map
	for(x1=0;x1<size_small_x;x1++){
		for(y1=0;y1<size_small_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}	
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for(y1=size_y - size_small_y+1;y1<size_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}	
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	for(x1=size_x - size_small_x+1;x1<size_x;x1++){
		for(y1=0;y1<size_small_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}	
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}

		for(y1=size_y - size_small_y+1;y1<size_y;y1++){
			for(z1=0;z1<size_small_z;z1++){

				// k1 is small momentum
				k1x = signed_long(x1,size_x-x1) * kpixX;
				k1y = signed_long(y1,size_y-y1) * kpixY;
				k1z = z1 * kpixZ;
				k1 = sqrt(k1x*k1x + k1y*k1y + k1z*k1z);
				p1 = x1*size_y*size_z + y1*size_z + z1;

				//Decide in which bin_small this small momentum pixel falls into
				for(b_small=0;b_small<smallbins;b_small++){
					if(k1>ksmall_edges[b_small] && k1<=ksmall_edges[b_small+1]){
						for(x2=0;x2<size_x;x2++){
							for(y2=0;y2<size_y;y2++){
								for(z2=0;z2<size_z;z2++){

									// k2 is large momentum
									k2x = signed_long(x2,size_x-x2) * kpixX;
									k2y = signed_long(y2,size_y-y2) * kpixY;
									k2z = z2 * kpixZ;
									k2 = sqrt(k2x*k2x + k2y*k2y + k2z*k2z);
									p2 = x2*size_y*size_z + y2*size_z + z2;
									// k3 = k1 + k2
									k3x = k1x+k2x;
									k3y = k1y+k2y;
									k3z = k1z+k2z;
									k3 = sqrt(k3x*k3x + k3y*k3y + k3z*k3z);

									if(odd_or_even == 1){
										if(k3x>=-kpixX*(size_x-2)/2 && k3x <= kpixX*size_x/2 && k3y>=-kpixY*(size_y-2)/2 && k3y<=kpixY*size_y/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}
												}
											}
										}
									}
									else{
										if(k3x>=-kpixX*(size_x-1)/2 && k3x <= kpixX*(size_x-1)/2 && k3y>=-kpixY*(size_y-1)/2 && k3y<=kpixY*(size_y-1)/2 && k3z<=kpixZ*(size_z-1)){
											if(k3x>=0){
												x3 = (long)(k3x/kpixX + 0.001); // Add or subtract 0.001 to avoid error due to finite precision of real numbers
											}
											else{
												x3 = size_x + (long)(k3x/kpixX - 0.001);
											}
											if(k3y>=0){
												y3 = (long)(k3y/kpixY + 0.001);
											}
											else{
												y3 = size_y + (long)(k3y/kpixY - 0.001);
											}
											z3 = (long)(k3z/kpixZ + 0.001);
											p3 = x3*size_y*size_z + y3*size_z + z3;
											for(b_large=0;b_large<largebins;b_large++){
												if(k2>klarge_edges[b_large] && k2<=klarge_edges[b_large+1]){
													cos_theta = (k1x*k2x + k1y*k2y + k1z*k2z)/(k1*k2);
													for(b_cosine=0;b_cosine<cosinebins;b_cosine++){
														if(cos_theta>cosines[b_cosine] && cos_theta<=cosines[b_cosine+1]){
															bispectrum_k[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine] += creal(ft_map1[p1]*ft_map2[p2]*conj(ft_map3[p3]));
															hits[b_small*largebins*cosinebins + b_large * cosinebins + b_cosine]++;
														}
													}	
												}
											}							
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//Return
	return 0;
}