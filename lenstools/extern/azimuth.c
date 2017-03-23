#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "coordinates.h"

/*Compute power spectral azimuthal averages of 2D real Fourier transforms of images*/
int azimuthal_rfft2(double _Complex *ft_map1,double _Complex *ft_map2,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *power_l,double *scale){

	//Define the pixel physical size in fourier space
	const double lpix = 360.0/map_angle_degrees;
	double lx,ly,l;

	//Take care of Fourier transforms normalization
	const double normalization = pow((map_angle_degrees * M_PI/180.0)/(size_x*size_x),2); 

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


	if(scale){

		///////////////////////////
		/*Apply scaling to pixels*/
		///////////////////////////

		//Cycle over pixels
		for(i=0;i<size_x;i++){

			lx = min_long(i,size_x-i) * lpix;

			for(j=0;j<size_y;j++){

				ly = j*lpix;
				l = sqrt(lx*lx + ly*ly);

				pixid = fourier_coordinate(i,j,size_x);

				//decide in which l bin this pixel falls into
				for(k=0;k<Nbins;k++){
				
					if(l>lvalues[k] && l<=lvalues[k+1]){

						power_l[k] += (creal(ft_map1[pixid])*creal(ft_map2[pixid]) + cimag(ft_map1[pixid])*cimag(ft_map2[pixid]))*scale[pixid];
						hits[k]++; 

					}

				}


			}

		}


	} else{

		//////////////////////////////////
		/*Do not apply scaling to pixels*/
		//////////////////////////////////

		//Cycle over pixels
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

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*General bispectrum calculator*/
int bispectrum(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *bispectrum_l,int (*k1tok2)(int,int,int*,int*,void*),void *args){

	//Define the pixel physical size in fourier space
	double lpix = 360.0/map_angle_degrees;
	double l;

	//Bispectrum normalization
	const double normalization = pow((map_angle_degrees*M_PI/180.0),4) / pow(size_x,6);

	//Placeholders
	double _Complex b1,b2,b3;
	int conjugate;

	//Allocate space for hits map, initialize to 0
	int *hits,k;
	int Nbins = Nvalues - 1;

	hits = (int *)malloc(sizeof(int)*Nbins);
	if(hits==NULL) return 1;

	for(k=0;k<Nbins;k++){
		hits[k] = 0;
	}

	//Multipoles
	int i,j,kx1,kx2,kx3,ky1,ky2,ky3;
	long n1,n2,n3;

	//Cycle over pixels in Fourier map
	for(i=0;i<size_x;i++){

		//Calculate integer wavenumber kx according to complex FFT frequencies
		if(i<size_x>>1){
			kx1 = i;
		} else{
			kx1 = i - size_x;
		}

		for(j=0;j<size_y;j++){

			//Compute array location of first leg
			n1 = size_y*i + j;
			b1 = ft_map1[n1];

			//Calculate integer wavenumber ky according to real FFT frequencies
			ky1 = j;

			//Calculate multipole
			l = sqrt(kx1*kx1 + ky1*ky1)*lpix;
			if(l>=lvalues[Nbins]){
				continue;
			}

			//Calcuate second leg
			if(k1tok2(kx1,ky1,&kx2,&ky2,args)){
				continue;
			}

			//Third leg closes the triangle
			kx3 = -kx1-kx2;
			ky3 = -ky1-ky2;

			//Compute locations in FFT arrays for 2nd and 3rd leg
			if(ky2<0){
				kx2 = -kx2;
				ky2 = -ky2;
				conjugate = 1;
			} else{
				conjugate = 0;
			}

			//Bounds check
			if(ky2>=size_y || kx2<=-size_y || kx2>=size_y){
				continue;
			}

			if(kx2<0){
				kx2 = size_x + kx2;
			}

			n2 = size_y*kx2 + ky2;
			b2 = ft_map2[n2];
			if(conjugate){
				b2 = conj(b2);
			}

			if(ky3<0){
				kx3 = -kx3;
				ky3 = -ky3;
				conjugate = 1;
			} else{
				conjugate = 0;
			}

			//Bounds check
			if(ky3>=size_y || kx3<=-size_y || kx3>=size_y){
				continue;
			}

			if(kx3<0){
				kx3 = size_x + kx3;
			} 

			n3 = size_y*kx3 + ky3;
			b3 = ft_map3[n3];
			if(conjugate){
				b3 = conj(b3);
			}


			//Find the correct bin
			for(k=0;k<Nbins;k++){
				
				if(l>lvalues[k] && l<=lvalues[k+1]){

					hits[k]++;
					bispectrum_l[k] += creal(b1)*creal(b2)*creal(b3) - creal(b1)*cimag(b2)*cimag(b3) - cimag(b1)*creal(b2)*cimag(b3) - cimag(b1)*cimag(b2)*creal(b3) ;


				}

			}

		
		}
	
	
	}

	//Normalize result
	for(k=0;k<Nbins;k++){
		if(hits[k]>0){
			bispectrum_l[k] = normalization * bispectrum_l[k]/hits[k];
		}
	}

	//Free hits map
	free(hits);

	//Return, no problem
	return 0;

}

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

/*Equilateral rotator*/
int k1tok2_equilateral(int kx1,int ky1,int *kx2,int *ky2,void *args){

	//sin, cos, tan of 120deg
	const double S120 = 0.8660254037844386;
	const double C120 = -0.5;
	const double T120 = -1.7320508075688772;

	//Skip to next if tan(angle)<tan(120)
	if(kx1<0 && ky1<kx1*T120){
		return 1;
	}

	*kx2 = (int)(C120*kx1 - S120*ky1);
	*ky2 = (int)(S120*kx1 + C120*ky1);

	return 0;

}

/*Equilateral bispectrum*/
int bispectrum_equilateral(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,int Nvalues,double *lvalues,double *bispectrum_l){

	return bispectrum(ft_map1,ft_map2,ft_map3,size_x,size_y,map_angle_degrees,Nvalues,lvalues,bispectrum_l,k1tok2_equilateral,NULL);

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*Folded rotator*/
int k1tok2_folded(int kx1,int ky1,int *kx2,int *ky2,void *args){

	double folding_ratio = *(double *)args;

	*kx2 = -(int)(folding_ratio*kx1);
	*ky2 = -(int)(folding_ratio*ky1);

	return 0;

}

/*Folded bispectrum*/
int bispectrum_folded(double _Complex *ft_map1,double _Complex *ft_map2,double _Complex *ft_map3,long size_x,long size_y,double map_angle_degrees,double folding_ratio,int Nvalues,double *lvalues,double *bispectrum_l){

	return bispectrum(ft_map1,ft_map2,ft_map3,size_x,size_y,map_angle_degrees,Nvalues,lvalues,bispectrum_l,k1tok2_folded,&folding_ratio);

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

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