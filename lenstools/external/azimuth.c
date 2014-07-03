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