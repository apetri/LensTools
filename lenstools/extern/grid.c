#include <stdlib.h>
#include <math.h>

#include "coordinates.h"
#include "grid.h"

#define WEIGHT_DEFAULT 1.0
#define CONCENTRATION_DEFAULT 1.0
#define NFW_CUT 0.1


//NFW density profile
double nfwKernel(double dsquared,double w,double rv,double c){

	double x = c * sqrt(dsquared) / rv;

	if(x<NFW_CUT){
		return w/(NFW_CUT*pow(1.0+NFW_CUT,2));
	} else{
		return w/(x*pow(1.0+x,2));
	}

}

//Two dimensional pixelization of a galaxy catalog
int grid2d(double *x,double *y,double *s,double *map,int Nobjects,int Npixel,double map_size){

	int ip,jp,p,n;
	double i,j;
	unsigned char *discovered;
	double resolution=map_size/Npixel;

	//First allocate an array that keeps track if a pixel is discovered (if there is at least one object that falls into it) or not
	if((discovered = (unsigned char *)malloc(sizeof(char)*Npixel*Npixel))==NULL){
		return 1;
	}

	//initialize to 0
	for(n=0;n<Npixel*Npixel;n++) discovered[n]=0;

	//cycle over objects and assign to the pixels in the grid
	for(n=0;n<Nobjects;n++){

		//Compute the position on the grid in the fastest way
		i = x[n] / resolution;
		j = y[n] / resolution;

		//If the particle lands on the grid, put it in the correct pixel
		if(i>=0 && i<Npixel && j>=0 && j<Npixel){

			ip = (int)i;
			jp = (int)j;
			p = ip*Npixel + jp;

			if(!discovered[p]){
				discovered[p]=1;
				map[p] = 0.0;
			}

			map[p] += s[n];

		}

	}

	//Release allocated memory and return to user
	free(discovered);
	return 0;
}


//Snap particles on a 3d regularly spaced grid
int grid3d(float *positions,float *weights,double *radius,double *concentration,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid,double(*kernel)(double,double,double,double)){

	int n;
	double i,j,k;

	//Cycle through the particles and for each one compute the position on the grid
	if(radius==NULL || concentration==NULL || kernel==NULL){
	
		if(weights==NULL){

			for(n=0;n<Npart;n++){

				//Compute the position on the grid in the fastest way
				i = (positions[3*n] - leftX)/sizeX;
				j = (positions[3*n + 1] - leftY)/sizeY;
				k = (positions[3*n + 2] - leftZ)/sizeZ;

				//If the particle lands on the grid, put it in the correct pixel
				if(i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz){

					grid[((int)i)*ny*nz + ((int)j)*nz + (int)k] += 1.0;

				}
			
			}

		} else{

			for(n=0;n<Npart;n++){

				//Compute the position on the grid in the fastest way
				i = (positions[3*n] - leftX)/sizeX;
				j = (positions[3*n + 1] - leftY)/sizeY;
				k = (positions[3*n + 2] - leftZ)/sizeZ;

				//If the particle lands on the grid, put it in the correct pixel
				if(i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz){

					grid[((int)i)*ny*nz + ((int)j)*nz + (int)k] += weights[n];

				}
			
			}

		}

	

	} else{

		int minI,maxI,minJ,maxJ,minK,maxK;
		int ii,jj,kk;
		double distanceSquared;

		for(n=0;n<Npart;n++){

			//Compute the position of the center on the grid in the fastest way
			i = (positions[3*n] - leftX)/sizeX;
			j = (positions[3*n + 1] - leftY)/sizeY;
			k = (positions[3*n + 2] - leftZ)/sizeZ;

			//Compute the extremes of the cloud corresponding to the particle
			minI = max_int((int)(i-radius[n]/sizeX),0);
			minI = min_int(minI,nx);
			maxI = min_int((int)(i+radius[n]/sizeX),nx);
			maxI = max_int(maxI,0);

			minJ = max_int((int)(j-radius[n]/sizeY),0);
			minJ = min_int(minJ,ny);
			maxJ = min_int((int)(j+radius[n]/sizeY),ny);
			maxJ = max_int(maxJ,0);

			minK = max_int((int)(k-radius[n]/sizeZ),0);
			minK = min_int(minK,nz);
			maxK = min_int((int)(k+radius[n]/sizeZ),nz);
			maxK = max_int(maxK,0);

			//Cycle over all the cloud pixels covered by the particle and assign the correct density
			for(ii=minI;ii<maxI;ii++){
				for(jj=minJ;jj<maxJ;jj++){
					for(kk=minK;kk<maxK;kk++){

						distanceSquared = pow(leftX + (ii+0.5)*sizeX - positions[3*n],2) + pow(leftY + (jj+0.5)*sizeY - positions[3*n+1],2) + pow(leftZ + (kk+0.5)*sizeZ - positions[3*n+2],2);
						if(distanceSquared<pow(radius[n],2)) grid[ii*ny*nz + jj*nz + kk] += (float)kernel(distanceSquared,(double)weights[n],radius[n],concentration[n]); 

					}
				}
			}



		}

	}

	return 0;


}


//adaptive smoothing
int adaptiveSmoothing(int NumPart,float *positions,float *weights,double *rp,double *concentration,double *binning0, double *binning1,double center,int direction0,int direction1,int normal,int size0,int size1,int projectAll,double *lensingPlane,double(*kernel)(double,double,double,double)){

	int i,j,p;
	float posNormal,posTransverse0,posTransverse1;
	double catchmentRadius,distanceSquared,w,c;
	int catchmentRadiusPixel,pos0Pixel,pos1Pixel,pixelLeft0,pixelRight0,pixelLeft1,pixelRight1;

	//Loop over particles
	for(p=0;p<NumPart;p++){

		//Set particle weight
		if(weights){
			w = (double)(weights[p]);
		} else{
			w = WEIGHT_DEFAULT;
		}

		//Set particle concentration
		if(concentration){
			c = concentration[p];
		} else{
			c = CONCENTRATION_DEFAULT;
		}


		//Compute transverse and longitudinal positions with respect to the plane
		posNormal = positions[3*p + normal];
		posTransverse0 = positions[3*p + direction0];
		posTransverse1 = positions[3*p + direction1];

		//If we don't want to collapse all the snapshot, and if the particle is too far, skip to the next
		if((!projectAll) && (fabs(posNormal-center)>rp[p])) continue;

		//Compute catchment radius
		if(projectAll){
			catchmentRadius = rp[p];
		} else{
			catchmentRadius = sqrt(pow(rp[p],2) - pow(posNormal-center,2));
		}
		
		catchmentRadiusPixel = (int)(catchmentRadius / (binning0[1] - binning0[0]));

		//Compute pixel extremes on the plane, enforcing the bounds
		pos0Pixel = (int)((posTransverse0 - binning0[0]) / (binning0[1] - binning0[0]));
		pixelLeft0 = max_int(pos0Pixel - catchmentRadiusPixel,0);
		pixelLeft0 = min_int(pixelLeft0,size0);
		pixelRight0 = min_int(pos0Pixel + catchmentRadiusPixel,size0);
		pixelRight0 = max_int(pixelRight0,0);

		pos1Pixel = (int)((posTransverse1 - binning1[0]) / (binning1[1] - binning1[0]));
		pixelLeft1 = max_int(pos1Pixel - catchmentRadiusPixel,0);
		pixelLeft1 = min_int(pixelLeft1,size1);
		pixelRight1 = min_int(pos1Pixel + catchmentRadiusPixel,size1);
		pixelRight1 = max_int(pixelRight1,0);

		//Snap the particles on the plane
		for(i=pixelLeft0;i<pixelRight0;i++){
			for(j=pixelLeft1;j<pixelRight1;j++){

				//Compute the distance between this plane pixel and the particle
				if(projectAll){
					distanceSquared = pow(posTransverse0 - 0.5*(binning0[i]+binning0[i+1]),2) + pow(posTransverse1 - 0.5*(binning1[j]+binning1[j+1]),2);
				} else{	
					distanceSquared = pow(posNormal-center,2) + pow(posTransverse0 - 0.5*(binning0[i]+binning0[i+1]),2) + pow(posTransverse1 - 0.5*(binning1[j]+binning1[j+1]),2);
				}

				//Add the corresponding contribution to the density
				if(distanceSquared<pow(rp[p],2)) lensingPlane[i*size0 + j] += kernel(distanceSquared,w,rp[p],c); 

			}
		}


	}

	return 0;

}