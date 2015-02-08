#include <math.h>

#include "coordinates.h"
#include "grid.h"


//Snap particles on a 3d regularly spaced grid
int grid3d(float *positions,int Npart,double leftX,double leftY,double leftZ,double sizeX,double sizeY,double sizeZ,int nx,int ny,int nz,float *grid){

	int n;
	double i,j,k;

	//Cycle through the particles and for each one compute the position on the grid
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


	return 0;


}


//adaptive smoothing
int adaptiveSmoothing(int NumPart,float *positions,double *rp,double *binning0, double *binning1,double center,int direction0,int direction1,int normal,int size0,int size1,int projectAll,double *lensingPlane){

	int i,j,p;
	float posNormal,posTransverse0,posTransverse1;
	double catchmentRadius,distanceSquared;
	int catchmentRadiusPixel,pos0Pixel,pos1Pixel,pixelLeft0,pixelRight0,pixelLeft1,pixelRight1;

	//Loop over particles
	for(p=0;p<NumPart;p++){

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
		pixelLeft0 = max_int(pos0Pixel - catchmentRadiusPixel/2,0);
		pixelLeft0 = min_int(pixelLeft0,size0);
		pixelRight0 = min_int(pos0Pixel + catchmentRadiusPixel/2,size0);
		pixelRight0 = max_int(pixelRight0,0);

		pos1Pixel = (int)((posTransverse1 - binning1[0]) / (binning1[1] - binning1[0]));
		pixelLeft1 = max_int(pos1Pixel - catchmentRadiusPixel/2,0);
		pixelLeft1 = min_int(pixelLeft1,size1);
		pixelRight1 = min_int(pos1Pixel + catchmentRadiusPixel/2,size1);
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
				if(distanceSquared<pow(rp[p],2)) lensingPlane[i*size0 + j] += quadraticKernel(distanceSquared,rp[p]); 

			}
		}


	}

	return 0;

}