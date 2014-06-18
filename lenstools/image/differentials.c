#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coordinates.h"

////////////Statistical operations on maps/////////////////////////

float average_map(float *map,long map_size){
	
	long k;
	float sum,average;
	
	sum=0.0;
	
	for(k=0;k<map_size*map_size;k++){
		
		sum+=map[k];
	}
	
	average=sum/(map_size*map_size);
	
	return average;
}

//subtract mean from the map
void average_subtract(float *map,long map_size){
	
	long k;
	float avg;
	
	avg=average_map(map,map_size);
	
	for(k=0;k<map_size*map_size;k++){
		
		map[k]=map[k]-avg;
	}

}

//variance of the map
float variance_map(float *map,long map_size){
	
	long k;
	float sum,sum2,average,average2,variance;
	
	sum=0.0;
	sum2=0.0;
	
	for(k=0;k<map_size*map_size;k++){
		
		sum+=map[k];
		sum2+=pow(map[k],2);
	}
	
	average=sum/(map_size*map_size);
	average2=sum2/(map_size*map_size);
	
	variance=average2-pow(average,2);
	
	return variance;
}

////////////Differential operations: finite difference/////////////

void gradient(float *map,float *grad_map,long map_size){
	
	long i,j;
	float grad_x,grad_y;
	
	for(i=0;i<map_size;i++){
		for(j=0;j<map_size;j++){
			
			grad_x=(map[coordinate(i+1,j,map_size)]-map[coordinate(i-1,j,map_size)])/2.0;
			grad_y=(map[coordinate(i,j+1,map_size)]-map[coordinate(i,j-1,map_size)])/2.0;
			
			grad_map[coordinate(i,j,map_size)]=sqrt(pow(grad_x,2)+pow(grad_y,2));
			
		}
	}
}

void gradient_xy(float *map,float *grad_map_x,float *grad_map_y,long map_size){
	
	long i,j;
	float grad_x,grad_y;
	
	for(i=0;i<map_size;i++){
		for(j=0;j<map_size;j++){
			
			
			grad_x=(map[coordinate(i+1,j,map_size)]-map[coordinate(i-1,j,map_size)])/2.0;
			grad_y=(map[coordinate(i,j+1,map_size)]-map[coordinate(i,j-1,map_size)])/2.0;
			
			grad_map_x[coordinate(i,j,map_size)]=grad_x;
			grad_map_y[coordinate(i,j,map_size)]=grad_y;
			
			
		}
	}
}

void hessian(float *map,float *hess_xx_map,float *hess_yy_map,float *hess_xy_map,long map_size){
	
	long i,j;
	float hessian_xx,hessian_yy,hessian_xy;
	
	for(i=0;i<map_size;i++){
		for(j=0;j<map_size;j++){
			
			hessian_xx=(map[coordinate(i+2,j,map_size)]+map[coordinate(i-2,j,map_size)]-2*map[coordinate(i,j,map_size)])/4.0;
			hessian_yy=(map[coordinate(i,j+2,map_size)]+map[coordinate(i,j-2,map_size)]-2*map[coordinate(i,j,map_size)])/4.0;
			hessian_xy=(map[coordinate(i+1,j+1,map_size)]+map[coordinate(i-1,j-1,map_size)]-map[coordinate(i-1,j+1,map_size)]-map[coordinate(i+1,j-1,map_size)])/4.0;
			
			hess_xx_map[coordinate(i,j,map_size)]=hessian_xx;
			hess_yy_map[coordinate(i,j,map_size)]=hessian_yy;
			hess_xy_map[coordinate(i,j,map_size)]=hessian_xy;
			
		}
	}
}