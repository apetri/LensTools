#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coordinates.h"

void gradient_xy(double *map,double *grad_map_x,double *grad_map_y,long map_size,int Npoints,int *x_points,int *y_points){
	
	long i,j;
	double grad_x,grad_y;

	if(Npoints<0){
	
		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){
			
			
				grad_x=(map[coordinate(i+1,j,map_size)]-map[coordinate(i-1,j,map_size)])/2.0;
				grad_y=(map[coordinate(i,j+1,map_size)]-map[coordinate(i,j-1,map_size)])/2.0;
			
				grad_map_x[coordinate(i,j,map_size)]=grad_x;
				grad_map_y[coordinate(i,j,map_size)]=grad_y;
			
			
			}
		}
	

	} else{

		for(i=0;i<Npoints;i++){

			grad_x=(map[coordinate(x_points[i]+1,y_points[i],map_size)]-map[coordinate(x_points[i]-1,y_points[i],map_size)])/2.0;
			grad_y=(map[coordinate(x_points[i],y_points[i]+1,map_size)]-map[coordinate(x_points[i],y_points[i]-1,map_size)])/2.0;
			
			grad_map_x[i]=grad_x;
			grad_map_y[i]=grad_y;

		}


	}
}

void hessian(double *map,double *hess_xx_map,double *hess_yy_map,double *hess_xy_map,long map_size,int Npoints, int *x_points,int *y_points){
	
	long i,j;
	double hessian_xx,hessian_yy,hessian_xy;
	

	if(Npoints<0){

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

	} else{


			for(i=0;i<Npoints;i++){
				
				hessian_xx=(map[coordinate(x_points[i]+2,y_points[i],map_size)]+map[coordinate(x_points[i]-2,y_points[i],map_size)]-2*map[coordinate(x_points[i],y_points[i],map_size)])/4.0;
				hessian_yy=(map[coordinate(x_points[i],y_points[i]+2,map_size)]+map[coordinate(x_points[i],y_points[i]-2,map_size)]-2*map[coordinate(x_points[i],y_points[i],map_size)])/4.0;
				hessian_xy=(map[coordinate(x_points[i]+1,y_points[i]+1,map_size)]+map[coordinate(x_points[i]-1,y_points[i]-1,map_size)]-map[coordinate(x_points[i]-1,y_points[i]+1,map_size)]-map[coordinate(x_points[i]+1,y_points[i]-1,map_size)])/4.0;
			
				hess_xx_map[i]=hessian_xx;
				hess_yy_map[i]=hessian_yy;
				hess_xy_map[i]=hessian_xy;
			
			}

	}


}


void gradLaplacian(double *map,double *grad_map_x,double *grad_map_y,long map_size,int Npoints,int *x_points,int *y_points){
	
	long i,j;
	double grad_x,grad_y;

	if(Npoints<0){
	
		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){
			
			
				grad_x = (map[coordinate(i+3,j,map_size)] + map[coordinate(i-1,j,map_size)] + map[coordinate(i+1,j+2,map_size)] + map[coordinate(i+1,j-2,map_size)] - 4*map[coordinate(i+1,j,map_size)])/8.0;
				grad_x -= (map[coordinate(i+1,j,map_size)] + map[coordinate(i-3,j,map_size)] + map[coordinate(i-1,j+2,map_size)] + map[coordinate(i-1,j-2,map_size)] - 4*map[coordinate(i-1,j,map_size)])/8.0;
				grad_y = (map[coordinate(i+2,j+1,map_size)] + map[coordinate(i-2,j+1,map_size)] + map[coordinate(i,j+3,map_size)] + map[coordinate(i,j-1,map_size)] - 4*map[coordinate(i,j+1,map_size)])/8.0;
				grad_y -= (map[coordinate(i+2,j-1,map_size)] + map[coordinate(i-2,j-1,map_size)] + map[coordinate(i,j+1,map_size)] + map[coordinate(i,j-3,map_size)] - 4*map[coordinate(i,j-1,map_size)])/8.0;
			
				grad_map_x[coordinate(i,j,map_size)]=grad_x;
				grad_map_y[coordinate(i,j,map_size)]=grad_y;
			
			
			}
		}
	

	} else{

		for(i=0;i<Npoints;i++){

			grad_x = (map[coordinate(x_points[i]+3,y_points[i],map_size)] + map[coordinate(x_points[i]-1,y_points[i],map_size)] + map[coordinate(x_points[i]+1,y_points[i]+2,map_size)] + map[coordinate(x_points[i]+1,y_points[i]-2,map_size)] - 4*map[coordinate(x_points[i]+1,y_points[i],map_size)])/8.0;
			grad_x -= (map[coordinate(x_points[i]+1,y_points[i],map_size)] + map[coordinate(x_points[i]-3,y_points[i],map_size)] + map[coordinate(x_points[i]-1,y_points[i]+2,map_size)] + map[coordinate(x_points[i]-1,y_points[i]-2,map_size)] - 4*map[coordinate(x_points[i]-1,y_points[i],map_size)])/8.0;
			grad_y = (map[coordinate(x_points[i]+2,y_points[i]+1,map_size)] + map[coordinate(x_points[i]-2,y_points[i]+1,map_size)] + map[coordinate(x_points[i],y_points[i]+3,map_size)] + map[coordinate(x_points[i],y_points[i]-1,map_size)] - 4*map[coordinate(x_points[i],y_points[i]+1,map_size)])/8.0;
			grad_y -= (map[coordinate(x_points[i]+2,y_points[i]-1,map_size)] + map[coordinate(x_points[i]-2,y_points[i]-1,map_size)] + map[coordinate(x_points[i],y_points[i]+1,map_size)] + map[coordinate(x_points[i],y_points[i]-3,map_size)] - 4*map[coordinate(x_points[i],y_points[i]-1,map_size)])/8.0;
			
			grad_map_x[i]=grad_x;
			grad_map_y[i]=grad_y;

		}


	}
}