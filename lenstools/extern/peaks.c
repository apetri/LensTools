#include "coordinates.h"

//decide if a map pixel corresponds to a peak looking at the values of its nearest neighbors
int is_peak(int i,int j,long map_size,double *map){
	
	float pixel_value = map[coordinate(i,j,map_size)];
	int condition;

	condition = pixel_value>map[coordinate(i+1,j,map_size)] && pixel_value>map[coordinate(i-1,j,map_size)];
	condition = condition && pixel_value>map[coordinate(i,j+1,map_size)] && pixel_value>map[coordinate(i,j-1,map_size)];
	condition = condition && pixel_value>map[coordinate(i+1,j+1,map_size)] && pixel_value>map[coordinate(i-1,j-1,map_size)];
	condition = condition && pixel_value>map[coordinate(i-1,j+1,map_size)] && pixel_value>map[coordinate(i+1,j-1,map_size)];
	
	if(condition){
		
		return 1;
		
	}
	else{
		
		return 0;
		
	}
	
}

//count the peaks in the map for varying threshold
void peak_count(double *map,unsigned char *mask,long map_size, double sigma, int Nthresh, double *thresholds, double *peaks){
	
	int Nbins=Nthresh-1,i,j,k;
	long l;

	if(mask){

		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){

				l = coordinate(i,j,map_size);

				//If the pixel is masked, skip it
				if(!mask[l]){
					continue;
				}
		
				if(is_peak(i,j,map_size,map)){
				
					for(k=0;k<Nbins;k++){
					
						if(map[l]>=thresholds[k]*sigma && map[l]<thresholds[k+1]*sigma){
							peaks[k]+= 1.0/(thresholds[k+1]-thresholds[k]);
						}
					
					}
				
				}			   
		
			}
		}

	}
	

	else{
	
		//If there is no mask available don't bother check if the pixel is masked or not
		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){

				l = coordinate(i,j,map_size);
		
				if(is_peak(i,j,map_size,map)){
				
					for(k=0;k<Nbins;k++){
					
						if(map[l]>=thresholds[k]*sigma && map[l]<thresholds[k+1]*sigma){
							peaks[k]+= 1.0/(thresholds[k+1]-thresholds[k]);
						}
					
					}
				
				}
			   
		
			}
		}
	}
}

//locate the peaks on the map
int peak_locations(double *map,unsigned char *mask,long map_size, double sigma, int Nthresh, double *thresholds,double *values,int *locations_x,int *locations_y){

	int Nbins=Nthresh-1,i,j,loc;
	long l;

	loc=0;

	if(mask){

		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){

				l = coordinate(i,j,map_size);

				//If the pixel is masked, skip it
				if(!mask[l]){
					continue;
				}
		
				if(is_peak(i,j,map_size,map)){
					
					if(map[l]>=thresholds[0]*sigma && map[l]<thresholds[Nbins]*sigma){
						values[loc] = map[l] / sigma;
						locations_x[loc] = i;
						locations_y[loc] = j;
						loc++;
					}
					
				}			   
		
			}
		}

	} else{
	
		//If there is no mask available don't bother check if the pixel is masked or not
		for(i=0;i<map_size;i++){
			for(j=0;j<map_size;j++){

				l = coordinate(i,j,map_size);
		
				if(is_peak(i,j,map_size,map)){
					
					if(map[l]>=thresholds[0]*sigma && map[l]<thresholds[Nbins]*sigma){
						values[loc] = map[l] / sigma;
						locations_x[loc] = i;
						locations_y[loc] = j;
						loc++;
					}
					
					
				
				}
			   
		
			}
		}
	}


	return loc;

}

		
			   