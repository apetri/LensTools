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
void peak_count(double *map,long map_size, double sigma, int Nthresh, double *thresholds, double *peaks){
	
	int Nbins=Nthresh-1,i,j,k;
	
	for(i=0;i<map_size;i++){
		for(j=0;j<map_size;j++){
		
			if(is_peak(i,j,map_size,map)){
				
				for(k=0;k<Nbins;k++){
					
					if(map[coordinate(i,j,map_size)]>=thresholds[k]*sigma && map[coordinate(i,j,map_size)]<thresholds[k+1]*sigma){
						peaks[k]+= 1.0/(thresholds[k+1]-thresholds[k]);
					}
					
				}
				
			}
			   
		
		}
	}

}

		
			   