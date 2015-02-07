#ifndef __COORDINATES_H
#define __COORDINATES_H

int min_int(int,int);
int max_int(int,int);
long min_long(long,long);

static inline long coordinate(long x,long y,long map_size){
	return ((y+map_size) % map_size)*map_size + ((x+map_size) % map_size);
}

static inline long fourier_coordinate(long x,long y,long map_size){
	return ((map_size/2+1)*x + y);
}

#endif