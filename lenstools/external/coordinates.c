int min_int(int i,int j){
	
	if(i<j){
		return i;
	}
	else{
		return j;
	}
}

long min_long(long i,long j){
	
	if(i<j){
		return i;
	}
	else{
		return j;
	}
}

/*From the linear coordinate of the array computes the pixel coordinate*/
long fourier_coordinate(long x,long y,long map_size){
	
	long c;
	
	c = ((map_size/2+1)*x + y);
	return c;
	
}

long coordinate(long x,long y,long map_size){
	
	long c;
	
	c = ((y+map_size) % map_size)*map_size + ((x+map_size) % map_size);
	return c;
	
}