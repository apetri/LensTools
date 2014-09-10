#include <stdio.h>
#include "read_gadget.h"

/* this routine loads particle positions and velocities from Gadget's default
binary file format*/
int getPosVel(FILE *fp,long offset,float *data,int Npart){

	int n;

	/*First offset the file pointer to go to where the first particle is*/
	if(fseek(fp,offset,SEEK_SET)) return -1;

	/*Next read the particle positions or velocities in the data array*/
	for(n=0;n<Npart;n++){
		if(fread(&data[3*n],sizeof(float),3,fp)!=3) return -1;
	}

	return 0;

}