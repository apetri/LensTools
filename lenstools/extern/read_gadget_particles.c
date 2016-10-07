#include <stdio.h>
#include <unistd.h>

#include "gadget2.h"

/*this routine loads particle positions and velocities from Gadget's default
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

int getPosVelFD(int fd,long offset,float *data,int Npart){

	int n;

	/*First offset the file pointer to go to where the first particle is*/
	if(lseek(fd,offset,SEEK_SET)!=offset) return -1;

	/*Read the particle positions or velocities in the data array*/
	for(n=0;n<Npart;n++){
		if(read(fd,&data[3*n],sizeof(float)*3)!=sizeof(float)*3) return -1;
	}

	return 0;

}

/*this routine loads in particle IDs (4 byte ints) from Gadget's default binary file format*/
int getID(FILE *fp,long offset,int *data,int Npart){

	int n;

	/*First offset the file pointer to go to where the first particle is*/
	if(fseek(fp,offset,SEEK_SET)) return -1;

	/*Next read the particle positions or velocities in the data array*/
	for(n=0;n<Npart;n++){
		if(fread(&data[n],sizeof(int),1,fp)!=1) return -1;
	}

	return 0;

}

int getIDFD(int fd,long offset,int *data,int Npart){

	int n;

	/*First offset the file pointer to go to where the first particle is*/
	if(lseek(fd,offset,SEEK_SET)!=offset) return -1;

	/*Read the particle positions or velocities in the data array*/
	for(n=0;n<Npart;n++){
		if(read(fd,&data[n],sizeof(int))!=sizeof(int)) return -1;
	}

	return 0;

}