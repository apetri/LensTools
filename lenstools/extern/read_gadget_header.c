#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "gadget2.h"

/*This routine reads the information from the snapshot header
and fills out a struct io_header1 instance*/
int getHeader(FILE *fp,struct io_header_1 *header){

	char buf[4];
	
	/*Skip the first 4 bytes of the snapshot*/
	if(fread(buf,sizeof(buf),1,fp)!=1) return -1;

	/*Read in the header*/
	if(fread(header,sizeof(struct io_header_1),1,fp)!=1) return -1;

	/*Perform an endianness check*/
	if(memcmp(buf,"\x00\x01\x00\x00",4)==0){	
		
		//little endian
		return 0;
	
	} else if(memcmp(buf,"\x00\x00\x01\x00",4)==0){
		
		//big endian
		return 1;
	
	} else{

		//not a valid gadget2 snapshot
		return -1;
	}

}

int getHeaderFD(int fd,struct io_header_1 *header){

	char buf[4];
	
	/*Skip the first 4 bytes of the snapshot*/
	if(read(fd,buf,4)!=4) return -1;

	/*Read in the header*/
	if(read(fd,header,sizeof(struct io_header_1))!=sizeof(struct io_header_1)) return -1;

	/*Perform an endianness check*/
	if(memcmp(buf,"\x00\x01\x00\x00",4)==0){	
		
		//little endian
		return 0;
	
	} else if(memcmp(buf,"\x00\x00\x01\x00",4)==0){
		
		//big endian
		return 1;
	
	} else{

		//not a valid gadget2 snapshot
		return -1;
	}

}