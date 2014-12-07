#include <stdio.h>
#include "gadget.h"

/*This routine reads the information from the snapshot header
and fills out a struct io_header1 instance*/
int getHeader(FILE *fp,struct io_header_1 *header){

	char buf[4];
	
	/*Skip the first 4 bytes of the snapshot*/
	if(fread(buf,sizeof(buf),1,fp)!=1) return -1;

	/*Read in the header*/
	if(fread(header,sizeof(struct io_header_1),1,fp)!=1) return -1;

	/*Perform an endianness check*/
	if(buf[1]==1){	
		
		//little endian
		return 0;
	
	} else if(buf[1]==0){
		
		//big endian
		return 1;
	
	} else{

		//not a valid gadget2 snapshot
		return -1;
	}

}