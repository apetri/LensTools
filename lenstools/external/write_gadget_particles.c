#include <stdio.h>
#include "read_gadget.h"

/*this routine writes particle positions and velocities to a snapshot in Gadget's default
binary file format*/
int writeSnapshot(FILE *fp,struct io_header_1 *header,float *positions,float *velocities,int firstID,int NumPart){

	int garbage=256;

	//the first 4 bytes are for endianness check
	if(fwrite(&garbage,sizeof(int),1,fp)!=1) return -1;

	//the header comes next
	if(fwrite(header,sizeof(struct io_header_1),1,fp)!=1) return -1;

	//the next 8 bytes are garbage
	if(fwrite(&garbage,sizeof(int),1,fp)!=1) return -1;
	if(fwrite(&garbage,sizeof(int),1,fp)!=1) return -1;


	return 0;
	
}