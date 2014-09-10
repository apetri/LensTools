#ifndef __READ_GADGET_H
#define __READ_GADGET_H

#include <stdio.h>

//Gadget2 snapshot header
struct io_header_1
{
  int npart[6];
  double mass[6];
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	/* fills to 256 Bytes */
};


//Methods
int getHeader(FILE *fp,struct io_header_1 *header);
int getPosVel(FILE *fp,int offset,float *data);

#endif