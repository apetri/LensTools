#ifndef __READ_GADGET2_H
#define __READ_GADGET2_H

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
  double w0;
  double wa;
  double comoving_distance;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int flag_entropy_instead_u;
  int nothing;
  char fill[32];	/* fills to 256 Bytes */
};


//Methods
int getHeader(FILE *fp,struct io_header_1 *header);
int getPosVel(FILE *fp,long offset,float *data,int Npart);
int getID(FILE *fp,long offset,int *data,int Npart);
int writeSnapshot(FILE *fp,struct io_header_1 *header,float *positions,float *velocities,int firstID,int NumPart,int writeVel);

int getHeaderFD(int fd,struct io_header_1 *header);
int getPosVelFD(int fd,long offset,float *data,int Npart);
int getIDFD(int fd,long offset,int *data,int Npart);
int writeSnapshotFD(int fd,struct io_header_1 *header,float *positions,float *velocities,int firstID,int NumPart,int writeVel);

#endif