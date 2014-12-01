#ifndef __NICAEA_WRAPPERS_H
#define __NICAEA_WRAPPERS_H

#include "lensing.h"
#include "errorlist.h"

static inline double xi_plus(cosmo_lens *model,double theta,int i,int j,error **err){
	return xi(model,1,theta,i,j,err);
}


static inline double xi_minus(cosmo_lens *model,double theta,int i,int j,error **err){
	return xi(model,-1,theta,i,j,err);
}

#endif