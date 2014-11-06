#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "design.h"

/*This function computes the p-distance between 2 points in D dimensions:
the exponent p is tunable*/

double distance(double *x,double *y,int D,double p){

	double dist = 0.0;
	int i;

	for(i=0;i<D;i++){
		dist += pow(fabs(x[i]-y[i]),p);
	}

	return pow(dist,1.0/p);

}

/*This is the cost function of the problem: it is a sum of all the reciprocal of 
the pairs distances, with some tunable exponent lambda; note that the final
1/lambda exponentiation is not performed here!*/

double cost(gsl_matrix *data,int Npoints,int D,double p,double lambda){

	double sum = 0.0;
	int i,j;

	for(i=0;i<Npoints;i++){
		for(j=i+1;j<Npoints;j++){

			//Add the contribution of pair (i,j) to the cost function
			sum += pow(pow(D,1.0/p)/distance(gsl_matrix_ptr(data,i,0),gsl_matrix_ptr(data,j,0),D,p),lambda); 

		}
	}

	return (2.0/(Npoints*(Npoints-1)))*sum;

}

/*This computes the cost function in the particular case in which all the points are 
equally spaced on the diagonal of the hypercube*/

double diagonalCost(int Npoints,double lambda){

	double sum = 0.0;
	int i,j;

	for(i=0;i<Npoints;i++){
		for(j=i+1;j<Npoints;j++){

			//Add the contribution of pair (i,j) to the cost function
			sum += pow((Npoints-1)*1.0/(j-i),lambda); 

		}
	}

	return (2.0/(Npoints*(Npoints-1)))*sum;

}

/*This function computes the variation of the cost when a pair of coordinates is exchanged;
this function also performs an in-place swap of the coordinates. More specifically:
this function swaps coordinate d of points i1 and i2 and returns the variation of the cost
function due to this swap*/

/**/

double swap(gsl_matrix *data,int Npoints,int D,double p,double lambda,int i1,int i2, int d){

	double costBefore,costAfter,temp;
	int i;

	//initialize to 0
	costBefore = costAfter = 0.0;

	/*compute the contribution of points i1 and i2 to the cost function, before swapping;
	sum over all the particles except i1 and i2*/
	for(i=0;i<Npoints;i++){
		
		if(i!=i1 && i!=i2){
			costBefore += pow(pow(D,1.0/p)/distance(gsl_matrix_ptr(data,i,0),gsl_matrix_ptr(data,i1,0),D,p),lambda) + pow(pow(D,1.0/p)/distance(gsl_matrix_ptr(data,i,0),gsl_matrix_ptr(data,i2,0),D,p),lambda);
		}
	
	}

	//perform the coordinate swap
	temp = gsl_matrix_get(data,i1,d);
	gsl_matrix_set(data,i1,d,gsl_matrix_get(data,i2,d));
	gsl_matrix_set(data,i2,d,temp);

	/*compute the contribution of points i1 and i2 to the cost function, after swapping;
	sum over all the particles except i1 and i2*/
	for(i=0;i<Npoints;i++){
		
		if(i!=i1 && i!=i2){
			costAfter += pow(pow(D,1.0/p)/distance(gsl_matrix_ptr(data,i,0),gsl_matrix_ptr(data,i1,0),D,p),lambda) + pow(pow(D,1.0/p)/distance(gsl_matrix_ptr(data,i,0),gsl_matrix_ptr(data,i2,0),D,p),lambda);
		}
	
	}

	//return the cost difference
	return (2.0/(Npoints*(Npoints-1)))*(costAfter - costBefore);

}

/*This function swaps the particle pair back: needed if the original swap didn't improve the cost function*/

void swapBack(gsl_matrix *data,int i1,int i2,int d){

	double temp;

	//perform the coordinate swap
	temp = gsl_matrix_get(data,i1,d);
	gsl_matrix_set(data,i1,d,gsl_matrix_get(data,i2,d));
	gsl_matrix_set(data,i2,d,temp);

}

/*Main design sampler*/

double sample(int Npoints,int D,double p,double lambda,int seed,int maxIterations,gsl_matrix *data,double *costValues){

	int i,d,i1,i2,iterCount;
	
	const gsl_rng_type *T;
	gsl_rng *r;
	gsl_permutation *perm = gsl_permutation_alloc(Npoints);

	double currentCost,deltaCost,deltaPerc;

	//Initialize random number generator
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	//Initialize permutation
	gsl_permutation_init(perm);

	//Initialize random number generator with provided seed
	gsl_rng_set(r,seed);

	//Initialize the point coordinates in data with random permutations of (1..Npoints) to enforce latin hypercube structure
	for(d=0;d<D;d++){

		//Shuffle the numbers
		gsl_ran_shuffle(r,perm->data,Npoints,sizeof(size_t));

		//Permute coordinates
		for(i=0;i<Npoints;i++){
			gsl_matrix_set(data,i,d,(double)perm->data[i]/(Npoints-1));
		}

	}

	/*The loop does the following: it swaps a random coordinate of a random pair,
	checks if the cost is lower. If so, it keeps the configuration, otherwise it
	reverses it and tries a new one.*/

	iterCount = 0;
	currentCost = cost(data,Npoints,D,p,lambda);
	deltaPerc = 0.0;

	while(1){

		//Decide which coordinate to swap of which pair

		i1 = gsl_rng_uniform_int(r,Npoints);
		while((i2=gsl_rng_uniform_int(r,Npoints))==i1);
		d = gsl_rng_uniform_int(r,D);

		//Compute the change in the cost function
		deltaCost = swap(data,Npoints,D,p,lambda,i1,i2,d);


		/*Check if gain in cost is positive or negative: if positive, revert the swap, if negative keep it;
		anyway, log the result*/
		if(deltaCost>=0){
			swapBack(data,i1,i2,d);
		} else{
			currentCost += deltaCost;
			deltaPerc = deltaCost/currentCost;
		}

		//Save the current cost to array
		costValues[iterCount] = currentCost;

		//Criterion to break the loop
		if(++iterCount == maxIterations) break;
	
	}

	//Release resources for random number generator and permutations
	gsl_rng_free(r);
	gsl_permutation_free(perm);

	//Return the relative cost change due to the last iteration
	return deltaPerc;

}