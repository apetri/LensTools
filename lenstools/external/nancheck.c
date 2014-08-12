#include <stdio.h>
#include <math.h>


int main(){

#ifdef NAN
	printf("nan is defined\n");
#endif

#ifdef INFINITY
	printf("infinity is defined\n");
#endif

	double inf = INFINITY;
	double NaN = NAN;

	double x = 2.0;
	double y = x - inf;

	printf("x=%lf y=%lf\n",x,y);

	printf("sizeof nan is %ld bytes\n",sizeof(inf));
	printf("size of infinity is %ld bytes\n",sizeof(NaN));

	return 0;

}