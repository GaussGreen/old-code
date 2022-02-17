#include <stdio.h>
#include <math.h>
#include <string.h>
#include "linalg.h" 

#define EPS 1.0e-9
#define INF exp(3000)



#include "xdfpmin.h"



double leastsq(ARM_Vector *p,
			   void **parameters,
			   ARM_Matrix *data,
			   T_FUNC func);

void searchq(int *pcnt,
				double *fnew,
				double oldx,
				ARM_Vector *matl,
				ARM_Vector *matx,
				double gdold,
				double *stepsize);	
				
double cubici2(double graddold,
			   ARM_Vector *matl,
			   ARM_Vector *matx);

void jacobian(ARM_Vector* x,
			  ARM_Matrix* grad,
			  void ** parameters,
			  ARM_Matrix* data,
			  T_FUNC func );

