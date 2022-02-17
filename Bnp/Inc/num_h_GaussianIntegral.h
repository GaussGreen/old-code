
#ifndef NUM_H_GAUSSIANINTEGRAL_H
#define NUM_H_GAUSSIANINTEGRAL_H

Err GaussianIntegral (double start,  // start < end
					  double end, 
					  int Is_start_Infinity,   //=1 if integral from -Infinity to b
					  int Is_end_Infinity,   //=1 if integral from a to Infinity
					  int n,		       //n has to be greater than or equal to 1
					  double *x,
					  double *w);

#endif