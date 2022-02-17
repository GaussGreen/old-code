/* --------------------------------------------------------
	FILENAME: 	num_h_levenberg.h
	
	FUNCTION:   levenberg_marquardt optimisation algorithm
   -------------------------------------------------------- */

#ifndef NUM_H_LEVENBERG_H
#define NUM_H_LEVENBERG_H

Err levenberg_marquardt(
		double *data,   /* From [1] to [ndata] */  
		double *target, /* From [1] to [ndata] */ 
		double *weight, /* From [1] to [ndata] */ 
		long ndata,
		double *param,  /* From [1] to [nparam] */ 
		long nparam,
		long niter,
		Err (*funcs)(double, double[], double*, double[], int),
		double *chisq);


Err levenberg_marquardt_select(
		double *data,   /* From [1] to [ndata] */  
		double *target, /* From [1] to [ndata] */ 
		double *weight, /* From [1] to [ndata] */ 
		long ndata,
		double *param,  /* From [1] to [nparam] */
		long *use_param,/* From [1] to [nparam] */
		long nparam,
		long niter,
		Err (*funcs)(double, double[], double*, double[], int),
		double *chisq);

#endif 
