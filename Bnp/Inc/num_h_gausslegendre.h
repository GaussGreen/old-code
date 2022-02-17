/* ======================================================
   FILENAME:  num_h_gausslegendre.h
   
   PURPOSE:   integration by Gauss Legendre (NUMC p152)
   ====================================================== */

#ifndef NUM_H_GAUSSLEGENDRE_H
#define NUM_H_GAUSSLEGENDRE_H

void gauleg(double x1,double x2,double x[], double w[],int n);

void GaussLeg(	double x1, 
				double x2, 
				double *x, 
				double *w, 
				int n);

#endif