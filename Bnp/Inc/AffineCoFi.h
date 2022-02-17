// AffineCoFi.h:  Cornish-Fisher expansion for affine models

#ifndef __AFFINECOFI_H__
#define __AFFINECOFI_H__

Err HestonDDensityFTExpansion(double sigma, double alpha, double rho, double lam, double mat,
							  int order, int nu, double *u_re, double *u_im, double *h_re, double *h_im);

Err HestonDCoFiMC(double sigma, double alpha, double rho, double lam, double mat, int order, int cf_ord,
				  long npaths, double fwd, int nK, double *K, double *pv, double *stddev);

#endif  // #ifndef __AFFINECOFI_H__