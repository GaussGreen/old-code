// FxSVOptMC.h : MC european option pricing in the FXSV model

#ifndef __FXSVOPTMC_H__
#define __FXSVOPTMC_H__

Err OptFXSVMC(char *fxundname, double gamma, double alpha, double **rho,
              long mat, double *K, int nK, long npaths, int nsteps, double *res,
              double *std);

#endif // #ifndef __FXSVOPTMC_H__