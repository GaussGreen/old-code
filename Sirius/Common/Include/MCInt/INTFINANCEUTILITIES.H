
#ifndef __INTFINANCEUTILITITES_H__
#define __INTFINANCEUTILITITES_H__

#include "INTSTL.h"
#include "INTUtilities.h"

//----------------------------------------------------------------------------------------------

// Black-Scholes formula for European Call/Put options
// Parameters type: call/put, S0: spot, K: strike, rho: interest rate, T: maturity, vol: volatility
 double BS(int type, double S0, double K, double rho, double T, double vol);

//----------------------------------------------------------------------------------------------

// Input utilities

// Parameters k: dimension and v: value(s) (redundant if to be generated uniformly at random)
// uni*(k): value(s) generated uniformly at random
// equal*(k, v): equal value in all dimensions, v^k
// user*(k, v): user-specified values

// Time steps t: 0 = t_0 <= t_1 <= ... <= t_k = T
 STLDoubleVector& unit		(STLDoubleVector& t, int k, double T);
 STLDoubleVector& equalt	(STLDoubleVector& t, int k, double T);
 STLDoubleVector& usert	(STLDoubleVector& t, int k, double T);

// Volatilities vol: (vol_1, ..., vol_k) in [vol_min, vol_max)^k
 STLDoubleVector& univol	(STLDoubleVector& vol, int k, double vol_min, double vol_max);
 STLDoubleVector& equalvol	(STLDoubleVector& vol, int k, double commonvol);
 STLDoubleVector& uservol	(STLDoubleVector& vol, int k);

// Correlations corr: kxk-matrix with corr_ii = 1, corr_ij = corr_ji in [corr_min, corr_max)
 STLDoubleVectorVector& unicorr	(STLDoubleVectorVector& corr, int k, double corr_min, double corr_max);
 STLDoubleVectorVector& equalcorr	(STLDoubleVectorVector& corr, int k, double commoncorr);
 STLDoubleVectorVector& usercorr	(STLDoubleVectorVector& corr, int k, ...);

// Stock prices S: (S_1, ..., S_k) in [S_min, S_max)^k
 STLDoubleVector& uniS		(STLDoubleVector& S, int k, double S_min, double S_max);
 STLDoubleVector& equalS	(STLDoubleVector& S, int k, double commonS);
 STLDoubleVector& userS	(STLDoubleVector& S, int k);

//----------------------------------------------------------------------------------------------


#endif
