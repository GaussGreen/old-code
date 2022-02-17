#include	"stdafx.h"


#include "MCInt\INTFinanceUtilities.h"



//----------------------------------------------------------------------------------------------

// Black-Scholes formula for European Call/Put options
// Parameters type: call/put, S0: spot, K: strike, rho: interest rate, T: maturity, vol: volatility
double BS(int type, double S0, double K, double rho, double T, double vol)
{
	double d1, bs;
	double volsqrtT;
	
	volsqrtT = vol*sqrt(T);
	
	d1 = (log(S0/K) + (rho + 0.5*vol*vol)*T)/volsqrtT;
	
	bs = S0*Ncdf(d1) - K*exp(-rho*T)*Ncdf(d1 - volsqrtT);

	/*
	switch(type)
    {
	}
	*/
	
	return bs;
}

//----------------------------------------------------------------------------------------------

// Time steps t: 0 = t_0 <= t_1 <= ... <= t_k = T
// Parameters T: maturity, t = [t_1, ..., t_{k-1}]: list of user defined intermediate/control times
STLDoubleVector& unit(STLDoubleVector& t, int k, double T)
{ 
  if (T < 0.0) ERROR("In unit: T must be non-negative. ")

  t.resize(k+1);
  
  t[0] = 0.0;
  for (int i = 1; i < k; i++)
      t[i] = randomU01();
  t[k] = T;

  sort(t.begin(), t.end());
  return t;
}

STLDoubleVector& equalt(STLDoubleVector& t, int k, double T)
{
  if (T < 0.0) ERROR("In equalt: T must be non-negative. ")

  t.resize(k+1);

  for (int i = 0; i <= k; i++)
    {
      t[i] = (double)i*T/(double)k;
    }

  return t;
}

STLDoubleVector& usert(STLDoubleVector& t, int k, double T)
{
  // Validate and sort user defined t
  if (T < 0.0) ERROR("In usert: T must be non-negative. ")
  if (t.size() != k-1) ERROR("In usert: wrong number of time steps. ")

  t.insert(t.begin(), 0);
  t.push_back(T);
  sort(t.begin(), t.end());

  if ((t[1] <= 0) || (t[k-1] >= T)) ERROR("In usert: t_i must be in (0, T) for 1 <= i <= k-1. ")

  return t;
}

//----------------------------------------------------------------------------------------------

// Volatilities vol: (vol_1, ..., vol_k) in [vol_min, vol_max)^k
// Parameters vol_min, vol_max: range of vol for uniform sampling resp. vol: desired value(s)
STLDoubleVector& univol(STLDoubleVector& vol, int k, double vol_min, double vol_max)
{ 
  if (vol_min < 0) ERROR("In univol: Lower bound of sample range must be non-negative. ")
  if (vol_max < vol_min) ERROR("In univol: vol_min must be less than vol_max. ")

  vol.resize(k);
  randomU01k(vol);
  vol *= vol_max - vol_min;
  vol += vol_min;

  return vol;
}

STLDoubleVector& equalvol(STLDoubleVector& vol, int k, double commonvol)
{ 
  if (commonvol < 0) ERROR("In equalvol: vol must be non-negative. ")

  vol.clear();
  vol.resize(k, commonvol);
  return vol;
}

STLDoubleVector& uservol(STLDoubleVector& vol, int k)
{ 
  // Validate user defined vol
  if (vol.size() != k) ERROR("In uservol: wrong number of volatilities. ")

  for (int i = 0; i < k; i++)
    {
      if (vol[i] < 0) ERROR("In uservol: volatilities must be non-negative. ")
    }

  return vol;
}

//----------------------------------------------------------------------------------------------

// Correlations corr: kxk-matrix with corr_ii = 1, corr_ij = corr_ji in [corr_min, corr_max)
// Parameters corr_min, corr_max: range of corr for uniform sampling resp. corr: desired value(s)
STLDoubleVectorVector& unicorr(STLDoubleVectorVector& corr, int k, double corr_min, double corr_max)
{ 
  if (corr_min < 0) ERROR("In unicorr: Lower bound of sample range must be non-negative. ")
  if (corr_max < corr_min) ERROR("In unicorr: corr_min must be less than corr_max. ")

  corr.resize(k, STLDoubleVector(k, 0.0)); 

  for (int r = 0; r < k; r++)
    {
      for (int c = 0; c < k; c++)
        {
          if (c < r) corr[r][c] = corr[c][r];
          else if (c == r) corr[r][c] = 1.0;
          else corr[r][c] = randomUab(corr_min, corr_max);
        }
    }

  return corr;
}

STLDoubleVectorVector& equalcorr(STLDoubleVectorVector& corr, int k, double commoncorr)
{ 
  if (commoncorr < 0) ERROR("In equalcorr: corr must be non-negative. ")

  corr.clear();
  corr.resize(k, STLDoubleVector(k, commoncorr));
  for (int r = 0; r < k; r++) corr[r][r] = 1;

  return corr;
}

// Parameter corr = [[corr_12, ..., corr_1k], ..., [corr_{k-1}k]] upper-diagonal entries
STLDoubleVectorVector& usercorr(STLDoubleVectorVector& corr, int k, ...)
{ 
  // Validate user defined corr and complete full matrix
  corr.resize(k, STLDoubleVector(k, 1.0));

  // Transcribe upper diagonal entries
  va_list ap;
  va_start(ap, k);

  int n = (k-1)*k/2;
  int r = 0, c = r;
  for (int i = 0; i < n; i++)
    {
      if (++c == k) { r++; c = r+1; }
      corr[r][c] = va_arg(ap, double);
      if (corr[r][c] < 0) ERROR("In usercorr: correlations must be non-negative. ")
    }
  va_end(ap);

  // Fill in entries on and below main diagonal
  for (r = 1; r < k; r++)
    {
      for (c = 0; c < r; c++) 
        {
          corr[r][c] = corr[c][r];
        }
    }

  return corr;
}

//----------------------------------------------------------------------------------------------

// Stock prices S: (S_1, ..., S_k) in [S_min, S_max)^k
// Parameters S_min, S_max: range of S for uniform sampling resp. S: desired value(s)
STLDoubleVector& uniS(STLDoubleVector& S, int k, double S_min, double S_max)
{ 
  S.resize(k);
  randomUabk(S, S_min, S_max);
  return S;
}

STLDoubleVector& equalS(STLDoubleVector& S, int k, double commonS)
{ 
  S.clear();
  S.resize(k, commonS);
  return S;
}

STLDoubleVector& userS(STLDoubleVector& S, int k)
{ 
  // Validate user defined S
  if (S.size() != k) ERROR("In userS: wrong number of stock prices. ")

  return S;
}

//----------------------------------------------------------------------------------------------

