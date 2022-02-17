
#ifndef __INTFINANCEINTEGRANDS_H__
#define __INZFINANCEINTEGRANDS_H__

#include "INTSTL.h"
#include "INTUtilities.h"
#include "INTIntegrand.h"
#include "INTFinanceUtilities.h"

//----------------------------------------------------------------------------------------------

// Plain vanilla European option computed over k time subintervals 0 = t_0 <= t_1 <= ... <= t_k = T
class  CEuropeanintegrand : public CIntegrand
{
 private:
  // Type (call or put), spot S0, strike K, interest rate rho, maturity T, volatility vol
  int _type;
  double _S0, _K, _rho, _T, _vol;

  // Time subintervals 0 = t_0 <= t_1 <= ... <= t_k = T
  STLDoubleVector _t;

  // Statistics from (t_i)
  STLDoubleVector _c;

 public:
  CEuropeanintegrand(const STLDoubleVector& t, int type = 1, double S0 = 100.0, 
                     double K = 110.0, double rho = 0.0, double T = 1.0, double vol = 0.2);

  double evaluate(const STLDoubleVector& x) const;
        
  double I(int k) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Asian option using geometric average of the spot price at control times 0 = t_0 <= t_1 <= ... <= t_k = T
class  CAsianintegrand : public CIntegrand
{
 private:
  // Type (call or put), spot S0, strike K, interest rate rho, maturity T, volatility vol
  int _type;
  double _S0, _K, _rho, _T, _vol;

  // Control times 0 = t_0 <= t_1 <= ... <= t_k = T
  STLDoubleVector _t;

  // Statistics from (t_i)
  double _c1;
  STLDoubleVector _c2;
  double _c3;

 public:
  CAsianintegrand(const STLDoubleVector& t, int type = 1, double S0 = 100.0, 
                  double K = 110.0, double rho = 0.0, double T = 1.0, double vol = 0.2);

  double evaluate(const STLDoubleVector& x) const;
        
  double I(int k) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Basket option for k assets with spots S0_i, volatilities vol_i and correlation STLDoubleVectorVector corr_ij
class  CBasketintegrand : public CIntegrand
{
 private:
  // Type (call or put), strike K, interest rate rho, maturity T
  int _type;
  double _K, _rho, _T;

  // Spots S0_i, volatilities vol_i, correlation matix corr_ij of k assets
  STLDoubleVector _S0;
  STLDoubleVector _vol;

  STLDoubleVectorVector _corr;

  // Statistics from (S0_i), (vol_i) and (corr_ij)
  double _c1, _c2;
  STLDoubleVector _c3;

 public:
  CBasketintegrand(const STLDoubleVector& S0, const STLDoubleVector& vol, const STLDoubleVectorVector& corr,
                   int type = 1, double K = 110.0, double rho = 0.0, double T = 1.0);

  double evaluate(const STLDoubleVector& x) const;
        
  double I(int k) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------


#endif
