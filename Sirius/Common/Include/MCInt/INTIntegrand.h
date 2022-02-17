
#ifndef __INTINTEGRAND_H__
#define __INTINTEGRAND_H__

#include "INTSTL.h"
#include "INTUtilities.h"
#include "INTDomain.h"

class CPeriodization;

//----------------------------------------------------------------------------------------------

// Base class for integrands
class  CIntegrand
{
 protected:
  // Acronym
  string _S;

 public:
  // Plot integrands of one or two variables after domain transformation (except rectilinear domain)
  static void plot(const vector<CIntegrand*>& f, const CDomain& D, 
                   const string& filename, int res = 20, int margin = 0, bool excel = true);

  // Plot integrands of one or two variables after domain transformation and periodization
  static void plot(const vector<CIntegrand*>& f, const CDomain& D, const CPeriodization& p, 
                   const string& filename, int res = 20, int margin = 0, bool excel = true);

  // Plot this integrand
  void plot(const CDomain& D, const string& filename = "", int res = 20, int margin = 0, bool excel = true);
  void plot(const CDomain& D, const CPeriodization& p, const string& filename = "", int res = 20, int margin = 0, bool excel = true);

  // Integrand function after domain transformation and possibly periodization
  // Note: This method will change the value of x.
  double evaluate(STLDoubleVector& x, const CDomain& D, const CPeriodization *p = NULL) const;

  // Integrand acronym
  const char *acronym(void) const;

  // Integrand function |R^k -> |R
  virtual double evaluate(const STLDoubleVector& x) const = 0;

  // Analytic integral of f over unit cube [0, 1)^k (if known)
  virtual double I(int k) const;

  // Report integrand
  inline virtual void report(void) const { std::cout << "other integrand"; }
};

//----------------------------------------------------------------------------------------------

// Integrand instantiated from a function |R^k -> |R or |R -> |R
class  Cf : public CIntegrand
{
 private:
  // Integrand function |R^k -> |R
  double (*_fk)(const STLDoubleVector&);

  // Integrand function |R -> |R
  double (*_f1)(double);

 public:
  Cf(double (*f)(const STLDoubleVector&)) : _fk(f), _f1(NULL) { _S = "f"; }
  Cf(double (*f)(double)) : _fk(NULL), _f1(f) { _S = "f"; }

  inline double evaluate(const STLDoubleVector& x) const { return (_f1 ? (*_f1)(x[0]) : (*_fk)(x)); }

  inline void report(void) const { std::cout << "integrand f"; }
};

//----------------------------------------------------------------------------------------------

// Integrand after domain transformation and (possibly) periodization
class  CTransformedintegrand : public CIntegrand
{
 private:
  const CIntegrand& _f;
  const CDomain& _D;
  const CPeriodization *_p;

 public:
  CTransformedintegrand(const CIntegrand& f, const CDomain& D, const CPeriodization *p = NULL);

  // Note: to ensure that x remains unaltered, we have to use a copy when calling _f.evaluate(x, D, p)
  inline double evaluate(const STLDoubleVector& x) const { STLDoubleVector y(x); return _f.evaluate(y, _D, _p); }

  // Note: periodizations preserve the integral over [0, 1)^k
  inline double I(int k) const { return _f.I(k); };

  inline void report(void) const { std::cout << "transformed "; _f.report(); }
};

//----------------------------------------------------------------------------------------------


#endif
