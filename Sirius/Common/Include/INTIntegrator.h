
#ifndef __INTINTEGRATOR_H__
#define __INTINTEGRATOR_H__

#include ".\MCInt\INTSTL.h"

#include ".\MCInt\INTAlpha.h"
#include ".\MCInt\INTWeights.h"
#include ".\MCInt\INTPeriodization.h"
#include ".\MCInt\INTIntegrand.h"
#include "smart.h"
#include "CMatrix.h"


//----------------------------------------------------------------------------------------------

// Base class for integrators
class  CIntegrator : public  RCObject
{
 protected:
  // Compute integral of f over domain D using N function samples
//  virtual double integrate(const CIntegrand& f, const CDomain& D, long N) = 0;
  virtual double integrate(const CIntegrand& f, const CDomain& D, long N)
  {assert(false);return -1e99;};
};

//----------------------------------------------------------------------------------------------

// Plain vanilla Monte Carlo integrator
class  CMCintegrator : public CIntegrator
{
 public:
  double integrate(const CIntegrand& f, const CDomain& D, long N);
};

//----------------------------------------------------------------------------------------------

// Parameters for Number Theoretic integration
class  CNTparameters : public  RCObject
{
  public:
   CAlpha* _a;
   CWeights* _w;
   CPeriodization* _p;

   CNTparameters(CAlpha* a = NULL, CWeights* w = NULL, CPeriodization* p = NULL) : _a(a), _w(w), _p(p) {}  
};

//----------------------------------------------------------------------------------------------

// Number Theoretic integrator
class  CNTintegrator : public CIntegrator
{
 private:
  // Default parameters
  CNTparameters _defaultpar;

 public:
  CNTintegrator() {}
  CNTintegrator(const CNTparameters& defaultpar) : _defaultpar(defaultpar) {}
    
  // Compute integral of f over domain D using N function samples and specified set of parameters
  double integrate(const CIntegrand& f, const CDomain& D, long N, const CNTparameters& par);

  void integrate(STLDoubleVector& integral,const vector<CIntegrand>& f, const CDomain& D, long N, const CNTparameters& par);

  void getSlicedIntegrationPoints(vector<double>& x,int islice, const CDomain& D, long N, const CNTparameters& par);

  void getIntegrationPoint( int m, vector<double>& x, const CDomain& D,  const CNTparameters& par);

  void getSlicedIntegrationPoints(CVector& x,int islice, const CDomain& D, long N, const CNTparameters& par);


  inline double integrate(const CIntegrand& f, const CDomain& D, long N) { return integrate(f, D, N, _defaultpar); }
};

//----------------------------------------------------------------------------------------------


#endif
