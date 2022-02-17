
#ifndef __INTPERIODIZATION_H__
#define __INTPERIODIZATION_H__

#include "INTSTL.h"
#include "INTUtilities.h"
#include "smart.h"


//----------------------------------------------------------------------------------------------

#ifdef _S
#undef _S
#endif

// Base class for periodizations
class  CPeriodization : public RCObject
{
 protected:
  // Acronym
  string _S;

  // First derivative of phi
//  virtual double phi1(double x) const = 0;
  virtual double phi1(double x) const 
  {assert(false);return -1e99;};

 public:

  // Change of variables mapping |R -> |R
//  virtual double phi(double x) const = 0;  
  virtual double phi(double x) const
  {assert(false);return -1e99;};


  // Plot phi(x), or phi1(x), or ... vs x for given periodizations
  static void plot(const vector<CPeriodization*>& p, int choice, const string& filename, int res = 20, int margin = 0);

  // Plot phi(x), or phi1(x), or ... vs x for this periodization
  void plot(int choice, const string& filename = "", int res = 20, int margin = 0);

  // Change of variables mapping |R^k -> |R^k
  STLDoubleVector& phi(STLDoubleVector& x) const;

  // Weights acronym
  inline const char *acronym(void) const { return _S.c_str(); }

  // Product Pi_{i=1}^k phi1(x_i)
  virtual double product(const STLDoubleVector& x) const;

  virtual int jacobianIsConstant(){return 0 ;};

  // Report periodization
  inline virtual void report(void) const {};
};

//----------------------------------------------------------------------------------------------

// Change of variables renders f continous on boundaries of periodic cell
class  CZerotypeperiodization : public CPeriodization
{
 protected:
  inline double phi1(double x) const { return 1.0; }

 public:
  inline double product(const STLDoubleVector& phi1x) { return 1.0; }
  virtual int jacobianIsConstant(){return 1 ;};

};

//----------------------------------------------------------------------------------------------

// Identity periodization id
class  CIdentityperiodization : public CZerotypeperiodization
{
 private:
  inline double phi(double x) const { return x; }

 public:
  CIdentityperiodization() { _S = "p-id"; }

  inline void report(void) const { std::cout << "Identity periodization"; }
};

//----------------------------------------------------------------------------------------------

// Baker periodization bak
class  CBakerperiodization : public CZerotypeperiodization
{
 private:
  inline double phi(double x) const { return 1.0 - fabs(2.0*x - 1.0); }

 public:
  CBakerperiodization() { _S = "p-bak"; } 

  inline void report(void) const { std::cout << "Baker periodization"; }
};

//----------------------------------------------------------------------------------------------

// Change of variables renders f continuous on boundaries of periodic cell and sets derivatives 
// f^(i) to zero for i in {1, ..., d}. Note: d = 0 has no effect on the integrand.
class  CFinitetypeperiodization : public CPeriodization
{
 protected:
  // Degree of symmetry d in {1, 2, ...}
  int _d;

  CFinitetypeperiodization(int d) : _d(d) {}
};

//----------------------------------------------------------------------------------------------

// Polynomial periodization poly [Korobov63]
class  CPolynomialperiodization : public CFinitetypeperiodization
{
 private:
  // Normalizing constant c
  double _c;

  // Weighting coefficients for summation in method phi(x)
  inline double A(int i) const { return ((i%2 == 0) ? 1.0 : -1.0)/(double)(_d + i + 1)*binomial(_d, i); }

  double phi(double x) const;
  
  inline double phi1(double x) const { return _c*pow(x*(1.0 - x), _d); }

 public:
  CPolynomialperiodization(int d);

  inline void report(void) const { std::cout << "Polynomial periodization: d = " << _d; }
};

//----------------------------------------------------------------------------------------------

// Trigonometric periodization trig [Sidi93]
// Note: For d = 1: 1/2*(1 - cos(Pi*x)), 
//           d = 2: x - 1/Pi*sin(Pi*x)*cos(Pi*x) = x - 1/(2*Pi)*sin(2*Pi*x)
class  CTrigonometricperiodization : public CFinitetypeperiodization
{
 private:
  // Normalizing constant c
  double _c;

  // Upper end of range for summation in method phi(x)
  int _i_max;

  double phi(double x) const;

  inline double phi1(double x) const { return _c*pow(sin(PI*x), _d); }

 public:
  CTrigonometricperiodization(int d);

  inline void report(void) const { std::cout << "Trigonometric periodization: d = " << _d; }
};

//----------------------------------------------------------------------------------------------

// Kress periodization kre
class  CKressperiodization : public CFinitetypeperiodization
{
 private:
  double phi(double x) const;

  double phi1(double x) const;

 public:
  CKressperiodization(int d) : CFinitetypeperiodization(d) { _S = string("p-kre(") + _d + ")"; }

  inline void report(void) const { std::cout << "Kress periodization: d = " << _d; }
};

//----------------------------------------------------------------------------------------------

// Change of variables renders f continuous on boundaries of periodic cell 
// and sets derivatives f^(i) to zero for all i in |N.
class  CInfinitetypeperiodization : public CPeriodization {};

//----------------------------------------------------------------------------------------------

// TANH periodization [SagSzekeres64])
// Note: with a = p = 1: exp(-1/x)/(exp(-1/x) + exp(-1/(1-x)))
class  CTANHperiodization : public CInfinitetypeperiodization
{
 private:
  // Parameters a in |R, p in |N
  double _a;
  int _p;

  double phi(double x) const;

  double phi1(double x) const;

 public:
  CTANHperiodization(double a, int p) : _a(a), _p(p) { _S = string("p-TANH(") + _a + ", " + _p + ")"; }

  inline void report(void) const { std::cout << "TANH periodization: a = " << _a << ", p = " << _p; }
};

//----------------------------------------------------------------------------------------------

// DE (Double Exponential) periodization
class  CDEperiodization : public CInfinitetypeperiodization
{
 private:
  // Parameter a in |R^+
  double _a;

  double phi(double x) const;

  double phi1(double x) const;

 public:
  CDEperiodization(double a) : _a(a) { _S = string("p-DE(") + _a + ")"; }

  inline void report(void) const { std::cout << "DE periodization: a = " << _a; }
};

//----------------------------------------------------------------------------------------------

// Gaussian periodization gauss
class  CGaussianperiodization : public CInfinitetypeperiodization
{
 private:
  // Parameter sigma in |R^+
  double _sigma;

  inline double phi(double x) const { return Ncdf(_sigma*(2.0*x - 1.0)); }

  inline double phi1(double x) const { return Ndf(_sigma*(2.0*x - 1.0))*2.0*_sigma; }

 public:
   CGaussianperiodization(double sigma) : _sigma(sigma) { _S = string("p-gauss(") + _sigma + ")"; }

  inline void report(void) const { std::cout << "Gaussian periodization: sigma = " << _sigma; }
};

//----------------------------------------------------------------------------------------------


#endif





