#ifndef __INTMISCINTEGRANDS_H__
#define __INTMISCINTEGRANDS_H__


#include "INTSTL.h"
#include "INTIntegrand.h"

//----------------------------------------------------------------------------------------------

// Identity integrand id
class  CIdentityintegrand : public CIntegrand
{
 public:
  CIdentityintegrand() { _S = "f-id"; }

  inline double evaluate(const STLDoubleVector& x) const { assert(x.size() == 1); return x[0]; };
  inline double I(int k) const { assert(k == 1); return 1.0; };

  inline void report(void) const { std::cout << "identity"; }
};

//----------------------------------------------------------------------------------------------

// Zaremba example zar1 [Zaremba70, p. 20]
class  CZaremba1integrand : public CIntegrand
{
 public:
  CZaremba1integrand() { _S = "f-zar1"; }

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { assert(k == 2); return 1.0; };

  inline void report(void) const { std::cout << "Zaremba example 1"; }
};

// Zaremba example zar2 [Zaremba70, p. 20]
class  CZaremba2integrand : public CIntegrand
{
 public:
  CZaremba2integrand() { _S = "f-zar2"; }

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { assert(k == 2); return 1.0/8.0*erf(1.0)*PI*sqrt(2.0)*erf(sqrt(2.0)); };

  inline void report(void) const { std::cout << "Zaremba example 2"; }
};

// Zaremba example zar3 [Zaremba70, p. 21]
class  CZaremba3integrand : public CIntegrand
{
 public:
  CZaremba3integrand() { _S = "f-zar3"; }

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { assert(k == 2); return 1.0; };

  inline void report(void) const { std::cout << "Zaremba example 3"; }
};

//----------------------------------------------------------------------------------------------

// Constant integrand const
class  CConstantintegrand : public CIntegrand
{
 private:
  // Constant value c of the integrand
  double _c;

 public:
  CConstantintegrand(double c) : _c(c) { _S = string("f-const(") + _c + ")"; }

  inline double evaluate(const STLDoubleVector& x) const { return _c; };
  inline double I(int k) const { return _c; };

  inline void report(void) const { std::cout << "constant integrand c = " << _c; }
};

//----------------------------------------------------------------------------------------------

// Haselgrove worst case integrands hwc [Haselgrove61], p.331
class  CHaselgroveWorstCaseintegrand : public CIntegrand
{
private:
  // Decay rate s in {2, 4} of Fourier coefficients
  int _s;

public:
  CHaselgroveWorstCaseintegrand(int s);

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { return pow(((_s == 2) ? 1.0/3.0 : 7.0/15.0), k); };

  inline void report(void) const { std::cout << "Haselgrove worst case integrand: s = " << _s; }
};

// Haselgrove example has [Haselgrove61], p. 337
class  CHaselgroveintegrand : public CIntegrand
{
 public:
  CHaselgroveintegrand() { _S = "f-has"; }

  double evaluate(const STLDoubleVector& x) const;

  // Note: analytic integral over [0, 1)^k only available for k = 5.
  inline double I(int k) const { assert(k == 5); return 0.97065719; };

  inline void report(void) const { std::cout << "Haselgrove integrand"; }
};

//----------------------------------------------------------------------------------------------

// Tsuda example tsu [Tsuda73], p. 384, in [SugiharaMurota82], p. 552
// Note: Strongly peaked integrand is designed to exhibit limitations of Haselgrove's method.
class  CTsudaintegrand : public CIntegrand
{
 private:
  // Parameter c in |R
  double _c;

 public:
  CTsudaintegrand(double c) : _c(c) { _S = string("f-tsu(") + _c + ")"; }

  double evaluate(const STLDoubleVector& x) const;    
  inline double I(int k) const { return 1.0; };

  inline void report(void) const { std::cout << "Tsuda integrand: c = " << _c; }
};

//----------------------------------------------------------------------------------------------

// Roos and Arnold example ra1 [RoosArnold63], in [DavisRabinowitz84], p. 405 and [Owen01], p. 20
class  CRoosArnold1integrand : public CIntegrand
{
 public:
  CRoosArnold1integrand() { _S = "f-ra1"; }

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { return 1.0; };

  inline void report(void) const { std::cout << "Roos and Arnold example 1"; }
};

// Roos and Arnold example ra2 [RoosArnold63], in [DavisRabinowitz84], p. 406, [Sobol01], and [Owen01], p. 20
class  CRoosArnold2integrand : public CIntegrand
{
 private:
  // Parameters a_1, ..., a_k in |R
  STLDoubleVector _a;
  
 public:
  CRoosArnold2integrand(const STLDoubleVector& a) : _a(a) {}

  double evaluate(const STLDoubleVector& x) const;    
  inline double I(int k) const { return 1.0; };
};

// Roos and Arnold example ra3 [RoosArnold63], in [DavisRabinowitz84], p. 405 and [Owen01], p. 20
class  CRoosArnold3integrand : public CIntegrand
{
public:
  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { return 1.0; };
};

//----------------------------------------------------------------------------------------------

class  CGenzintegrand : public CIntegrand
{
 private:
  // Parameters a_1, ..., a_k in |R
  STLDoubleVector _a;

  // Parameters u_1, ..., u_k in [0, 1)
  STLDoubleVector _u;
  
 public:
  CGenzintegrand(const STLDoubleVector& a, const STLDoubleVector& u) : _a(a), _u(u) {}
  
  // Note: for analytical integrals, refer to the MAPLE-version
  inline double I(int k) const { assert(false); return 1.0;};
};

// Genz example gen1: product peak [Genz84], in [Owen01], p. 22
class  CGenz1integrand : public CGenzintegrand
{
public:
  CGenz1integrand(const STLDoubleVector& a, const STLDoubleVector& u) : CGenzintegrand(a, u) {}

  double evaluate(const STLDoubleVector& x) const;
};

// Genz example gen2: product peak [Genz84], in [Owen01], p. 22
class  CGenz2integrand : public CGenzintegrand
{
public:
  CGenz2integrand(const STLDoubleVector& a, const STLDoubleVector& u) : CGenzintegrand(a, u) {}

  double evaluate(const STLDoubleVector& x) const;
};

// Genz example gen3: product peak [Genz84], in [Owen01], p. 22
class  CGenz3integrand : public CGenzintegrand
{
public:
  CGenz3integrand(const STLDoubleVector& a, const STLDoubleVector& u) : CGenzintegrand(a, u) {}

  double evaluate(const STLDoubleVector& x) const;
};

// Genz example gen4: product peak [Genz84], in [Owen01], p. 22
class  CGenz4integrand : public CGenzintegrand
{
public:
  CGenz4integrand(const STLDoubleVector& a, const STLDoubleVector& u) : CGenzintegrand(a, u) {}

  double evaluate(const STLDoubleVector& x) const;
};

//----------------------------------------------------------------------------------------------

// Schmitzberger example sch [Schmitzberger94], p. 23
class  CSchmitzbergerintegrand : public CIntegrand
{
public:
  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { return 1.0; };
};

//----------------------------------------------------------------------------------------------

// Capstick and Keister example ck [CapstickKeister96]
class  CCapstickKeisterintegrand : public CIntegrand
{
public:
	double evaluate(const STLDoubleVector& x) const;

  // Note: analytic I over [0, 1)^k only available for k = 1, 25.
  double I(int k) const;
};

//----------------------------------------------------------------------------------------------

// Hellekalek example hel [Hellekalek98], in [Owen97], [Owen98] (with alpha=1) and [Owen01], p. 17
// This integrand is of full dimension k in the sense of Owen: QMC integration should perform badly.
class  CHellekalekintegrand : public CIntegrand
{
 private:
  // Parameter alpha in {1,2,3}.
  int _alpha;

 public:
   CHellekalekintegrand(int alpha) : _alpha(alpha) {}

  double evaluate(const STLDoubleVector& x) const;
  inline double I(int k) const { return 0.0; };
};

//----------------------------------------------------------------------------------------------


#endif
