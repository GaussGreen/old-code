
#ifndef __INTALPHA_H__
#define __INTALPHA_H__



#include "INTSTL.h"
#include "INTUtilities.h"
#include "INTSequence.h"
#include "smart.h"


//typedef CC_STL_MAP( int, STLDoubleVector ) AlphaDB;
typedef std::map <int, STLDoubleVector > AlphaDB;

//----------------------------------------------------------------------------------------------


#ifdef _S
#undef _S
#endif

// Base class for alphas
class  CAlpha : public RCObject
{
protected:
	// Buffered vectors alpha for various dimensions k
	AlphaDB _alphaDB;
	
	// Acronym
	string _S;
	
	// Compute and buffer alpha for new dimension
	//  virtual STLDoubleVector& compute(int k) = 0;
	virtual STLDoubleVector& compute(int k)
	{
		assert(false);
		static STLDoubleVector x; 
		return x;
	}
	
	// Update buffered alpha (or insert if no alpha buffered for dimension k)
	virtual STLDoubleVector& update(const STLDoubleVector& alpha, int k = 0);
	
public:
	// Retrieve buffered alpha for dimension k (compute if not available)
	STLDoubleVector& get(STLDoubleVector& alpha, int k);
	
	// Plot distribution of sample points {m*alpha} in [0,1]^k (for k in {1, 2, 3})
	void plotdistribution(int k, int N, const string& filename = "");
	
	// Alpha acronym
	inline const char *acronym(void) const { return _S.c_str(); }
	
	// Report alpha
	virtual void report(int k = 0);
};

//----------------------------------------------------------------------------------------------

// User-defined alpha
class  CUseralpha : public CAlpha
{
 private:
  // Read in alpha for new dimension
  STLDoubleVector& compute(int k);

 public:
  CUseralpha(const STLDoubleVector& alpha) { update(alpha); _S = "a-user"; }

  // Only for user-defined alphas is updating publicly accessible
  STLDoubleVector& update(const STLDoubleVector& alpha, int k = 0) { return CAlpha::update(alpha, k); };

  // Remove buffered alpha for dimension k
  inline void remove(int k) { assert(_alphaDB.erase(k)); };

  void report(int k = 0) { std::cout << "User-defined alpha"; CAlpha::report(k); }
};

//----------------------------------------------------------------------------------------------

// Haselgrove alpha has [Haselgrove61], pp. 333--334
class  CHaselgrovealpha : public CAlpha
{
 private:
  // Decay rate s in {2, 4} of the integrand's Fourier coefficients
  int _s;

  STLDoubleVector& compute(int k);

 public:
  CHaselgrovealpha(int s = 2);

  inline void report(int k = 0) { std::cout << "Haselgrove alpha (s = " << _s << ") "; CAlpha::report(k); }
};

//----------------------------------------------------------------------------------------------

// Cyclotomic alpha cyc: alpha_i = 2*cos(2*Pi*i/(2*k + 3))
class  CCyclotomicalpha : public CAlpha
{
 private:
  STLDoubleVector& compute(int k);

 public:
   CCyclotomicalpha() { _S = "a-cyc"; }

  inline void report(int k = 0) { std::cout << "Cyclotomic alpha"; CAlpha::report(k); };
};

//----------------------------------------------------------------------------------------------

// Baker alpha bak: alpha_i = exp(m*r_i)
class  CBakeralpha : public CAlpha
{
 private:
  // Largest exponent m >= 2 in the sequence, as (r_i) enumerate the rationals in [0,1)
  double _m;

  // Sequence of rationals (r_i) enumerating |Q intersect [0,1)
  CRationalsequence _r;

  STLDoubleVector& compute(int k);

 public:
  CBakeralpha(double m = 3.0) : _m(m) { _S = string("a-bak(") + _m + ")"; }

  inline void report(int k = 0) { std::cout << "Baker alpha (m = " << _m << ")"; CAlpha::report(k); };
};

//----------------------------------------------------------------------------------------------

// Niederreiter alpha nie: alpha_i = 2^(i/(k+1))
class  CNiederreiteralpha : public CAlpha
{
 private:
  STLDoubleVector& compute(int k);

 public:
  CNiederreiteralpha() { _S = "a-nie"; } 

  inline void report(int k = 0) { std::cout << "Niederreiter alpha"; CAlpha::report(k); };
};

//----------------------------------------------------------------------------------------------

// Prime root alpha pr: alpha_i = sqrt(p_i)
class  CPrimerootalpha : public CAlpha
{
 private:
  // Sequence of primes (p_i)
  CPrimesequence _p;

  STLDoubleVector& compute(int k);

 public:
  CPrimerootalpha() { _S = "a-pr"; }

  inline void report(int k = 0) { std::cout << "Prime root alpha"; CAlpha::report(k); };
};

//----------------------------------------------------------------------------------------------


#endif
