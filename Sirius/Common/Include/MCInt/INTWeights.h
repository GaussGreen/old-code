#ifndef __INTWEIGHTS_H__
#define __INTWEIGHTS_H__


#include "INTSTL.h"
#include "INTUtilities.h"
#include "smart.h"


//----------------------------------------------------------------------------------------------

#ifdef _S
#undef _S
#endif

// Base class for weights
class  CWeights : public RCObject
{
public:
	// Rate of convergence O(1/n^r) with r in |N
	// Note: r redundant for arithmetic mean weights, and in {1, 2, 3, 4} for Haselgrove weights
	int _r;
	
	// Acronym
	string _S;
	
	CWeights(int r) : _r(r) {}
	
	// Plot weight distributions (m_i, C*A(m_i)) vs i
	static void plotdistribution(const vector < CWeights*>& w, long N, const string& filename, int res = 20);
	
	// Plot this weight distribution
	void plotdistribution(long N, const string& filename = "", int res = 20);
	
	// Weights acronym
	inline const char *acronym(void) const { return _S.c_str(); }
	
	// Overall constant factor
	//  virtual double C(long N) const = 0;
	virtual double C(long N) const {assert(false);return -1e99;};
	
	// Weighting coefficient for f(m*alpha)
	//  virtual double A(long N, long m) const = 0;      
	virtual double A(long N, long m) const {assert(false);return -1e99;};      
	
	virtual int AisConstant(){return 0;};
	
	// Rescaled number of sampling points (to ensure f is evaluated N times independently of weights)
	// virtual long Neff(long N) const = 0;
	virtual long Neff(long N) const  
	{
		assert(false);
		return 0L;
	}
	
	// Range {N_min, ..., N_max} of m, where Nrange is the pair (N_min, N_max)
	//  virtual EQSP_CC_STL_PAIR(long, long)& Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const = 0;
	virtual EQSP_CC_STL_PAIR(long, long)& Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const 
	{
		assert(false);
		static EQSP_CC_STL_PAIR(long, long) x;
		return x;
	}
	
	// Shift applied to components of x = m*alpha (0 except for Haselgrove weights)
	inline virtual double xshift(void) const { return 0; }
	
	// Report weights
	virtual void report(void) const { std::cout << ": r = " << _r; }
};

//----------------------------------------------------------------------------------------------

class  CArithmeticmeanweights : public CWeights
{
 public:
  CArithmeticmeanweights() : CWeights(1) { _S = "w-am"; }

  inline double C(long N) const { return 1.0/(double)N; }

  inline double A(long N, long m) const { return 1.0; }
  int AisConstant(){return 1;};


  inline long Neff(long N) const { return N; }

  inline EQSP_CC_STL_PAIR(long, long)& Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const{ Nrange = make_pair((long)1, N); return Nrange; }
  
  inline void report(void) const { std::cout << "arithmetic mean weights"; CWeights::report(); }
};

//----------------------------------------------------------------------------------------------

class  CHaselgroveweights : public CWeights
{
 private:
  inline double Ar(long N, long m) const { return binomial(N - (long)fabs(m) + (long)(_r - 1), (long)(_r - 1)); }

 public:
  CHaselgroveweights(int r);

  double C(long N) const;

  double A(long N, long m) const;

  long Neff(long N) const;

  EQSP_CC_STL_PAIR(long, long)& Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const;
  inline void report(void) const { std::cout << "Haselgrove weights"; CWeights::report(); }
};

//----------------------------------------------------------------------------------------------

class  CNiederreiterweights : public CWeights
{
 public:
  CNiederreiterweights(int r);

  inline double C(long N) const { return 1.0/pow(N + 1, _r); }

  double A(long N, long m) const;

  inline long Neff(long N) const { return (double)(N - 1)/(double)_r; }

  EQSP_CC_STL_PAIR(long, long)& Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const;
  inline void report(void) const { std::cout << "Niederreiter weights"; CWeights::report(); }
};

//----------------------------------------------------------------------------------------------

class  CSugiharaMurotaweights : public CWeights
{
 private:
  double _c;

  inline double w(double x) const { return _c * pow(x * (1.0 - x), _r - 1); }

 public:
  CSugiharaMurotaweights(int r);
       
  inline double C(long N) const { return 1.0/(double)(N + 1); }

  inline double A(long N, long m) const { return w((double)m/(double)(N + 1)); }


  inline long Neff(long N) const { return N - 1; }

  inline EQSP_CC_STL_PAIR(long, long) & Nrange(EQSP_CC_STL_PAIR(long, long)& Nrange, long N) const 
		{ Nrange = make_pair((long)0, N); return Nrange; }

  inline void report(void) const { std::cout << "Sugihara and Murota weights"; CWeights::report(); }
};

//----------------------------------------------------------------------------------------------


#endif
