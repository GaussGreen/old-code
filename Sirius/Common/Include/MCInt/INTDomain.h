
#ifndef __INTDOMAIN_H__
#define __INTDOMAIN_H__


#include "INTSTL.h"

#include "INTUtilities.h"
#include "smart.h"

//----------------------------------------------------------------------------------------------

#ifdef _S
#undef _S
#endif

// Base class for domains
class  CDomain : public RCObject
{
protected:
	// Acronym
	string _S;
	
public:
	// Dimension k
	int _k;
	
	CDomain(int k) : _k(k) {};
	
	// Domain transformation
	//  virtual STLDoubleVector& psi(STLDoubleVector& x) const = 0;
	virtual STLDoubleVector& psi(STLDoubleVector& x) const
	{
		assert(false);
		static STLDoubleVector y;
		return y;
	}
	
	// Jacobian | det( dpsi_i/dx_j (x) ) |
	//  virtual double J(const STLDoubleVector& x) const = 0;
	virtual double J(const STLDoubleVector& x) const
	{ assert(false);return -1e99;}
	
	// Indicates if the geometry of the domain D requires to plot f(psi(x))*J(x) vs x in [0, 1)^k 
	// (rather than f(x) vs x in D). True by default except for rectilinear domains.
	inline virtual bool transform_f_to_plot(void) const { return true; }
	
	// Domain acronym
	inline const char *acronym(void) const { return _S.c_str(); };
	
	// Report domain
	// virtual void report(void) const = 0; 
	virtual void report(void) const {}; 

};

//----------------------------------------------------------------------------------------------

class  CLineardomain : public CDomain
{
 protected:
  // Constant Jacobian J
  double _J;

 public:
  CLineardomain(int k) : CDomain(k) {};

  inline double J(const STLDoubleVector& x) const { return _J; };
};

//----------------------------------------------------------------------------------------------

class  CParallelepiped : public CLineardomain
{
 private:
  // Linear transformation matrix M
  STLDoubleVectorVector _M;

 public:
  CParallelepiped(int k, const STLDoubleVectorVector& M);

  STLDoubleVector& psi(STLDoubleVector& x) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

class  CRectilineardomain : public CLineardomain
{
 public:
  CRectilineardomain(int k) : CLineardomain(k) {};

  inline double J(const STLDoubleVector& x) const { return _J; };

  bool plottransformed(void) const { return false; }
};

//----------------------------------------------------------------------------------------------

#ifdef _L
#undef _L
#endif

// General cube [L_11, L_12) x ... x [L_k1, L_k2)
class  CGeneralcube : public CRectilineardomain
{
 private:
  // Limits L
  vector< EQSP_CC_STL_PAIR(double, double) > _L;

 public:
  // Constructor parameters k: dimension, L_11, L_12, ..., L_k1, L_k2: list of limits
  CGeneralcube(int k, ...);

  STLDoubleVector& psi(STLDoubleVector& x) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Cube [L_1, L_2)^k
class  CCube : public CRectilineardomain
{
 private:
  // Limits L
  EQSP_CC_STL_PAIR(double,double) _L;

 public:
  // Constructor parameters k: dimension, a, b: common limits
  CCube(int k, double a, double b);

  STLDoubleVector& psi(STLDoubleVector& x) const;

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Unit cube: [0, 1)^k
class  CUnitcube : public CRectilineardomain
{
 public:
  CUnitcube(int k) : CRectilineardomain(k) { _J = 1.0; _S = string("D-[0,1)^") + _k; };

  inline STLDoubleVector& psi(STLDoubleVector& x) const { return x; };

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

#ifdef _C
#undef _C
#endif

class  CEllipsoidaldomain : public CDomain
{
 protected:
  // Center C
  STLDoubleVector& _C;

 public:
  CEllipsoidaldomain(int k, STLDoubleVector& C) : CDomain(k), _C(C) {}

  STLDoubleVector& psi(STLDoubleVector& x) const;

  double J(const STLDoubleVector& x) const;

  // Radii
  virtual double getr(int i) const = 0;
};

//----------------------------------------------------------------------------------------------

// Ellipsoid (x_1 - C_1)^2/r_1^2 + ... + (x_k - C_k)^2/r_k^2 <= 1
class  CEllipsoid : public CEllipsoidaldomain
{
 protected:
  // Radii r_i
  STLDoubleVector& _r;

 public:
  CEllipsoid(int k, STLDoubleVector& r, STLDoubleVector& C);

  inline double getr(int i) const { return _r[i-1]; };

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Sphere (x_1 - C_1)^2 + ... + (x_k - C_k)^2 <= r^2
class  CSphere : public CEllipsoidaldomain
{
 protected:
  // Radius r
  double _r;

 public:
  CSphere(int k, double r, STLDoubleVector& C);

  inline double getr(int i) const { return _r; };

  void report(void) const;
};

//----------------------------------------------------------------------------------------------

// Unit sphere x_1^2 + ... + x_k^2 <= 1
class  CUnitsphere : public CEllipsoidaldomain
{
 public:
  CUnitsphere(int k, STLDoubleVector& C);

  inline double getr(int i) const { return 1.0; };

  void report(void) const;
};

//----------------------------------------------------------------------------------------------



#endif
