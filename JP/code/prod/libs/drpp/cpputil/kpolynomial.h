// Polynomial.h: interface for the Polynomial class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POLYNOMIAL_H__C2FDE8CE_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
/**@#-*/
#define AFX_POLYNOMIAL_H__C2FDE8CE_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_
/**@#+*/

#include "kplatdep.h"
#include <vector>

// Simple polynomial class
/** Simple polynomial class */
template <class T>
class KPolynomial {	
public:
	/** Creates a polynomial from a vector of coefficients */
	KPolynomial (const KVector(double)& c) : _c(c) {};
	// if poly is c0 + c1 x + c2 x^2 ..., enter coeffs as c0 c1 c2 ...

	/** Evaluates the polynomial for a given value */
	T operator()(const T& a) const
	{
		T ans = 0;
		for (int i = _c.size()-1; i>=0; i--) {
			ans *= a;
			ans += _c[i];
		}
		return ans;
	}
private:
	KVector(double) _c;
};


#endif // !defined(AFX_POLYNOMIAL_H__C2FDE8CE_69C3_11D2_97E7_00C04FD8EB9A__INCLUDED_)
