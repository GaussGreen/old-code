//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids
//
//   Filename    : QuadraticProg.hpp
//
//   Description : Wrapper around IMSL_D_QUADRATIC_PROG
//
//   Date        : 4 Oct 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/AtomicArray.hpp"
#include "edginc/IMSLException.hpp"


#ifndef EDR_QUADRATIC_PROG_HPP
#define EDR_QUADRATIC_PROG_HPP

DRLIB_BEGIN_NAMESPACE


class UTIL_DLL QuadraticProg
{
public:


	/* returns array with n components containing the solution */
	static DoubleArraySP minimize(
		int numLinConstraints,				// Total number of linear constraints (m)
		int numVariables,					// number of variables (n)
		int numLinEqualityConstraints,		// number of linear equality constraints (meq)
		DoubleArray& constraintFuncs,		// flattened array of size m x n containing equality constraints in first meq rows followed by inequality constraints
		DoubleArray& constraintVals,		// array with m components containing RHS of linear constraints (b)
		DoubleArray& linCoeffs,				// array with n components containing the coeffs of linear term of objective function (g)
		DoubleArray& hessian				// array (m x n) containing Hessian matix of objective function (h)
		);	
 



};



DRLIB_END_NAMESPACE

#endif

