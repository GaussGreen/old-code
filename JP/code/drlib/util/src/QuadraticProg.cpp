//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids
//
//   Filename    : QuadraticProg.cpp
//
//   Description : Wrapper around IMSL_D_QUADRATIC_PROG
//
//   Date        : 4 Oct 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QuadraticProg.hpp"


DRLIB_BEGIN_NAMESPACE

/** minimization using imsl function. Throws IMSLException on failure */
DoubleArraySP QuadraticProg::minimize(
		int numLinConstraints,				// Total number of linear constraints (m)
		int numVariables,					// number of variables (n)
		int numLinEqualityConstraints,		// number of linear equality constraints (meq)
		DoubleArray& constraintFuncs,		// flattened array of size m x n containing equality constraints in first meq rows followed by inequality constraints
		DoubleArray& constraintVals,		// array with m components containing RHS of linear constraints (b)
		DoubleArray& linCoeffs,				// array with n components containing the coeffs of linear term of objective function (g)
		DoubleArray& hessian				// array (m x n) containing Hessian matix of objective function (h)
		) 
{
	static const string method = "QuadraticProg::minimize";

	
	if(numVariables<1) throw ModelException("numVariables < 1");
    
/*  
 *  Adrian: this check should be taken out
 *  int numGridPoints;
    if (numVariables % 2 != 0)
	{
		throw ModelException("numVariables is not of form 2N-2", method);
	}
    else
	{
        numGridPoints = numVariables / 2 + 1;
	}
*/
	if(numLinEqualityConstraints<1)
	{
		throw ModelException("numLinEqualityConstraints < 1",method);
	}
	if(numLinConstraints<numLinEqualityConstraints)
	{
		throw ModelException("numLinConstraints < numLinEqualityConstraints",method);
	}
	if(constraintFuncs.size() != numLinConstraints*numVariables) 
	{
		throw ModelException("constraintFuncs size is not numLinConstraints*numVariables",method);
	}
	if(constraintVals.size() != numLinConstraints)
	{
		throw ModelException("constraintVals size not same as numLinConstraints",method);
	}
	if(linCoeffs.size() != numVariables)
	{
		throw ModelException("linCoeffs size is not numVariables",method);
	}
	if(hessian.size() !=numVariables*numVariables)
	{
		throw ModelException("hessian is not size numVariables*numVariables",method);
	}
	

	/** result array */
	DoubleArraySP outx(new DoubleArray(numVariables));

	imsl_d_quadratic_prog(
		numLinConstraints,
		numVariables,
		numLinEqualityConstraints,
		&constraintFuncs[0], 
		&constraintVals[0], 
		&linCoeffs[0], 
		&hessian[0],
		IMSL_RETURN_USER, &(*outx)[0],
		0);

	// throw IMSLException if quadratic prog fails
	IMSLException::throwIMSLExceptionIfError();

	return outx;
		
}							


DRLIB_END_NAMESPACE
