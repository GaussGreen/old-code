#ifndef _BASIS_FUNCTION_CONTAINER_HPP
#define _BASIS_FUNCTION_CONTAINER_HPP

#include "edginc/config.hpp"
#include "edginc/IFunctionBasis.hpp"
#include <cassert>
DRLIB_BEGIN_NAMESPACE


// a class that describes a basis function used in the fitting procedure. 
// The basis function container may consists of several basis functions.
class   RADAR_DLL PolynomialBasis : public IFunctionBasis
{
	int numVars, maxDeg;
	int basisSize;
public: 
	PolynomialBasis(size_t _numVars, size_t _maxDeg); 
    //PolynomialBasis(int _numVars, int _maxDeg);

	// input:	a vector of (transformed) fitting variables
	// output:	a vector of corresponding values of the basis functions
	virtual BasisValArray operator()( const TransformedArray& x ) const;

	// similar to operator() but returns derivatives
	// returns df[i] / dx[j] for all basis functions and fitting var.
	// output: M[i][j], so each line corresponds to a basis function and each column to a variable
	virtual BasisDerivativeMatrix getDerivative( const TransformedArray& x ) const; 
	// returns the number of basis functions
	virtual size_t getNumBasisFunctions() const; 

	// returns the number of variables (i.e. the required size of TransformedArray)
	virtual size_t getNumVariables() const;

private:

	/* Algorithm:
	   1. To generate all monomials of degree <= maxDeg on N variables is the same as to generate 
	   all monomials of degre == maxDeg with N+1 variables and then "erase" the N+1-th one.
	   2. We put the extra variable at the position "-1" and start from the last variable.
	   3. The fill() then goes over all available degrees for the current last variable and having fixed it, 
          recursively generates values of all subsequent monomials.
	   4. Each thus obtained value is scaled by the current power and appended to the final list.
	*/
	BasisValArray fill(const TransformedArray& x, size_t active, size_t maxDeg) const;
};

DECLARE_REF_COUNT(PolynomialBasis);

DRLIB_END_NAMESPACE

#endif
