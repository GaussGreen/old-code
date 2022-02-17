#include "edginc/config.hpp"
#include "edginc/PolynomialBasis.hpp"

DRLIB_BEGIN_NAMESPACE
static long Combinations(long n, long k)
{
    ASSERT(k >= 0);
    return (k==0) ? 1 : (Combinations(n, k-1) * (n+1-k)/k);
}

//PolynomialBasis::PolynomialBasis(int _numVars, int _maxDeg) :
PolynomialBasis::PolynomialBasis(size_t _numVars, size_t _maxDeg) :
    numVars(_numVars), maxDeg(_maxDeg) 
{
    basisSize = Combinations(numVars+maxDeg, numVars); // proof: exercise
}

// input:	a vector of (transformed) fitting variables
// output:	a vector of corresponding values of the basis functions
BasisValArray PolynomialBasis::operator()( const TransformedArray& x ) const
{
    return fill(x, x.size(), maxDeg);
}

// similar to operator() but returns derivatives
// returns df[i] / dx[j] for all basis functions and fitting var.
// output: M[i][j], so each line corresponds to a basis function and each column to a variable
BasisDerivativeMatrix PolynomialBasis::getDerivative( const TransformedArray& x ) const 
{
    throw "Forgot to implement";
}

// returns the number of basis functions
size_t PolynomialBasis::getNumBasisFunctions() const {return basisSize;}

// returns the number of variables (i.e. the required size of TransformedArray)
size_t PolynomialBasis::getNumVariables() const { return numVars;}

    /* Algorithm:
    1. To generate all monomials of degree <= maxDeg on N variables is the same as to generate 
    all monomials of degre == maxDeg with N+1 variables and then "erase" the N+1-th one.
    2. We put the extra variable at the position "-1" and start from the last variable.
    3. The fill() then goes over all available degrees for the current last variable and having fixed it, 
    recursively generates values of all subsequent monomials.
    4. Each thus obtained value is scaled by the current power and appended to the final list.
    */

BasisValArray PolynomialBasis::fill(const TransformedArray& x, size_t active, size_t maxDeg) const
{
    if (active == 0 || maxDeg == 0)
            return BasisValArray(1, 1.);
    BasisValArray result;
    double power = 1.0;

    for(size_t i=0; i <= maxDeg; ++i)
    {
        BasisValArray reduced = fill(x, active -1, maxDeg - i);
        // Now we need to multiply by x[active]^i;
        for(size_t j=0; j < reduced.size(); ++j)
            reduced[j] *= power;
        result.insert(result.end(), reduced.begin(), reduced.end());
        power *= x[active-1];
     }
     return result;
}

DRLIB_END_NAMESPACE
