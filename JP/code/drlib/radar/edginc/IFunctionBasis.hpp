#ifndef _IBASIS_FUNCTION_CONTAINER_HPP
#define _IBASIS_FUNCTION_CONTAINER_HPP

#include "edginc/config.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/RadarRepUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// a class that describes a basis function used in the fitting procedure. 
// The basis function container may consists of several basis functions.
class  RADAR_DLL IFunctionBasis {
public:
    // input:	a vector of (transformed) fitting variables
    // output:	a vector of corresponding values of the basis functions
    virtual BasisValArray operator()( const TransformedArray& x ) const = 0;

    // similar to operator() but returns derivatives
    // returns df[i] / dx[j] for all basis functions and fitting var.
    // output: M[i][j], so each line corresponds to a basis function and each column to a variable
    virtual BasisDerivativeMatrix getDerivative( const TransformedArray& x ) const = 0; 
    // returns the number of basis functions
    virtual size_t getNumBasisFunctions() const = 0;

    // returns the number of variables (i.e. the required size of TransformedArray)
    virtual size_t getNumVariables() const = 0;
    virtual ~IFunctionBasis() {}
};

DECLARE_REF_COUNT( IFunctionBasis);

DRLIB_END_NAMESPACE
#endif
