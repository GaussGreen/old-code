//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseLinearMappingFunction.hpp
//
//   Description : Piecewise Linear Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_LINEAR_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_LINEAR_MAPPING_FUNCTION_HPP

#include "edginc/PiecewiseMappingFunction.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise linear mapping function */
class MARKET_DLL PiecewiseLinearMappingFunction :
    public PiecewiseMappingFunction
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Constructor from the fields */
    PiecewiseLinearMappingFunction(CDoubleArraySP x, CDoubleArraySP y);

    /** Destructor */
    virtual ~PiecewiseLinearMappingFunction();

    /** The mapping function */
    virtual double map(double x) const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    PiecewiseLinearMappingFunction();
        
    /** Default constructor */
    static IObject* defaultConstructor();
};

typedef smartPtr<PiecewiseLinearMappingFunction> PiecewiseLinearMappingFunctionSP;
typedef const smartPtr<PiecewiseLinearMappingFunction> PiecewiseLinearMappingFunctionConstSP;


DRLIB_END_NAMESPACE

#endif

