//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseFlatMappingFunction.hpp
//
//   Description : Piecewise Flat Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_FLAT_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_FLAT_MAPPING_FUNCTION_HPP

#include "edginc/PiecewiseMappingFunction.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise Flat mapping function */
class MARKET_DLL PiecewiseFlatMappingFunction :
    public PiecewiseMappingFunction
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Constructor from the fields */
    PiecewiseFlatMappingFunction(CDoubleArraySP x, CDoubleArraySP y);

    /** Destructor */
    virtual ~PiecewiseFlatMappingFunction();

    /** The mapping function */
    virtual double map(double x) const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Only build instances of that class using reflection */
    PiecewiseFlatMappingFunction();
    
    /** Default constructor */
    static IObject* defaultConstructor();
};

typedef smartPtr<PiecewiseFlatMappingFunction> PiecewiseFlatMappingFunctionSP;
typedef const smartPtr<PiecewiseFlatMappingFunction> PiecewiseFlatMappingFunctionConstSP;


DRLIB_END_NAMESPACE

#endif

