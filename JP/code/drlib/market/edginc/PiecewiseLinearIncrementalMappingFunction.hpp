//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseLinearIncrementalMappingFunction.hpp
//
//   Description : Piecewise Linear Incremental Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_LINEAR_INCREMENTAL_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_LINEAR_INCREMENTAL_MAPPING_FUNCTION_HPP

#include "edginc/PiecewiseIncrementalMappingFunction.hpp"
#include "edginc/PiecewiseLinearMappingFunction.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise linear mapping function */
class MARKET_DLL PiecewiseLinearIncrementalMappingFunction :
    public PiecewiseIncrementalMappingFunction
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
        
    /** Destructor */
    virtual ~PiecewiseLinearIncrementalMappingFunction();

    /** The mapping function */
    virtual double map(double x) const;

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Called immediately after object modified */
    virtual void fieldsUpdated(const CFieldArray& fields);
    
    /* Utilities to get the input arrays and fields */
    virtual CDoubleArrayConstSP getX() const;
    virtual CDoubleArrayConstSP getY() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Only build instances of that class using reflection */
    PiecewiseLinearIncrementalMappingFunction();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // -------------
    // HIDDEN FIELDS
    // -------------
    /** corresponding mapping function */
    PiecewiseLinearMappingFunctionSP mappingFunction;
};

typedef smartPtr<PiecewiseLinearIncrementalMappingFunction> PiecewiseLinearIncrementalMappingFunctionSP;
typedef const smartPtr<PiecewiseLinearIncrementalMappingFunction> PiecewiseLinearIncrementalMappingFunctionConstSP;


DRLIB_END_NAMESPACE

#endif

