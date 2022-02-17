//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseFlatIncrementalMappingFunction.hpp
//
//   Description : Piecewise Flat Incremental Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_FLAT_INCREMENTAL_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_FLAT_INCREMENTAL_MAPPING_FUNCTION_HPP

#include "edginc/PiecewiseIncrementalMappingFunction.hpp"
#include "edginc/PiecewiseFlatMappingFunction.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise Flat mapping function */
class MARKET_DLL PiecewiseFlatIncrementalMappingFunction :
    public PiecewiseIncrementalMappingFunction
{
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
        
    /** Destructor */
    virtual ~PiecewiseFlatIncrementalMappingFunction();

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
    PiecewiseFlatIncrementalMappingFunction();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // -------------
    // HIDDEN FIELDS
    // -------------
    /** corresponding mapping function */
    PiecewiseFlatMappingFunctionSP mappingFunction;

};

typedef smartPtr<PiecewiseFlatIncrementalMappingFunction> PiecewiseFlatIncrementalMappingFunctionSP;
typedef const smartPtr<PiecewiseFlatIncrementalMappingFunction> PiecewiseFlatIncrementalMappingFunctionConstSP;


DRLIB_END_NAMESPACE

#endif

