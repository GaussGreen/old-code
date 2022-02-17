//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseIncrementalMappingFunction.hpp
//
//   Description : Piecewise Incremental Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_INCREMENTAL_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_INCREMENTAL_MAPPING_FUNCTION_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/PiecewiseMappingFunction.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/SqueezeParallelTweak.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise Mapping function interface */
class MARKET_DLL PiecewiseIncrementalMappingFunction :
    public PiecewiseMappingFunction,
    public virtual TweakableWith<SqueezeParallelTwk>,
    public virtual Calibrator::IAdjustable
{

public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~PiecewiseIncrementalMappingFunction();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    // Squeeze parallel tweak (sensitivity) support
    virtual string sensName(SqueezeParallelTwk* shift) const;
    // Squeeze parallel tweak (sensitivity) support
    virtual bool sensShift(SqueezeParallelTwk* shift);

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const;

    double getBasePoint() const;


private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    /** Only build instances of that class using reflection */
    PiecewiseIncrementalMappingFunction(CClassConstSP clazz = TYPE);

    /**
     * Utility function to compute the xPoints and yPoints of a PiecewiseMappingFunction from
     * the base points, increments and values of a PiecewiseIncrementalMappingFunction
     * */
    void computeXandY(CDoubleArraySP x, CDoubleArraySP y) const;

    // ----------------
    // MANDATORY FIELDS
    // ----------------

    /** initial x point */
    double basePoint;
};

typedef smartPtr<PiecewiseIncrementalMappingFunction> PiecewiseIncrementalMappingFunctionSP;

DRLIB_END_NAMESPACE

#endif

