//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseMappingFunction.hpp
//
//   Description : Piecewise Mapping function interface
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PIECEWISE_MAPPING_FUNCTION_HPP
#define EDR_PIECEWISE_MAPPING_FUNCTION_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/SqueezeParallelTweak.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Piecewise Mapping function interface */
class MARKET_DLL PiecewiseMappingFunction :
    public MarketObject,
    public virtual IMappingFunction,
    public virtual TweakableWith<SqueezeParallelTwk>,
    public virtual QuasiContractualBaseCorrelation::IShift,
    public virtual Calibrator::IAdjustable
{

public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~PiecewiseMappingFunction();

    /** Set a name for this function (useful for tweaking) */
    virtual void setName(string name);

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    // Squeeze parallel tweak (sensitivity) support
    virtual string sensName(SqueezeParallelTwk* shift) const;
    // Squeeze parallel tweak (sensitivity) support
    virtual bool sensShift(SqueezeParallelTwk* shift);

    // Support to flatten the yPoints
    virtual bool sensShift(QuasiContractualBaseCorrelation* shift);

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const;

    /* Utilities to get the input arrays */
    virtual CDoubleArrayConstSP getX() const;
    virtual CDoubleArrayConstSP getY() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

protected:
    /** Only build instances of that class using reflection */
    PiecewiseMappingFunction(CClassConstSP clazz = TYPE);
        
    // ----------------
    // MANDATORY FIELDS
    // ----------------

    /** Array of 'x' points */
    CDoubleArraySP xPoints;
    
    /** Array of 'y'=map(x) points */
    CDoubleArraySP yPoints;
    
    // -------------
    // HIDDEN FIELDS
    // -------------
    
    /** Name (useful for tweaking)*/
    string name;
};

typedef smartPtr<PiecewiseMappingFunction> PiecewiseMappingFunctionSP;
typedef smartConstPtr<PiecewiseMappingFunction> PiecewiseMappingFunctionConstSP;

DRLIB_END_NAMESPACE

#endif

