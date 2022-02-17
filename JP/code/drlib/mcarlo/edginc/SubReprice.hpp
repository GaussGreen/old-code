//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SubReprice.hpp
//
//   Description : Captures whether paths in monte carlo simulation can be 
//                 skipped when calculating sensitivities
//
//   Author      : Mark A Robson
//
//   Date        : 29 November 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SUB_REPRICE_HPP
#define EDR_SUB_REPRICE_HPP

#include "edginc/smartPtr.hpp"
#include "edginc/IndexedPerfList.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/ScalarShift.hpp"

DRLIB_BEGIN_NAMESPACE
/** Captures whether paths in monte carlo simulation can be skipped
    when calculating sensitivities. Each ISubReprice object is describing
    a function f:R^n -> R. */
class MCARLO_DLL ISubReprice{
public:
    /** Returns a deep copy of this object. Although there should be no
        need to deep copy read only data (eg saved by store methods). 
        So clients must not clone this object and try to store different
        data in each */
    virtual ISubReprice* clone() const = 0;

    /** Returns true if the first derivative is zero for the specified path.
        Here the derivative is the mathematical definition. */
    virtual bool firstDerivZero(int path) const = 0;

    /** Returns true if the first 'numerical derivative' is zero for
        the specified path. Here 'numerical derivative' means
        calculating a derivative using a Taylor expansion using some
        previously specified shift. */
    virtual bool firstNumericalDerivZero(int path) const = 0;

    /** Returns true if the 'cross' 'numerical derivative' is zero for
        the specified path. Here 'numerical derivative' has the same
        meaning as in firstNumericalDerivZero(). 'Cross' derivative
        means the second order derivatives wrt to different
        variables. */
    virtual bool crossNumericalDerivZero(int path) const = 0;
    
    virtual ~ISubReprice();

    /*** Captures parameters needed for setForOneSidedGreek
         method. This allows derived classes to use different data whilst the
         'set' method is still accessed through a generic method */
    class MCARLO_DLL IOneSidedGreekParamsSet{
    public:
        virtual ~IOneSidedGreekParamsSet();
    };
    /** Sets object in mode ready to pricing for a [1st order] greek
        which is done by a one sided taylor expansion */
    virtual void setForOneSidedGreek(IOneSidedGreekParamsSet* params)= 0;

    /*** Captures parameters needed for setForTwoSidedGreek
         method. This allows derived classes to use different data whilst
         the 'set' method is still accessed through a generic method */
    class MCARLO_DLL ITwoSidedGreekParamsSet{
    public:
        virtual ~ITwoSidedGreekParamsSet();
    };
    /** Sets object in mode ready to pricing for a [1st or 2nd order] greek
        which is done by a two sided taylor expansion */
    virtual void setForTwoSidedGreek(ITwoSidedGreekParamsSet* params)= 0;

    /*** Captures parameters needed for setForXGamma method. This
         allows derived classes to use different data whilst the 'set'
         method is still accessed through a generic method */
    class MCARLO_DLL IXGammaParamsSet{
    public:
        virtual ~IXGammaParamsSet();
    };
    /** Sets object in mode ready to pricing for cross gamma */
    virtual void setForXGamma(IXGammaParamsSet* params)= 0;

    /** If the class supports 'chaining' then this returns the
        repriceForCmpt as supplied in that class's constructor (may be
        null) */
    virtual refCountPtr<ISubReprice> getRepriceForCmpt() = 0;

    class Rainbow;
}; 

typedef refCountPtr<ISubReprice> ISubRepriceSP;

class MCPathGenerator;
class IMultiFactors;
class IMCProduct;


/** For describing rainbow type functions ie sorted sum */
class MCARLO_DLL ISubReprice::Rainbow: public virtual ISubReprice{
private:
    bool                isStoringMode; // true: record values
    DoubleArrayArraySP  perfDiff;
    int                 asset1; // index into perfDiff
    int                 asset2; // index into perfDiff
    double              perf1Tol;  // tolerance for asset1
    double              perf2Tol;  // tolerance for asset1

public:
    /** Constructor. mode is the bit-wise int passed to createOrigPrices()
        method */
    Rainbow(int mode, int numPaths);
            
    /** Returns a deep copy of this object. Storage data is shallow copied */
    virtual ISubReprice* clone() const;

    /** Sets object in mode ready to pricing for a [1st order] greek
        which is done by a one sided taylor expansion. Use NULL for params */
    virtual void setForOneSidedGreek(IOneSidedGreekParamsSet* params);

    /** Sets object in mode ready to pricing for a [1st or 2nd order] greek
        which is done by a two sided taylor expansion. Use NULL for params */
    virtual void setForTwoSidedGreek(ITwoSidedGreekParamsSet* params);

    class MCARLO_DLL XGammaParamsSet: virtual public IXGammaParamsSet{
    public:
        //// Primitive approach
        XGammaParamsSet(int             asset1,
                        double          perf1Tol,
                        int             asset2,
                        double          perf2Tol);
        //// Does the work for you
        XGammaParamsSet(
            const DoubleArray&               maxCmptScalingFactor,
            const MCPathGenerator*           futurePathGen,
            const IMultiFactors*             assets, 
            const ScalarShiftArray&          sens);
        
         ~XGammaParamsSet();
        friend class Rainbow;
    private:
        int                 asset1; // index into perfDiff
        int                 asset2; // index into perfDiff
        double              perf1Tol;  // tolerance for asset1
        double              perf2Tol;  // tolerance for asset1
    };
    /** Sets object in mode ready to pricing for cross gamma. Use
        XGammaParamsSet for an instance of IXGammaParamsSet */
    virtual void setForXGamma(IXGammaParamsSet* params);

    /** Stores the necessary data for the current path on initial
        pricing run and thus allows methods below to be
        implemented. The method must be invoked for each simulated
        path and called in order. The IndexedPerfList must have already
        been sorted */
    void store(const IndexedPerfList&           indexPerfs,
               const IMCProduct*                 product,
               const MCPathGenerator* pathGen);

    /** For best or worst of payoffs. Uses the supplied information to
        store data for the current path on initial pricing run and this
        allows methods below to be implemented. The method must be invoked
        for each simulated path and called in order. The performances must not
        be sorted */
    void store(const DoubleArray&               performances,
               int                              pickedAsset,
               const IMCProduct*                 product,
               const MCPathGenerator* pathGen);

    /** Returns false */
    virtual bool firstDerivZero(int path) const;

    /** Returns false */
    virtual bool firstNumericalDerivZero(int path) const;

    /** Returns true if order of assets in rainbow will not change under
        tweak */
    virtual bool crossNumericalDerivZero(int path) const;
        
    /** Returns the largest absolute value in a double array */
    static double largestAbsWeight(const DoubleArray& dbles);
 
    /** Returns null SP */
    virtual refCountPtr<ISubReprice> getRepriceForCmpt();

    /** Calculates tolerances needed for Rainbow under cross gamma. It
        returns the id's of the two assets being shifted together with
        the maximum size that they can move (scaled by
        maxCmptScalingFactor eg 'participation' but note rainbow
        weights should not be included in maxCmptScalingFactor) */
    static void crossGammaTolerances(
        const DoubleArray&               maxCmptScalingFactor,
        const MCPathGenerator*           futurePathGen,
        const IMultiFactors*             assets,    
        const ScalarShiftArray&          sens,
        int                   assetIDs[2],    // which assets are being tweaked
        double                tolerances[2]); // tolerances for these assets

#if 0
    /* Removed for now as we need a path generator specific method to
       do this calculation. The implementation here is ok for log normal
       but not implied */
    /** Calculates the largest 'derivative' of the current simulated
        path wrt to spot price (of specified asset) across the
        different simulation dates. If it returns a value of lambda
        then the value of spot on a simulated path cannot move more
        than lambda * delta S when the initial spot S is moved by
        delta S. */
    static double pathWiseMaxDriftProduct(
        const MCPathGenerator* futurePathGen,
        int                              iAsset,
        int                              iPath);
#endif
private:
#if 0
    static void divideByPathWiseMaxDriftProduct(
        DoubleArray&                     toScale,
        const MCPathGenerator* futurePathGen,
        int                              iPath);
#endif

};
typedef refCountPtr<ISubReprice::Rainbow> SubRepriceRainbowSP;

DRLIB_END_NAMESPACE
#endif
