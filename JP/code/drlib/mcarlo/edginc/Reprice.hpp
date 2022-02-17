//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Reprice.hpp
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

#ifndef EDR_REPRICE_HPP
#define EDR_REPRICE_HPP
#include "edginc/SubReprice.hpp"

DRLIB_BEGIN_NAMESPACE
/** Extends ISubReprice interface with additional information needed by the
    IMCPrices object when skipping a path */
class MCARLO_DLL IReprice: public virtual ISubReprice{
public:
    /** Returns the original value (ie on pricing run) for specified path.
        Note this is the real price (ie after any MAX etc) */
    virtual double originalPrice(int path) const = 0;

    virtual ~IReprice();

    class Vanilla;
    class Spread;

    /** Calculates threshold for quick greeks when doing a two sided
        shift (eg delta/gamma/cross gamma). If 
        B = sumProduct(maxDeltaScalingFactor, assets) then this returns
        the largest absolute amount that B will move by under
        the specified shifts */
    static double baskTwoFactorGreekTolerance(
        const DoubleArray&                  maxDeltaScalingFactor,
        const IMultiFactors*                assets,
        const MCPathGenerator*    futurePathGen,
        const ScalarShiftArray&             sens);

    /** Given the a ScalarShift (which contains an OutputName). This
        returns (|shift size|)/(1-|shift size|) *(spot of asset
        i)*(scaleFactors[i]), where asset i is the one which is being
        shifted because it is sensitive to OutputName. (The division by
        1-|shiftSize| is to allow for the supplied market data being
        tweaked already - it is conservative).If there is more than one
        asset sensitive to the OutputName, the sum of the above is used. If
        there are no assets sensitive, 0 is returned.  The SensControl
        must implement the IAssetSpotGreek interface */
    static double calcMaxAssetShiftSize(const IMultiFactors* assets,
                                        const ScalarShift*   delta,
                                        const DoubleArray&   scaleFactors);
    
}; 

typedef refCountPtr<IReprice> IRepriceSP; 

/** For describing payoffs of form MAX(f-k,0). f can be optionally be
    described by a IReprice object */
class MCARLO_DLL IReprice::Vanilla: public virtual IReprice{
private:
    bool                    isStoringMode;
    double                  notional;
    double                  tolerance;
    ISubRepriceSP           repriceForCmpt;
    DoubleArraySP           untweakedPrices;
public:
    ~Vanilla();

    /** Constructor. 'store' should be set to true if the object needs
        to record values */
    Vanilla(int mode, int numPaths, double notional);

    /** Same as Constructor above but treats reprice as a reprices object
        for the f in MAX(f-k,0) (otherwise f is assumed to have zero
        cross gamma eg sum of spots). Note that only a reference is
        taken to repriceForCmpt. The caller must ensure that the right
        methods (eg setForQuickGreeks() etc) are called at the right
        time for this object */
    Vanilla(const ISubRepriceSP& repriceForCmpt, int mode, int numPaths,
            double notional);

    /** Returns a deep copy of this object. Storage data is shallow copied */
    virtual ISubReprice* clone() const;

    /** Sets object in mode ready to pricing for a [1st order] greek
        which is done by a one sided taylor expansion. Use NULL for params */
    virtual void setForOneSidedGreek(IOneSidedGreekParamsSet* params);

    class MCARLO_DLL TwoSidedGreekParamsSet: virtual public ITwoSidedGreekParamsSet,
                                  virtual public IXGammaParamsSet{
    public:
        //// Primitive approach
        TwoSidedGreekParamsSet(double tolerance);
        //// Does the work for you
        TwoSidedGreekParamsSet(
            const DoubleArray&               maxCmptScalingFactor,
            const MCPathGenerator*           futurePathGen,
            const IMultiFactors*             assets, 
            const ScalarShiftArray&          sens);
        
         ~TwoSidedGreekParamsSet();
        friend class Vanilla;
        friend class Spread;
    private:
        double              tolerance;
    };
    /** Sets object in mode ready to pricing for a [1st or 2nd order]
        greek which is done by a two sided taylor expansion. Use
        TwoSidedGreekParamsSet for an instance of
        ITwoSidedGreekParamsSet */
    virtual void setForTwoSidedGreek(ITwoSidedGreekParamsSet* params);

    /** Sets object in mode ready to pricing for cross gamma.  Use
        TwoSidedGreekParamsSet for an instance of
        ITwoSidedGreekParamsSet */
    virtual void setForXGamma(IXGammaParamsSet* params);

    /** Stores the supplied value for the current path on initial
        pricing run and this allows methods below to be
        implemented. The method must be invoked for each simulated
        path and called in order. The price should be before the
        MAX with zero is applied and without any scaling by notional */
    void store(double preMaxPrice);

    /** Returns the repriceForCmpt as supplied in the constructor (may
        be null) */
    virtual ISubRepriceSP getRepriceForCmpt();

    /** Returns true if original value < 0.0 */
    virtual bool firstDerivZero(int path) const;

    /** Returns true if original value < -tolerance. */
    virtual bool firstNumericalDerivZero(int path) const;
   
    /** Returns true if original value < -tolerance or 
        (original value > tolerance and 
        repriceForCmpt.crossNumericalDerivZero() is true) */
    virtual bool crossNumericalDerivZero(int path) const;

    /** Returns the original value (ie on pricing run) for specified path */
    virtual double originalPrice(int path) const;
};

typedef refCountPtr<IReprice::Vanilla> RepriceVanillaSP;

/** For describing payoffs of form MAX(f-k1,0) - MAX(f-k2). f can be
    optionally be described by a IReprice object.
    Note: need to review this class in light of the 'modifier' approach to
    mc. Do we need to support a single strike? Would it be easier not to
    have all the instrument type data below? */
class MCARLO_DLL IReprice::Spread: public virtual IReprice{
private:
    bool                    isStoringMode;
    bool                    isCall;
    double                  notional;
    double                  loStrike;
    double                  hiStrike;
    double                  tolerance;
    DoubleArraySP           basketPrices;
    ISubRepriceSP           repriceForCmpt;
public:
    ~Spread();

    /** Constructor. mode is as per IQuickGreeksCommon::createOrigPrices().
        payoff=(isCall?1: -1)*notional*(MAX(P1,0)-MAX(P2, 0)) with 
        Pi = (isCall?1: -1)*(f-ki). k1 = loStrike, k2=hiStrike. If not doing
        crossGamma then the important part is that the optionality
        lies around k1 and k2 (the overall price being
        irrelevant). The repriceForCmpt is optional and can be
        null. If present it is treated as a reprices object for the f
        in MAX(f-k,0) (otherwise f is assumed to have zero cross gamma
        eg sum of spots). Note that only a reference is taken to
        repriceForCmpt. The caller must ensure that the right methods
        (eg setForQuickGreeks() etc) are called at the right time for
        this object */
    Spread(const ISubRepriceSP& repriceForCmpt, // optional
           int                  mode, 
           int                  numPaths,
           bool                 isCall, 
           double               notional,
           double               loStrike,
           double               hiStrike);

    /** Returns a deep copy of this object. Storage data is shallow copied */
    virtual ISubReprice* clone() const;

    /** Sets object in mode ready to pricing for a [1st order] greek
        which is done by a one sided taylor expansion. Use NULL for params */
    virtual void setForOneSidedGreek(IOneSidedGreekParamsSet* params);

    /** Sets object in mode ready to pricing for a [1st or 2nd order]
        greek which is done by a two sided taylor expansion. Use
        Vanilla::TwoSidedGreekParamsSet for an instance of
        ITwoSidedGreekParamsSet */
    virtual void setForTwoSidedGreek(ITwoSidedGreekParamsSet* params);

    /** Sets object in mode ready to pricing for cross gamma.  Use
        Vanilla::TwoSidedGreekParamsSet for an instance of
        ITwoSidedGreekParamsSet */
    virtual void setForXGamma(IXGammaParamsSet* params);

    /** Stores the supplied value for the current path on initial
        pricing run and this allows methods below to be
        implemented. The method must be invoked for each simulated
        path and called in order. The basketValue supplied should be as
        described in the constructor (ie the variable f) */
    void store(double basketPrice);

    /** Returns the repriceForCmpt as supplied in the constructor (may
        be null) */
    virtual ISubRepriceSP getRepriceForCmpt();

    /** Returns true if path lies out of the money or is above the
        hiStrike */
    virtual bool firstDerivZero(int path) const;

    /** Returns true if path is below loStrike-tolerance or is above
        the hiStrike+tolerance */
    virtual bool firstNumericalDerivZero(int path) const;
   
    /** Returns true if path lies away from loStrike and hiStrike (using
        tolerance) */
    virtual bool crossNumericalDerivZero(int path) const;

    /** Returns the original value (ie on pricing run) for specified path.
        Uses formula specified in constructor */
    virtual double originalPrice(int path) const;
};

typedef refCountPtr<IReprice::Spread> RepriceSpreadSP;


DRLIB_END_NAMESPACE
#endif
