#ifndef EDR_MCPRICES_HPP
#define EDR_MCPRICES_HPP

#include "edginc/config.hpp"
#include "edginc/Reprice.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/*
**** Background info on design of quick greeks/ quick X gamma ****

  Need to allow IMCPrices classes to work together for the purpose of quick
  greeks. Eg think of a MAX(.,0) or a spread at the option level. And
  quick x gamma for a rainbow.
  Clearly need only one place to store actual prices.

  Need methods to reflect the different properties we need to know about:
  We view the RePrices object as describing a mapping f(R^n) -> R.
  (a) is f constant (wrt some predefined variable) ie first deriv is zero.
  The 'derivatives' below are those calculated numerically using a shift of
  size  delta.
  (b) is the first 'derivative' zero (wrt predefined variable)
  (c) is the cross 'derivative' (wrt some pair of predefined variables) zero

  These three rules drive the current functionality. (a) is used for 1st order
  quick greeks where we skip paths that are constant (even though the shifted
  path might occasionally change in value). (b) is used for quick greeks where
  the greek is calculated via a 2 sided shift (eg delta). (c) is used for
  quick x gamma where we skip paths where the cross gamma is zero.
  For (a) and (b) we need to know what the value of (f) is if the condition is
  true. For (c) we can escape this requirement since we are currently reporting
  the change in value.
  This suggests three methods:
  virtual bool firstDerivZero(int path, double& constVal) const;
  virtual bool firstNumericalDerivZero(int path, double& constVal) const;
  virtual bool crossNumericalDerivZero(int path) const;

  Three implementations to consider:
  (1) MAX(S-k,0). Needs to store S for each path. For (a) true if S < k.
  (b) true if S < (k-tolerance). (c) true if S < (k-tolerance) or
  S > (k+tolerance)
  (2) Spread ie MAX(S-k1,0) - MAX(S-k2,0). For (a) true if S < k1 or S > k2
  (b) true if S < (k1-tolerance) or S > (k2+-tolerance)
  (c) true if S < (k1-tolerance) or (S > (k1+tolerance) and S < k2-tolerance)
  or S > (k2+tolerance)
  (3) Rainbow operation: (a) false always (b) false always (c) true if order
  of assets does not change.

  How would we use these instances? Think we need to be able to combine them
  so that for example for SuperRainbow our reprices object for MAX(.,0) would
  hold a reprice object for the rainbow effect. In this scenario we have to
  reconsider the implementations of (a) to (c). (a) and (b) remain unaffected
  via chain rule. For (c) extra condition to still be true:
  the component reprices object returns true for this condition.

  Finally, have a IMCPrices object which contains a std simple prices object and a
  reprices object. Have a setForQuickGreeks(bool for2ndOrderGreek) and
  setForQuickXGamma(void) methods as well as getReprices().

 */

class MCARLO_DLL IGreeks
{
public:
    virtual void reset() = 0; // Seems to be QuickGreeks specific?
    /** Support product-level caching across tweaks */
    virtual int storagePerPath(IMCProduct* product) const = 0;
    virtual void configureCache(const IntArray& changedAssets) = 0;
    /** On pricing run returns MAX(x, 0.0). It should be used for applying
    the 'final point of optionality' within the product. This allows
    QuickGreeks type of IMCPrices not to apply the max when doing first
    order greeks (this may sound strange but see the docs for why) */
    virtual double maxWithZero(double x) const { return Maths::max(x,0.0); }
	virtual ~IGreeks() {}
};

/** base class for IMCPricesGreek -- pricing with Quick Greeks */
class MCARLO_DLL IQuickGreeks : public virtual IGreeks
{
public:
    /** Used to control how the reprices object is used. */
    enum QuickGreekType{
        NONE = 0, // repriceForGreek returns true
        FIRST_ORDER,  // firstDerivZero() on reprice object is used
        SECOND_ORDER, // eg Delta/Gamma. firstNumericalDerivZero() is used
        CROSS_GAMMA   // crossNumericalDerivZero() is used
    };

    /** Returns true if the path, identified by pathIdx, should be
        repriced. If false is returned then there will be no add()
        method called for this path and the IMCPrices object must
        take any appropriate action */
    virtual bool repriceForGreek(int pathIdx) = 0;

    /** Reset this object so that it can be used for the same operation
        again. Normally, a new IMCPrices object is created for each
        pricing run. However, for quick x gamma, it is important to use
        the same one for each pair of assets */
    virtual void reset() = 0;


    /** Returns the suggested mode for handling specified sensitivity.
        This is either NONE, FIRST_ORDER, or SECOND_ORDER only */
    static QuickGreekType defaultMode(const Sensitivity* sens);

    /** Returns true if 'quick greeks' is applicable for this
        sensitivity */
    static bool doQuickGreeks(const Sensitivity* sens);
    


};
DECLARE_REF_COUNT(IQuickGreeks);


/** base class for IMCPrices */
// Should support resetting itself to the initial state
class MCARLO_DLL IMCPrices : public virtual IGreeks
{
public:
        /** Returns a deep copy of this object */
    virtual IMCPrices* clone() const = 0;

    /** adds supplied price to this set of IMCPrices */
    virtual void add(double price) = 0; // FIXME: not a good place
        
    /** Returns the last price stored. Undefined behaviour if no
        prices yet stored */
//    virtual double lastPrice() const = 0;


    /** Returns the averaged result, together with the standard error */
/*    virtual void getResult(double* result,
                           double* resultStdErr) const = 0;*/
    
    virtual ~IMCPrices();

protected:

    /** Ease cloning */
    virtual IMCPrices* emptyConstructor() const = 0;

};
DECLARE_REF_COUNT(IMCPrices);


class IMCProduct;

class MCARLO_DLL IMCPricesSimple : public IMCPrices
{
    public:
    
    /** adds supplied price to this set of IMCPrices. Running totals of
        sum and sum^2 (wrt sample) are kept. */
    virtual void add(double price) = 0;

    /** Returns the averaged result, together with the standard error */
    virtual void getResult(double* result,
                           double* resultStdErr) const = 0;

    /** Returns the last price stored. Undefined behaviour if no
        prices yet stored */
    virtual double lastPrice() const = 0;
    
};

DECLARE_REF_COUNT(IMCPricesSimple);

/** Simplest implementation - repriceForGreek always returns true */
class MCARLO_DLL MCPricesSimple:    public IMCPricesSimple
{
public:
    MCPricesSimple( int             NbIter,
                    int             NbSubSamples);

    /** Returns a deep copy of this object */
    IMCPrices* clone() const;

    /** adds supplied price to this set of IMCPrices. Running totals of
        sum and sum^2 (wrt sample) are kept. */
    virtual void add(double price);

    /** Returns the averaged result, together with the standard error */
    virtual void getResult(double* result,
                           double* resultStdErr) const;

    /** Returns the last price stored. Undefined behaviour if no
        prices yet stored */
    virtual double lastPrice() const;
    
    /** Clears out SumSubPrices and resets iSubSample */
    virtual void reset();
    
    /** Support product-level caching across tweaks */
    virtual int storagePerPath(IMCProduct* product) const;
    virtual void configureCache(const IntArray& changedAssets);
    virtual double maxWithZero(double x) const;
    
    virtual ~MCPricesSimple();
protected:
    virtual IMCPrices* emptyConstructor() const;

    // fields
    int           NbIter;
    int           NbSubSamples;
    int           subSampPos;//counts how far through a sub sample we are
    int           numPerSubSample; // NbIter/NbSubSamples
    double        sampleSum; // current sum for this sample
    double        sumSoFar;
    double        sumSqrSoFar;
    double        lastPriceStored;

};

DECLARE_REF_COUNT(MCPricesSimple);




/** Implementation of IMCPrices interface that supports a Reprice object
    to control the skipping of paths */
class MCARLO_DLL MCPricesGeneral:   public MCPricesSimple,
                                    public virtual IQuickGreeks
{
public:
    /** Constructor: First two params as per MCPricesSimple. The reprice
        parameter is used to control which paths are skipped. Note that
        a reference is taken to the reprice object. The caller must ensure
        that the right data is passed to it and that it is updated
        correctly for each greek. Additionally the caller must ensure
        that setMode() is called appropriately on this object */
    MCPricesGeneral(int               nbIter,
                  int               nbSubSamples,
                  const IRepriceSP& reprice);

    /** Returns a deep copy of this object */
    IMCPrices* clone() const;

    /** Sets what type of sensitivity is being calculated. (Perhaps
        could consider getting MonteCarlo to call this method?) */
    void setMode(QuickGreekType mode);

    /** Returns the reprice object */
    IRepriceSP getReprice();

    /** Uses internal Reprice object to determine whether path can
        be skipped or not. If the path is to be skipped add() is called
        with the relevant value */
    virtual bool repriceForGreek(int pathIdx);

    /** adds supplied price to this set of IMCPrices. Running totals of
        sum and sum^2 (wrt sample) are kept. Note that this must be the
        real price (ie after any MAX etc) */
    virtual void add(double price);

    /** Resets internal variables ready for another pricing run */
    virtual void reset();

    /**  returns MAX(x, 0.0) unless doing first order quick greeks
         in which case x is returned */
    virtual double maxWithZero(double x) const;

    virtual ~MCPricesGeneral();
    
protected:
    virtual IMCPrices* emptyConstructor() const;

private:
    QuickGreekType mode;
    int            path;
    IRepriceSP     reprice;
};
DECLARE_REF_COUNT(MCPricesGeneral);

//////////////////////////////////////////////////////////////////
/** Class to mark abstract information we want to collect from payoff
 * IPayoffEvent will be the mean of communications between a product and customized IMCPricesGeneric class
*/

class MCARLO_DLL IPayoffEvent 
{
public:
    virtual ~IPayoffEvent() {}
};
DECLARE_REF_COUNT(IPayoffEvent);

// class MCARLO_DLL IGenericPricesResult {
// public:
//     virtual ~IGenericPricesResult() {};    
// };
// DECLARE_REF_COUNT(IGenericPricesResult);


/** Provides mechanism to accumulate results that are not a single "double"
 * Each path can be weighted (by PathConfig or product)
 * Weight 1.0 means the usual MC. (weight = probability of strata).
 */

class MCARLO_DLL IGenericPrices 
{
public:
    virtual void add(const IPayoffEvent&, double weight ) = 0;
    /*  Converts accumulated information in a form suitable for Qlib.*/
    virtual IObjectSP getResult() const = 0;
    virtual ~IGenericPrices() {}
};
DECLARE_REF_COUNT(IGenericPrices);


/** Abstract class to handle arbitrary (non-scalar) payoffs
 * One needs to implement IGreeks, IMCPrices and IGenericPrices
 */

class MCARLO_DLL IMCPricesGeneric : public IMCPrices,
                                    public virtual IGenericPrices
{
public:
    virtual void add(double price) 
    {
        QLIB_VERIFY(0, "IMCPricesGeneric doesn't support void add(double price)");
    }
    virtual void add(const IPayoffEvent&, double weight ) = 0; // pass abstract payment

};

DECLARE_REF_COUNT(IMCPricesGeneric);

DRLIB_END_NAMESPACE
#endif // EDR_MCPRICES_HPP
