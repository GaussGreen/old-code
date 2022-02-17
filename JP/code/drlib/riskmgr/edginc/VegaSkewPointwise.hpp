//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewPointwise.cpp
//
//   Description : Controls calculation of Vega pointwise
//
//   Author      : Mark A Robson
//
//   Date        : 5 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VEGA_SKEW_POINTWISE_HPP
#define EDG_VEGA_SKEW_POINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

class VegaSkewProxyPointwise;

/** Sens Control for vega skew pointwise */
class RISKMGR_DLL VegaSkewPointwise: public VectorShift,
                         public Additive {
public:
    friend class VegaSkewPointwiseHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class VegaSkewPointwiseHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(VegaSkewPointwise* shift) const = 0;

        /** Name matching mechanism with default single name matching on sensName */
        virtual bool sensNameMatches(VegaSkewPointwise* shift, const OutputName& name) const;

        /** Name adding mechanism with default single name adding of sensName */
        virtual void sensAppendName(VegaSkewPointwise* shift, OutputNameArray& namesList) const;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(
            VegaSkewPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaSkewPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class VegaSkewPointwiseHelper;
        static CClassConstSP const TYPE;
        IRestorableShift();
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaSkewPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    VegaSkewPointwise(double     shiftSize);

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    CClassConstSP restorableShiftInterface() const;


    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     VegaSkewPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing vega pointwise and this array is returned. The supplied
        object must implement the VegaSkewPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the VegaSkewPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** Returns spot value of associated asset. Fails if spot has not been
        set */
    double getSpot() const;

    /** Stores the spot value of the associated asset. Allows vol to 
        implement skew tweak */
    void setSpot(double spotValue);

    /** return skew shift for given strike i.e.
        -shift * log(strike/spot)/log(SKEW_NORMALISATION) */
    double skewShift(double strike) const;

    static VegaSkewPointwise* fromProxy(VegaSkewProxyPointwise* proxy);

private:
    /** for reflection */
    VegaSkewPointwise();
    VegaSkewPointwise(const VegaSkewPointwise &rhs);
    VegaSkewPointwise& operator=(const VegaSkewPointwise& rhs);
    //// fields ////
    double spot;  // note: populated by asset during tweak
    bool   spotSet; // indicates if spot value has been set
};

typedef smartConstPtr<VegaSkewPointwise> VegaSkewPointwiseConstSP;
typedef smartPtr<VegaSkewPointwise> VegaSkewPointwiseSP;


DRLIB_END_NAMESPACE

#endif
