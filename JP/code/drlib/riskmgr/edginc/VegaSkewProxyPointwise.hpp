//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewProxyPointwise.hpp
//
//   Description : Controls calculation of fund proxy vega skew
//
//   Author      : Andrew J Swain
//
//   Date        : 12 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef _VEGASKEWPROXYPOINTWISE_HPP
#define _VEGASKEWPROXYPOINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for proxy vega skew pointwise */
class RISKMGR_DLL VegaSkewProxyPointwise: public VectorShift,
                              public Additive {
public:
    friend class VegaSkewProxyPointwiseHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for VEGA_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class VegaSkewProxyPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns true if this object matches the supplied name with 
            respect to the sensitivity */
        virtual bool sensNameMatches(VegaSkewProxyPointwise* shift,
                                     const OutputName&       name) const = 0;

        /** Appends the name(s) of this object with respect to
            the sensitivity to the supplied list. */
        virtual void sensAppendName(VegaSkewProxyPointwise* shift,
                                    OutputNameArray&        namesList) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this vol */
        virtual ExpiryArrayConstSP sensExpiries(
            VegaSkewProxyPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaSkewProxyPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class VegaSkewProxyPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaSkewProxyPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    VegaSkewProxyPointwise(double     shiftSize);

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
     VegaSkewProxyPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing vega pointwise and this array is returned. The supplied
        object must implement the VegaSkewProxyPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the VegaSkewProxyPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    VegaSkewProxyPointwise();
    VegaSkewProxyPointwise(const VegaSkewProxyPointwise &rhs);
    VegaSkewProxyPointwise& operator=(const VegaSkewProxyPointwise& rhs);
};

typedef smartConstPtr<VegaSkewProxyPointwise> VegaSkewProxyPointwiseConstSP;
typedef smartPtr<VegaSkewProxyPointwise> VegaSkewProxyPointwiseSP;


DRLIB_END_NAMESPACE

#endif
