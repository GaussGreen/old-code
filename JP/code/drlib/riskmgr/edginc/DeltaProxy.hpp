//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaProxy.hpp
//
//   Description : Fund proxy delta sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 7 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef _DELTAPROXY_HPP
#define _DELTAPROXY_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/Delta.hpp"
#include "edginc/Additive.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for proxy delta - a scalar shift where the derivative is
    calculated via a two sided tweak operation */
class RISKMGR_DLL DeltaProxy: public ScalarShift,
                  virtual public Additive,
                  virtual public IAssetSpotGreek,
                  virtual public ITwoSidedDeriv {
public:
    friend class DeltaProxyHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static string SECOND_ORDER_NAME;
    const static double DEFAULT_SHIFT;
    const static double MINIMUM_SHIFT;

    /** What an object must implement to be tweakable for DeltaProxy */
    class RISKMGR_DLL Shift{
    public:
        friend class DeltaProxyHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();

        /** Returns true if this object matches the supplied name with 
            respect to the DeltaProxy sensitivity */
        virtual bool sensNameMatches(DeltaProxy*       shift,
                                     const OutputName& name) const = 0;

        /** Appends the name(s) of this object with respect to
            the DeltaProxy sensitivity to the supplied list. */
        virtual void sensAppendName(DeltaProxy*      shift,
                                    OutputNameArray& namesList) const = 0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(DeltaProxy* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for DELTA. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class DeltaProxyHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(DeltaProxy* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** calculates delta via 2 sided tweak */
    virtual void calculate(TweakGroup*      tweakGroup,
                           CResults*        results);

    /** Scales delta and gamma numbers in results object */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const;

    /** Adds delta and gamma numbers in results object */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

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
     VegaParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the DeltaProxy.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** constructor with explicit shift size */
    DeltaProxy(double shiftSize);

    /** make a delta shift using own shift or Delta's shift (if available) */
    DeltaSP makeDeltaShift() const;

    /** Returns the shift which has been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const;

private:
    /** for reflection */
    DeltaProxy();
    DeltaProxy(const DeltaProxy &rhs);
    DeltaProxy& operator=(const DeltaProxy& rhs);
};

typedef DeltaProxy CDeltaProxy;
typedef smartConstPtr<DeltaProxy> DeltaProxyConstSP;
typedef smartPtr<DeltaProxy> DeltaProxySP;
typedef DeltaProxyConstSP CDeltaProxyConstSP;
typedef DeltaProxySP CDeltaProxySP;


DRLIB_END_NAMESPACE

#endif
