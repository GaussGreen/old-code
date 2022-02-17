//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootTimeVegaProxy.hpp
//
//   Description : Controls calculation of proxy Root TIme Vega
//
//   Author      : Andrew J Swain
//
//   Date        : 13 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef _ROOTTIMEVEGAPROXY_HPP
#define _ROOTTIMEVEGAPROXY_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for proxy root time vega */
class RISKMGR_DLL RootTimeVegaProxy: public ScalarShift,
                         public Additive {
public:
    friend class RootTimeVegaProxyHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for VEGA_PARALLEL */
    class RISKMGR_DLL IShift{
    public:
        friend class RootTimeVegaProxyHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns true if this object matches the supplied name with 
            respect to the sensitivity */
        virtual bool sensNameMatches(RootTimeVegaProxy* shift,
                                     const OutputName&  name) const = 0;

        /** Appends the name(s) of this object with respect to
            the sensitivity to the supplied list. */
        virtual void sensAppendName(RootTimeVegaProxy* shift,
                                    OutputNameArray&   namesList) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(RootTimeVegaProxy* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class RootTimeVegaProxyHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(RootTimeVegaProxy* shift) = 0;
    };

    /** constructor with explicit shift size */
    RootTimeVegaProxy(double     shiftSize);

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
     RootTimeVegaProxy.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** @param param1 The object to shift. The
     object must implement the RootTimeVegaProxy.RestorableShift interface */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    RootTimeVegaProxy();
    RootTimeVegaProxy(const RootTimeVegaProxy &rhs);
    RootTimeVegaProxy& operator=(const RootTimeVegaProxy& rhs);
};

typedef smartConstPtr<RootTimeVegaProxy> RootTimeVegaProxyConstSP;
typedef smartPtr<RootTimeVegaProxy> RootTimeVegaProxySP;

DRLIB_END_NAMESPACE

#endif
