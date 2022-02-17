//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaProxyParallel.hpp
//
//   Description : Fund proxy vega sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 8 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef _VEGAPROXYPARALLEL_HPP
#define _VEGAPROXYPARALLEL_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for proxy vega parallel */
class RISKMGR_DLL VegaProxyParallel: public ScalarShift,
                         public virtual Additive {
public:
    friend class VegaProxyParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for VEGA_PARALLEL */
    class RISKMGR_DLL Shift{
    public:
        friend class VegaProxyParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~Shift();

        /** Returns true if this object matches the supplied name with 
            respect to the  VegaProxyParallelsensitivity */
        virtual bool sensNameMatches(VegaProxyParallel* shift,
                                     const OutputName&  name) const = 0;

        /** Appends the name(s) of this object with respect to
            the DeltaProxy sensitivity to the supplied list. */
        virtual void sensAppendName(VegaProxyParallel* shift,
                                    OutputNameArray&   namesList) const = 0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaProxyParallel* shift) = 0;
    };
    typedef Shift IShift;

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL RestorableShift: public virtual Shift{
    public:
        friend class VegaProxyParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~RestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaProxyParallel* shift) = 0;
    };
    typedef RestorableShift IRestorableShift;

    /** constructor with explicit shift size */
    VegaProxyParallel(double shiftSize);

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
        for this sensitivity.  The object must implement the
        VegaProxyParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     VegaProxyParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the VegaProxyParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    VegaProxyParallel();
    VegaProxyParallel(const VegaProxyParallel &rhs);
    VegaProxyParallel& operator=(const VegaProxyParallel& rhs);
};

typedef VegaProxyParallel CVegaProxyParallel;
typedef smartConstPtr<VegaProxyParallel> VegaProxyParallelConstSP;
typedef smartPtr<VegaProxyParallel> VegaProxyParallelSP;
typedef VegaProxyParallelConstSP CVegaProxyParallelConstSP;
typedef VegaProxyParallelSP CVegaProxyParallelSP;

DRLIB_END_NAMESPACE

#endif
