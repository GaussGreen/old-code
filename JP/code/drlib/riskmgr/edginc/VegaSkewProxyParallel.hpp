//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaSkewProxyParallel.hpp
//
//   Description : Controls calculation of fund proxy vega skew
//
//   Author      : Andrew J Swain
//
//   Date        : 12 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef _VEGASKEWPROXYPARALLEL_HPP
#define _VEGASKEWPROXYPARALLEL_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for proxy vega skew parallel */
class RISKMGR_DLL VegaSkewProxyParallel: public ScalarShift,
                             public Additive {
public:
    friend class VegaSkewProxyParallelHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for VEGA_SKEW_PARALLEL */
    class RISKMGR_DLL IShift{
    public:
        friend class VegaSkewProxyParallelHelper;
        static CClassConstSP const TYPE;

        virtual ~IShift();
        /** Returns true if this object matches the supplied name with 
            respect to the sensitivity */
        virtual bool sensNameMatches(VegaSkewProxyParallel* shift,
                                     const OutputName&      name) const = 0;

        /** Appends the name(s) of this object with respect to
            the sensitivity to the supplied list. */
        virtual void sensAppendName(VegaSkewProxyParallel* shift,
                                    OutputNameArray&       namesList) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VegaSkewProxyParallel* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class VegaSkewProxyParallelHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(VegaSkewProxyParallel* shift) = 0;
    };

    /** constructor with explicit shift size */
    VegaSkewProxyParallel(double     shiftSize);

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
     VegaSkewProxyParallel.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the VegaSkewProxyParallel.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    VegaSkewProxyParallel();
    VegaSkewProxyParallel(const VegaSkewProxyParallel &rhs);
    VegaSkewProxyParallel& operator=(const VegaSkewProxyParallel& rhs);
};

typedef smartConstPtr<VegaSkewProxyParallel> VegaSkewProxyParallelConstSP;
typedef smartPtr<VegaSkewProxyParallel> VegaSkewProxyParallelSP;

DRLIB_END_NAMESPACE

#endif
