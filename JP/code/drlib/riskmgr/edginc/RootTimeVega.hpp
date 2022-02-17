//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootTimeVega.hpp
//
//   Description : Controls calculation of Root TIme Vega
//
//   Author      : Mark A Robson
//
//   Date        : 7 March 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ROOT_TIME_VEGA_HPP
#define EDG_ROOT_TIME_VEGA_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for root time vega */
class RISKMGR_DLL RootTimeVega: public ScalarShift,
                    public Additive {
public:
    friend class RootTimeVegaHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_PARALLEL */
    class RISKMGR_DLL IShift{
    public:
        friend class RootTimeVegaHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(RootTimeVega* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(RootTimeVega* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class RootTimeVegaHelper;
        static CClassConstSP const TYPE;
        IRestorableShift();
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(RootTimeVega* shift) = 0;
    };

    /** constructor with explicit shift size */
    RootTimeVega(double     shiftSize);

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
     RootTimeVega.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param obj The object to query
     * @return Whether the given object implements RootTimeVega.RestorableShift
     */
    virtual bool restorableShift(IObjectConstSP& obj) const;

    /**
     *
     * @param param1 The object to shift. The
     object must implement the RootTimeVega.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** return shift/sqrt(time) */
    double rtVegaShift(double years) const;

private:
    /** for reflection */
    RootTimeVega();
    RootTimeVega(const RootTimeVega &rhs);
    RootTimeVega& operator=(const RootTimeVega& rhs);
};

typedef smartConstPtr<RootTimeVega> RootTimeVegaConstSP;
typedef smartPtr<RootTimeVega> RootTimeVegaSP;

DRLIB_END_NAMESPACE

#endif
