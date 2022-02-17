//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FXVega.cpp
//
//   Description : Controls calculation of FX_VEGA
//
//   Author      : Mark A Robson
//
//   Date        : 15 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FX_VEGA_HPP
#define EDG_FX_VEGA_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for FX_VEGA - Controls calculation of FX_VEGA */
class RISKMGR_DLL FXVega: public ScalarShift,
              public virtual Additive {
public:
    friend class FXVegaHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for VEGA_PARALLEL */
    class RISKMGR_DLL IShift{
    public:
        friend class FXVegaHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the vol - used to determine whether to tweak
            the object */
        virtual string sensName(FXVega* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(FXVega* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for VEGA_PARALLEL. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class FXVegaHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(FXVega* shift) = 0;
    };

    /** constructor with explicit shift size */
    FXVega(double     shiftSize);

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
        FXVega.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP          obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP          obj);

    /**
     * @param obj The object to shift. The object must implement the
     FXVega.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /**
     * @param param1 The object to shift. The
     object must implement the FXVega.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

protected:
    FXVega(CClassConstSP clazz, const string& sensName);

private:
    /** for reflection */
    FXVega();
    FXVega(const FXVega &rhs);
    FXVega& operator=(const FXVega& rhs);
};

typedef smartConstPtr<FXVega> FXVegaConstSP;
typedef smartPtr<FXVega> FXVegaSP;


DRLIB_END_NAMESPACE

#endif
