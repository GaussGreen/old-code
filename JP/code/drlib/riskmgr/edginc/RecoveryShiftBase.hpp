//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RecoveryShiftBase.hpp
//
//   Description : Sensitivity to Recovery
//
//   Author      : Stephen Hope
//
//   Date        : 19 August 2002
//
//----------------------------------------------------------------------------

#ifndef EDG_RECOVERYSHIFTBASE_HPP
#define EDG_RECOVERYSHIFTBASE_HPP

#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/RecoveryPerturb.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for RecoveryShiftBase */
class RISKMGR_DLL RecoveryShiftBase: public ScalarShift,
                         public virtual Additive,
                         public virtual IRecoveryPerturb{
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for types derived from
        RecoveryShiftBase. */
    class RISKMGR_DLL IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the par curve - used to determine
            whether to tweak the object */
        virtual string sensName(IRecoveryPerturb* shift) const = 0;
        
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(IRecoveryPerturb* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for RECOVERY_TWEAK. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(IRecoveryPerturb* shift) = 0;
    };

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
        RecoveryShiftBase.Shift interface */
    virtual bool nameMatches(const OutputName& name,
                             IObjectConstSP    obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&  namesList,
                            IObjectConstSP    obj);

     
    /**
     * @param obj The object to shift. The object must implement the
     RecoveryShiftBase.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** @param param1 The object to shift. The
        object must implement the RecoveryShiftBase.RestorableShift interface */
    virtual void restore(IObjectSP obj);

protected:
    /** for derived classes */
    RecoveryShiftBase(CClassConstSP clazz,
                      const string& outputName,
                      double        shiftSize);

    RecoveryShiftBase(CClassConstSP clazz,
                      const string& outputName);

    /** Will change the sign of the shift */
    // may be called from applyShift methods of concrete instances
    void negateShift();

private:
    static void load(CClassSP& clazz);
};


DRLIB_END_NAMESPACE

#endif
