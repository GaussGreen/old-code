//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ICreditDefaultSens.hpp
//
//   Description : Sensitivities corresponding to 'shifting' credit default, 
//                 with the possibility to override recovery :
//                 - CreditDefaultSens : use current recovery
//                 - CreditDefaultSensWithZeroRecovery : use recovery = 0
//                 - CreditDefaultSensWithSpecifiedRecovery : use recovery = user specified recovery
//
//   Author      : Antoine Gregoire
//
//   Date        : February 2005
//
//----------------------------------------------------------------------------

#ifndef EDR_CREDIT_DEFAULT_SENS_BASE_HPP
#define EDR_CREDIT_DEFAULT_SENS_BASE_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Additive.hpp"
#include "edginc/SensControlPerName.hpp"


DRLIB_BEGIN_NAMESPACE

// -------------------
// ABSTRACT BASE CLASS
// -------------------
class RISKMGR_DLL CreditDefaultSensBase: public SensControlPerName,
                                         public virtual Additive {
public:
    static CClassConstSP const TYPE;
    
    static double const DIVISOR;

    virtual ~CreditDefaultSensBase();
    
    /** Return the overriden recovery given the current recovery */    
    virtual double getOverriddenRecovery(double currentRecovery) = 0;
    
    /** What an object must implement to be 'tweakable' for Credit Default */
    class RISKMGR_DLL IShift {
    public:
        static CClassConstSP const TYPE;
        virtual ~IShift();
        
        /** Used to determine whether to tweak the object */
        virtual string sensName(CreditDefaultSensBase* sens) const = 0;
        
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditDefaultSensBase* sens) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        'tweak' for Credit Default. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift {
    public:
        friend class CreditDefaultSensHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        
        /** Restores the object to its original form */
        virtual void sensRestore(CreditDefaultSensBase* sens) = 0;
    };
    
    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. */
    virtual bool discreteShift() const;

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    virtual double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    virtual CClassConstSP shiftInterface() const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    virtual CClassConstSP restorableShiftInterface() const;

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity.  The object must implement the
        VegaParallel.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /** Shifts the object (which supports being tweaked
        by this type of sens control) using given shift. The return value
        indicates whether or not components of this object need to be
        tweaked ie true: infrastructure should continue to recurse through
        components tweaking them; false: the infrastructure shouldn't
        touch any components within this object */
    virtual bool shift(IObjectSP obj);

    /** Restores the object (which supports being tweaked
        by this type of sens control) to its original form */
    virtual void restore(IObjectSP obj);
    
    /** Returns the date of the credit event */
    DateTime getCreditEventDate(const DateTime& valueDate) const;

protected:
    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*     tweakGroup,
                           CResults*       results);
                           
    CreditDefaultSensBase(
        const CClassConstSP& clazz,
        const string& outputName);
                               
private:
    static void load(CClassSP& clazz);

    // Fields
    DateTime defaultDate; // optional
    
    // not implemented !
    CreditDefaultSensBase(const CreditDefaultSensBase &rhs);
    CreditDefaultSensBase& operator=(const CreditDefaultSensBase& rhs);
};

typedef smartConstPtr<CreditDefaultSensBase> CreditDefaultSensBaseConstSP;
typedef smartPtr<CreditDefaultSensBase> CreditDefaultSensBaseSP;


DRLIB_END_NAMESPACE

#endif

