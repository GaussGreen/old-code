//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CCMBetaSens.hpp
//
//   Description : Sensitivity corresponding to shifting CCM Beta parameter
//
//   Author      : Antoine Gregoire
//
//   Date        : October 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_BETA_HPP
#define EDR_BETA_HPP

#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sensitivity corresponding to shifting CCM Beta parameter */
class RISKMGR_DLL CCMBetaSens: public ScalarShift,
                   public virtual Additive /* you can add results together */ {
public:
    friend class CCMBetaSensHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;


    /** What an object must implement to be tweakable for BETA */
    class RISKMGR_DLL IShift{
    public:
        friend class CCMBetaSensHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();

        /** Used to determine whether to tweak the object */
        virtual string sensName(CCMBetaSens* shift) const = 0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CCMBetaSens* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for BETA. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class CCMBetaSensHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(CCMBetaSens* shift) = 0;
    };

    virtual ~CCMBetaSens();

    /** constructor with explicit shift size */
    CCMBetaSens(double     shiftSize);

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
        CCMBetaSens.Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /**
     * @param obj The object to shift. The object must implement the
     CCMBetaSens.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    

    /**
     * @param param1 The object to shift. The
     object must implement the CCMBetaSens.RestorableShift interface
    */
    virtual void restore(IObjectSP obj);

private:
    /** for reflection */
    CCMBetaSens();
    CCMBetaSens(const CCMBetaSens &rhs);
    CCMBetaSens& operator=(const CCMBetaSens& rhs);
};

typedef smartConstPtr<CCMBetaSens> CCMBetaSensConstSP;
typedef smartPtr<CCMBetaSens> CCMBetaSensSP;

DRLIB_END_NAMESPACE

#endif
