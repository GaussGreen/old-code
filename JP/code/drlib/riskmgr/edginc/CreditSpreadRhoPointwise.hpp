//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadRhoPointwise.hpp
//
//   Description : Controls calculation of spread rho pointwise
//
//   Author      : Andr� Segger
//
//   Date        : 16 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef CS_RHOPOINTWISE_HPP
#define CS_RHOPOINTWISE_HPP
#include "edginc/VectorShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for rho pointwise */
class RISKMGR_DLL CreditSpreadRhoPointwise: public VectorShift,
                                public virtual Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** What an object must implement to be tweakable for SPREAD_RHO_POINTWISE */
    class RISKMGR_DLL IShift{
    public:
        friend class CreditSpreadRhoPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the yield curve - used to determine whether 
            to tweak the object */
        virtual string sensName(CreditSpreadRhoPointwise* shift) const = 0;

        /** Return the array of expiries (ie maturities/benchmark dates) that
            need to be tweaked for this  yield curve */
        virtual ExpiryArrayConstSP sensExpiries(CreditSpreadRhoPointwise* shift) const =0;

        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditSpreadRhoPointwise* shift) = 0;
    };

    /** What an object must implement to be able to perform a restorable
        tweak for RHO_POINTWISE. This is used by {@link SensMgr SensMgr} to
        determine whether it should try and tweak a particular object. */
    class RISKMGR_DLL IRestorableShift: public virtual IShift{
    public:
        friend class CreditSpreadRhoPointwiseHelper;
        static CClassConstSP const TYPE;
        virtual ~IRestorableShift();
        /** Restores the object to its original form */
        virtual void sensRestore(CreditSpreadRhoPointwise* shift) = 0;
    };

    /** constructor with explicit shift size */
    CreditSpreadRhoPointwise(double shiftSize);

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
     RhoPointwise.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
    
    /** The supplied object is queried for the expiries array needed
        for doing rho pointwise and this array is returned. The supplied
        object must implement the RhoPointwise.Shift interface */
    virtual IObjectConstSP qualifier(IObjectConstSP obj);

    /**
     * @param param1 The object to shift. The
     object must implement the RhoPointwise.RestorableShift interface
     */
    virtual void restore(IObjectSP obj);

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

private:
    friend class CreditSpreadRhoPointwiseHelper;
    /** for reflection */
    CreditSpreadRhoPointwise();
    CreditSpreadRhoPointwise(const CreditSpreadRhoPointwise &rhs);
    CreditSpreadRhoPointwise& operator=(const CreditSpreadRhoPointwise& rhs);
};

typedef smartConstPtr<CreditSpreadRhoPointwise> CreditSpreadRhoPointwiseConstSP;
typedef smartPtr<CreditSpreadRhoPointwise> CreditSpreadRhoPointwiseSP;


DRLIB_END_NAMESPACE

#endif
