//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadPropShift.hpp
//
//   Description : Credit Spread Proportional Shift scenario - apply proportional shift to CS curve
//
//   Author      : Andrew McCleery
//
//   Date        : 24 February 2004
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_SPREAD_PROP_SHIFT_HPP
#define CREDIT_SPREAD_PROP_SHIFT_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for cs shift scenario - apply proportional shift to CS curve */
class RISKMGR_DLL CreditSpreadPropShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support CreditSpreadPropShift */
    class RISKMGR_DLL IShift{
    public:
        friend class CreditSpreadPropShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(CreditSpreadPropShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditSpreadPropShift* shift) = 0;
    };

    /** constructor with relative shift amount */
    CreditSpreadPropShift(double shift);

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents */
    CClassConstSP shiftInterface() const;

    // methods to shift given objects

    /** Returns true if the supplied object matches the supplied name
        for this sensitivity. The object must implement this sensitivity's
        Shift interface */
    virtual bool nameMatches(const OutputName&         name,
                             IObjectConstSP            obj);

    /** Appends the name(s) of the supplied object with respect to
        this sensitivity to the supplied list. The object must
        implement this sensitivity's Shift interface */
    virtual void appendName(OutputNameArray&          namesList,
                            IObjectConstSP            obj);

    /**
     * @param obj The object to shift. The object must implement the
     CreditSpreadPropShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

private:
    friend class CreditSpreadPropShiftHelper;
    /** for reflection */
    CreditSpreadPropShift();
    CreditSpreadPropShift(const CreditSpreadPropShift &rhs);
    CreditSpreadPropShift& operator=(const CreditSpreadPropShift& rhs);
};


typedef smartConstPtr<CreditSpreadPropShift> CreditSpreadPropShiftConstSP;
typedef smartPtr<CreditSpreadPropShift> CreditSpreadPropShiftSP;

DRLIB_END_NAMESPACE

#endif
