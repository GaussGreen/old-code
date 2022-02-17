//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditSpreadParallelShift.hpp
//
//   Description : Credit Spread Parallel Shift Scenario
//                 - apply parallel shift to credit spread curve
//
//   Author      : Andrew McCleery
//
//   Date        : 25 March 2004
//
//
//----------------------------------------------------------------------------

#ifndef CREDIT_SPREAD_PARALLEL_SHIFT_HPP
#define CREDIT_SPREAD_PARALLEL_SHIFT_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

// Sens Control for CS shift scenario - apply parallel shift to credit spread curve
class RISKMGR_DLL CreditSpreadParallelShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support CreditSpreadParallelShift */
    class RISKMGR_DLL IShift{
    public:
        friend class CreditSpreadParallelShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the CS curve - used to determine
            whether to shift the object */
        virtual string sensName(CreditSpreadParallelShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(CreditSpreadParallelShift* shift) = 0;
    };

    /** constructor with explicit shift size */
    CreditSpreadParallelShift(double shiftSize);

    // pretty meaningless for a scenario
    double divisor() const;

    /** returns the interface identifying what an object has to do in order
        to support the tweak that this object represents */
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
     CreditSpreadParallelShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

private:
    friend class CreditSpreadParallelShiftHelper;
    /** for reflection */
    CreditSpreadParallelShift();
    CreditSpreadParallelShift(const CreditSpreadParallelShift &rhs);
    CreditSpreadParallelShift& operator=(const CreditSpreadParallelShift& rhs);
};

typedef smartConstPtr<CreditSpreadParallelShift> CreditSpreadParallelShiftConstSP;
typedef smartPtr<CreditSpreadParallelShift> CreditSpreadParallelShiftSP;

DRLIB_END_NAMESPACE

#endif
