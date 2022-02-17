//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadPropShift.cpp
//
//   Description : CDS Par Spread Proportional Shift scenario 
//                 - apply proportional shift to CDS par spread curve
//
//   Author      : Andrew McCleery
//
//   Date        : 16 March 2004
//
//
//----------------------------------------------------------------------------

#ifndef PAR_SPREAD_PROP_SHIFT_HPP
#define PAR_SPREAD_PROP_SHIFT_HPP
#include "edginc/ScalarPerturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for ps shift scenario - apply proportional shift to CDS spread curve */
class RISKMGR_DLL ParSpreadPropShift: public ScalarPerturbation {
public:
    static CClassConstSP const TYPE;

    /** What an object must implement to support ParSpreadPropShift */
    class RISKMGR_DLL IShift{
    public:
        friend class ParSpreadPropShiftHelper;
        static CClassConstSP const TYPE;
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(ParSpreadPropShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(ParSpreadPropShift* shift) = 0;
    };

    /** constructor with relative shift amount */
    ParSpreadPropShift(double shift);

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
     ParSpreadPropShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);

private:
    friend class ParSpreadPropShiftHelper;
    /** for reflection */
    ParSpreadPropShift();
    ParSpreadPropShift(const ParSpreadPropShift &rhs);
    ParSpreadPropShift& operator=(const ParSpreadPropShift& rhs);
};


typedef smartConstPtr<ParSpreadPropShift> ParSpreadPropShiftConstSP;
typedef smartPtr<ParSpreadPropShift> ParSpreadPropShiftSP;

DRLIB_END_NAMESPACE

#endif
