//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolAbsoluteShift.hpp
//
//   Description : Scenario shift where an absolute shift (e.g. a 10% shift
//                 to 20% vol gives 30% vol) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 9 May 2003
//
//
//----------------------------------------------------------------------------

#ifndef _VOLABSOLUTESHIFT__HPP
#define _VOLABSOLUTESHIFT__HPP
#include "edginc/MultiExpiryStepShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Scenario shift where an absolute shift (e.g. a 10% shift
//  to 20% vol gives 30% vol) is applied for each benchmark.
//  Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//  1Y -> 5Y, 5Y ->  30Y etc).
*/

class RISKMGR_DLL VolAbsoluteShift: public MultiExpiryStepShift {
public:    
    static CClassConstSP const TYPE;

    /** What an object must implement to support VolAbsoluteShift */
    class RISKMGR_DLL IShift{
    public:
        friend class VolAbsoluteShiftHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(VolAbsoluteShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(VolAbsoluteShift* shift) = 0;
    };

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
     VolAbsoluteShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
 

    virtual void validatePop2Object();

private:
    friend class VolAbsoluteShiftHelper;
    // for reflection
    VolAbsoluteShift();
};


typedef smartConstPtr<VolAbsoluteShift> VolAbsoluteShiftConstSP;
typedef smartPtr<VolAbsoluteShift> VolAbsoluteShiftSP;

DRLIB_END_NAMESPACE

#endif
