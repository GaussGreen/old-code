//----------------------------------------------------------------------------
//
//   Group       : EDG Quantitative Research
//
//   Filename    : ContangoTermShift.hpp
//
//   Description : Scenario shift where an absolute shift (e.g. a 10% shift
//                 to 20% rate gives 30% rate) is applied for each benchmark.
//                 Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//                 1Y -> 5Y, 5Y ->  30Y etc).
//
//   Author      : Andrew J Swain
//
//   Date        : 19 April 2006
//
//
//----------------------------------------------------------------------------

#ifndef _CONTANGOTERMSHIFT__HPP
#define _CONTANGOTERMSHIFT__HPP
#include "edginc/MultiExpiryStepShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Scenario shift where an absolute shift (e.g. a 10% shift
//  to 20% rate gives 30% rate) is applied for each benchmark.
//  Shifts are defined for ranges of benchmarks (e.g. <= 1Y,
//  1Y -> 5Y, 5Y ->  30Y etc).
*/

class RISKMGR_DLL ContangoTermShift: public MultiExpiryStepShift {
public:    
    static CClassConstSP const TYPE;

    /** What an object must implement to support ContangoTermShift */
    class RISKMGR_DLL IShift{
    public:
        friend class ContangoTermShiftHelper;
        static CClassConstSP const TYPE;
        IShift();
        virtual ~IShift();
        /** Returns the name of the stock/asset - used to determine
            whether to shift the object */
        virtual string sensName(ContangoTermShift* shift) const = 0;
        /** Shifts the object using given shift. Return true to make the
            infrastructure keep tweaking the components within the object
            which implements this interface */
        virtual bool sensShift(ContangoTermShift* shift) = 0;
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
     ContangoTermShift.Shift interface. The return value indicates whether
     or not components of this object need to be tweaked ie true:
     infrastructure should continue to recurse through components
     tweaking them; false: the infrastructure shouldn't touch any
     components within this object */
    virtual bool shift(IObjectSP obj);
 

    virtual void validatePop2Object();

private:
    friend class ContangoTermShiftHelper;
    // for reflection
    ContangoTermShift();
};


typedef smartConstPtr<ContangoTermShift> ContangoTermShiftConstSP;
typedef smartPtr<ContangoTermShift> ContangoTermShiftSP;

DRLIB_END_NAMESPACE

#endif
