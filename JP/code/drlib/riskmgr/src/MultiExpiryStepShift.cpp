//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiExpiryStepShift.cpp
//
//   Description : Specialised Perturbation where shifts for a range of 
//                 expiries are defined together in buckets. 
//                 Really meant as a scenario not as a greek.
//                 Shifts are defined for dates <= corresponding expiry
//                 in a step-wise manner.
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MultiExpiryStepShift.hpp"

DRLIB_BEGIN_NAMESPACE

MultiExpiryStepShift::~MultiExpiryStepShift(){}

// what's the shift for a given date ? */
double MultiExpiryStepShift::shiftSize(
    const DateTime& today,     // to anchor expiries
    const DateTime& shiftDate) const {
    static const string method("MultiExpiryStepShift::shiftSize");
    try {
        double shift;
        bool   found = false;
        for (int i = 0; i < expiries.size() && !found; i++) {
            if (shiftDate <= expiries[i]->toDate(today)) {
                shift = shifts[i];
                found = true;
            }
        }

        // if not found, we after last expiry, so use upper bound
        if (!found) {
            shift = upperShift;
        }
                         
        return shift;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

MultiExpiryStepShift::MultiExpiryStepShift(const CClassConstSP& clazz):
    MultiExpiryShift(clazz){}


class MultiExpiryStepShiftHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(MultiExpiryStepShift, clazz);
        SUPERCLASS(MultiExpiryShift);
        FIELD(upperShift, "used after last expiry");
    }
};

CClassConstSP const MultiExpiryStepShift::TYPE = CClass::registerClassLoadMethod(
    "MultiExpiryStepShift", typeid(MultiExpiryStepShift), 
    MultiExpiryStepShiftHelper::load);


DRLIB_END_NAMESPACE

