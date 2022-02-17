//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RateShift.cpp
//
//   Description : Scenario shift: Yield Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Yield Curve benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//                 ** NOTE **
//                 ==========
//                 This class is to be retired and should not be used.
//                 Use class YCWeightedAdditiveShift instead.
//
//   Author      : Stephen Hope
//
//   Date        : 7 May 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/RateShift.hpp"

DRLIB_BEGIN_NAMESPACE

RateShift::RateShift() : YCWeightedAdditiveShift(TYPE) {}

class RateShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute shifts to yield curve by benchmark");
        REGISTER(RateShift, clazz);
        SUPERCLASS(YCWeightedAdditiveShift);
        EMPTY_SHELL_METHOD(defaultRateShift);
    }

    static IObject* defaultRateShift(){
        return new RateShift();
    }
};

CClassConstSP const RateShift::TYPE = CClass::registerClassLoadMethod(
    "RateShift", typeid(RateShift), RateShiftHelper::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool RateShiftLinkIn() {
    return RateShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
