//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YCWeightedAdditiveShift.cpp
//
//   Description : Scenario shift: Yield Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Yield Curve benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/YCWeightedAdditiveShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** for reflection */
YCWeightedAdditiveShift::YCWeightedAdditiveShift(): 
    YCWeightedShift(TYPE){}

YCWeightedAdditiveShift::YCWeightedAdditiveShift(const CClassConstSP& clazz):
    YCWeightedShift(clazz){}

/* Additive shift to the DoubleArray passed as an argument */
void YCWeightedAdditiveShift::shiftArray(ExpiryArraySP rate_expiries, 
                                         DoubleArray* rates,
                                         DateTime& today)
{
    /* Call parent indicating we want an additive shift */
    YCWeightedShift::shiftArray(true /*additive*/, rate_expiries, rates, today);
}


class YCWeightedAdditiveShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute additive shifts to yield curve by benchmark");
        REGISTER(YCWeightedAdditiveShift, clazz);
        SUPERCLASS(YCWeightedShift);
        EMPTY_SHELL_METHOD(defaultYCWeightedAdditiveShift);
        /* Note that we do not register how to build this sensitivity because
         * this is not really a sensitivity but a scenario  */
    }

    static IObject* defaultYCWeightedAdditiveShift(){
        return new YCWeightedAdditiveShift();
    }
};


CClassConstSP const YCWeightedAdditiveShift::TYPE = 
    CClass::registerClassLoadMethod(
                                    "YCWeightedAdditiveShift", 
                                    typeid(YCWeightedAdditiveShift), 
                                    YCWeightedAdditiveShiftHelper::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool YCWeightedAdditiveShiftLinkIn() {
    return YCWeightedAdditiveShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
