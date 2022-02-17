//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YCWeightedMultiplicativeShift.cpp
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
#include "edginc/YCWeightedMultiplicativeShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** for reflection */
YCWeightedMultiplicativeShift::YCWeightedMultiplicativeShift(): 
    YCWeightedShift(TYPE){}


/* Additive shift to the DoubleArray passed as an argument */
void YCWeightedMultiplicativeShift::shiftArray(ExpiryArraySP rateExpiries, 
                                               DoubleArray* rates,
                                               DateTime& today)
{
    /* Call parent indicating we want an additive shift */
    YCWeightedShift::shiftArray(false /*NOT additive but multiplicative*/, 
                                rateExpiries, rates, today);
}


class YCWeightedMultiplicativeShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute multiplicative shifts to yield curve by benchmark");
        REGISTER(YCWeightedMultiplicativeShift, clazz);
        SUPERCLASS(YCWeightedShift);
        EMPTY_SHELL_METHOD(defaultYCWeightedMultiplicativeShift);
        /* Note that we do not register how to build this sensitivity because
         * this is not really a sensitivity but a scenario  */
    }

    static IObject* defaultYCWeightedMultiplicativeShift(){
        return new YCWeightedMultiplicativeShift();
    }
};

CClassConstSP const YCWeightedMultiplicativeShift::TYPE = 
    CClass::registerClassLoadMethod(
           "YCWeightedMultiplicativeShift", 
           typeid(YCWeightedMultiplicativeShift), 
           YCWeightedMultiplicativeShiftHelper::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool YCWeightedMultiplicativeShiftLinkIn() {
    return YCWeightedMultiplicativeShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
