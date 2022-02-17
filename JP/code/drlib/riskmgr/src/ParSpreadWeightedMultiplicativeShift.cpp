//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadWeightedMultiplicativeShift.cpp
//
//   Description : Scenario shift: Par Spread Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Par Spread benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ParSpreadWeightedMultiplicativeShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** for reflection */
ParSpreadWeightedMultiplicativeShift::ParSpreadWeightedMultiplicativeShift(): 
    ParSpreadWeightedShift(TYPE){}


/* Additive shift to the DoubleArray passed as an argument */
void ParSpreadWeightedMultiplicativeShift::shiftArray(
    ExpiryArraySP rateExpiries, 
    DoubleArray* rates,
    DateTime& today){
    /* Call parent indicating we want an additive shift */
    ParSpreadWeightedShift::shiftArray(false/*NOT additive but multiplicative*/,
                                       rateExpiries, rates, today);
}


class ParSpreadWeightedMultiplicativeShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute multiplicative shifts to Par Spread curves by benchmark");
        REGISTER(ParSpreadWeightedMultiplicativeShift, clazz);
        SUPERCLASS(ParSpreadWeightedShift);
        EMPTY_SHELL_METHOD(defaultParSpreadWeightedMultiplicativeShift);
        /* Note that we do not register how to build this sensitivity 
         * because this is not really a sensitivity but a scenario  */
    }

    static IObject* defaultParSpreadWeightedMultiplicativeShift(){
        return new ParSpreadWeightedMultiplicativeShift();
    }
};

CClassConstSP const ParSpreadWeightedMultiplicativeShift::TYPE = 
    CClass::registerClassLoadMethod(
           "ParSpreadWeightedMultiplicativeShift", 
           typeid(ParSpreadWeightedMultiplicativeShift), 
           ParSpreadWeightedMultiplicativeShiftHelper::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool ParSpreadWeightedMultiplicativeShiftLinkIn() {
    return ParSpreadWeightedMultiplicativeShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
