//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadWeightedAdditiveShift.cpp
//
//   Description : Scenario shift: Par Spread curves benchmarks can be shifted by
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
#include "edginc/ParSpreadWeightedAdditiveShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** for reflection */
ParSpreadWeightedAdditiveShift::ParSpreadWeightedAdditiveShift(): 
    ParSpreadWeightedShift(TYPE){}


/* Additive shift to the DoubleArray passed as an argument */
void ParSpreadWeightedAdditiveShift::shiftArray(ExpiryArraySP rateExpiries, 
                                                DoubleArray* rates,
                                                DateTime& today)
{
    /* Call parent indicating we want an additive shift */
    ParSpreadWeightedShift::shiftArray(true /*additive*/, rateExpiries, rates, today);
}


class ParSpreadWeightedAdditiveShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute additive shifts to Par Spread curve by benchmark");
        REGISTER(ParSpreadWeightedAdditiveShift, clazz);
        SUPERCLASS(ParSpreadWeightedShift);
        EMPTY_SHELL_METHOD(defaultParSpreadWeightedAdditiveShift);
        /* Note that we do not register how to build this sensitivity because
         * this is not really a sensitivity but a scenario  */
    }

    static IObject* defaultParSpreadWeightedAdditiveShift(){
        return new ParSpreadWeightedAdditiveShift();
    }
};


CClassConstSP const ParSpreadWeightedAdditiveShift::TYPE = 
    CClass::registerClassLoadMethod(
                                    "ParSpreadWeightedAdditiveShift", 
                                    typeid(ParSpreadWeightedAdditiveShift), 
                                    ParSpreadWeightedAdditiveShiftHelper::load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe.
 */
bool ParSpreadWeightedAdditiveShiftLinkIn() {
    return ParSpreadWeightedAdditiveShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
