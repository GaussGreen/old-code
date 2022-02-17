//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingSet.cpp
//
//   Description : Tweak to set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/BCStrikeMappingSet.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
BCStrikeMappingSet::BCStrikeMappingSet() :
    BCStrikeMappingTweakBase(TYPE)
{}


/** Returns the StrikeMapping defined as shift size */
double BCStrikeMappingSet::applyShift(double& unadjStrikeMapping) {
    double strikeMapping = getShiftSize();

    if ((strikeMapping < 0.0) || (strikeMapping > 1.0)) {
        throw ModelException("BCStrikeMappingSet::applyShift",
                             "strike mapping is " +
                             Format::toString(strikeMapping) +
                             ". It must lie in [0,1].");
    }
    return strikeMapping;
}


IObject* BCStrikeMappingSet::defaultConstructor() {
    return new BCStrikeMappingSet();
}


/** Invoked when class is 'loaded' */
void BCStrikeMappingSet::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BCStrikeMappingSet, clazz);
    SUPERCLASS(BCStrikeMappingTweakBase);
    EMPTY_SHELL_METHOD(defaultConstructor);
};


CClassConstSP const BCStrikeMappingSet::TYPE = 
    CClass::registerClassLoadMethod("BCStrikeMappingSet", 
                                    typeid(BCStrikeMappingSet), 
                                    load);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe
 */
bool BCStrikeMappingSetLinkIn() {
    return BCStrikeMappingSet::TYPE != NULL;
}

DRLIB_END_NAMESPACE
