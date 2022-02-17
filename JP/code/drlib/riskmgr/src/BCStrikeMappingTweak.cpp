//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingTweak.cpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/BCStrikeMappingTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/** constructor */
BCStrikeMappingTweak::BCStrikeMappingTweak() :
    BCStrikeMappingTweakBase(TYPE), tweakType(IPerturbation::SET)
{}


/** Shifts the original StrikeMapping according to the tweakType and
    returns the new one */
double BCStrikeMappingTweak::applyShift(double& unadjStrikeMapping) {
    static const string method = "BCStrikeMappingTweak::applyShift";
    double shiftSize = getShiftSize();
    double strikeMapping;

    // make adjustment
    if (tweakType == IPerturbation::SET) {
        if (shiftSize < 0.0 || shiftSize > 1.0) {
            throw ModelException("BCStrikeMappingTweak::applyShift",
                                 "strike mapping is " +
                                 Format::toString(shiftSize) +
                                 ". It must lie in [0,1].");
        }
        // Set the strikeMapping
        strikeMapping = shiftSize;
    }
    else {
        if (tweakType == IPerturbation::ABSOLUTE) {
            // Add the shift size as an absolut value
            strikeMapping = unadjStrikeMapping + shiftSize;
        }
        else if (tweakType == IPerturbation::RELATIVE) {
            // Add the shift size as a relative value
            strikeMapping = unadjStrikeMapping * (1 + shiftSize);
        }
        else {
            throw ModelException(method,
                                 "Invalid tweakType: '" + tweakType + 
                                 "'. Valid values are '" +
                                 IPerturbation::SET + "', '" +
                                 IPerturbation::ABSOLUTE + "' and '" +
                                 IPerturbation::RELATIVE + "'.");
        }

        //Cap and floor tweaked value
        if (strikeMapping < 0.0) {
            strikeMapping = 0.0;
        }
        if (strikeMapping > 1.0) {
            strikeMapping = 1.0;
        }
    }
    return strikeMapping;
}


IObject* BCStrikeMappingTweak::defaultConstructor() {
    return new BCStrikeMappingTweak();
}


/** Invoked when class is 'loaded' */
void BCStrikeMappingTweak::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BCStrikeMappingTweak, clazz);
    SUPERCLASS(BCStrikeMappingTweakBase);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(tweakType, "Indicates whether to '" + 
                 IPerturbation::SET + "' (Set) the shiftsize (default), '" + 
                 IPerturbation::ABSOLUTE + "' (Add) it in Absolute terms "
                 "(ie, += shiftsize), or add it in '" +
                 IPerturbation::RELATIVE + "' (Relative) terms "
                 "(ie, *= [1+shiftsize]).");

    FIELD_MAKE_OPTIONAL(tweakType);
};


CClassConstSP const BCStrikeMappingTweak::TYPE = 
    CClass::registerClassLoadMethod("BCStrikeMappingTweak", 
                                    typeid(BCStrikeMappingTweak), 
                                    load);


/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe
 */
bool BCStrikeMappingTweakLinkIn() {
    return BCStrikeMappingTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE
