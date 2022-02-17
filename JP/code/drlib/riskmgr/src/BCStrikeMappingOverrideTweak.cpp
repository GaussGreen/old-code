//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingOverrideTweak.cpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//                 Override parameter
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/GenericAllNameScalarShift.hpp"
#include "edginc/BCStrikeMappingOverrideTweak.hpp"

DRLIB_BEGIN_NAMESPACE

BCStrikeMappingOverrideTwk::~BCStrikeMappingOverrideTwk()
{}

template<> CClassConstSP const TweakableWith<BCStrikeMappingOverrideTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<BCStrikeMappingOverrideTwk>",
        typeid(TweakableWith<BCStrikeMappingOverrideTwk>),
        TweakableWith<BCStrikeMappingOverrideTwk>::load);

template<> CClassConstSP const RestorableWith<BCStrikeMappingOverrideTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<BCStrikeMappingOverrideTwk>",
        typeid(RestorableWith<BCStrikeMappingOverrideTwk>),
        RestorableWith<BCStrikeMappingOverrideTwk>::load);

typedef GenericAllNameScalarShift<BCStrikeMappingOverrideTwk, true> BCStrikeMappingOverrideTweakBase;



class BCStrikeMappingOverrideTweak : public BCStrikeMappingOverrideTweakBase {
public:
    /** Reflection type of this class */
    static CClassConstSP const TYPE;

    /** A StrikeMappingOverrideTweak with shift size = DEFAULT_SHIFT */
    BCStrikeMappingOverrideTweak() :
        BCStrikeMappingOverrideTweakBase(TYPE, NAME, DEFAULT_SHIFT),
        tweakType(IPerturbation::SET)
    {}

    /** A StrikeMappingOverrideTweak with given shift size */
    BCStrikeMappingOverrideTweak(double shift) :
        BCStrikeMappingOverrideTweakBase(TYPE, NAME, shift),
        tweakType(IPerturbation::SET)
    {}

    /** Shifts the original StrikeMappingOverride according to the tweakType 
     * and returns the new one */
    virtual double applyShift(double& unadjStrikeMappingOverride) {
        static const string method = "BCStrikeMappingOverrideTweak::applyShift";
        double shiftSize = getShiftSize();
        double strikeMappingOverride;

        // make adjustment
        if (tweakType == IPerturbation::SET) {
            if (shiftSize < 0.0 || shiftSize > 1.0) {
                throw ModelException("BCStrikeMappingOverrideTweak::applyShift",
                                     "strike mapping override is " +
                                     Format::toString(shiftSize) +
                                     ". It must lie in [0,1].");
            }
            // Set the strikeMapping
            strikeMappingOverride = shiftSize;
        }
        else {
            if (tweakType == IPerturbation::ABSOLUTE) {
                // Add the shift size as an absolut value
                strikeMappingOverride = unadjStrikeMappingOverride + shiftSize;
            }
            else if (tweakType == IPerturbation::RELATIVE) {
                // Add the shift size as a relative value
                strikeMappingOverride = unadjStrikeMappingOverride * (1 + shiftSize);
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
            if (strikeMappingOverride < 0.0) {
                strikeMappingOverride = 0.0;
            }
            if (strikeMappingOverride > 1.0) {
                strikeMappingOverride = 1.0;
            }
        }
        return strikeMappingOverride;
    }

    /** Output name for first derivative */
    const static string NAME;

    /** Shift size to use if none provided */
    const static double DEFAULT_SHIFT;

private:
    string tweakType; // optional - default to "S"et (as opposed to 
                      // "A"bsolute or "R"elative shift)

    /** for reflection */
    static IObject* defaultBCStrikeMappingOverrideTweak(){
        return new BCStrikeMappingOverrideTweak();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BCStrikeMappingOverrideTweak, clazz);
        SUPERCLASS(BCStrikeMappingOverrideTweakBase);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultBCStrikeMappingOverrideTweak);
        FIELD(tweakType, "Indicates whether to '" + 
                     IPerturbation::SET + "' (Set) the shiftsize (default), '" + 
                     IPerturbation::ABSOLUTE + "' (Add) it in Absolute terms "
                     "(ie, += shiftsize), or add it in '" +
                     IPerturbation::RELATIVE + "' (Relative) terms "
                     "(ie, *= [1+shiftsize]).");
        FIELD_MAKE_OPTIONAL(tweakType);
    }
};

template<> CClassConstSP const GenericAllNameScalarShift<BCStrikeMappingOverrideTwk, true>::TYPE =
    CClass::registerClassLoadMethod(
        "GenericAllNameScalarShift<BCStrikeMappingOverrideTwk, true>",
        typeid(GenericAllNameScalarShift<BCStrikeMappingOverrideTwk, true>),
        GenericAllNameScalarShift<BCStrikeMappingOverrideTwk, true>::load);

CClassConstSP const BCStrikeMappingOverrideTweak::TYPE =
CClass::registerClassLoadMethod("BCStrikeMappingOverrideTweak",
                                typeid(BCStrikeMappingOverrideTweak),
                                load);

const string BCStrikeMappingOverrideTweak::NAME = "BC_STRIKE_MAPPING_OVERRIDE";

const double BCStrikeMappingOverrideTweak::DEFAULT_SHIFT = 0.1;



/*******************************************************************************
 * Need to create a different tweak that just sets the StrikeMappingOverride
 * (rather than using BCStrikeMappingOverrideTweak with a tweakType of "S"et) 
 * because, at the moment, Archimedes does not handle strings well - to ease 
 * integration this new (redundant) tweak has been added.
 *******************************************************************************/

class BCStrikeMappingOverrideSet : public BCStrikeMappingOverrideTweakBase {
public:
    /** Reflection type of this class */
    static CClassConstSP const TYPE;

    /** A StrikeMappingOverrideSet with shift size = DEFAULT_SHIFT */
    BCStrikeMappingOverrideSet():
        BCStrikeMappingOverrideTweakBase(TYPE, NAME, DEFAULT_SHIFT) 
    {}

    /** A StrikeMappingOverrideSet with given shift size */
    BCStrikeMappingOverrideSet(double shift):
        BCStrikeMappingOverrideTweakBase(TYPE, NAME, shift) 
    {}

    /** Shifts the original StrikeMappingOverride according to the tweakType 
     * and returns the new one */
    virtual double applyShift(double& unadjStrikeMappingOverride) {
        double strikeMappingOverride = getShiftSize();

        if ((strikeMappingOverride < 0.0) || (strikeMappingOverride > 1.0)) {
            throw ModelException("BCStrikeMappingOverrideSet::applyShift",
                                 "strike mapping override is " +
                                 Format::toString(strikeMappingOverride) +
                                 ". It must lie in [0,1]");
        }
        return strikeMappingOverride;
    }

    /** Output name for first derivative */
    const static string NAME;

    /** Shift size to use if none provided */
    const static double DEFAULT_SHIFT;

private:
    /** for reflection */
    static IObject* defaultBCStrikeMappingOverrideSet(){
        return new BCStrikeMappingOverrideSet();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BCStrikeMappingOverrideSet, clazz);
        SUPERCLASS(BCStrikeMappingOverrideTweakBase);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultBCStrikeMappingOverrideSet);
    }
};

CClassConstSP const BCStrikeMappingOverrideSet::TYPE =
CClass::registerClassLoadMethod("BCStrikeMappingOverrideSet",
                                typeid(BCStrikeMappingOverrideSet),
                                load);

const string BCStrikeMappingOverrideSet::NAME = "BC_STRIKE_MAPPING_OVERRIDE_SET";

const double BCStrikeMappingOverrideSet::DEFAULT_SHIFT = 0.1;




/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe
 */
bool BCStrikeMappingOverrideTweakLinkIn() {
    return BCStrikeMappingOverrideTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE
