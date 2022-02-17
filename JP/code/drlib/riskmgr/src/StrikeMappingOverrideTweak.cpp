//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : StrikeMappingOverrideTweak.cpp
//
//   Description : Tweaks the overridden strike mapping parameter in base correlation models
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericAllNameScalarShift.hpp"
//#include "edginc/GenericScalarOneSidedShift.hpp"
#include "edginc/StrikeMappingOverrideTweak.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<StrikeMappingOverrideTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<StrikeMappingOverrideTwk>",
        typeid(TweakableWith<StrikeMappingOverrideTwk>),
        TweakableWith<StrikeMappingOverrideTwk>::load);

template<> CClassConstSP const RestorableWith<StrikeMappingOverrideTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<StrikeMappingOverrideTwk>",
        typeid(RestorableWith<StrikeMappingOverrideTwk>),
        RestorableWith<StrikeMappingOverrideTwk>::load);

typedef GenericAllNameScalarShift<StrikeMappingOverrideTwk, true> StrikeMappingOverrideTweakBase;

/** Derive and override getPacketName() method */
class StrikeMappingOverrideTweak:
    public StrikeMappingOverrideTweakBase {
public:
    /** Reflection type of this class */
    static CClassConstSP const TYPE;

    /** Note that we just use NAME and DEFAULT_SHIFT from our parent class */

    /** A StrikeMappingOverrideTweak with shift size = DEFAULT_SHIFT */
    StrikeMappingOverrideTweak():
        StrikeMappingOverrideTweakBase(TYPE, NAME, DEFAULT_SHIFT) {}

    /** A StrikeMappingOverrideTweak with given shift size */
    StrikeMappingOverrideTweak(double shift):
        StrikeMappingOverrideTweakBase(TYPE, NAME, shift) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(StrikeMappingOverrideTweak, clazz);
        SUPERCLASS(StrikeMappingOverrideTweakBase);
        IMPLEMENTS(Additive);
        EMPTY_SHELL_METHOD(defaultStrikeMappingOverrideTweak);
        // no fields
        // register how to build our sensitivity
//        SensitivityFactory::addSens(
//            NAME,
//            new GenericSensitivityFactory<StrikeMappingOverrideTweak>(),
//            new StrikeMappingOverrideTweak(DEFAULT_SHIFT),
//            Tweakable::TYPE);
    }

    /**
     * Output name for first derivative
     */

    const static string NAME;

    /**
     * Shift size to use if none provided
     */

    const static double DEFAULT_SHIFT;

private:

    static IObject* defaultStrikeMappingOverrideTweak(){
        return new StrikeMappingOverrideTweak();
    }
};

template<> CClassConstSP const GenericAllNameScalarShift<StrikeMappingOverrideTwk, true>::TYPE =
    CClass::registerClassLoadMethod(
        "GenericAllNameScalarShift<StrikeMappingOverrideTwk, true>",
        typeid(GenericAllNameScalarShift<StrikeMappingOverrideTwk, true>),
        GenericAllNameScalarShift<StrikeMappingOverrideTwk, true>::load);

CClassConstSP const StrikeMappingOverrideTweak::TYPE =
CClass::registerClassLoadMethod("StrikeMappingOverrideTweak",
                                typeid(StrikeMappingOverrideTweak),
                                load);

const string StrikeMappingOverrideTweak::NAME = "STRIKE_MAPPING_OVERRIDE";

const double StrikeMappingOverrideTweak::DEFAULT_SHIFT = 0.1;

/**
 * Included in RiskMgrLib::linkInClasses() to force StrikeMappingOverrideTweak
 * to get linked into the Windows exe.
 */
bool StrikeMappingOverrideTweakLinkIn() {
    return StrikeMappingOverrideTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

