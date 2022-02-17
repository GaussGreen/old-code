//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : StrikeMappingTweak.cpp
//
//   Description : Strike mapping tweak for base correlation
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericScalarOneSidedShift.hpp"
#include "edginc/StrikeMappingTweak.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<StrikeMappingTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<StrikeMappingTwk>",
        typeid(TweakableWith<StrikeMappingTwk>),
        TweakableWith<StrikeMappingTwk>::load);

template<> CClassConstSP const RestorableWith<StrikeMappingTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<StrikeMappingTwk>",
        typeid(RestorableWith<StrikeMappingTwk>),
        RestorableWith<StrikeMappingTwk>::load);

typedef GenericScalarOneSidedShift<
    StrikeMappingTwk, true /* constant divisor */> StrikeMappingTweak;

template<> CClassConstSP const StrikeMappingTweak::TYPE =
    CClass::registerClassLoadMethod("StrikeMappingTweak",
                                    typeid(StrikeMappingTweak),
                                    StrikeMappingTweak::load);

template<> const string StrikeMappingTweak::NAME = "STRIKE_MAPPING";

template<> const double StrikeMappingTweak::DEFAULT_SHIFT = 0.1;
template<> const double StrikeMappingTweak::SENSITIVITY_UNIT = 1.0;

/**
 * Included in RiskMgrLib::linkInClasses() to force StrikeMappingTweak
 * to get linked into the Windows exe.
 */
bool StrikeMappingTweakLinkIn() {
    return StrikeMappingTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

