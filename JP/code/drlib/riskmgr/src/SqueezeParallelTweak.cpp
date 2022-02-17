//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : SqueezeParallelTweak.cpp
//
//   Description : Multiplicative parallel tweak for squeeze functions
//
//   Author      : Antoine Gregoire
//
//   Date        : June 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericScalarOneSidedShift.hpp"
#include "edginc/SqueezeParallelTweak.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<SqueezeParallelTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<SqueezeParallelTwk>",
        typeid(TweakableWith<SqueezeParallelTwk>),
        TweakableWith<SqueezeParallelTwk>::load);

template<> CClassConstSP const RestorableWith<SqueezeParallelTwk>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<SqueezeParallelTwk>",
        typeid(RestorableWith<SqueezeParallelTwk>),
        RestorableWith<SqueezeParallelTwk>::load);

typedef GenericScalarOneSidedShift<
    SqueezeParallelTwk, true /* constant divisor */> SqueezeParallelTweak;

template<> CClassConstSP const SqueezeParallelTweak::TYPE =
    CClass::registerClassLoadMethod("SqueezeParallelTweak",
                                    typeid(SqueezeParallelTweak),
                                    SqueezeParallelTweak::load);

template<> const string SqueezeParallelTweak::NAME = "SQUEEZE_PARALLEL";

template<> const double SqueezeParallelTweak::DEFAULT_SHIFT = 1.1;
template<> const double SqueezeParallelTweak::SENSITIVITY_UNIT = 1.0;

/**
 * Included in RiskMgrLib::linkInClasses() to force SqueezeParallelTweak
 * to get linked into the Windows exe.
 */
bool SqueezeParallelTweakLinkIn() {
    return SqueezeParallelTweak::TYPE != NULL;
}

DRLIB_END_NAMESPACE

