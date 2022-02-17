//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskyCurve.cpp
//
//   Description : Interface for risky curves
//
//   Author      : Andrew J Swain
//
//   Date        : 17 November 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/RiskyCurve.hpp"

DRLIB_BEGIN_NAMESPACE

#if 0
// apparently this is a sensible default behaviour regardless of the value
// of useAssetRecovery
double IRiskyCurve::riskyPV(const DateTime& lodate,
                            const DateTime& hidate,
                            double          cashFlow,
                            double          recoveryNotional,
                            bool            useAssetRecovery,
                            double          assetRecovery) const {
    return riskyPV(lodate, hidate, cashFlow, recoveryNotional);
}
#endif

IRiskyCurve::IRiskyCurve() {}
IRiskyCurve::~IRiskyCurve() {}

CClassConstSP const IRiskyCurve::TYPE = CClass::registerInterfaceLoadMethod(
    "IRiskyCurve", typeid(IRiskyCurve), 0);

DRLIB_END_NAMESPACE
