/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/LegalBasisMultiplierRecoveryTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<LegalBasisMultiplierRecoveryTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<LegalBasisMultiplierRecoveryTweak>",
        typeid(TweakableWith<LegalBasisMultiplierRecoveryTweak>),
        TweakableWith<LegalBasisMultiplierRecoveryTweak>::load);

template<> CClassConstSP const RestorableWith<LegalBasisMultiplierRecoveryTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<LegalBasisMultiplierRecoveryTweak>",
        typeid(RestorableWith<LegalBasisMultiplierRecoveryTweak>),
        RestorableWith<LegalBasisMultiplierRecoveryTweak>::load);

DRLIB_END_NAMESPACE
