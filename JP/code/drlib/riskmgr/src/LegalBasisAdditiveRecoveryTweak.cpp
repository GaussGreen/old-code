/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/LegalBasisAdditiveRecoveryTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<LegalBasisAdditiveRecoveryTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<LegalBasisAdditiveRecoveryTweak>",
        typeid(TweakableWith<LegalBasisAdditiveRecoveryTweak>),
        TweakableWith<LegalBasisAdditiveRecoveryTweak>::load);

template<> CClassConstSP const RestorableWith<LegalBasisAdditiveRecoveryTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<LegalBasisAdditiveRecoveryTweak>",
        typeid(RestorableWith<LegalBasisAdditiveRecoveryTweak>),
        RestorableWith<LegalBasisAdditiveRecoveryTweak>::load);

DRLIB_END_NAMESPACE
