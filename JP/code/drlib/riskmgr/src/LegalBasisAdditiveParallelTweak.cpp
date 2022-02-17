/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/LegalBasisAdditiveParallelTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<LegalBasisAdditiveParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<LegalBasisAdditiveParallelTweak>",
        typeid(TweakableWith<LegalBasisAdditiveParallelTweak>),
        TweakableWith<LegalBasisAdditiveParallelTweak>::load);

template<> CClassConstSP const RestorableWith<LegalBasisAdditiveParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<LegalBasisAdditiveParallelTweak>",
        typeid(RestorableWith<LegalBasisAdditiveParallelTweak>),
        RestorableWith<LegalBasisAdditiveParallelTweak>::load);

DRLIB_END_NAMESPACE
