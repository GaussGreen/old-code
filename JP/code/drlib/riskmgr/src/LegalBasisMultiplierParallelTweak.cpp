/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/LegalBasisMultiplierParallelTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<LegalBasisMultiplierParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<LegalBasisMultiplierParallelTweak>",
        typeid(TweakableWith<LegalBasisMultiplierParallelTweak>),
        TweakableWith<LegalBasisMultiplierParallelTweak>::load);

template<> CClassConstSP const RestorableWith<LegalBasisMultiplierParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<LegalBasisMultiplierParallelTweak>",
        typeid(RestorableWith<LegalBasisMultiplierParallelTweak>),
        RestorableWith<LegalBasisMultiplierParallelTweak>::load);

DRLIB_END_NAMESPACE
