/**
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/CurrencyBasisRhoParallelTweak.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<CurrencyBasisRhoParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<CurrencyBasisRhoParallelTweak>",
        typeid(TweakableWith<CurrencyBasisRhoParallelTweak>),
        TweakableWith<CurrencyBasisRhoParallelTweak>::load);

template<> CClassConstSP const RestorableWith<CurrencyBasisRhoParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<CurrencyBasisRhoParallelTweak>",
        typeid(RestorableWith<CurrencyBasisRhoParallelTweak>),
        RestorableWith<CurrencyBasisRhoParallelTweak>::load);

DRLIB_END_NAMESPACE
