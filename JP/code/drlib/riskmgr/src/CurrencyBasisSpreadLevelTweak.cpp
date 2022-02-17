/**
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/CurrencyBasisSpreadLevelTweak.hpp"

DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<CurrencyBasisSpreadLevelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<CurrencyBasisSpreadLevelTweak>",
        typeid(TweakableWith<CurrencyBasisSpreadLevelTweak>),
        TweakableWith<CurrencyBasisSpreadLevelTweak>::load);

template<> CClassConstSP const RestorableWith<CurrencyBasisSpreadLevelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<CurrencyBasisSpreadLevelTweak>",
        typeid(RestorableWith<CurrencyBasisSpreadLevelTweak>),
        RestorableWith<CurrencyBasisSpreadLevelTweak>::load);

DRLIB_END_NAMESPACE
