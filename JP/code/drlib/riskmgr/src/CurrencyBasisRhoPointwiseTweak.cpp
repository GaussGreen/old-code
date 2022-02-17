/**
 *
 *
 */

#include "edginc/config.hpp"
#include "edginc/PointwiseRestorableWith.hpp"
#include "edginc/CurrencyBasisRhoPointwiseTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>",
        typeid(PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>),
        PointwiseTweakableWith<CurrencyBasisRhoPointwiseTweak>::load);

template<> CClassConstSP const PointwiseRestorableWith<CurrencyBasisRhoPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseRestorableWith<CurrencyBasisRhoPointwiseTweak>",
        typeid(PointwiseRestorableWith<CurrencyBasisRhoPointwiseTweak>),
        PointwiseRestorableWith<CurrencyBasisRhoPointwiseTweak>::load);

DRLIB_END_NAMESPACE
