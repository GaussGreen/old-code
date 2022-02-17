/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/SpotVegaParallelTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<SpotVegaParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<SpotVegaParallelTweak>",
        typeid(TweakableWith<SpotVegaParallelTweak>),
        TweakableWith<SpotVegaParallelTweak>::load);

template<> CClassConstSP const RestorableWith<SpotVegaParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<SpotVegaParallelTweak>",
        typeid(RestorableWith<SpotVegaParallelTweak>),
        RestorableWith<SpotVegaParallelTweak>::load);

DRLIB_END_NAMESPACE
