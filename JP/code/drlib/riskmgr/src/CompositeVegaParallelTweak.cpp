/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/RestorableWith.hpp"
#include "edginc/CompositeVegaParallelTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const TweakableWith<CompositeVegaParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "TweakableWith<CompositeVegaParallelTweak>",
        typeid(TweakableWith<CompositeVegaParallelTweak>),
        TweakableWith<CompositeVegaParallelTweak>::load);

template<> CClassConstSP const RestorableWith<CompositeVegaParallelTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "RestorableWith<CompositeVegaParallelTweak>",
        typeid(RestorableWith<CompositeVegaParallelTweak>),
        RestorableWith<CompositeVegaParallelTweak>::load);

DRLIB_END_NAMESPACE
