/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/PointwiseRestorableWith.hpp"
#include "edginc/CompositeVegaPointwiseTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const PointwiseTweakableWith<CompositeVegaPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseTweakableWith<CompositeVegaPointwiseTweak>",
        typeid(PointwiseTweakableWith<CompositeVegaPointwiseTweak>),
        PointwiseTweakableWith<CompositeVegaPointwiseTweak>::load);

template<> CClassConstSP const PointwiseRestorableWith<CompositeVegaPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseRestorableWith<CompositeVegaPointwiseTweak>",
        typeid(PointwiseRestorableWith<CompositeVegaPointwiseTweak>),
        PointwiseRestorableWith<CompositeVegaPointwiseTweak>::load);

DRLIB_END_NAMESPACE
