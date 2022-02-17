/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/PointwiseRestorableWith.hpp"
#include "edginc/SpotVegaPointwiseTweak.hpp"


DRLIB_BEGIN_NAMESPACE

template<> CClassConstSP const PointwiseTweakableWith<SpotVegaPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseTweakableWith<SpotVegaPointwiseTweak>",
        typeid(PointwiseTweakableWith<SpotVegaPointwiseTweak>),
        PointwiseTweakableWith<SpotVegaPointwiseTweak>::load);

template<> CClassConstSP const PointwiseRestorableWith<SpotVegaPointwiseTweak>::TYPE =
    CClass::registerInterfaceLoadMethod(
        "PointwiseRestorableWith<SpotVegaPointwiseTweak>",
        typeid(PointwiseRestorableWith<SpotVegaPointwiseTweak>),
        PointwiseRestorableWith<SpotVegaPointwiseTweak>::load);

DRLIB_END_NAMESPACE
