/**
 *
 */

#include "edginc/config.hpp"
#include "edginc/GenericVectorOneSidedShift.hpp"
#include "edginc/SpotVegaPointwiseTweak.hpp"


DRLIB_BEGIN_NAMESPACE

typedef GenericVectorOneSidedShift<SpotVegaPointwiseTweak> SpotVegaPointwise;

template<> CClassConstSP const SpotVegaPointwise::TYPE =
    CClass::registerClassLoadMethod("SpotVegaPointwise",
                                    typeid(SpotVegaPointwise),
                                    load);

template<> const string SpotVegaPointwise::NAME = "SPOT_VEGA_POINTWISE";

const double ONE_BASIS_POINT = 0.0001;
template<> const double SpotVegaPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;
template<> const double SpotVegaPointwise::SENSITIVITY_UNIT = ONE_BASIS_POINT;

bool SpotVegaPointwiseLinkIn() {
    return SpotVegaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE
