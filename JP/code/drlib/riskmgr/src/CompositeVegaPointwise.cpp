/**
 *
 */


#include "edginc/config.hpp"
#include "edginc/GenericVectorOneSidedShift.hpp"
#include "edginc/CompositeVegaPointwiseTweak.hpp"


DRLIB_BEGIN_NAMESPACE

typedef GenericVectorOneSidedShift<
    CompositeVegaPointwiseTweak> CompositeVegaPointwise;

template<> CClassConstSP const CompositeVegaPointwise::TYPE =
    CClass::registerClassLoadMethod("CompositeVegaPointwise",
                                    typeid(CompositeVegaPointwise),
                                    load);

template<> const string CompositeVegaPointwise::NAME = "COMP_VEGA_POINTWISE";

const double ONE_BASIS_POINT = 0.0001;
template<> const double CompositeVegaPointwise::DEFAULT_SHIFT = ONE_BASIS_POINT;
template<> const double CompositeVegaPointwise::SENSITIVITY_UNIT = ONE_BASIS_POINT;


bool CompositeVegaPointwiseLinkIn() {
    return CompositeVegaPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE
