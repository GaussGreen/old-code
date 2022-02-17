/**
 * @file RiskProperty_TYPES.hpp
 */

#ifndef QLIB_RiskProperty_TYPES_H
#define QLIB_RiskProperty_TYPES_H

#include "edginc/MultiTweakGroup.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/PropertyPerturbation.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

/**
 * Defines ::TYPE descriptors for all template classes associated with a "property tag".
 *
 * See VolPointwise for an example.
 */

#define RiskProperty_TYPES(TAG) \
DEFINE_TEMPLATE_TYPE(RiskProperty<TAG>) \
DEFINE_TEMPLATE_TYPE(PropertyRiskAxis<TAG>) \
DEFINE_TEMPLATE_TYPE(PropertyPerturbation<TAG>) \
DEFINE_TEMPLATE_TYPE(PropertyTweakHypothesis<TAG>) \
DEFINE_INTERFACE_TEMPLATE_TYPE(ITweakableWithRespectTo<TAG>) \
DEFINE_INTERFACE_TEMPLATE_TYPE(IRestorableWithRespectTo<TAG>)

#endif
