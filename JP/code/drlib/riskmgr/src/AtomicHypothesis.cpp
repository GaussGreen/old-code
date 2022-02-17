/**
 * @file AtomicHypothesis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/AtomicHypothesis.hpp"

DRLIB_BEGIN_NAMESPACE

AtomicHypothesis::AtomicHypothesis(CClassConstSP type):
    CObject(type),
	approxOrder(1)
{}

AtomicHypothesis::~AtomicHypothesis() {}

int AtomicHypothesis::numAtomics() const {
    return 1;
}

AtomicHypothesisConstSP AtomicHypothesis::atomic(int i) const {
    ASSERT(i == 0);
    return AtomicHypothesisConstSP::attachToRef(this);
}

IHypothesis::IDistanceMetricConstSP AtomicHypothesis::distanceMetric() const {
    return IHypothesis::IDistanceMetric::last();
}

void AtomicHypothesis::load(CClassSP& clazz) {
    REGISTER(AtomicHypothesis, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IHypothesis);
	FIELD(approxOrder, "approxOrder");
	FIELD_MAKE_OPTIONAL(approxOrder);
}

int AtomicHypothesis::getApproxOrder() const {
    return approxOrder;
}

void AtomicHypothesis::setApproxOrder(int order){
    approxOrder = order;
}


CClassConstSP const AtomicHypothesis::TYPE = CClass::registerClassLoadMethod(
    "AtomicHypothesis", typeid(AtomicHypothesis), load);

DEFINE_TEMPLATE_TYPE(AtomicHypothesisArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
