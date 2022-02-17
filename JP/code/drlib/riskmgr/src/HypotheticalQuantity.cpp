/**
 * @file HypotheticalQuantity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/IHypothesis.hpp"
#include "edginc/HypotheticalQuantity.hpp"

DRLIB_BEGIN_NAMESPACE

HypotheticalQuantity::HypotheticalQuantity(
        IHypothesisConstSP hypothesis,
        IResultsFunctionConstSP quantity):
    CObject(TYPE),
    _hypothesis(hypothesis),
    _quantity(quantity)
{}

HypotheticalQuantitySP HypotheticalQuantity::SP(
        IHypothesisConstSP hypothesis,
        IResultsFunctionConstSP quantity) {
    return HypotheticalQuantitySP(
         new HypotheticalQuantity(hypothesis, quantity));
}

HypotheticalQuantity::~HypotheticalQuantity() {}

IHypothesisConstSP HypotheticalQuantity::hypothesis() const {
    return _hypothesis;
}

IResultsFunctionConstSP HypotheticalQuantity::quantity() const {
    return _quantity;
}

IObject* defaultHypotheticalQuantity() {
    return new HypotheticalQuantity(IHypothesisConstSP(),
                                    IResultsFunctionConstSP());
}

void HypotheticalQuantity::load(CClassSP& clazz) {
    REGISTER(HypotheticalQuantity, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultHypotheticalQuantity);
    FIELD(_hypothesis, "hypothesis");
    FIELD(_quantity, "quantity");
}

CClassConstSP const HypotheticalQuantity::TYPE =
    CClass::registerClassLoadMethod(
        "HypotheticalQuantity", typeid(HypotheticalQuantity), load);

DEFINE_TEMPLATE_TYPE(HypotheticalQuantityArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
