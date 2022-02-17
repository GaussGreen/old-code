/**
 * @file VegaParallel2Sided.cpp
 */

#include "edginc/config.hpp"
#define QLIB_VEGAPARALLEL2SIDED_CPP
#include "edginc/VegaParallel2Sided.hpp"
#include "edginc/VegaParallel.hpp"

DRLIB_BEGIN_NAMESPACE

template <> const string GenericScalarTwoSidedShift<VolParallel, false, true>::NAME = "VEGA2SIDED";
template <> const string GenericScalarTwoSidedShift<VolParallel, false, true>::SECOND_ORDER_NAME = "VOLGAMMA";
template <> const double GenericScalarTwoSidedShift<VolParallel, false, true>::DEFAULT_SHIFT = 0.0001;
template <> const double GenericScalarTwoSidedShift<VolParallel, false, true>::SENSITIVITY_UNIT = VegaParallel::SENSITIVITY_UNIT;
template <> CClassConstSP const GenericScalarTwoSidedShift<VolParallel, false, true>::TYPE = CClass::registerClassLoadMethod(
    "GenericScalarTwoSidedShift<VolParallel, false, true>", typeid(GenericScalarTwoSidedShift<VolParallel, false, true>), load);

VegaParallel2Sided::VegaParallel2Sided(double shiftSize):
    Super(TYPE, NAME, shiftSize)
{}

VegaParallel2Sided::VegaParallel2Sided(double shiftSize,
                                       IModel* model, CControl* control):
    Super(TYPE, NAME, shiftSize)
{
    this->algorithm = model;
    this->control = control;
}

IObject* VegaParallel2Sided::defaultConstructor() { 
    return new VegaParallel2Sided(DEFAULT_SHIFT);
}

void VegaParallel2Sided::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaParallel2Sided, clazz);
    SUPERCLASS(Super);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(VegaPar, "contains a VegaParallel object");
    FIELD_MAKE_OPTIONAL(VegaPar);

    // Some numpty made VegaPar non-transient at some point, so it's in
    // many of our tests and maybe production too, who knows.  Hence I've
    // kept it for back-compat.

    //FIELD_MAKE_TRANSIENT(VegaPar);

    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<VegaParallel2Sided>(), 
                                new VegaParallel(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<VolParallel>::TYPE);
}

CClassConstSP const VegaParallel2Sided::TYPE = CClass::registerClassLoadMethod(
    "VegaParallel2Sided", typeid(VegaParallel2Sided), load);

bool VegaParallel2SidedLinkIn() {
    return VegaParallel2Sided::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
