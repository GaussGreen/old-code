//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : FeeLegCouponShift.cpp
//
//   Description : Sensitivity to a shift in CDO fee leg coupon
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CreditFeeLegCouponShift.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE
const string CreditFeeLegCouponShift::NAME = "CREDIT_FEE_LEG_COUPON_SHIFT";
const double CreditFeeLegCouponShift::DEFAULT_SHIFT = 0.001;
const double CreditFeeLegCouponShift::SENSITIVITY_UNIT = 1;

IHypothesis::AlternateWorldSP CreditFeeLegCouponShift::appliedTo(
        OutputNameConstSP name,
        double shiftSize,
        IObjectSP world) {
    try {
        return property()->axisFor(name, VoidSP())->hypothesis(shiftSize)->
                   appliedTo(world);
    }
    catch (exception& e) {
        throw ModelException(e, "CreditFeeLegCouponShift::appliedTo()");
    }
}

double CreditFeeLegCouponShift::getShiftSize() const {
    return shiftSize;
}

OutputNameArrayConstSP CreditFeeLegCouponShift::allNames(const IObject* object) const {
    return OutputNameArrayConstSP(new OutputNameArray());//Empty
}

CreditFeeLegCouponShift::CreditFeeLegCouponShift(double shiftSize):
    AllNamesRiskPropertySensitivity(TYPE, NAME, shiftSize)
{}

AllNamesRiskPropertySensitivity::Deriv CreditFeeLegCouponShift::deriv() const {
    return Deriv(IResultsFunction::price(),
                 RiskProperty<CreditFeeLegCoupon>::SP(),
                 IScalarDerivative::oneSided(),
                 SENSITIVITY_UNIT);
}

CreditFeeLegCouponShift::~CreditFeeLegCouponShift() {}

void CreditFeeLegCouponShift::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(CreditFeeLegCouponShift, clazz);
    SUPERCLASS(AllNamesRiskPropertySensitivity);
    IMPLEMENTS(IScalarPerNameShift);
    EMPTY_SHELL_METHOD(&DefaultConstructor<CreditFeeLegCouponShift>::iObject);
    SensitivityFactory::addSens(NAME,
                                new GenericSensitivityFactory<CreditFeeLegCouponShift>(), 
                                new CreditFeeLegCouponShift(DEFAULT_SHIFT),
                                ITweakableWithRespectTo<CreditFeeLegCoupon>::TYPE);
}

CClassConstSP const CreditFeeLegCouponShift::TYPE = CClass::registerClassLoadMethod(
    "CreditFeeLegCouponShift", typeid(CreditFeeLegCouponShift), load);

bool CreditFeeLegCouponShiftLinkIn() {
    return CreditFeeLegCouponShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
