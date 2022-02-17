//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CreditFeeLegCoupon.cpp
//
//   Description : Tag class for shift in CDO fee leg coupon
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskProperty_TYPES.hpp"
#include "edginc/CreditFeeLegCoupon.hpp"

DRLIB_BEGIN_NAMESPACE

CreditFeeLegCoupon::CreditFeeLegCoupon(): CObject(TYPE) {}
CreditFeeLegCoupon::~CreditFeeLegCoupon() {}

static void CreditFeeLegCoupon_load(CClassSP& clazz) {
    REGISTER(CreditFeeLegCoupon, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(DefaultConstructor<CreditFeeLegCoupon>::iObject);
}

CClassConstSP const CreditFeeLegCoupon::TYPE = CClass::registerClassLoadMethod("CreditFeeLegCoupon", typeid(CreditFeeLegCoupon), CreditFeeLegCoupon_load);

RiskProperty_TYPES(CreditFeeLegCoupon)

DRLIB_END_NAMESPACE
