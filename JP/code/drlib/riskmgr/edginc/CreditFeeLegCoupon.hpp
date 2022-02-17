//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CreditFeeLegCoupon.hpp
//
//   Description : Tag class for shift in CDO fee leg coupon
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#ifndef DRLIB_FEE_LEG_COUPON_H
#define DRLIB_FEE_LEG_COUPON_H

#include "edginc/Void.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

struct RISKMGR_DLL CreditFeeLegCoupon: CObject {
    static CClassConstSP const TYPE;
    CreditFeeLegCoupon(); ~CreditFeeLegCoupon();
    typedef Void Qualifier;
    //Is continuous, i.e. tweaks to it can be made arbitrarily small
    enum { discrete = 0};
};

DRLIB_END_NAMESPACE

#endif
