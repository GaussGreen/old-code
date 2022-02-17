//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CreditFeeLegCouponShift.hpp
//
//   Description : Sensitivity to a shift in CDO fee leg coupon
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CREDITF_EE_LEG_COUPON_SHIFT_H
#define QLIB_CREDITF_EE_LEG_COUPON_SHIFT_H

/* A sensitivity to a shift in fee leg coupon, used by managed trades
 * through ImpliedScalarShiftMulti.
 * Inheriting from AllNamesRiskPropertySensitivity, so it's a declarative greek.
 * Implementing IScalarPerNameShift, so it can be used for solving.
 */

#include "edginc/CreditFeeLegCoupon.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"


DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL CreditFeeLegCouponShift : public AllNamesRiskPropertySensitivity, 
                                            virtual public IScalarPerNameShift {
 
 public:
    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    CreditFeeLegCouponShift(double shiftSize = DEFAULT_SHIFT);
    
    //// To implement IScalarPerNameShift
    IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                            double shiftSize,
                                            IObjectSP world);
   
    OutputNameArrayConstSP allNames(const IObject* object) const;
    
 private:
     ~CreditFeeLegCouponShift();
    CreditFeeLegCouponShift(const CreditFeeLegCouponShift& rhs); 
    CreditFeeLegCouponShift& operator=(const CreditFeeLegCouponShift& rhs);
    static void load(CClassSP& clazz);
    double getShiftSize() const;
    static const double SENSITIVITY_UNIT;
    static const string NAME;
    Deriv deriv() const;
};

FORWARD_DECLARE(CreditFeeLegCouponShift)

DRLIB_END_NAMESPACE

#endif
