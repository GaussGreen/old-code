//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : NullDecretionCurve.hpp
//
//   Description : A null decretion curve which always returns principal balance as 1
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#ifndef NULLDECRETIONCURVE_HPP
#define NULLDECRETIONCURVE_HPP

#include "edginc/config.hpp"
#include "edginc/DecretionCurve.hpp"

DRLIB_BEGIN_NAMESPACE
class MARKET_DLL NullDecretionCurve : public DecretionCurve {
public:
    static CClassConstSP const TYPE;

    NullDecretionCurve(const string& name);
    ~NullDecretionCurve();

    //------------------------------------------
    // override IDecretionCurve methods
    //------------------------------------------

    /** pv, which actually is balance here, can be relative to a start date */
    virtual double pv(const DateTime& startDate,
                      const DateTime& endDate) const;
    
    /** pv (balance) relative to initial balance */
    virtual double pv(const DateTime& endDate) const;

    /** return decretion speed on a date */
    virtual double getDecretionSpeed(const DateTime& date) const;

    /** return if balances are stepwise or continuous */ 
    virtual bool isStepBalances() const;

    /** Returns a reference to the value date */
    virtual double getFactor(const DateTime&) const;

    /**/
    virtual DateTimeArraySP getStepDates() const;
    virtual DateTimeArraySP getAlternativeStepDates() const;

    /** Returns settlement rule for prepay */
    virtual SettlementConstSP getSettlement() const;

protected:
    NullDecretionCurve(CClassConstSP clazz = TYPE);
    static void load(CClassSP& clazz);
    static IObject* defaultNullDecretionCurve();
    SettlementSP prepaySettle;
};

typedef smartConstPtr<NullDecretionCurve> NullDecretionCurveConstSP;
typedef smartPtr<NullDecretionCurve>      NullDecretionCurveSP;
typedef MarketWrapper<NullDecretionCurve> NullDecretionCurveWrapper;

DRLIB_END_NAMESPACE
#endif
