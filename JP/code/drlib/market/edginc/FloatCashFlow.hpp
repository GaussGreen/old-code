//-------------------------------------------------------------
//
//   Group       : Credit Hybrids QR&D
//
//   Filename    : FloatCashFlow.hpp
//
//   Description : A floating rate coupon that pays over a specified
//                 accrual period
//
//   Author      : Gordon Stephens
//
//   Date        : 8 June 2005
//-------------------------------------------------------------

#ifndef FLOATCASHFLOW_HPP
#define FLOATCASHFLOW_HPP

#include "edginc/RiskFreeCashFlow.hpp"
#include "edginc/DayCountConvention.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CreditFeeLegCoupon.hpp"
#include "edginc/MQQuasiIRVolCMS.hpp"
#include "edginc/MQCMSPricer.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL FloatCashFlow : public RiskFreeCashFlow,
                                 public virtual IRestorableWithRespectTo<CreditFeeLegCoupon>,
                                 public virtual Theta::Shift,
                                 public virtual MQCMSPricer::IIntoProduct,
                                 public virtual ClosedFormForwardRatePricer::IIntoProduct

{
public:
    static CClassConstSP const TYPE;

    virtual bool                 isZero() const;
    virtual double               getAmount(const ZeroCurve& curve, 
                                           const IRVolBase* volModelIR) const;
    virtual void                 getCriticalDates(DateTimeArray& dates) const;
    virtual const DateTime       getMaturityDate() const;

	virtual double               getAmount(IForwardRatePricerSP model) const;
    virtual double               getNotional() const;
    virtual void                 setNotional(double newNotional);
    virtual double               getRate() const; //implicitly a fixed rate
    virtual void                 setRate(double newRate); //implicitly a fixed rate

    virtual DateTime             getAccrueStart() const;
    virtual DateTime             getAccrueEnd() const;
    virtual DayCountConventionConstSP getAccrualDcc() const;
    virtual bool                 isAdhoc() const;    

    /** Returns cash flow amount using given rate override (if relevant) */
    virtual double getAmount(double rateOverride) const;

    FloatCashFlow();
    FloatCashFlow(  DateTime             valueDate,
                    DateTime             payDate,
                    double               notional,
                    DateTime             accrueStart,
                    DateTime             accrueEnd,
                    DayCountConventionSP accrualDcc,
                    ExpirySP             rateTenor,
                    double               spread,
                    bool                 isCMS,
                    DateTime             refixDate,
                    DayCountConventionSP rateDcc,
                    BadDayConventionSP   rateBdc,
                    YieldCurveWrapper    rateCurrency,
                    double               weight,
                    double*              fixing,
                    int                  applyAdjustments);

    FloatCashFlow(
        const DateTime&           payDate, 
        double                    notional, 
        double                    spread, 
        bool                      additive,
        const DayCountConvention* dcc, 
        const DateTime&           accrueStartDate, 
        const DateTime&           accrueEndDate, 
        const DateTime&           rateStartDate, 
        const DateTime&           rateEndDate);

    FloatCashFlow(  DateTime                valueDate,
                    DateTime                payDate,
                    double                  notional,
                    DateTime                accrueStart,
                    DateTime                accrueEnd,
                    DayCountConventionSP    accrualDcc,
                    ExpirySP                rateTenor,
                    double                  spread,
                    bool                    isCMS,
                    DateTime                refixDate,
                    DayCountConventionSP    rateDcc,
                    BadDayConventionSP      rateBdc,
                    YieldCurveWrapper       rateCurrency,
                    double                  weight,
                    double*                 fixing,
                    int                     applyAdjustments,
                    bool                    isAdditive,
                    MaturityPeriodSP        rateFreq,
                    double                  capLmt,
                    double                  floorLmt
                    );

    ~FloatCashFlow();

    /** When to stop tweaking for Yield Curve type tweaks */
    virtual DateTime lastYCSensDate() const;
    //Theta::Shift methods
    bool sensShift(Theta* shift);

    //// Implementing tweakable with respect to CreditFeelegCoupon 
    //// - a shift in spread
    TweakOutcome sensShift(const PropertyTweak<CreditFeeLegCoupon>& shift);
    string sensName(const CreditFeeLegCoupon* shift) const;
    void sensRestore(const PropertyTweak<CreditFeeLegCoupon>& shift);

    //CObject methods
    virtual void validatePop2Object();
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the accrual period (start and end dates) corresponding to 
     *  this cash flow. */
    virtual AccrualPeriodSP getAccrualPeriod() const;
    
    //MQCMSPricer::IIntoProduct method
    virtual MQCMSPricer::IProduct* createProduct(MQCMSPricer* model) const;

    //ClosedFormForwardRatePricer::IIntoProduct method
    virtual ClosedFormForwardRatePricer::IProduct* createProduct(ClosedFormForwardRatePricer* model) const;

    //forward method
    double fwdRateMQ(MQCMSPricer* model) const;
    double fwdRateAlib(ClosedFormForwardRatePricer* model) const;

private:
    /** Updates fixings for theta shift (if needed) */
    void setFixingforThetaShift(int rollDate);

    //infrastructure support
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //fields
    double                  notional;         // the size of the cashflow
    DateTime                accrueStart;      // the start of the accrual period
    DateTime                accrueEnd;        // the end of the accrual period
    string                  accrualDcc;       // the accrual day count convention
    string                  rateTenor;        // the tenor of the floating rate
    double                  spread;           // additive to percentage of calculated rate
    bool                    isAdditive;       // true [default] if spread is added to rate
    CBoolSP                 isCMS;            // false=>Simple, true=>Compounding (at underlying ccy swap frequency)
    DateTime                refixDate;        // when the rate is observed
    string                  rateDcc;          // the rate day count convention
    string                  rateBdc;          // the rate bad day convention 
    YieldCurveWrapper       rateCurrency;     // the curve used to calculate the forward
    double                  weight;           // a scaling factor for the coupon amount
    double*                 fixing;           // if refixDate is in the past, the observed rate (spread will be added)
    int                     applyAdjustments; // if floating rate adjustments are required (>0)
    MaturityPeriodSP        rateFreq;         // rate reset frequency.  Default to rateTenor
    double                  capLmt;           // cap
    double                  floorLmt;         // floor

    //transient fields
    DayCountConventionConstSP _accrualDcc;    // the accrual day count convention
    ExpirySP                  _rateTenor;     // the tenor of the floating rate
    DayCountConventionConstSP _rateDcc;       // the rate day count convention
    BadDayConventionSP        _rateBdc;       // the rate bad day convention
    bool                      _isCMS;         // if not supplied, we look up _rateTenor in rateCurrency
                                              // and if MM then false, if SW then true, else failure

    DateTime                rateStart;
    DateTime                rateEnd;
};

typedef smartPtr<FloatCashFlow> FloatCashFlowSP;
typedef array<FloatCashFlowSP, FloatCashFlow>    FloatCashFlowArray;

DRLIB_END_NAMESPACE

#endif
