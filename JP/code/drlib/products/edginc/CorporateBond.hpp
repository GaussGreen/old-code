//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CorporateBond.hpp
//
//   Description : Corporate Bond
//
//   Author      : André Segger
//
//   Date        : 09 September 2003
//
//
//----------------------------------------------------------------------------

#ifndef CORPORATE_BOND_HPP
#define CORPORATE_BOND_HPP
#include "edginc/Instrument.hpp"
#include "edginc/Generic1FactorCredit.hpp"
#include "edginc/ClosedForm.hpp"
#include "edginc/ClosedFormCDSPS.hpp"
#include "edginc/ClosedFormFA.hpp"
#include "edginc/ClosedFormCDSPSandFA.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/DeltaToCredit.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/BootstrappedYieldCurve.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/Bond.hpp"
#include "edginc/FD1FGeneric.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/InstrumentAsAsset.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL CorporateBond: public Generic1FactorCredit,
                     virtual public ClosedFormCDSPS::IIntoProduct,
                     virtual public ClosedFormFA::IIntoProduct,
                     virtual public ClosedFormCDSPSandFA::IIntoProduct,
                     virtual public FD1FGeneric::IIntoProduct,
                     virtual public LastSensDate,
                     virtual public DeltaToCredit::IShift,
                     virtual public ObjectIteration::IOverride, 
                     virtual public IGetMarket,
                     virtual public IInstrumentAsAsset{
public:
    static CClassConstSP const TYPE;
    friend class CorporateBondHelper;
    friend class CorporateBondClosedFormCDSPS;
    friend class CorporateBondClosedFormFA;
    friend class CorporateBondClosedFormCDSPSandFA;
    friend class BondProd;
    friend class BondFDProd;

    static const double CDS_ADJUSTMENT;

    virtual void validatePop2Object(); 

    virtual void Validate();

    /** Returns all known cash flows */
    CashFlowArraySP knownCashFlows()const;

    /** when do payments occur ? */
    DateTimeArraySP paymentDates() const;

    /** used to override CREDIT_SPREAD_RHO on cdsParSpreads
        ObjectIteration::IOverride interface */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** instrument specific market data handling */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Implementation of ClosedFormCDSPS::IntoProduct interface */
    virtual ClosedFormCDSPS::IProduct* createProduct(ClosedFormCDSPS* model) const;

    /** Implementation of ClosedFormFA::IntoProduct interface */
    virtual ClosedFormFA::IProduct* createProduct(ClosedFormFA* model) const;

    /** Implementation of ClosedFormCDSPSandFA::IntoProduct interface */
    virtual ClosedFormCDSPSandFA::IProduct* createProduct(ClosedFormCDSPSandFA* model) const;

    /** Implementation of FD1FGeneric::IntoProduct interface */
    virtual FD1FGeneric::IProduct* createProduct(FD1FGeneric* model) const;

    /** what's today ? */
    virtual DateTime getValueDate() const;

    /** when to stop tweaking */
    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(DeltaToCredit* shift) const;

    /** Shifts the object using given shift */
    virtual bool sensShift(DeltaToCredit* shift);

    double  calcThetaAccruedInterest() const;

    bool    getPutLevel(const DateTime& putDate,
                        double&         putLevel) const;

    bool    getCallLevel(const DateTime& callDate,
                         double&         callLevel) const;

    /** Part of IInstrumentAsAsset interface.
        Returns a date after which the instrument can no longer be used as
        an asset  */
    virtual DateTime maturityDate() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the 'coupons' or payments that this instrument will make during
        its lifetime. This can include historic payments. */
    virtual CashFlowArraySP getCoupons() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the dates on which the instrument has to be held in order to
        hold the right to the corresponding coupon as returned by 
        getCoupons() */
    virtual DateTimeArraySP getExCouponDates() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the yield curve used for discounting */
    virtual YieldCurveConstSP getDiscount() const;

    /** Part of IInstrumentAsAsset interface.
        Returns the accured interest (if any) to date */
    virtual double getAccrued() const;

private:
    CorporateBond();
    CorporateBond(const CorporateBond& rhs);
    CorporateBond& operator=(const CorporateBond& rhs);
    void priceParSpreads(CResults* results, Control* control) const;
    void priceFirmAsset(CResults* results, Control* control) const;
    void addRequests(Control*                    control, 
                     CResults*                   results, 
                     CashFlowArraySP             cleanSpreadCurve,
                     IObjectSP                   currentSpread) const;

    double parSpreadsValue(bool           priceE2C,
                           bool           generateRiskyStream,
                           ObjectArraySP& riskyStreamDetails)  const;

    double  priceCashflowStream(const DateTime&                       maturity,
                                CashFlowArrayConstSP                  cashFlows,
                                DefaultRatesSP                        defRates,
                                const auto_ptr<YieldCurve::IKey>&     discFactorKey) const;

    void calculateDefaultPayments(const DateTime&                       startDate,
                                  const DateTime&                       maturity,
                                  DateTimeArraySP                       critDates,
                                  CashFlowArrayConstSP                  cashFlows,
                                  DefaultRatesSP                        defRates,
                                  const auto_ptr<YieldCurve::IKey>&     discFactorKey,
                                  double&                               defaultPayment,
                                  double&                               accruedPayment) const;

    double getBondRecovery() const;

    YieldCurveSP getRiskyCurve() const;

    // Inputs
    ScheduleSP              callSchedule;
    bool                    callAdjustForAccrued;
    ScheduleSP              putSchedule;
    bool                    putAdjustForAccrued;
    double                  bondRecovery;
    bool                    useBondRecovery;
    bool                    payAccruedUponDefault;
    BondSP                  bond;
    bool                    koBeforeIssueDate;

    bool                    inDefault;
    DateTime                defaultDate;

    string                  dcc; // $unregistered
    string                  bdc;

    // Transient 
    mutable DayCountConventionSP bondAccrualDCC;          // used if payAccruedFee = true 
    mutable BadDayConventionSP   BDC;                    // bad day conv for BM adjustment 
    mutable bool                 createE2Csensitivities; // whether to create debt/equity sensitivities
    mutable bool                 e2cBasePriceCalculated;
    mutable double               e2cBasePrice;
};

typedef smartConstPtr<CorporateBond>  CorporateBondConstSP;
typedef smartPtr<CorporateBond>       CorporateBondSP;


DRLIB_END_NAMESPACE
#endif
