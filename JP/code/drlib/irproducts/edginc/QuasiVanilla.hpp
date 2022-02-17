//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : QuasiVanilla.hpp
//
//   Description : First attempt at implementing quasi-vanilla in qlib
//                 subject to change in the future
//
//----------------------------------------------------------------------------

#ifndef QUASI_VANILLA_HPP
#define QUASI_VANILLA_HPP
#include "edginc/Model.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CouponSched.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/QuasiMQ.hpp"
#include "edginc/KComponent.hpp"

DRLIB_BEGIN_NAMESPACE


/***************************************************************************************************************
*
* Define Univariate Interest Rate Coupon object.
*
* This defines the payment details of an instrument with a payoff contingent on one interest rate. The
* instrument makes a certain payout (defined by a payoff enum) at intervals defined by the coupon
* frequency. This instrument does NOT support a RIB schedule.
*
* Example : Cap on CMS rate. CMS floating leg.
*
* Samy A Mohammed (22 August 2006)
***************************************************************************************************************/
class UnivarIRCpn : public KComponent,
                    virtual public QuasiMQ::IIntoProduct
                    // virtual public anyModel::IIntoProduct
{
public:

    struct UnivarIRCpnProdMQ : virtual QuasiMQ::Product
    {
        UnivarIRCpn const *inst;
        virtual void price(Control*   control,
                           CResults*  results);
        UnivarIRCpnProdMQ(QuasiMQ *model, UnivarIRCpn const *inst)
            : QuasiMQ::Product(model), inst(inst) {}
    };
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultConstructor(void) { return new UnivarIRCpn(); }

    struct PayoffType {
        enum Enum {CALL, PUT, BINARY_CALL, BINARY_PUT, INSIDE_RIB, OUTSIDE_RIB};
        // ANNUITY_CALL, ANNUITY_PUT, TEC_CALL, TEC_PUT);
    };
    typedef BoxedEnum<PayoffType::Enum> PayoffTypeBoxedEnum;

    UnivarIRCpn() : KComponent(TYPE) {}

    // CInstrument interface
    virtual void validatePop2Object();
    virtual void setup(const IModel* model, const MarketData *market)
    {
        try
        {
            rate->setup(model, market);
        }
        catch (exception& e)
        {
            throw ModelException(e, __FUNCTION__);
        }
    }

    // QuasiMQ::IIntoProduct interface
    virtual QuasiMQ::ProductSP createProduct(QuasiMQ* model) const;
    //UnivarIRCpn::UnivariateIRCouponProd* createProduct(QuasiMQ* model) const;

    ///////// Exported Fields /////////

    DateTime                start;              // Accrual Start date
    DateTime                end;                // Accrual End date
    MaturityPeriodSP        couponFreq;         // accrual period e.g. 1Y
    DayCountConventionSP    accrualDCC;         // dcc for coupon accruals NOT for the underlying payment rate CMS
    DayCountConventionSP    expiryDCC;          // dcc for coupon accruals NOT for the underlying payment rate CMS
    BadDayConventionSP      resetBadDayConv;    // BDC for rate reset date
    BadDayConventionSP      accrualBadDayConv;  // BDC for rate accrual st/end dates
    BadDayConventionSP      paymentBadDayConv;  // BDC for coupon payment date
    IndexSpecIRSP           rate;               // Define rate to observe - include rate start/end
    bool                    arrears;            // Defines if coupon rate resets in arrears
    PayoffType::Enum        payoffType;         // Defines payoff
    DoubleArray             payoffParams;       // e.g. leverage, spread, etc.
    StringArray             payoffParamsLabel;  // ???
    MaturityPeriodSP        paymentDelay;       // delay from each internal coupon reset date
    string                  settlementType;     // Cash/Physical
    string                  stub;               // 'F'ront or 'Back'
    DateTimeArray           amortisationDates;  // array of step up dates
    DoubleArray             notionals;          // array of notional steps
    DateTimeArray           strikeStepUpDates;  // array of strike step up dates
    DoubleArray             strikes;            // array of strike steps
    CouponSchedDatesSP      sched;              //Accrual st, end, pay
};

typedef smartPtr<UnivarIRCpn> UnivarIRCpnSP;
typedef smartConstPtr<UnivarIRCpn> UnivarIRCpnConstSP;



/***************************************************************************************************************
*
* Define Univariate Interest Rate RIB Coupon object.
*
* This defines the payment details of an instrument with fixed RIB coupon payment (observe on
* floating rate, pay fixed rate).
*
* Example : Cap on CMS rate.
*
* Samy A Mohammed (22 August 2006)
*************************************************************************************************************
class UnivarIRRibCpn :  public KComponent,
                        virtual public QuasiMQ::IIntoProduct
{
public:

    struct UnivarIRRibCpnProd2Q : virtual QuasiMQ::Product
    {
        UnivarIRRibCpn const * inst;

        virtual void price(Control*   control,
                           CResults*  results);
        UnivarIRRibCpnProd2Q(QuasiMQ *model, UnivarIRRibCpn const * inst)
            : QuasiMQ::Product(model), inst(inst) {}
    };
    static CClassConstSP const TYPE;
    static void load(CClassSP& );
    static IObject* defaultConstructor(void) { return new UnivarIRRibCpn(); }

    struct PayoffType {
        enum Enum {CALL, PUT, BINARY_CALL, BINARY_PUT, INSIDE_RIB,
                   OUTSIDE_RIB};  // ANNUITY_CALL, ANNUITY_PUT, TEC_CALL, TEC_PUT};
    };

    UnivarIRRibCpn() : KComponent(TYPE) {}

    // CInstrument interface
    virtual void validatePop2Object();
    virtual void Validate() {}

	// QuasiMQ::IIntoProduct interface
    virtual QuasiMQ::ProductSP createProduct(QuasiMQ* model) const;
    //UnivarIRRibCpn::UnivariateIRCouponProd* createProduct(QuasiMQ* model) const;

    ///////// Exported Fields /////////

    DateTime                start;              // Accrual Start date
    DateTime                end;                // Accrual End date
    MaturityPeriodSP        couponFreq;         // accrual period e.g. 1Y
    MaturityPeriodSP        obsFreq;            // Float idx observation freq, e.g 1M
    DayCountConventionSP    accrualDCC;         // dcc for coupon accruals NOT for the underlying payment rate CMS
    DayCountConventionSP    expiryDCC;          // dcc for coupon accruals NOT for the underlying payment rate CMS
    BadDayConventionSP      resetBadDayConv;    // BDC for rate reset date
    BadDayConventionSP      accrualBadDayConv;  // BDC for rate accrual st/end dates
    BadDayConventionSP      paymentBadDayConv;  // BDC for coupon payment date
    IndexSpecIRSP           rate;               // Define rate to observe - include rate start/end
    bool                    arrears;            // Defines if coupon rate resets in arrears
    PayoffTypeEnum          payoffType;         // Defines payoff
    DoubleArray             payoffParams;       // e.g. leverage, spread, etc.
    StringArray             payoffParamsLabel;  // ???
    MaturityPeriodSP        paymentDelay;       // delay from each internal coupon reset date
    string                  settlementType;     // Cash/Physical
    string                  stub;               // 'F'ront or 'Back'
    DateTimeArray           amortisationDates;  // array of step up dates
    DoubleArray             notionals;          // array of notional steps
    DateTimeArray           strikeStepUpDates;  // array of strike step up dates
    DoubleArray             strikes;            // array of strike steps
};

typedef smartPtr<UnivarIRRibCpn> UnivarIRRibCpnSP;
typedef smartConstPtr<UnivarIRRibCpn> UnivarIRRibCpnConstSP;
*/
DRLIB_END_NAMESPACE

#endif
