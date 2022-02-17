//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SwapLegIntFace.hpp
//
//   Description   Contains LiborLeg and Fixed Leg for cash flows.
//
//
//----------------------------------------------------------------------------

#ifndef EDG_LIBOR_LEG_HPP
#define EDG_LIBOR_LEG_HPP

#include "edginc/Class.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/FloatRate.hpp"
#include "edginc/Theta.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"
#include "edginc/StateVariable.hpp"
#include "edginc/StateVariableCollector.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/PayStream.hpp" // formerly just a class declaration - why?

DRLIB_BEGIN_NAMESPACE

const double tweak_spread_size = 0.0050;   //+0.50% for tweak.

class IBarrierPay;
class LiborLeg;
#ifndef QLIB_SWAPLEGINTFACE_CPP
EXTERN_TEMPLATE(class PRODUCTS_DLL_SP smartPtr<LiborLeg>);
#else
INSTANTIATE_TEMPLATE(class PRODUCTS_DLL smartPtr<LiborLeg>);
#endif

// libor leg
class PRODUCTS_DLL LiborLeg  : public CObject,
                  virtual public IGetMarket,
                  virtual public Theta::Shift
{
public:
    friend class LiborLegHelper;
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    //for PayStream
    DateTimeArray           RefixDates;
    DateTimeArray           AccrualDates;
    DateTimeArray           PayDates;
    CDoubleArray            Fixings;
    CDoubleArray            Spreads;
    CDoubleArray            Weights;
    CDoubleArray            Notionals;
    CBoolArray              Compounding;
    string                  PayDCC;
    bool                    CompoundFlat;

    //for FloatRate
    string                  RateDCC;
    string                  RateType;      //e.g. 3M or 6M....
    string                  BadDayConvention;

    //to get market data
    YieldCurveWrapper  couponCurve;
    DateTime           valueDate;

    LiborLeg();

    //** overrides default */
    virtual void validatePop2Object();

    int getSize();

    // get the last coupon payment date
    DateTime getLastPayDate();

    FloatRateSP getFloatRate(const YieldCurve* discount);


    CashFlowArrayConstSP getCashFlowArray();

    CashFlowArrayConstSP getCashFlowArray(const DateTime& valueDate,
                                          const YieldCurve* yield);

    CashFlowArraySP getCashFlowArray(const DateTime& valueDate,
                                     const DateTime& fromDate,
                                     const YieldCurve* yield);

    double getPV(const DateTime& valueDate, const YieldCurve* discount);

    // get pv of cashflows paying between fromDate & toDate as of fromDate
    double getPV(const DateTime&   valueDate,
                 const DateTime&   fromDate,
                 const DateTime&   toDate,
                 const YieldCurve* discount);

    // compute accrued interest as of given date
    double getAccrued(const DateTime&   valueDate,
                      const DateTime&   when,
                      const YieldCurve* yield);

    // next coupon payment date
    bool getNextPayDate(const DateTime&   when,
                        DateTime& nextPayDate);

    void setFixingforThetaShift(const DateTime& valueDate,
                                const YieldCurve* discount,
                                const DateTime& rollDate);

    // populate couponCurve from product, if it doesn't have own Coupon Curve
    void setCouponCurve(const YieldCurveWrapper instCouponCurve);

    // populate from market cache
    void getMarket(const IModel* model, const MarketData* market);

    // check market data is already set.
    void checkMarket();

    // prepare KNOWN_CASH_FLOW
    void makeKnownCashFlow();

    CashFlowArraySP getKnownCashFlows();

    // tweak spreads for client valuation purpose
    // isTweak = true, shift the spreads.  false recover the original values.
    void tweakSpreads(bool isTweak);

    // Implementation of the Theta shift interface
    bool sensShift(Theta* shift);

    // create a PayStream
    PayStream* makePayStream(const YieldCurve* discount);

    // support barrier pay
    refCountPtr<IBarrierPay> createBarrierPay(bool isOut, const string& koStubRule);

    /** Interface for the state variable for LiborLeg. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getLiborLegSV() method below */
    class PRODUCTS_DLL LiborLegSV: public virtual IStateVariable{
    public:
        LiborLegSV(smartPtr<LiborLeg>            libor, // would like const but getFloatRate & makePayStream stop me
                   PayStreamSP                   paystream,
                   YieldCurveConstSP             discount,
                   const DateTimeArray&          matDates,  // convenient
                   vector<SVGenExpectedDiscFactorSP> expDfGens,
                   SVGenDiscFactorSP                paymentDfGen,
                   IStateVariableGen::IStateGen* pathGen);

        virtual ~LiborLegSV();

        virtual bool doingPast() const;

        void update(IStateVariableGen::IStateGen* pathGen);

        // TARN calls this from constructor
        PayStreamSP getPayStream() const;

        // reference only - rates updated each iter
        const CashFlowArray* getCashFlowArray(const DateTime& valueDate) const;

        // present value of libor leg
        double getPV(const DateTime& valueDate) const;

        // PV factor from cf pay date to today (valueDate?)
        double payDatePV(int idx) const;

        // return libor leg
        smartPtr<LiborLeg> getLibor() const;

#if 0
        consider later
        setFixingforThetaShift; // called from sensShift(Theta* shift)
#endif
    private:
        LiborLegSV(const LiborLegSV& rhs);
        LiborLegSV& operator=(const LiborLegSV& rhs);
        class Imp; // hide implementation
        auto_ptr<Imp> me;
    };
    typedef smartPtr<LiborLegSV> LiborLegSVSP;

    /** A Generator of MC LiborLeg Variables  */
    class PRODUCTS_DLL LiborLegSVGen: public virtual IStateVariableGen,
        public virtual IStateVariableClient,
        public virtual VirtualDestructorBase
    {
    public:

        // Not quite sure what is needed here...
        LiborLegSVGen(smartPtr<LiborLeg> liborLeg,
                      YieldCurveConstSP  discount);

        virtual ~LiborLegSVGen();

        IStateVariableSP create(IStateVariableSP             oldStateVar,
                                IStateVariableGen::IStateGen* pathGen) const;

        /** Returns a LiborLeg state variable which then provides
            access to various rate-based methods. This is the method that
            products should call to get a LiborLeg::ILiborLeg. The
            previous IRefLevel::ILiborLegSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar.*/
        virtual LiborLegSVSP getLiborLegSV(LiborLegSVSP               oldStateVar,
                                           IStateVariableGen::IStateGen* pathGen) const;

        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. Implementations typically call
            IStateVariableCollector::append */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const;

    private:
        LiborLegSVGen(const LiborLegSVGen& rhs);
        LiborLegSVGen& operator=(const LiborLegSVGen& rhs);
        class Imp; // hide implementation
        auto_ptr<Imp> me;
    };
    typedef smartPtr<LiborLegSVGen> LiborLegSVGenSP;

    /** How to get a IStateVariableGen from a LiborLeg */
    virtual LiborLegSVGen* createLiborLegSVGen(YieldCurveConstSP discount) const;

private:
    CashFlowArraySP     knownCFL;          // knownCFL $unregistered
//    KOSettleSP          koSettle;
    bool                hasKOSettle;        // whether LiborLeg has KOSettle Rule or not. $unregistered

    CDoubleArray        orgSpreads;    // to keep the original spreads.
};

typedef smartPtr<LiborLeg> LiborLegSP;

// fixed leg
class PRODUCTS_DLL FixedLeg  : public CObject
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    friend class FixedLegBarrierPay;
    friend class FixedLegHelper;

    //field (better to move to Private.....)
    DateTimeArray       AccrueStartDates;  //FixedCoupon start date array
    DateTimeArray       AccrueEndDates;    //FixedCoupon end date array
    DateTimeArray       PaymentDatesArray; //FixedCoupon payment dates array
    CDoubleArray        CouponAmounts;     //fixed coupon amounts

    FixedLeg();

    //** overrides default */
    virtual void validatePop2Object();

    int getSize();

    // get the last coupon payment date
    DateTime getLastPayDate();

    CashFlowArrayConstSP getCashFlowArray();

    //retern Fixed Cash Flow payment after baseDate //
    CashFlowArraySP getCashFlowArray(const DateTime& baseDate);

    double getPV(const DateTime& valueDate, const YieldCurve* discount);

    // get pv of cashflows paying between fromDate & toDate as of fromDate
    double getPV(const DateTime&   valueDate,
                 const DateTime&   fromDate,
                 const DateTime&   toDate,
                 const YieldCurve* discount);

    // compute accrued fixed coupon as of given date
    double getAccrued(const DateTime&   valueDate,
                      const DateTime&   when,
                      const YieldCurve* discount);

    // compute accrued coupon amount (no PV)
    double getAccrued(const DateTime&   when);

    // next coupon payment date
    bool getNextPayDate(const DateTime&   when,
                        DateTime& nextPayDate);

    bool isDCC();

    // support barrier pay
    refCountPtr<IBarrierPay> createBarrierPay(bool isOut, const string& koStubRule);

private:
    CashFlowArraySP     knownCFL;          // knownCFL $unregistered
    bool                hasKOSettle;        // whether LiborLeg has KOSettle Rule or not. $unregistered

public:
    //field
    string              dcc;               //fixed couopn day count
};

typedef smartPtr<FixedLeg> FixedLegSP;
#ifndef QLIB_SWAPLEGINTFACE_CPP
EXTERN_TEMPLATE(class PRODUCTS_DLL_SP smartPtr<FixedLeg>);
#else
INSTANTIATE_TEMPLATE(class PRODUCTS_DLL smartPtr<FixedLeg>);
#endif

DRLIB_END_NAMESPACE
#endif
