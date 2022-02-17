//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Swaption.cpp
//
//   Description : Swaption
//
//   Author      : Andrew J Swain
//
//   Date        : 1 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericSimpleIR.hpp"
#include "edginc/ClosedFormLN.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/Black.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VolProcessedBSIR.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/SensitiveIRVolPoints.hpp"
#include "edginc/SwapTool.hpp"
#include "edginc/B30360.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/SVGenIRSwap.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/ClosedFormIRLN.hpp"
#include "edginc/DateTime.hpp"
DRLIB_BEGIN_NAMESPACE

class Swaption: public GenericSimpleIR,
                public virtual CClosedFormLN::IIntoProduct,
                public virtual ClosedFormIRLN::IIntoProduct,
                public virtual ISensitiveIRVolPoints,
                public virtual LastSensDate,
                public virtual IMCIntoProduct{
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method = "Swaption::validatePop2Object";
    }

    virtual void Validate() {
        static const string method = "Swaption::Validate";
        try {
            if (!vol.get()) {
                throw ModelException(method, "no IR vol supplied");
            }

            if (swapStartDate.isLess(valueDate)) {
                throw ModelException(method,
                                     "swap start (" + swapStartDate.toString()+
                                     ") is before today ("+valueDate.toString()+
                                     ")");
            }

            if (optExpDate.isGreater(swapStartDate)) {
                throw ModelException(method,
                                     "option expiry (" + optExpDate.toString()+
                                     ") is after swap start date (" +
                                     swapStartDate.toString() + ")");
            }

            if (swapMatDate.isLess(swapStartDate)) {
                throw ModelException(method,
                                     "swapStartDate ("+swapStartDate.toString()+
                                     ") is > or = to swapMatDate (" +
                                     swapMatDate.toString() + ")");
            }

            if (Maths::isNegative(strikeRate)) {
                throw ModelException(method,
                                     "negative strike rate (" +
                                     Format::toString(strikeRate) + ")");
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** MonteCarlo::IntoProduct interface */
    IMCProduct* createProduct(const MonteCarlo* model) const;
    /** Swaption::IntoProduct interface */
    CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;
    ClosedFormIRLN::IProduct* createProduct(ClosedFormIRLN* model) const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Swaption, clazz);
        SUPERCLASS(GenericSimpleIR);
        IMPLEMENTS(CClosedFormLN::IIntoProduct);
        IMPLEMENTS(ClosedFormIRLN::IIntoProduct);
        IMPLEMENTS(ISensitiveIRVolPoints);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultSwaption);
        FIELD(swapStartDate, "start date of underlying swap");
        FIELD(swapMatDate, "maturity date of underlying swap");
        FIELD(fixedPayInterval, "payment interval of fixed leg");
        FIELD(fixedDCC, "day count convention of fixed leg");
        FIELD(stubType, "type of front or back stub");
        FIELD(stubAtEnd, "only matters if stubType != N");
        FIELD(strikeRate, "fixed interest strike level");

        FIELD(accrueBadDayConv, "Accrual date bad day convention");
        FIELD(payBadDayConv, "Payment date bad day convention");
        FIELD(resetBadDayConv, "Reset date bad day convention");
        FIELD_MAKE_OPTIONAL(accrueBadDayConv);
        FIELD_MAKE_OPTIONAL(payBadDayConv);
        FIELD_MAKE_OPTIONAL(resetBadDayConv);

        FIELD(optExpDate, "when option expires");
        FIELD(IsCall, "TRUE = receive fixed leg (call option)");
        FIELD(cashPhys, "T=PV wrt ytm, F=PV wrt zcurve");
        FIELD(swapNotional, "notional");
    }

    static IObject* defaultSwaption(){
        return new Swaption();
    }

private:
    friend class SwaptionClosedForm;

    Swaption():GenericSimpleIR(TYPE), accrueBadDayConv("NONE"),
                                      payBadDayConv("NONE"),
                                      resetBadDayConv("NONE") {};
    Swaption(const Swaption& rhs);
    Swaption& operator=(const Swaption& rhs);

    double annuityScale(double parFixedRate,
                        double fixedAnnuity) const {
        static const string method = "Swaption::annuityScale";
        try  {
            // The annuity we multiply by depends on whether it's
            // cash-settled or not. Currently, we use yield to maturity
            // for cash-settled deals, and zero curve for physical maturity
            // settled deals
            if (!cashPhys) {
                return fixedAnnuity;
            }
            else {
                // cash settled
                // Compute ytmAnnuity discounted to startDate using simple bond
                // maths. Note that we ignore holidays because in practice an
                // HP-12 calculator is used which also ignores holidays.
                // Hence we need to re-calculate the cash flows at this point.
                B30360 dcc;

                CashFlowArray cfl = SwapTool::cashflows(swapStartDate,
                                                        swapMatDate,
                                                        false, // stub @ start
                                                        1.0,
                                                        1,
                                                        fixedPayInterval,
                                                        &dcc);
                // remove principal
                if (cfl.empty()) {
                    throw ModelException(method,
                                         "no swap cashflows");
                }
                cfl[cfl.getLength()-1].amount -= 1.0;

                // get discount factor @ swap start date
                double startZeroPrice = discount->pv(swapStartDate);

                double         ytmAnnuity = 0.0;
                MaturityPeriod interval(1, fixedPayInterval);
                int            basis = (int)floor(0.5 + 1.0/interval.toYears());

                for (int i = 0; i < cfl.size(); i++) {
                    double years = dcc.years(swapStartDate, cfl[i].date);

                    double pv = RateConversion::rateToDiscountYearFrac(parFixedRate,
                                                                       years,
                                                                       basis);

                    ytmAnnuity += cfl[i].amount * pv;
                }

                return ytmAnnuity * startZeroPrice;
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    void price(ClosedFormIRLN* model, Control* control, CResults* results)const{
        static const string method = "Swaption::price";

        try  {
            double value = 0.0;
            DayCountConventionSP dcc(DayCountConventionFactory::make(fixedDCC));
            StubSP               stub(StubFactory::make(stubType));
            BadDayConventionSP   accrueBDC(BadDayConventionFactory::make(accrueBadDayConv));
            BadDayConventionSP   payBDC(BadDayConventionFactory::make(payBadDayConv));
            MaturityPeriod       interval(1, fixedPayInterval);
            HolidaySP            hols(Holiday::noHolidays());

            if (optExpDate.isGreaterOrEqual(valueDate)) {
                // Make cash flow list for fixed side of swap, with coupon = 1
                CashFlowArray fixedCFL = SwapTool::cashflows(swapStartDate,
                                                             swapMatDate,
                                                             stub.get(),
                                                             stubAtEnd,
                                                             accrueBDC.get(),
                                                             payBDC.get(),
                                                             hols.get(),
                                                             false, // dump start date
                                                             false, // add principal
                                                             1.0,
                                                             1,
                                                             fixedPayInterval,
                                                             dcc.get());

                // now pv it
                int i;
                double fixedAnnuity = 0.0;
                for (i = 0; i < fixedCFL.size(); i++) {
                    fixedAnnuity += discount->pv(fixedCFL[i].date)*fixedCFL[i].amount;
                }

                // paranoia
                if (Maths::isZero(fixedAnnuity)) {
                    throw ModelException(method, "fixed cashflows have no value");
                }

                // get par fixed rate
                double parFixedRate = discount->couponRate(swapStartDate,
                                                           swapMatDate,
                                                           interval,
                                                           stubAtEnd,
                                                           dcc.get());

                // compute option price. Note that a call is the option to receive
                // fixed. This is a *put* on the parfixed rate
                double pv = annuityScale(parFixedRate, fixedAnnuity);

                MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(swapMatDate,
                                                                    swapStartDate));

                CVolRequestSP volRequest(new SwapMaturityVolRequest(tenor.get()));
                CVolProcessedSP volCurve(vol->getProcessedVol(volRequest.get(), 0));
                CVolProcessedBS& volBS = dynamic_cast<CVolProcessedBS&>(*volCurve);
                double variance = volBS.CalcVar(valueDate, optExpDate);

                if (model) {
                    value = model->black(!IsCall,
                                         parFixedRate,
                                         strikeRate,
                                         pv,
                                         variance,
                                         vol.get());
                }
                else {
                    value = Black::price(!IsCall,
                                         parFixedRate,
                                         strikeRate,
                                         pv,
                                         variance);
                }

                value *= swapNotional;

                if (control && control->isPricing()) {
                    OutputRequest* request =
                        control->requestsOutput(OutputRequest::IND_VOL);
                    if (request) {
                        const double IMP_VOL_TOL = 0.0001;
                        double tradTime = volBS.calcTradingTime(valueDate, optExpDate);
                        double varTol = IMP_VOL_TOL*IMP_VOL_TOL*tradTime;

                        // handle case of using 2Q smile by backing out vol from price
                        double impVar;

                        if (Black::impliedVariance(!IsCall,
                                                   parFixedRate,
                                                   strikeRate,
                                                   pv,
                                                   variance,
                                                   value/swapNotional,
                                                   varTol,
                                                   impVar)) {

                            if (!Maths::isZero(tradTime)) {
                                double ivol = sqrt(impVar/tradTime);
                                results->storeRequestResult(request, ivol);
                            }
                        }
                    }

                    request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                    if (request) {
                        // know payment if cash settled & today's the last day
                        if (optExpDate.equals(valueDate) && cashPhys) {
                            CashFlow cf(optExpDate, value);
                            CashFlowArray cfl(1, cf);
                            OutputRequestUtil::recordKnownCashflows(control,
                                                                    results,
                                                                    discount->getCcy(),
                                                                    &cfl);
                        }
                    }

                    request = control->requestsOutput(OutputRequest::FWD_AT_MAT);
                    if (request) {
                        // Report the forward swap rate against the yield curve name
                        results->storeRequestResult(request, parFixedRate, discount->getName());
                    }
                }
            }

            results->storePrice(value, discount->getCcy());

            if (control && control->isPricing()) {
                OutputRequest* request =
                    control->requestsOutput(OutputRequest::PAYMENT_DATES);
                if (request) {
                    DateTimeArray paydate(1, optExpDate);
                    OutputRequestUtil::recordPaymentDates(control,
                                                          results,&paydate);
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    DateTime endDate(const Sensitivity* sensControl) const {
        return swapMatDate;
    }

    virtual IRGridPointAbsArraySP getSensitiveIRVolPoints(
        OutputNameConstSP outputName,
        const IModel*      model) const {
        static const string method = "Swaption::getSensitiveIRVolPoints";
        try {
            MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(swapMatDate,
                                                                swapStartDate));

            SwapMaturityVolRequest volRequest(tenor.get());
            DateTimeArray dates(1, optExpDate);
            return VolProcessedBSIR::sensitiveIRVolPoints(vol.get(),
                                                          discount.get(),
                                                          &volRequest, dates);
#if 0
            IRGridPointArraySP points(new IRGridPointArray(0));
            MaturityPeriodSP tenor(MaturityPeriod::dateSubtract(swapMatDate,
                                                                swapStartDate));

            CVolRequestSP volRequest(new SwapMaturityVolRequest(tenor.get()));
            DateTimeArray dates(1);
            dates[0] = optExpDate;

            vol->sensitiveIRVolPoints(volRequest.get(),
                                      outputName,
                                      dates,
                                      points);
            return points;
#endif
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

private:
    class MC;
    friend class MC;
    DateTime swapStartDate;       // start date of underlying swap
    DateTime swapMatDate;         // maturity date of underlying swap
    string   fixedPayInterval;    // payment interval of fixed leg
    string   fixedDCC;            // day count convention of fixed leg
    string   stubType;            // type of front or back stub
    bool     stubAtEnd;           // only matters if stubType != NONE
    double   strikeRate;          // fixed interest strike level
    string   accrueBadDayConv;
    string   payBadDayConv;
    string   resetBadDayConv;
    DateTime optExpDate;          // when option expires
    bool     IsCall;              // TRUE = receive fixed leg (call option)
                                  // FALSE = pay fixed leg (put option)
    bool     cashPhys;            // T=PV wrt ytm, F=PV wrt zcurve
    double   swapNotional;        // notional
};

/******************************************************************************/
/********************* end of Swaption Class declaration **********************/
/******************************************************************************/

CClassConstSP const Swaption::TYPE = CClass::registerClassLoadMethod(
    "Swaption", typeid(Swaption), Swaption::load
);


/** private class */
class SwaptionClosedForm: public virtual CClosedFormLN::IProduct,
                          public virtual ClosedFormIRLN::IProduct {
private:
    const Swaption &opt;

public:
    SwaptionClosedForm(const Swaption &opt): opt(opt){}

    void price(CClosedFormLN* model,
               Control*       control,
               CResults*      results) const{
        opt.price(0, control, results);
    }

    void price(ClosedFormIRLN* model,
               Control*        control,
               CResults*       results) const{
        opt.price(model, control, results);
    }
};

/** Implementation of ClosedFormIRLN::IntoProduct interface */
CClosedFormLN::IProduct* Swaption::createProduct(CClosedFormLN* model) const {
    return new SwaptionClosedForm(*this);
}

ClosedFormIRLN::IProduct* Swaption::createProduct(ClosedFormIRLN* model) const {
    return new SwaptionClosedForm(*this);
}

///////////// Experimental, FIXME: delete
class SwaptionPathPayoff : public IPayoffEvent {
    double  val;
    public:
        SwaptionPathPayoff(double p) : val(p) {}
        double  getVal() const {return val;}
};
///////////// End of experimental


/* MC product class for Swaption */
class Swaption::MC : public MCProductClient,
                     public IMCStatelessProductClient {
private:
    SVGenIRSwap::IStateVarSP     swapSV;  // swap state variable
    SVDiscFactorSP dfSV;    // df state variable
    double                    strike;  // from instrument
    double                    notional; // from instrument
    bool                      isCall;   // from instrument
    SVGenIRSwapSP                swapGen; // generator for swap
    SVGenDiscFactorSP            dfGen;   // generator for discount factors
    DateTime                  optExpDate; // For new payoff interface
protected:
    /** Override default method on IMCProduct. This method is called every time
        the path generator is changed (which is, at the moment, when the
        past path generator is created, and then when the future path
        generator is created  */
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        swapSV = swapGen->getIRSwapSV(newPathGen);
        dfSV = dfGen->getSVDiscFactor(newPathGen);
    };
public:
    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        svCollector->append(swapGen.get()); // ask for a swap one
        svCollector->append(dfGen.get()); // and a DiscFactor one
    }

    /** equivalent to InstIntoMCProduct. Need to call parent's constructor */
    MC(const Swaption*          inst,
       const SimSeriesSP&       simSeries,
       InstrumentSettlementSP   instSettle):
        MCProductClient(IMultiMarketFactors::asMulti(inst->coupon).get(),
                        inst->valueDate,
                        inst->discount.get(),
                        IRefLevelSP(IRefLevel::Util::
                                    makeZero(inst->valueDate)), // fix!
                        simSeries, // fix!
                        IPastValuesSP(IPastValues::Util::
                                      makeTrivial(inst->valueDate, 0.0)), // fix
                        instSettle.get(),
                        inst->optExpDate),
        strike(inst->strikeRate), notional(inst->swapNotional),
        isCall(inst->IsCall),
        swapGen(new SVGenIRSwap(inst->coupon,
                             // hard coded to use (Z0-Zn)/A
                             YieldCurveConstSP(), // check for ccy basis?
                             inst->swapStartDate, inst->swapMatDate,
                             inst->fixedPayInterval, inst->fixedDCC,
                             inst->stubType, inst->stubAtEnd,
                             inst->accrueBadDayConv, inst->payBadDayConv,
                             HolidaySP(Holiday::noHolidays()), inst->cashPhys)),
        dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                               inst->optExpDate)),
        optExpDate(inst->optExpDate) { }

    /** Called within the simulation loop */
    virtual void payoff(const IPathGenerator*  pathGen,
                        IMCPrices&                prices)
    {
        payoffHelper( prices );
    }

    virtual double pvFromPaymentDate() const{
        return 1.0; // already applied pv factor in payoff
    }

    // IMCStatelessProductClient

    virtual IHistoricalContextSP createHistoricalContext()
    {
        // No path dependency
        return IHistoricalContextSP(   );
    }

    virtual IHistoricalContextSP getInitialHistoricalContext()
    {
        return IHistoricalContextSP(   );
    }

    virtual DateTimeArray getPastDates()
    {
        return DateTimeArray();
    }

    virtual vector<int> finalize( const DateTimeArray& simDates )
    {
        const DateTimeArray& D = simDates;
        vector<int> keyDate( 1 );
        array<DateTime>::const_iterator I;
        I = std::find( D.begin(), D.end(), optExpDate );
        if ( I == D.end() ) {
            throw ModelException("Could not find date '"
                + optExpDate.toString() + "' in the timeline");
        }
        keyDate[0] = I - D.begin();
        return keyDate;
    }

    virtual void statelessPayOff(
        int currentDateIdx,
        IHistoricalContextSP history,
        IMCPrices& prices )
    {
        if ( currentDateIdx != 0 ) {
            throw ModelException("Date index '"
                + Format::toString( currentDateIdx )
                + "' is invalid.  Swaption should only have 1 date");
        }
        payoffHelper( prices );
    }

private:
    void payoffHelper( IMCPrices& prices ) const
    {
        ASSERT( optExpDate == paymentDate );
        double parYield = 0, annuity = 0;
        // calculate rate for a par swap and what coupons are worth for a rate
        // of 1.0
        swapSV->parYield(parYield, annuity);
        double myPayoff = Maths::max(0.0, isCall?
                                     (strike - parYield): (parYield - strike));
        // then multiple reduced rate by annuity
        myPayoff *= annuity;
        myPayoff *= dfSV->firstDF(); // and then discount to today
        prices.add(notional * myPayoff); // finally scale by notional

    }

    int maturityIdx;
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* Swaption::createProduct(const MonteCarlo* model) const {
    // the simSeries passed to the IMCProduct is redundant
    // well not quite. The MC wants to know the last sim date
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(DateTimeArray(1, optExpDate));
    InstrumentSettlementSP instSettle(new CashSettlePeriod(0));
    return new MC(this, simSeries, instSettle);
}

// for class loading
bool SwaptionLoad() {
    return (Swaption::TYPE != 0);
}

    

DRLIB_END_NAMESPACE
