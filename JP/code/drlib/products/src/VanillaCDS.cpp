//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : VanillaCDS.cpp
//
//   Author      : Doris Morris
//
//   Date        : February 2006
//
//
//----------------------------------------------------------------------------
#define QLIB_VANILLACDS_CPP

#include "edginc/config.hpp"
#include "edginc/VanillaCDS.hpp"
#include "edginc/StubFactory.hpp"
#include "edginc/VanillaCreditFeeLeg.hpp"
#include "edginc/VanillaCreditContingentLeg.hpp"
#include "edginc/Results.hpp"
#include "edginc/Settlement.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/IForwardRatePricer.hpp"

//#include "edginc/BadDayConventionFactory.hpp"
//#include "edginc/DayCountConventionFactory.hpp"
//#include "edginc/YieldCurve.hpp"
//#include "edginc/ICDSParSpreads.hpp"
//#include "edginc/ICreditFeeLeg.hpp"
//#include "edginc/IRiskyContingentLeg.hpp"

DRLIB_BEGIN_NAMESPACE

// Closed form product
class VanillaCDSClosedFormCDSPS: public ClosedFormCDSPS::IProduct
{
private:
    const VanillaCDS* cfCDS;

public:
    VanillaCDSClosedFormCDSPS(const VanillaCDS* cfCDS): cfCDS(cfCDS)
    {}

    void price(ClosedFormCDSPS* model,
               Control*         control,
               CResults*        results) const
    {
        // retrieve the model to be used in the calculation
        // of fee leg forward rates
        IHasForwardRatePricer* ihfrp = dynamic_cast<IHasForwardRatePricer*>(model);
        if (!ihfrp)
        {
            throw ModelException(
                "VanillaCDSClosedFormCDSPS::price",
                "Model must implement IHasForwardRateModel");
        }
        IForwardRatePricerSP frModel = ihfrp->getForwardRatePricer();

        //now price
        cfCDS->price(results, control, frModel);
    }
};

//-------------------
// VanillaCDS methods
//-------------------

VanillaCDS::~VanillaCDS(){}

/** Price the instrument
 *  Price (dirty) = - IRiskyFeeLeg.getFeeLegPV(DateTime valuationDate, IDiscountCurveRisky)
                    + IRiskyContingentLeg.getContingentLegPV(DateTime valuationDate, IDiscountCurveRisky)
                    since positive notional = long protection
 *  The notional for each fee interval can be varying.
 *  However varying notional has yet to be implemented.
 *  The Vanilla CDS will retrieve its appropriate IDiscountCurve risky
 *  The price of the instrument is stored in CResults to fit into the QLib framework.
 */
void VanillaCDS::price(CResults* results, Control* control, IForwardRatePricerSP model) const
{
    static const string method = "VanillaCDS::price";
    double price = 0.0;
    double defaultedPrice = 0.0;
    double ai = 0.0;

    try {

        /* see if there is a default */
        if(riskyCurve->defaulted())
        {
            const DateTime& defaultDate = riskyCurve->getDefaultDate();
            defaultedPrice = getPVGivenDefault(defaultDate, defaultDate, *(riskyCurve.getSP()), model);
            results->storePrice(defaultedPrice, discount->getCcy());        
        }
        else
        {
            OutputRequest* request;

            // retrieve the pre-payment curve
            IDecretionCurveConstSP prepay = riskyCurve->getPrepayCurve();

            /* Find the settlement date */        
            DateTime settlementDate;
            settlementDate = settlementOffset->settles(valueDate);

            /* Adjust the settlement date for holidays */
            DateTime adjSettlementDate;
            adjSettlementDate = pmtBdc->adjust(settlementDate, settlementOffset->getMarketHolidays().get());

            /* Calculate the pv as of valueDate conditional on on default to adjSettlementDate */
            price = getPV(valueDate, adjSettlementDate, *(riskyCurve.getSP()), prepay, model, 
                          IBadDayAdjusterConstSP::attachToRef(this));
            results->storePrice(price, discount->getCcy());

            // ACCRUED_INTEREST
            request = control->requestsOutput(OutputRequest::ACCRUED_INTEREST);
            if(request)
            {
                ai = getAccruedInterest(adjSettlementDate, model);

                results->storeRequestResult(request, ai);
            }
            
            // If pricing, produce the output requests too (ie, not if we
            // are tweaking)
            if (control && control->isPricing()) 
            {
                addOutputRequests(control, results, prepay, model);
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//--------------------
// CInstrument methods
//--------------------

/**
 * Called immediately after object constructed
 */
void VanillaCDS::validatePop2Object()
{
    static const string method("VanillaCDS::validatePop2Object");

    if(effectiveDate > maturityDate)
    {
        throw ModelException(method, "protection start date (" +
                             effectiveDate.toString() + 
                             ") should be before protectionEndDate (" +
                             maturityDate.toString() + ").");
    }
}

/*  Allow the instrument to retrieve its market data
    Called by riskmgr before validate() */
void VanillaCDS::GetMarket(const IModel* model,
                           const CMarketDataSP market)
{
    const static string method = "VanillaCDS::getMarket";

    market->GetReferenceDate(valueDate);

    /* Get the interest rate discount curve */
    if(!discount.isEmpty())
    {
        discount.getData(model, market);
    }

    /* Get the credit risky curve */
    if(!riskyCurve.isEmpty())
    {
        riskyCurve.getData(model, market);
    }

    if(!pmtHol.isEmpty())
    {
        pmtHol.getData(model, market);
    }
    else
    {
        // default it to weekends only
        pmtHol.setObject(MarketObjectSP(Holiday::weekendsOnly()));
    }

    if(!accrualHol.isEmpty())
    {
        accrualHol.getData(model, market);
    }
    else
    {
        // default it to weekends only
        accrualHol.setObject(MarketObjectSP(Holiday::weekendsOnly()));
    }

}

/** Called after market data has been retrieved */
void VanillaCDS::Validate()
{
    static const string method = "VanillaCDS::validate";


   /* Set the optional fields */
   if(!accrualDcc)
   {
       accrualDcc = DayCountConventionSP(pmtDcc);
   }

   if(!accrualBdc)
   {
       accrualBdc = BadDayConventionSP(pmtBdc);
   }

   if(!valuationBdc)
   {
       valuationBdc = BadDayConventionSP(pmtBdc);
   }
   
   /* make Stub payment type */
   StubSP stub(StubFactory::make("Simple"));
   stubPaymentType = stub;

   /* check delay, if not set*/

   /* create a fee leg*/
   VanillaCreditFeeLegSP newFeeLeg(new VanillaCreditFeeLeg(baseNotional,
                                                           valueDate,
                                                           effectiveDate,
                                                           maturityDate,
                                                           feeRate,
                                                           paymentFreq,
                                                           fcd,
                                                           lcd,
                                                           pmtDcc,
                                                           accrualDcc,
                                                           pmtBdc,
                                                           accrualBdc,
                                                           valuationBdc,
                                                           stubPaymentType,
                                                           pmtHol,
                                                           accrualHol,
                                                           discount,
                                                           payAccruedFee,
                                                           delay));

   feeLeg = newFeeLeg;

   /* find the protection end date */   
   DateTime protectionEndDate = feeLeg->getRiskyCashFlowDates()->back();

   /* create a contingent leg, whose protection end date is the last date of the fee payments */
   VanillaCreditContingentLegSP newContingentLeg (new VanillaCreditContingentLeg(valueDate,
                                                                                 effectiveDate,
                                                                                 protectionEndDate,
                                                                                 baseNotional,
                                                                                 recoveryRate,
                                                                                 delay,
                                                                                 discount,
                                                                                 pmtBdc));

  
   contingentLeg = newContingentLeg;



}

DateTime VanillaCDS::getValueDate() const
{
    return valueDate;
}

/** Returns the name of the instrument's discount currency. */
string VanillaCDS::discountYieldCurveName() const
{
    return discount.getName();
}

//-------------
// ICDS methods
//-------------

/** Return the fee leg */
ICreditFeeLegSP VanillaCDS::getFeeLeg() const
{
    return feeLeg;
}

/** Return the contingent leg */
ICreditContingentLegSP VanillaCDS::getContingentLeg() const
{
    return contingentLeg;
}

/**
 * Compute PV of contingent leg at valuation date, with instrument default
 * settlement and payment date behaviour.
 */
double VanillaCDS::getContingentLegPV(const DateTime&              valuationDate,
                                      const IDiscountCurveRisky&   crv) const
{
    return contingentLeg->getContingentLegPV(valuationDate, crv,
                                             IBadDayAdjusterConstSP::attachToRef(this));
}


/**
 * Compute PV of contingent leg at valuationDate, but for unconditional payment
 * at paymentDate, conditional on no new defaults before valuationDate.
 */
double VanillaCDS::getContingentLegPV(const DateTime&              valuationDate,
                                      const DateTime&              paymentDate,
                                      const IDiscountCurveRisky&   crv) const
{
    return contingentLeg->getContingentLegPV(valuationDate, paymentDate, crv, 
                                             IBadDayAdjusterConstSP::attachToRef(this));
}

double VanillaCDS::getAccruedInterest(const DateTime& settlementDate,
                                      IForwardRatePricerSP model) const
{
    double ai = feeLeg->getAccruedInterest(settlementDate);

    /* positive notional = long protection */
    return -ai;
}

/** Calculates the PV of an instrument on defaultDate, given that default has just occured.
 *  Note that the curve is needed in case there is a delay between default and recovery
 *  payment. This method assumes a complete default, so should be interpreted with care
 *  when the instrument/curve has multiple names.
 */
double VanillaCDS::getPVGivenDefault(const DateTime&            valuationDate,
                                     const DateTime&            defaultDate,
                                     const IDiscountCurveRisky& crv,
                                     IForwardRatePricerSP       model) const
{
    const static string method = "VanillaCDS::getPVGivenDefault";

    /* default date is in the future */
    if(defaultDate > valueDate)
    {
        throw ModelException(method, "Default date is in the future " 
                             "(Default date is " + defaultDate.toString() +
                             ", Value date is " + valueDate.toString() + ")");
    }

    double feeLegPVGD = baseNotional * (payAccruedFee ? feeLeg->getAccruedInterest(defaultDate): 0.0);

    double ctgLegPVGD = baseNotional * (1 - recoveryRate);

    double df = crv.risklessPV(valuationDate, defaultDate);

    /* positive notional = long protection */
    return (df * ctgLegPVGD) - feeLegPVGD;
}

/** A VanillaCDS is NOT perpetual. */
bool VanillaCDS::hasFiniteMaturity() const
{
    return true;
}

DateTime VanillaCDS::getMaturity() const
{
    return maturityDate;
}

/**
 * Returns the earliest of the first accrual period, start of the contingent leg,
 * the first cash flow, etc.
 */
DateTime VanillaCDS::getStartDate() const
{
    return effectiveDate;
}

double VanillaCDS::getNotional() const
{
    return baseNotional;
}

void VanillaCDS::setNotional(double newNotional)
{
    const static string method = "VanillaCDS::setNotional";

    if(baseNotional == 0.0)
    {
        throw ModelException(method, "You cannot reset the notional of a "
                             "VanillaCDS which has zero notional!");
    }

    baseNotional = newNotional;
    feeLeg->setNotional(newNotional);
    contingentLeg->setNotional(newNotional);
}

/** Returns the day count convention used for accruals */
DayCountConventionSP VanillaCDS::getAccrualDcc() const
{
    return accrualDcc;
}

/* return only the fee leg cashflow */
CashFlowArraySP VanillaCDS::getInstrumentCashFlows(IForwardRatePricerSP model) const
{
    AbstractCashFlowArrayConstSP feeCfls = feeLeg->getCashFlows(model);
    //convert to simple cashflows
    return AbstractCashFlow::asCashFlows(feeCfls, model);
}    

/** Get the wrapper describing the discount curve */
YieldCurveWrapper VanillaCDS::getYieldCurveWrapper() const
{
    return YieldCurveWrapper(discount.getName());
}

/** Returns the wrapper describing the credit curve. */
ICDSParSpreadsWrapper VanillaCDS::getParSpreadsWrapper() const
{
    return ICDSParSpreadsWrapper(riskyCurve.getName());
}

double VanillaCDS::getPV(const DateTime&              valuationDate,
                         const IDiscountCurveRisky&   crv,
                         const IDecretionCurveConstSP prepay,
                         IForwardRatePricerSP         model,
                         IBadDayAdjusterConstSP       bda) const
{
    double pv = 0.0;

    /* positive notional = long protection */
    DateTime earliestRiskyDate = contingentLeg->firstObservationStartDate();
    DateTime latestRiskyDate = contingentLeg->lastObservationEndDate();
    pv = contingentLeg->getContingentLegPV(valuationDate, crv, bda) -
         feeLeg->getFeeLegPV(valuationDate, earliestRiskyDate, 
                             latestRiskyDate, crv, crv, prepay, 
                             true, getAccrualDcc(), model);

    return pv;
}

double VanillaCDS::getPV(const DateTime&              valuationDate,
                         const DateTime&              paymentDate,
                         const IDiscountCurveRisky&   crv,
                         const IDecretionCurveConstSP prepay,
                         IForwardRatePricerSP         model,
                         IBadDayAdjusterConstSP       bda) const
{
    double pv = 0.0;
 
    /* positive notional = long protection */
    DateTime earliestRiskyDate = contingentLeg->firstObservationStartDate();
    DateTime latestRiskyDate = contingentLeg->lastObservationEndDate();
    pv = contingentLeg->getContingentLegPV(valuationDate, paymentDate, crv, bda)-
         feeLeg->getFeeLegPV(valuationDate, paymentDate, earliestRiskyDate, latestRiskyDate, crv, crv, prepay, true, getAccrualDcc(), model);

    return pv;
}

//-----------------------
// ICDSConvention methods
//-----------------------

ICDSSP VanillaCDS::generateCDS(const DateTime&              startDate,
                               const DateTime&              endDate,
                               double                       feeRate) const
{
    const static string method = "VanillaCDS::generateCDS";
   
    /* checks */
    if (startDate>=endDate) throw ModelException(method,"startDate must be before endDate.");

    /* deep copy */
    VanillaCDSSP newCDS(dynamic_cast<VanillaCDS*>(this->clone()));

    if (!newCDS) 
    {
        throw ModelException(method, "Generating a VanillaCDS failed.");
    }

    /* reset the protection start and end date of the cds */
    newCDS->effectiveDate = startDate;
    newCDS->maturityDate  = endDate;

    /* reset the fee and contingent legs */
    VanillaCreditFeeLegSP newFeeLeg (dynamic_cast<VanillaCreditFeeLeg*>(feeLeg->clone()));
    newCDS->feeLeg = VanillaCreditFeeLegSP::dynamicCast(
        newFeeLeg->generateCreditFeeLeg(startDate, endDate, feeRate));

    VanillaCreditContingentLegSP newContingentLeg (dynamic_cast<VanillaCreditContingentLeg*>(contingentLeg->clone()));
    newCDS->contingentLeg = VanillaCreditContingentLegSP::dynamicCast(
        newContingentLeg->generateCreditContingentLeg(startDate, endDate));

    return newCDS;
}

//++++++++++++++++++++++++++++++++++++++++
//  IBadDayAdjuster methods
//
DateTime VanillaCDS::badDayAdjust(const DateTime& date) const {
    //no adjustment
    return date;
}

/** Add a number of business days to a date */
DateTime VanillaCDS::addBusinessDays(const DateTime& from, int busDays) const {
    //no adjustment
    return from.rollDate(busDays);
}
//---------------------
// LastSensDate methods
//---------------------

/**
  * When to stop tweaking.
  * Implementation is the same as that from CredDefSwap.
  * When re-routing of endDate through Model has changed to only be via
  * the instrument, this should be updated accordingly.
  */
DateTime VanillaCDS::endDate(const Sensitivity* sensControl) const
{
    const DateTime& lastSensDate = feeLeg->getCashFlowDates()->back();

    if(!sensControl)
    {
        return lastSensDate;
    }

    /** Mess brought over whereby endDate is passed through the Model 
     *  instead of via the Instrument.
     */
    IModel* model = sensControl->getModel();
    ClosedFormCDSPS* cdsModel = dynamic_cast<ClosedFormCDSPS*>(model);
    return cdsModel ? cdsModel->endDate(sensControl, this, lastSensDate) : 
                      lastSensDate;
}

//--------------------------------------
// ClosedFormCDSPS::IIntoProduct methods
//--------------------------------------

ClosedFormCDSPS::IProduct* VanillaCDS::createProduct(ClosedFormCDSPS* model) const
{
    return new VanillaCDSClosedFormCDSPS(this);
}

//---------------------
// Theta::Shift methods
//---------------------

bool VanillaCDS::sensShift(Theta* shift)
{
    valueDate = shift->rollDate(valueDate);
    return true;
}

//-----------------------------
// VanillaCDS methods (private)
//-----------------------------

void VanillaCDS::addOutputRequests(Control* control,
                                   Results* results,
                                   const IDecretionCurveConstSP prepay,
                                   IForwardRatePricerSP model) const
{
    static const string method = "VanillaCDS::addOutputRequests";
    double sp            = 1.0;
    double defaultValue  = 0.0;
    double parSpread     = 0.0;
    double riskyDuration = 0.0;

    try
    {
        OutputRequest* request;
        
        // IMPLIED_DEFAULT_PROBABILITY at maturity
        request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
        sp = riskyCurve->survivalProb(valueDate, maturityDate);
        if(request)
        {
            results->storeRequestResult(request, 1-sp);
        }
              
        // RECOVERY_VALUE
        request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
        if(request)
        {
             defaultValue = recoveryRate * baseNotional;
             results->storeRequestResult(request, defaultValue);
        }

        // CURRENT_SPREAD
        request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
        if(request)
        {
            results->storeRequestResult(request, feeRate);
        }

        request = control->requestsOutput(OutputRequest::IMPLIED_CDS_SPREAD);
        if(request)
        {
            DateTime earliestRiskyDate = valueDate;
            DateTime latestRiskyDate = feeLeg->getLastPayDate();
            parSpread = getContingentLegPV(valueDate, *(riskyCurve.getSP())) / 
                        feeLeg->getFeeLegPV(valueDate, earliestRiskyDate, latestRiskyDate, *(riskyCurve.getSP()), *(riskyCurve.getSP()), prepay, true, getAccrualDcc(), model);

            // quote in rate terms
            parSpread /= 100;

            // store the implied par spread
            results->storeRequestResult(request, parSpread);
        }

        // RISKY_DURATION - isn't picked up in addOutputRequest
        request = control->requestsOutput(OutputRequest::CDS_RISKY_DURATION);
        if(request)
        {
            DateTime earliestRiskyDate = valueDate;
            DateTime latestRiskyDate = feeLeg->getLastPayDate();
            riskyDuration = feeLeg->getFeeLegPV(valueDate, earliestRiskyDate, latestRiskyDate, *(riskyCurve.getSP()), *(riskyCurve.getSP()), prepay, true, getAccrualDcc(), model);

            /* unadjust notional, adjust in rate terms */
            riskyDuration /= baseNotional * .01;
            
            /* positive notional = long protection */
            /* adjust risky duration accordingly such that long */
            /* protection, pay fee, we are negative risky duration. */
            riskyDuration *= -1;

            results->storeRequestResult(request, riskyDuration);
        }

        // PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(paymentDates());
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }

        // KNOWN_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request) 
        {
            CashFlowArraySP cfl(getInstrumentCashFlows(model));
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get()); 

        }
        
    }
    catch(exception e)
    {
        throw ModelException(&e, method);
    }
}

/** when do payments occur ? */
DateTimeArraySP VanillaCDS::paymentDates() const {
    return feeLeg->getCashFlowDates();
}

/**Constructor*/
VanillaCDS::VanillaCDS(DateTime valueDate,
		       DateTime effectiveDate,     //protection start date 
		       DateTime maturityDate,      //maturity date i.e. protection end date 
		       bool koBeforeIssueDate, //knock out cds if the underlying credit defaults before issue date. 
		       YieldCurveWrapper discount,          //discount curve 
		       ICDSParSpreadsWrapper riskyCurve,        //risky curve 
		       bool payAccruedFee,     //pay accrued fee on default or not 
		       double feeRate,           //fee 
		       MaturityPeriodSP paymentFreq,       //payment frequency 
		       DayCountConventionSP pmtDcc,            //payment day count convention 
		       DayCountConventionSP accrualDcc,        //accrual day count convention 
		       BadDayConventionSP pmtBdc,            //payment bad day convention 
		       BadDayConventionSP accrualBdc,        //accrual bad day convention 
		       BadDayConventionSP valuationBdc,      //valuation bad day convention 
		       HolidayWrapper pmtHol,            //payment holiday(s) 
		       HolidayWrapper accrualHol,        //accrual holiday(s) 
		       string frontStubType,     //sf (short front stub), lf(long front stub) 
		       double recoveryRate,      //recovery rate 
		       double delay,             //delay between default date and settlement upon default 
		       SettlementSP settlementOffset)  //offset from a given trade date for settlement
  :CInstrument(TYPE),
   baseNotional(1.0),
   valueDate(valueDate),
   effectiveDate(effectiveDate),
   maturityDate(maturityDate),
   koBeforeIssueDate(koBeforeIssueDate),
   discount(discount),
   riskyCurve(riskyCurve),
   payAccruedFee(payAccruedFee),
   feeRate(feeRate),
   paymentFreq(paymentFreq),
   pmtDcc(pmtDcc),
   accrualDcc(accrualDcc),
   pmtBdc(pmtBdc),
   accrualBdc(accrualBdc),
   valuationBdc(valuationBdc),
   pmtHol(pmtHol),
   accrualHol(accrualHol),
   frontStubType(frontStubType),
   recoveryRate(recoveryRate),
   delay(delay),
   settlementOffset(settlementOffset)
{
  /* A VanillaCDS will always be 1-R */
    recoveryTypeString == "RECOVER_1_MINUS_R";
    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;

    /* Check inputs*/
    validatePop2Object();

    /* Defaults defined. Fee leg and contingent leg constructed*/
    Validate();
}

/**Default Constructor*/	     
VanillaCDS::VanillaCDS() : CInstrument(TYPE),
                           baseNotional(1.0),
                           payAccruedFee(true),
                           feeRate(1.0)
{
    /* A VanillaCDS will always be 1-R */
    recoveryTypeString == "RECOVER_1_MINUS_R";
    recovery = IDiscountCurveRisky::RECOVER_1_MINUS_R;
    delay = 0.0;
}

void VanillaCDS::load(CClassSP& clazz){
    clazz->setPublic();                // make visible to EAS/spreadsheet
    REGISTER(VanillaCDS, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(ICDS);
    IMPLEMENTS(ICDSConvention);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(name,               "name");
    FIELD(baseNotional,       "initial notional");
    FIELD(valueDate,          "valuation date");
    FIELD(effectiveDate,      "protection start date");
    FIELD(maturityDate,       "maturity date i.e. protection end date");
    FIELD(koBeforeIssueDate,  "knock out cds if the underlying credit defaults before issue date.");
    FIELD(discount,           "discount curve");
    FIELD(riskyCurve,         "risky curve");
    FIELD(payAccruedFee,      "pay accrued fee on default or not");
    FIELD(feeRate,            "fee");
    FIELD(paymentFreq,        "payment frequency");
    FIELD(fcd,                "first coupon date");
    FIELD(lcd,                "last coupon date");
    FIELD(pmtDcc,             "payment day count convention");
    FIELD(accrualDcc,         "accrual day count convention");
    FIELD(pmtBdc,             "payment bad day convention");
    FIELD(accrualBdc,         "accrual bad day convention");
    FIELD(valuationBdc,       "valuation bad day convention");
    FIELD(pmtHol,             "payment holiday(s)");
    FIELD(accrualHol,         "accrual holiday(s)");
    FIELD(frontStubType,      "sf (short front stub), lf(long front stub)");
    FIELD(backStubType,       "sb (short back stub), lb(long back stub)");
    FIELD(recoveryTypeString, "recovery type");
    FIELD(recoveryRate,       "recovery rate");
    FIELD(delay,              "delay between default date and settlement upon default");
    FIELD(settlementOffset,   "offset from a given trade date for settlement");
    FIELD(settlementType,     "type of settlement, physical or cash delivery");
    FIELD(feeLeg,             "fee leg");
    FIELD(contingentLeg,      "contingent leg");

    FIELD_MAKE_TRANSIENT(contingentLeg);
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD_MAKE_OPTIONAL(fcd);
    FIELD_MAKE_OPTIONAL(lcd);
    FIELD_MAKE_OPTIONAL(accrualDcc);
    FIELD_MAKE_OPTIONAL(accrualBdc);
    FIELD_MAKE_OPTIONAL(valuationBdc);
    FIELD_MAKE_OPTIONAL(accrualHol);
    FIELD_MAKE_OPTIONAL(backStubType);
    FIELD_MAKE_OPTIONAL(delay);

    FIELD_MAKE_TRANSIENT(recoveryTypeString);
    FIELD_MAKE_TRANSIENT(feeLeg);
}

IObject* VanillaCDS::defaultConstructor()
{
     return new VanillaCDS();
}

CClassConstSP const VanillaCDS::TYPE =
CClass::registerClassLoadMethod("VanillaCDS", typeid(VanillaCDS), load);

DEFINE_TEMPLATE_TYPE(VanillaCDSWrapper);

/* for class loading (avoid having header file) */
bool VanillaCDSLoad() {
    return (VanillaCDS::TYPE != 0);
}

//
///************** ICreditVanillaInstrument*************************/
//
//
//
///***  The model implementations and product definitions ***/
//
///* private class */
//
//
//
///******* End Implementation of ClosedFormCDSPS::IntoProduct interface ********/
//
//
///************************** Start Stuff for Reflection ********************/
//
//
//
//
//


///**************** Fee Leg ***************************************/
//
///**
// * Returns the fee rate on the IRiskyFeeLeg.
// */
//double VanillaCDS::getRate() const
//{
//    return (*feeLeg).getRate();
//}
//
///**
// * Change the fee rate on the fee leg.
// */
//void VanillaCDS::setRate(double newRate)
//{
//      feeRate = newRate;
//    (*feeLeg).setRate(newRate);
//}
//
//
///**
// * Compute PV of fee leg at valuation date, with instrument default settlement
// * and payment date behaviour.
// */
//double VanillaCDS::getFeeLegPV(const DateTime&              valuationDate,
//                               const IDiscountCurveRisky&   crv) const
//{
//    const static string method("VanillaCDS::getFeeLegPV");
//
//    return feeLeg->getPV(valuationDate, crv);
//}
//
///**
// * Compute PV of fee leg at valuationDate, but for unconditional payment at
// * paymentDate, conditional on no new defaults before valuationDate.
// */
//double VanillaCDS::getFeeLegPV(const DateTime&              valuationDate,
//                               const DateTime&              paymentDate,
//                               const IDiscountCurveRisky&   crv) const
//{
//    const static string method("VanillaCDS::getFeeLegPV");
//
//    return feeLeg->getPV(valuationDate, paymentDate, crv);
//}
//
//
//double VanillaCDS::getFeeLegPV(const DateTime&                 valuationDate,
//                               const DateTime&                 paymentDate,
//                               const IDiscountCurveRisky&      crv,
//                               bool                            defaultValueOnly) const
//{
//    const static string method("VanillaCDS::getFeeLegPV");
//
//    return feeLeg->getFeeLegPV(valuationDate, paymentDate, crv, defaultValueOnly);
//}
//
//

DRLIB_END_NAMESPACE
