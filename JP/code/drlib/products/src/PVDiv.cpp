//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PVDiv.cpp
//
//   Description   value the synthetic dividend (or PVDiv) instrument
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PVDiv.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CashSettlePeriod.hpp"
#include "edginc/CashSettleDate.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/AssetUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
void CPVDiv::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CPVDiv, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(MuParallel::IShift);
    IMPLEMENTS(MuPointwise::IShift);
    IMPLEMENTS(MuSpecial::IShift);
    EMPTY_SHELL_METHOD(defaultPVDiv);
    FIELD(valueDate,        "valuation Date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(isFixedNotional,  "Is contract schedule for fixed notional "
                 "(otherwise fixed num of shares) ?");
    FIELD(isCumulative,     "Is it a PVDiv of isCumulative style?");
    FIELD(contractSched,           "Contract schedule - sample date, "
          "size, hist spot");
    FIELD(equityAsset,      "equityAsset");
    FIELD_MAKE_OPTIONAL(equityAsset);
    FIELD(discount,      "Discount curve");
    FIELD(ccyTreatment,  "Currency Treatment");
    FIELD(instSettle,           "Instrument settlement at maturity");
    FIELD(averagingType, "Averaging Type: "
                 "either InvestOnContractDates or InvestOnWeightedAverage");
    FIELD_MAKE_OPTIONAL(averagingType);
    FIELD(averageSched, "Averaging schedule");
    FIELD_MAKE_OPTIONAL(averageSched);
    FIELD(useSynthetic, "if true, use synthetic dividend");
    FIELD_MAKE_OPTIONAL(useSynthetic);
    FIELD(synthDivs, "synthetic dividend input");
    FIELD_MAKE_OPTIONAL(synthDivs);
    FIELD(lastPayDate, "Cached last pay date");
    FIELD_MAKE_TRANSIENT(lastPayDate);
}

CClassConstSP const CPVDiv::TYPE = CClass::registerClassLoadMethod(
    "PVDiv", typeid(CPVDiv), load);

bool  CPVDivLoad() {
    return (CPVDiv::TYPE != 0);
   }


// constructor
CPVDiv::CPVDiv(): CInstrument(TYPE), averagingType("NA"), useSynthetic(false) {
    synthDivs = DividendListSP(   );
}

void CPVDiv::Validate()
{
    static const string method = "PVDiv::Validate";
    // validate against currency protected type
    if (CString::equalsIgnoreCase(ccyTreatment,
                                  CAsset::CCY_TREATMENT_PROTECTED))
    {
        throw ModelException(method, 
                             "PVDiv does not support ccy protected yet!");
    }

    DateTime contrStartDate, endDate;
    contractSched->getBoundingDates(contrStartDate,endDate);
    // general validation
    AssetUtil::assetCrossValidate(equityAsset.get(),
                                  contrStartDate.isGreater(valueDate),
                                  contrStartDate,
                                  valueDate,
                                  discount,
                                  this);


    //for fixedNotional, "InvestOnContractDates" pastValues must be defined and they must be positive  
    if (isFixedNotional && averagingType == "InvestOnContractDates")
    {  
        DateTime contractDate = valueDate;
        DateTime previousDate;
        double   size, pastSpot;
        while (contractSched->getPreviousSample(contractDate,previousDate,
                                                pastSpot,size))
        {
            if (!Maths::isPositive(pastSpot) )
            {
                throw ModelException(method, "Historical spot values "
                                     "must be positive");
            }
            contractDate = previousDate.rollDate(-1);
        }
    }

    // getAllDivsBetweenDates(divStartDate,endDate) and check them out
    DateTime divStartDate = contrStartDate;
    if (!isCumulative && valueDate >= contrStartDate) {
        divStartDate = valueDate;
    }

    // validate settlement
    if ((instSettle->isPhysical())||(instSettle->isMargin()))
    {// check that settle is not physical or margin
        throw ModelException(method, "PVDiv must be cash settled.\n");
    }
    else if (CashSettleDate::TYPE->isInstance(instSettle.get()) && 
             !isCumulative )
    {// if instSettle is CashSettleDate, must check that isCumulative 
        throw ModelException(method, "If PVDiv settled as "
                             "CashSettleDate then it must be cumulative");
    }
    else if (isCumulative)
    {
        if (CashSettlePeriod::TYPE->isInstance(instSettle.get()))
        {// if instSettle is CashSettlePeriod, it must not be cumulative
            throw ModelException(method, "PVDiv in cumulative mode "
                                 "must settle on a fixed date");
        }
        
    }
    if (isFixedNotional)
    {
        if (averagingType == "NA")
        {// if instSettle is CashSettlePeriod, it must not be cumulative
            throw ModelException(method, "fixed notional type requires "
                                 "averagingType string to be specified");
        }

        if (isCcyStruck())
        {// if instSettle is CashSettlePeriod, it must not be cumulative
            throw ModelException(method, "currency struck for fixed "
                                 "notional type is not allowed");
        }

        if (!(averagingType == "InvestOnContractDates" || 
              averagingType == "InvestOnWeightedAverage"))
        {// if instSettle is CashSettlePeriod, it must not be cumulative
            throw ModelException(method, "averagingType must be InvestOn"
                                 "ContractDates or InvestOnWeightedAverage");
        }
    }

    
    // validation for InvestOnWeightedAverage
    if (isFixedNotional && averagingType == "InvestOnWeightedAverage")
    {
        if (!averageSched || averageSched->getDates().size()==0)
        {
            throw ModelException(method, "no averaging schedule for "
                                 "type InvestOnWeightedAverage");
        }

        if (!isCumulative)
        {   // check that average ends before contract schedule begins
            // now only require this validation in non cum mode 
            if (averageSched->getLastDate() > contractSched->getFirstDate())
            {
                throw ModelException(
                    method, "For non cumulative mode, "
                    "the last averaging date (" +
                    averageSched->getLastDate().toString() + 
                    ") cannot be after the first contract date (" + 
                    contractSched->getFirstDate().toString() + ")");
            }
        }
    }
    
    if (isCumulative)
    {
        DateTime theCumPayDate = getCumPayDate();
        const DateTime& theLastContractDate = contractSched->getLastDate();
        if (theCumPayDate < theLastContractDate)
        {
            throw ModelException(
                method, "Cumulative mode: The cum pay date ("+
                theCumPayDate.toString() + 
                ") cannot lie before the last contract date (" +
                theLastContractDate.toString() + ").");
        } 
    }
    // validate synthetic divs
    if (useSynthetic && (!synthDivs))
    {
        throw ModelException(method, "If use synthetic dividend flag = true, "
                               "then synthetic dividends must be supplied.");
    } 
}

// initiate GetMarket 
void CPVDiv::GetMarket(const IModel* model, const CMarketDataSP market)
{
    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                               discount, equityAsset);
    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
   
    // since the equityAssets have changed, we have to validate the instrument
    validatePop2Object();

    // equityAsset and discount curve must not be null 
    if (!equityAsset)
    {
        throw ModelException("CPVDiv::GetMarket","Asset is NULL.");
    }

    if (!discount)
    {
        throw ModelException("CPVDiv::GetMarket","Discount curve is NULL.");
    }
}

/** private class */
class CPVDivClosedFormProd: public CClosedFormLN::IProduct{
private:
    const CPVDiv*  instrPVDiv; // a reference

    //private methods
    double getCashAmount(Dividend& div, double fwdValue) const;
    
    /** compute expected value of 1/(sum of w_i * S_i)) */
    double  calcExpectedReciprocalAve() const;
public:
    CPVDivClosedFormProd(const CPVDiv* instr): instrPVDiv(instr){}

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results) const;
};

/** compute expected value of 1/(sum of w_i * S_i)) */
double CPVDivClosedFormProd::calcExpectedReciprocalAve() const
{
    try {
        double aveFwd = instrPVDiv->
            averageSched->expectedAverage(instrPVDiv->equityAsset.get(),
                                          instrPVDiv->valueDate);

        // choose atm vol 
        ATMVolRequestSP volRequest(new ATMVolRequest());
        
        // interpolate the vol using our LN request
        CVolProcessedBSSP procVol(instrPVDiv->equityAsset->
                                  getProcessedVol(volRequest.get()));
        // calculate the variance
        double aveVar = instrPVDiv->averageSched->
            averageVariance(procVol.get(), instrPVDiv->valueDate, true);
        
        double result = exp(aveVar)/aveFwd;

        return result;
    }
    catch (exception& e) {
        throw ModelException(e, "CPVDivClosedFormProd::"
                             "calcExpectedReciprocalAve");
    }
}

/** class for adjusting dividend pay dates - could put method on CPVDiv
    but tidier to make it internal */
class CPVDiv::DivAdjuster: public virtual DividendList::IDivAdjuster{
public:
    // fields
    double         sizeAdjust;
    DateTime       cumPayDate;
    const CPVDiv*  pvDiv;
public:
    // simple constructor
    DivAdjuster(double          sizeAdjust,
                const DateTime& cumPayDate,
                const CPVDiv*   pvDiv):
        sizeAdjust(sizeAdjust), cumPayDate(cumPayDate), 
        pvDiv(pvDiv) {}

    /** Returns the new pay date for the specified dividend */
    virtual void adjustDividend(Dividend& div){
        pvDiv->adjustDividend(sizeAdjust, cumPayDate, div);
    }
    virtual ~DivAdjuster(){}
};

/** Adjusts the pay date for the specified dividend and also adjusts its
    size. Note dividend is in its underlying currency. Any yield divs are
    converted after this method is called. Additionally any struck divs are
    scaled using the div's pay date */
void CPVDiv::adjustDividend(double          sizeAdjust,
                            const DateTime& cumPayDate,
                            Dividend&       div) const {
    static const string method("CPVDiv::adjustDividend");
    try{
        const DateTime& exDate = div.getExDate();
        DateTime divStartDate(exDate.getDate()-1, DateTime::END_OF_DAY_TIME);
        // we must have valid sizes from contract schedule
        DateTime      contractDate;
        double        pastSpot;
        double        size = 1.0;
        if (!contractSched->getPreviousSample(divStartDate, 
                                              contractDate, pastSpot, size)){
            throw ModelException(method, "unable to determine contract "
                                 "size from schedule");
        }
        //different treatment for fixed nominal and fixed notional list
        double cashAmount = div.getDivAmount(); //either $ or %
        if (!isFixedNotional) 
        {   /** for nominal fixing, there is no size adjustment */
            sizeAdjust = 1.0;
        }
        else
        {   /** for notional fixing compute sizeAdjust to get
                equivalent nominal */
            if (valueDate >= exDate) 
            {
                if (averagingType == "InvestOnContractDates"){
                    sizeAdjust = 1.0/pastSpot;
                }
            }
            else
            {   // to do: we can allow yield dividend here - require 
                // calling average ratio pricer
                if ((averagingType == "InvestOnWeightedAverage") &&
                    (div.getDivType() == Dividend::PERCENT) )
                {
                    // don't need a convexity adjustment for yield
                    // dividends in type 2.  but type 2 will
                    // adjust all divs automatically.  see below
                    cashAmount /= sizeAdjust;
                    const DateTime& firstSchedDate =
                        averageSched->getFirstDate();
                    if (valueDate >= firstSchedDate){
                        cashAmount /= averageSched->averageToDate(valueDate);
                    } else {
                        cashAmount /= equityAsset->fwdValue(firstSchedDate);
                    }
                }

                // ordinary fixedNotional case 
                if (averagingType == "InvestOnContractDates")
                {
                    if (contractDate > valueDate)
                    {
                        if (div.getDivType() == Dividend::PERCENT)
                        {
                            sizeAdjust = 
                                1.0/equityAsset->fwdValue(contractDate);
                        }
                        else
                        {   //div is dollar type
                            sizeAdjust = equityAsset->expNumberShares(
                                valueDate, contractDate, true);
                        }
                    }
                    else
                    {
                        sizeAdjust = 1.0/pastSpot;
                    }
                }
            } 
        }  
        double newDivAmount = size * sizeAdjust * cashAmount;
        div.setDivAmount(newDivAmount);

        // then sort out div pay date
        DateTime payDate;
        if (isCumulative){
            payDate = cumPayDate;
        } else {
            payDate = instSettle->settles(div.getPayDate(), equityAsset.get());
        }
        div.setPayDate(payDate);
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

void CPVDivClosedFormProd::price(CClosedFormLN*  model,
                                 Control*        control, 
                                 CResults*       results) const
{
    static const string method = "CPVDivClosedForm::price";
    try {
        //use local objects instead of going through instrument each time
        DateTime                 valueDate        = instrPVDiv->valueDate;    
        bool                     isCumulative     = instrPVDiv->isCumulative;
        CAssetWrapper            equityAsset      = instrPVDiv->equityAsset;
        YieldCurveWrapper        discount         = instrPVDiv->discount;
    
        if (control && control->isPricing()){
            instrPVDiv->lastPayDate = valueDate; // in case no divs etc
        }
        /** Find dividend startDate:  */        
        DateTime divStartDate, contrStartDate, endDate;
        instrPVDiv->contractSched->getBoundingDates(contrStartDate, endDate);
        divStartDate = contrStartDate;

        // just for getting cum pay date for now, to be done in better way
        // to allow theta to work if cumpay date=value date
        DateTime cumPayDate(instrPVDiv->getCumPayDate());        
                
        // for InvestOnWeightedAverage, pre-compute E[1/sum(w_i * S_i)]
        double  sizeAdjust=1.0;   
        if (instrPVDiv->isFixedNotional && 
            instrPVDiv->averagingType == "InvestOnWeightedAverage")
        {
            sizeAdjust = calcExpectedReciprocalAve();
        }
        
        /** Get hold of the dividends - use DividendCollector to manage this
            process. Will call adjustDividend() for each div */
        // First set up our call back to override dividends
        OutputRequest*   cfRequest = 0;
        if (control && control->isPricing()){
            // find out whether to keep track of known cash flows
            cfRequest = 
                control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        }
        CPVDiv::DivAdjuster divAdjuster(sizeAdjust, cumPayDate, instrPVDiv);
        // Then initialise our DividendCollector
        OutputRequestUtil::KnownCashFlows knownCFs; // for KNOWN_CASHFLOWS
        DividendCollector collector(&divAdjuster, valueDate,
                                    divStartDate, endDate,
                                    cfRequest? &knownCFs: 0);
          // process the divs
        if (instrPVDiv->useSynthetic) {
            collector.go(equityAsset.get(), instrPVDiv->synthDivs);
        } else {
            collector.go(equityAsset.get());
        }

        // some validation which is easiest to defer until pricing
        if (instrPVDiv->isFixedNotional && 
            instrPVDiv->averagingType == "InvestOnWeightedAverage" &&
            collector.hasYieldDividends() && 
            instrPVDiv->averageSched->getDates().size() > 1){
            throw ModelException(method, "yield dividend not supported for "
                                 "InvestOnWeightedAverage if there is more "
                                 "than one averaging date");
        }
        /* get hold of divs (these come back in discount ccy, with yield divs
           having been converted) */
        DividendListSP divsInUse = collector.getDividends();
        const DividendArray& divArray = divsInUse->getArray();
        // just need to pv them
        DateTimeArray payDates; // for payment dates
        double  pvDividend = 0.0;
        for (int divIdx = 0; divIdx < divArray.size(); divIdx++){
            /*** NOTE: this block of code must only do pv'ing. Any other 
                 scaling done here will mess up the KNOWN_CASHFLOWS code ***/
            double amount = divArray[divIdx].getDivAmount();
            const DateTime& payDate = divArray[divIdx].getPayDate();
            if (payDate <= valueDate){
                // historic
                if (!isCumulative){
                    payDates.push_back(payDate);
                }
            } else if (isCumulative){
                pvDividend += amount;
            } else {
                // pay date set to settlement date in divAdjuster
                pvDividend += amount * discount->pv(payDate);
                payDates.push_back(payDate);
            }
            if (control && control->isPricing() && 
                (divIdx == 0 || payDate.isGreater(instrPVDiv->lastPayDate))){
                // cache lastPayDate (this is invariant under tweaks)
                instrPVDiv->lastPayDate = payDate;
            }
        }
        if (isCumulative) {
            payDates.push_back(cumPayDate);
            if (cumPayDate > valueDate){
                pvDividend *= discount->pv(valueDate, cumPayDate);
            }
        }

        //Store result
        results->storePrice(pvDividend, discount->getCcy());

        if (control && control->isPricing()) {
           
            // do payment dates/known cashflows
            OutputRequest* request =
                control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request){
                DateTime::removeDuplicates(payDates, false);
                OutputRequestUtil::recordPaymentDates(control,results, 
                                                      &payDates);
            }
            if (cfRequest){
                knownCFs.recordKnownCashFlows(control, results);
            }
        }
    }   
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

DateTime CPVDiv::getValueDate() const
{
    return valueDate;
}

// when to stop tweaking
DateTime CPVDiv::endDate(const Sensitivity* sensControl) const
{
    if (lastPayDate.empty()) 
    {  // lastPayDate is set in price, but credit calls endDate before price, so
       // need to compute on the fly here.
        DateTime lastPayDate2(0, 0); 
        try
        {
            double  sizeAdjust=1.0;   
            /** Find dividend startDate:  */        
            DateTime divStartDate, contrStartDate, endDate;
            contractSched->getBoundingDates(contrStartDate, endDate);
            divStartDate = contrStartDate;
            
            // default to start date in case no divs below
            lastPayDate2 = divStartDate;

            // just for getting cum pay date for now, to be done in better way
            // to allow theta to work if cumpay date=value date
            DateTime cumPayDate(getCumPayDate());        
            
            CPVDiv::DivAdjuster divAdjuster(sizeAdjust, cumPayDate, this);
            // Then initialise our DividendCollector
            DividendCollector collector(&divAdjuster, valueDate,
                divStartDate, endDate, 0);

            if (useSynthetic) {
                collector.go(equityAsset.get(), synthDivs);  // process the divs
            } else {
                collector.go(equityAsset.get());  // process the divs
            }
            
            DividendListSP divsInUse = collector.getDividends();
            const DividendArray& divArray = divsInUse->getArray();
            
            for (int divIdx = 0; divIdx < divArray.size(); divIdx++)
            {
                const DateTime& payDate = divArray[divIdx].getPayDate();
                if (divIdx == 0 || payDate.isGreater(lastPayDate2))
                {
                    lastPayDate2 = payDate;
                }
            }
        }
        catch (exception& e) 
        {
            throw ModelException(e, "CPVDiv::endDate",
                "Internal eror - couldn't compute last pay date");
        }
        
        return lastPayDate2;
    }
    else
    {
        return lastPayDate;
    }
}


/** Rolls the value date and sets sample if needed  */
bool CPVDiv::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    // this is for theta at spot
    CAssetConstSP ptr(equityAsset.getSP());
    contractSched->roll(ptr.get(), valueDate, newDate, !shift->useAssetFwds());

    if (!!averageSched){
		// if there isn't an averageSched, but should be one, validate will catch it
        averageSched->roll(ptr.get(), valueDate, newDate, 
                           !shift->useAssetFwds());
    }
    // validate against yield div between value date and roll date in synthetic divs
    if (useSynthetic && synthDivs.get()) {
        DividendListSP div(synthDivs->getAllDivsBetweenDates(valueDate, newDate));
        if (div->hasYieldDividend()) {
            throw ModelException("CPVDiv::sensShift",
                                 "cannot roll over yield dividend for synthetic div when calculating theta.");
        }
    }
    // roll today 
    valueDate = newDate;
    return true;
}

//** Methods for synthetic dividend tweaking: */

/** Returns name identifying vol for MU_PARALLEL */
string CPVDiv::sensName(MuParallel* shift) const{
    return equityAsset->getTrueName();
}

/** Returns name identifying vol for MU_POINTWISE */
string CPVDiv::sensName(MuPointwise* shift) const{
    return equityAsset->getTrueName();
}

/** Returns name identifying vol for MU_SPECIAL */
string CPVDiv::sensName(MuSpecial* shift) const{
    return equityAsset->getTrueName();
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_PARALLEL shift method */
bool CPVDiv::sensShift(MuParallel* shift)
{
    if (useSynthetic) {
        synthDivs->sensShift(shift, valueDate);
    }

    return true;
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_POINTWISE shift method */
bool CPVDiv::sensShift(MuPointwise* shift)
{
    if (useSynthetic) {
        bool result = false;
        ExpiryArrayConstSP expiries = shift->getExpiries();

        DateTime bucketStartDate, bucketEndDate;
        shift->getBucketDates(bucketStartDate, bucketEndDate);

        result = synthDivs->sensShift(shift, valueDate,
                                              bucketStartDate, bucketEndDate);
    }
    return true;
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_S shift method */
bool CPVDiv::sensShift(MuSpecial* shift)
{
    if (useSynthetic) {
        DateTime expiry;

        // get the current expiry date being tweaked and ensure that it
        // is not before today
        expiry = shift->getExpiry()->toDate(expiry);

        if (!expiry.isLess(valueDate)) {
            // note if it's struck we want to make sure we shift wrt underlying spot price (no FX)
            double spot = equityAsset->getSpot();
            if (isCcyStruck()) {
                Asset::IStruck* struckAsset =  &dynamic_cast<Asset::IStruck&>(*equityAsset.get());             
                spot /= struckAsset->getFXSpot();
            }
            synthDivs->sensShift(shift, spot, expiry);
        }
    }
    return true;
}

bool CPVDiv::isCcyStruck() const
{
    return (CString::equalsIgnoreCase(ccyTreatment,
                                      CAsset::CCY_TREATMENT_STRUCK) ||
            StruckEquity::TYPE->isInstance(equityAsset.get()));
}

// returns true if PVDiv is fixed notional, invest on wgted avg with 1 avg date, which is
// the edg definition of fwd starting.
bool CPVDiv::isLegacyFwdStarting() const
{
    return isFixedNotional && averagingType == "InvestOnWeightedAverage" && 
				!!averageSched && averageSched->getDates().size() == 1;

}

DateTime CPVDiv::getCumPayDate() const
{
    DateTime      dummy(0,0);
    return instSettle->settles(dummy, equityAsset.get());
}

IObject* CPVDiv::defaultPVDiv(){
    return new CPVDiv();
}

CreditSupportSP CPVDiv::createCreditSupport(CMarketDataSP market){
    return CreditSupportSP(new PVDivCreditSupport(this, market));
}


/** Returns the name of the instrument's discount currency. */
string CPVDiv::discountYieldCurveName() const {
    return discount.getName();
}


/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CPVDiv::createProduct(
    CClosedFormLN* model) const{
    return new CPVDivClosedFormProd(this);
}

DRLIB_END_NAMESPACE
