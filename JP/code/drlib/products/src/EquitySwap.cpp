//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquitySwap.cpp
//
//   Description   equity swap model
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EquitySwap.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Business252.hpp"

#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/EndDateCollector.hpp"

DRLIB_BEGIN_NAMESPACE

// helpers
void CEquitySwap::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(CEquitySwap, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(MuParallel::IShift);
    IMPLEMENTS(MuPointwise::IShift);
    IMPLEMENTS(MuSpecial::IShift);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(LastSensDate);
    EMPTY_SHELL_METHOD(defaultEquitySwap);
    FIELD(valueDate,        "valuation Date");
    FIELD_MAKE_OPTIONAL(valueDate);
    FIELD(excludeCashFlowsToday, "no longer used");
    FIELD_MAKE_OPTIONAL(excludeCashFlowsToday);
    FIELD(isCallable, "true, swap will be termininated on call date");
    FIELD(callDate, "swap call date");
    FIELD_MAKE_OPTIONAL(callDate);
    FIELD(asset, "swap underlying asset");
    FIELD(discount, "Discount ccy");
    FIELD(isFixedNotional, "if it is fixed notional type, otherwise fixed num of shares");
    FIELD(eqLeg, "equity leg handle");
    FIELD(liborLeg, "libor leg handle");
    FIELD(divLeg, "dividend leg handle");
    FIELD_MAKE_OPTIONAL(divLeg);
    FIELD(cashflowLeg, "cashflow leg handle");
    FIELD_MAKE_OPTIONAL(cashflowLeg);
    FIELD(validateSamples, "whether to validate historical samples");
    FIELD_MAKE_TRANSIENT(validateSamples);

    FIELD(swapType, " STANDARD: standard swap , BRZ: brazil eq swap, DCC : Bus/252");
    FIELD_MAKE_OPTIONAL(swapType);
    
}

CClassConstSP const CEquitySwap::TYPE = CClass::registerClassLoadMethod(
    "EquitySwap", typeid(CEquitySwap), load);

bool  CEquitySwapLoad() {
    return (CEquitySwap::TYPE != 0);
}



void CEquitySwap::validatePop2Object()
{
    // validate against call date in average out schedule here.  We push up call date to tomorrow
    // in ::Init, and at that point we check this condition again, but don't fail if the call date
    // lands in avg out, just turn off callable
    if (isCallable && !!eqLeg->avOut && callDate > eqLeg->avOut->sampleDates[0])
    {
        throw ModelException("CEquitySwap::validatePop2Object",
                             "Call date must lie on or before first averaging out date.");
    }

    // sets liborLeg->ESWAveType to AVE_IN, AVE_OUT, or HYBRID.
    setLibAvType();

    liborLeg->needSpread = (liborLeg->isFloating
                            || (!!divLeg && divLeg->accrueDivType == 2)
                            || eqLeg->accrueEqType == 2
                            || (!!eqLeg->avOut && eqLeg->avOut->accrueAvOutType == 2));

    // validation below:

    // validate agains unsupported, unsigned off functionality
//    if (!!divLeg && divLeg->accrueDivType != 0 ||
//        eqLeg->accrueEqType != 0 ||
//        !!eqLeg->avOut && eqLeg->avOut->accrueAvOutType != 0)
//    {
//        throw ModelException("Equity Swap","Reinvestment is currently not supported.");
//    }

}

void CEquitySwap::Validate()
{
    static const string method = "EquitySwap::Validate";

//      if (isCallable && (!!eqLeg->avIn || !!eqLeg->avOut))
//    {
//        throw ModelException("CEquitySwap::Validate",
//                             "Support for Callable feature in averaging swaps is pending DR signoff.");
//    }

    validateHistoricalSamples();

    // validate Equity leg
    eqLeg->Validate(liborLeg->accrueDates, liborLeg->useExplicitLinking, liborLeg->isFloating,
                    callDate, isCallable, valueDate);

    // validate Libor leg
    liborLeg->Validate(eqLeg,valueDate, eqLeg->eqRefixDates[0], callDate);// are these dates correct ?

    if(!!cashflowLeg)
    {
        // validate CashFlow leg
        cashflowLeg->Validate();
    }

    if(!!divLeg)
    {
        // validate Dividend leg
        divLeg->Validate();
    }

    // validate fx and ccy's.

    IObjectSP imnt(IObjectSP::attachToRef(this));
    CFXRateCollector::validateAllFXRates(imnt);

    AssetCcyCollector::validateAllCurrencies(imnt, discount.get());

    // no longer validate against callDate before start of swap
    // EAS request, due to increases/unwinds
    // get swap starting halfway through for example
    // should have same call date as orig swap
    // We set now call date to tomorrow in this case as well, in EquitySwap::Init

//    if (isCallable && callDate <= eqLeg->eqRefixDates[0])
//    {
//        throw ModelException("CEquitySwap::Validate",
//                             "call date cannot be before first equity start date.");
//    }


//    if (
//        !!eqLeg->avIn &&
//        !!eqLeg->avOut &&
//        !eqLeg->avIn->isConsolidated)
//    {
//        throw ModelException("CEquitySwap::Validate",
//                             "Hybrid swaps: averaging in must be consolidated.");
//    }
//

    if (swapType == "BRZ" ){ //Brazil Swap, :have to been Bus/252 as DCC
        if (!dynamic_cast<const Business252*>(liborLeg->getpayDCC().get())){
            throw ModelException("CEquitySwap::Init",
                                 "Brazil swaps should have Bus/252 dcc>");
        }

        //should check yield curve's DCC as well???  To review


        if (isCallable || (!!eqLeg->avIn || !!eqLeg->avOut) || !!divLeg){
            throw ModelException("CEquitySwap::Init",
                                 "Callable, averaging or div Leg are not available for Brazil equity swap.");
        }
    }

    //DCC Bus/252 is only available for Brazil swap 
    if(dynamic_cast<const Business252*>(liborLeg->getpayDCC().get()) && swapType != "BRZ"){
        throw ModelException("CEquitySwap::Init",
                             "Bus/252 dcc should only be specified for Brazil swaps");
    }
}

void CEquitySwap::validateHistoricalSamples()
{
    static const string method = "CEquitySwap::validateHistoricalSamples";
    try {

        if (validateSamples) {
            // Validate no. of past fixings and set future samples to zero

            eqLeg->validateHistoricalSamples(valueDate);
            liborLeg->validateHistoricalSamples(valueDate);

            if (!!eqLeg->avIn)
            {
                eqLeg->avIn->ccyTreatment = eqLeg->ccyTreatment;
                eqLeg->avIn->validateHistoricalSamples(valueDate, liborLeg->isFloating);
                CEquitySwap::zeroFutureSamples(eqLeg->avIn->sampleDates, eqLeg->avIn->sampleLevels, valueDate, 0);
                if (eqLeg->ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
                {
                    // to do: remove repeat resize in init, if any
                    eqLeg->avIn->fxFixing.resize(eqLeg->avIn->sampleDates.size());
                    CEquitySwap::zeroFutureSamples(eqLeg->avIn->sampleDates, eqLeg->avIn->fxFixing, valueDate, 0);
                }
            }
    
            if (!!eqLeg->avOut)
            {
                eqLeg->avOut->ccyTreatment = eqLeg->ccyTreatment;
                eqLeg->avOut->validateHistoricalSamples(valueDate);
                CEquitySwap::zeroFutureSamples(eqLeg->avOut->sampleDates, eqLeg->avOut->sampleLevels, valueDate, 0);
                if (eqLeg->ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
                {
                    // to do: remove repeat resize in init, if any
                    eqLeg->avOut->fxFixing.resize(eqLeg->avOut->sampleDates.size());
                    CEquitySwap::zeroFutureSamples(eqLeg->avOut->sampleDates, eqLeg->avOut->fxFixing, valueDate, 0);
                }
            }

            CEquitySwap::zeroFutureSamples(eqLeg->eqRefixDates,eqLeg->eqStartLevel, valueDate, 0);
            CEquitySwap::zeroFutureSamples(eqLeg->eqRefixDates,eqLeg->eqEndLevel, valueDate, 1);
            if (eqLeg->ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
            {
                // to do: remove repeat resize in init, if any
                eqLeg->fxFixing.resize(eqLeg->eqRefixDates.size());
                CEquitySwap::zeroFutureSamples(eqLeg->eqRefixDates,eqLeg->fxFixing, valueDate, 0);
            }
            validateSamples = false;
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}



// initiate GetMarket
void CEquitySwap::GetMarket(const IModel* model, const CMarketDataSP market)
{
    market->GetReferenceDate(valueDate);
    CAsset::getAssetMarketData(model, market.get(), eqLeg->ccyTreatment, discount, asset);
    discount.getData(model, market);

    if (!!cashflowLeg)
    {
        cashflowLeg->discount.getData(model, market); // gets its own curve
    }

    liborLeg->GetMarket(model, market);

    // since the assets have changed, we have to validate the instrument
    //validatePop2Object();

    // asset and discount curve must not be null
    if (!asset)
    {
        throw ModelException("CEquitySwap::GetMarket",
                             "Asset is NULL.");
    }

    if (!discount)
    {
        throw ModelException("CEquitySwap::GetMarket",
                             "Discount curve is NULL.");
    }

}


/** returns current value date */
DateTime CEquitySwap::getValueDate() const
{
    return valueDate;
}

// constructor
CEquitySwap::CEquitySwap(): CInstrument(TYPE), validateSamples(true) 
{
    swapType = "STANDARD";  //STANDARD: DEFAULT, ie: no brazil swap
}

// initialisation
// This method calls init() in each leg.  Callable swaps are curtailed in each
// leg's init() method.  The rules for curtailing legs are provided here:
// (courtesy of Mark Robson)
//
// * equity final valuation date = call date
// * fx final valuation date = call date
// * equity payment date = equity settlement date for the call date
// * IR last accrual date = equity settlement date for the call date
// * IR last reset date if the reset is in arrears = call date
// * IR payment date = equity settlement date for the call date
// * Dividend final cumulative date  = call date
// * Dividend payment date = equity settlement date for the call date
//
// subject to the following extra rules
// 1.  Any other equity payment dates after the settlement date (as defined above)
//         are set to the settlement date
// 2.  If the accrue end date of the final IR period is before/on the settlement date,
//         no action is taken (except under rules 3 and 4). Otherwise, all IR periods whose
//         accrue end date is before/on the call date are left unaffected (subject to rule 3). Note that if an accrue end date
//         of a period occurs between the call date and its settlement date,
//         This results in the period being extended rather than having a stub for
//     the time between the period accrual end date and the settlement date.
// 3.  Any other IR payment dates after the settlement date (as defined above) are set to the settlement date
// 4.  All IR refix dates after the call date are set to the call date.

void CEquitySwap::Init() {


    // if call date in past, or, for forward starting swaps, if
    // call date is before start of swap, set call date to next
    // business day after max (today, start of swap).

    if (isCallable && (valueDate >= callDate || eqLeg->eqRefixDates[0] >= callDate))
    {
        HolidayConstSP hol = AssetUtil::getHoliday(asset.get());
        // ideally the above inequalities should both be >=
        // would make the code cleaner!
        // (should have DateTime FirstAllowableCallableDate = max(valueDate, startDate)
        // then wouldn't need to treat condition to reset and how to reset differently

        if (valueDate >= eqLeg->eqRefixDates[0])
        {   // don't set to tomorrow for forward starting swaps!
            callDate = hol->addBusinessDays(valueDate,1);
        }
        else
        {   // it's a forward starting swap, so call date is set to one bus day after swap begin
            callDate = hol->addBusinessDays(eqLeg->eqRefixDates[0],1);
        }
    }


    // since we moved up the call date above, validate an early call date here.
    if (isCallable && !!eqLeg->avIn && callDate < eqLeg->avIn->sampleDates[eqLeg->avIn->sampleDates.size() - 1])
    {
        throw ModelException("CEquitySwap::Init",
                             "Call date must lie on or after last averaging in date.");
    }

    // if the call date was pushed up into the averaging out period
    if (isCallable && !!eqLeg->avOut && callDate > eqLeg->avOut->sampleDates[0])
    {
        // turn it off
        isCallable = false;
    }

    // if callable and average out, simply delete the averaging out schedule,
    // as we have validated that the call date will lie before the avg out schedule.
    if (isCallable && !!eqLeg->avOut)
    {
        eqLeg->avOut = ESWAvOutSP(   );
        // now since we no longer have av out, need to reset libor leg's av out type:
        setLibAvType();
    }

    //get payDCC from LiborLeg, and set holidays. 
    DayCountConventionSP& payDCC =liborLeg->getpayDCC(); 
    if(dynamic_cast<const Business252*>(payDCC.get())){
        HolidayConstSP  hol = liborLeg->getHolidays();
        dynamic_cast< Business252*>(payDCC.get())->setHoliday(hol);
    }

    //set holidays to rateDCC in LiborLeg 
    DayCountConventionSP& rateDCC =liborLeg->getrateDCC(); 
    if(dynamic_cast<const Business252*>(rateDCC.get())){
        HolidayConstSP hol = liborLeg->getHolidays();
        dynamic_cast< Business252*>(rateDCC.get())->setHoliday(hol);
    }

    // temp until full support of Brazil YieldCurve
//    if (swapType == "B" ){ //Brazil Swap, :have to been Bus/252 as DC
//        discount = ESWBrazil::replaceCurve(discount.getSP());
//    }

    eqLeg->init(asset, discount, valueDate, callDate, isCallable, isFixedNotional, 
            liborLeg->getpayDCC(),  /*get the payDCC from liborLeg*/ swapType);
    
    liborLeg->init(eqLeg, valueDate, callDate, eqLeg->EqCallSettleDate, isCallable,
                   YieldCurveSP::dynamicCast((IObjectSP)discount.getSP()), swapType);

    liborLeg->eqAccrualFactors(eqLeg, valueDate, discount.getSP());

    if (!!(divLeg))
    {
        divLeg->init(eqLeg, liborLeg, discount.getSP(), callDate,
                     isCallable, valueDate);
    }
}

/** private class */
class CEquitySwapClosedFormProd: public CClosedFormLN::IProduct{
private:
    // bad but true, ESW changes instrument definition !!! - callDate,
    // instrumenet deep copied
    CEquitySwapSP instrESW;
    //const CEquitySwap*  instrESW; // a reference

public:
    CEquitySwapClosedFormProd(const CEquitySwap* instr):
        instrESW(CEquitySwapSP(dynamic_cast<CEquitySwap*>(instr->clone()))){}

    void price(CClosedFormLN*   model,
               Control*        control,
               CResults*       results) const;

    void writeLegRequestedOutput(
        Control*            control,
        CResults*           results,
        double              legPrice,
        const string&       ccyISOCode,
        OutputNameConstSP&  output) const;
};

// pricing formula
void CEquitySwapClosedFormProd::price(CClosedFormLN*   model,
                                      Control*        control,
                                      CResults*       results) const
{
    static const string method = "CCEquitySwapClosedForm::price";
    try {
        double      value;         // the fair value
        double          eqPrice = 0.0, divPrice=0.0, cfPrice=0.0, irPrice=0.0;

        // need name for outputting
        string      ccyName = instrESW->discount->getCcy();

        // should not be needed as it's called by model - but need to move
        // some codes from each legs validate() into init() as Validate() should only be called once
        instrESW->Validate();

        // initialize knownCFs.  Each leg's pricing call will update.
        // Note: do not currently separate these as we already are returning all cf's known and unknown.
        OutputRequestUtil::KnownCashFlows knownCFs; // for KNOWN_CASHFLOWS
        DateTimeArray payment_dates; // for PAYMENT_DATES

        // Init() is not a const function - it changes various attributes
        // This method isn't const either - have to cheat for now
        CEquitySwap* nonConstESW = const_cast<CEquitySwap*>(instrESW.get());
        nonConstESW->Init();

        eqPrice = instrESW->eqLeg->priceEq(instrESW->discount, control, results, knownCFs, payment_dates);
        OutputNameConstSP eqOutput(new OutputName(instrESW->eqLeg->LegName));
        //results->storeGreek(CDoubleSP(CDouble::create(eqPrice)),"PRICE",eqOutput);

        writeLegRequestedOutput(control, results, eqPrice, ccyName, eqOutput);

        irPrice = instrESW->liborLeg->priceIR(control, results, knownCFs, payment_dates);
        OutputNameConstSP irOutput(new OutputName(instrESW->liborLeg->LegName));
        //results->storeGreek(CDoubleSP(CDouble::create(irPrice)),"PRICE",irOutput);

        writeLegRequestedOutput(control, results, irPrice, ccyName, irOutput);

        // cashflow and div legs may not be present
        if (!!(instrESW->cashflowLeg))
        {
            cfPrice = instrESW->cashflowLeg->priceCF(instrESW->valueDate,
                                                     instrESW->callDate,
                                                     instrESW->eqLeg->EqCallSettleDate,
                                                     instrESW->isCallable,
                                                     control,
                                                     results,
                                                     knownCFs,
                                                     payment_dates);
            OutputNameConstSP cfOutput(new OutputName(instrESW->cashflowLeg->LegName));
            //results->storeGreek(CDoubleSP(CDouble::create(cfPrice)),"PRICE",cfOutput);

            writeLegRequestedOutput(control, results, cfPrice, ccyName, cfOutput);
        }

        if (!!(instrESW->divLeg))
        {
            divPrice = instrESW->divLeg->priceDiv(instrESW->discount,
                                                  instrESW->liborLeg,
                                                  instrESW->valueDate,
                                                  instrESW->eqLeg->asset,
                                                  instrESW->eqLeg->ccyTreatment,
                                                  instrESW->callDate,
                                                  instrESW->eqLeg->EqCallSettleDate,
                                                  instrESW->isCallable,
                                                  control,
                                                  results,
                                                  knownCFs,
                                                  payment_dates);

            OutputNameConstSP divOutput(new OutputName(instrESW->divLeg->LegName));
            //results->storeGreek(CDoubleSP(CDouble::create(divPrice)),"PRICE",divOutput);

            writeLegRequestedOutput(control, results, divPrice, ccyName, divOutput);
        }

        // aggregate result
        value = eqPrice + divPrice + cfPrice + irPrice;
        results->storePrice(value, ccyName);

        // FWD_AT_MAT
        InstrumentUtil::recordFwdAtMat(control,
            results,
            instrESW->eqLeg->eqRefixDates[instrESW->eqLeg->eqRefixDates.size()-1],
            instrESW->valueDate,
            instrESW->eqLeg->asset.get());

        //KNOWN_CASHFLOWS
        knownCFs.recordKnownCashFlows(control, results);

        //PAYMENT_DATES
        OutputRequestUtil::recordPaymentDates(control, results, &payment_dates);
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void CEquitySwapClosedFormProd::writeLegRequestedOutput(
    Control*            control,
    CResults*           results,
    double              legPrice,
    const string&       ccyISOCode,
    OutputNameConstSP&  outputName) const
{
    OutputRequest * request = 0;
    if (control && control->isPricing() && control->requestsOutput(OutputRequest::ESW_LEG_PRICE, request))
    {
        results->storeRequestResult(request,
                                    legPrice,
                                    ccyISOCode,
                                    outputName);
    }
}

DateTime CEquitySwap::endDate(const Sensitivity* sensControl) const
{
    // to do: this should be max of settlement of all pay dates.
    DateTime maxDate = DateTime(0,0);
    int i;

    // eq leg pay dates
    for (i = 0; i < eqLeg->payDates.size(); i++)
    {
        if (eqLeg->payDates[i] > maxDate)
        {
            maxDate = eqLeg->payDates[i];
        }
    }

    // lib leg pay dates
    for (i = 0; i < liborLeg->payDates.size(); i++)
    {
        if (liborLeg->payDates[i] > maxDate)
        {
            maxDate = liborLeg->payDates[i];
        }
    }

    // div leg pay dates (for cum case)
    if (!!divLeg)
    {
        for (i = 0; i < divLeg->payDates.size(); i++)
        {
            if (divLeg->payDates[i] > maxDate)
            {
                maxDate = divLeg->payDates[i];
            }
        }
    }

    // cf leg
    if (!!cashflowLeg)
    {
        for (i = 0; i < cashflowLeg->PayDates.size(); i++)
        {
            if (cashflowLeg->PayDates[i] > maxDate)
            {
                maxDate = cashflowLeg->PayDates[i];
            }
        }
    }

    // Check any components of the asset implementing the LastSensDate interface
    EndDateCollector collector(maxDate, sensControl);

    maxDate = collector.getMaxEndDate(IObjectSP(asset.getMO()));

    // add settlement and 1 day, if necessary.
    DateTime settleMaxDate = asset->settleDate(maxDate).rollDate(1);
    if (maxDate < settleMaxDate)
    {
        maxDate = settleMaxDate;
    }

    // now max with refix dates
    if (liborLeg->isFloating)
    {
        for (i = 0; i < liborLeg->refixDates.size(); i++)
        {
            MaturityPeriodSP rateMP;
            if (i==0 && liborLeg->frontStubType != "")
            {
                rateMP = MaturityPeriodSP(new MaturityPeriod(liborLeg->frontStubType));
            }
            else if (i==liborLeg->refixDates.size() - 1 && liborLeg->backStubType != "")
            {
                rateMP = MaturityPeriodSP(new MaturityPeriod(liborLeg->backStubType));
            }
            else
            {
                rateMP = MaturityPeriodSP(new MaturityPeriod(liborLeg->rateType));
            }

            DateTime refixEndDate = rateMP->toDate(liborLeg->refixDates[i]);

            if (refixEndDate > maxDate)
            {
                maxDate = refixEndDate;
            }
        }
    }

    // to do: avging out pay dates, avg in refix dates.
    return maxDate;
}

/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* CEquitySwap::createProduct(
    CClosedFormLN* model) const{
    return new CEquitySwapClosedFormProd(this);
}

bool CEquitySwap::sensShift(Theta* shift)
{
    DateTime newDate = shift->rollDate(valueDate);

    // Validate samples up to valueDate if not done yet
    validateHistoricalSamples();    

    // set fixings on valueDate
    eqLeg->setFixing(valueDate, newDate, asset, shift);

    liborLeg->setFixing(valueDate, newDate, eqLeg->avIn, discount.getSP());

    if (!!divLeg)
        divLeg->setFixing(valueDate, newDate, liborLeg);

    // validate against yield div between value date and roll date
    if (!!divLeg && !! (divLeg->synthDivs) && divLeg->useSynthetic)
    {
        DividendListSP div(divLeg->synthDivs->getAllDivsBetweenDates(valueDate, newDate));
        if (div->hasYieldDividend())
            throw ModelException("CEquitySwap::sensShift",
                                 "cannot roll over yield dividend for synthetic div when calc theta.");
    }

    // roll today
    valueDate = newDate;
    return true;
};

//** Methods for synthetic dividend tweaking: */

/** Returns name identifying vol for MU_PARALLEL */
string CEquitySwap::sensName(MuParallel* shift) const{
    return asset->getTrueName();
}

/** Returns name identifying vol for MU_POINTWISE */
string CEquitySwap::sensName(MuPointwise* shift) const{
    return asset->getTrueName();
}

/** Returns name identifying vol for MU_SPECIAL */
string CEquitySwap::sensName(MuSpecial* shift) const{
    return asset->getTrueName();
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_PARALLEL shift method */
bool CEquitySwap::sensShift(MuParallel* shift)
{
    if (!!divLeg && !! (divLeg->synthDivs))
    {
        divLeg->synthDivs->sensShift(shift, valueDate);
    }

    return true;
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_POINTWISE shift method */
bool CEquitySwap::sensShift(MuPointwise* shift)
{
    if (!!divLeg && !! (divLeg->synthDivs))
    {
        bool result = false;
        ExpiryArrayConstSP expiries = shift->getExpiries();

        DateTime bucketStartDate, bucketEndDate;
        shift->getBucketDates(bucketStartDate, bucketEndDate);

        result = divLeg->synthDivs->sensShift(shift, valueDate,
                                              bucketStartDate, bucketEndDate);
    }
    return true;
}

/** Shifts the object using given shift. This is a wrapper for the
    DividendList MU_S shift method */
bool CEquitySwap::sensShift(MuSpecial* shift)
{
    if (!!divLeg && !! (divLeg->synthDivs))
    {
        DateTime expiry;

        // get the current expiry date being tweaked and ensure that it
        // is not before today
        expiry = shift->getExpiry()->toDate(expiry);

        if (!expiry.isLess(valueDate))
        {
            divLeg->synthDivs->sensShift(shift, asset->getSpot(), expiry);
        }
    }
    return true; // CEquitySwap contains an asset which is sensitive to MuSpecial
}

// sets all future sampleLevels to 0.
// offset is used to associate sampleLevels[i] to sampleDates[i + offset]
// e.g. offset will be 1 for end levels, 0 for start levels
void CEquitySwap::zeroFutureSamples(const DateTimeArray&        sampleDates,
                                    DoubleArray&         sampleLevels,
                                    const DateTime&      valueDate,
                                    int offset)
{
    const DateTime tomorrow = valueDate.rollDate(1);
    if (sampleDates.size() < sampleLevels.size() + offset)
    {
        throw ModelException("EquitySwap::zeroFutureSamples", "Not enough sample dates!");
    }
    for (int iLevel=0; iLevel < sampleLevels.size(); iLevel++)
    {
        if (sampleDates[iLevel + offset].getDate() >= tomorrow.getDate())
        {
            sampleLevels[iLevel] = 0.0;
        }
    }
}


void CEquitySwap::setLibAvType()
{
    // set some variables
    bool isIn = false;
    bool isOut = false;

    liborLeg->ESWAveType = "NO_AVE";
    if(!!eqLeg->avIn && !eqLeg->avIn->isConsolidated)
    {
        isIn = true;
    }
    if(!!eqLeg->avOut && !eqLeg->avOut->accrueIRonFullNotional)
    {
        isOut = true;
    }

    if (isIn && !isOut) { liborLeg->ESWAveType = "AVE_IN"; }
    if (!isIn && isOut) { liborLeg->ESWAveType = "AVE_OUT"; }
    if (isIn && isOut)  { liborLeg->ESWAveType = "HYBRID"; }
}

/** Returns the name of the instrument's discount currency. */
string CEquitySwap::discountYieldCurveName() const {
    return discount.getName();
}


DRLIB_END_NAMESPACE
