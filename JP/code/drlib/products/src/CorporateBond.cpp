//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CorporateBond.cpp
//
//   Description : Corporate Bond Model
//
//   Author      : André Segger
//
//   Date        : 11 September 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CorporateBond.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/FD1FE2C.hpp"
#include "edginc/RiskProperty.hpp"

DRLIB_BEGIN_NAMESPACE

static CFieldConstSP cdsParSpreadsField;
static CFieldConstSP assetField;

bool CorporateBond::recurse(const CFieldConstSP& field,
                            const CClassConstSP& targetClass) const
{
    // this gets called as part of the tweaking and allows us to specify 
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted
    if (field == assetField ) {
        if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE                    ||
            targetClass == IRestorableWithRespectTo<VolParallel>::TYPE          ||
            targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE    ||
            targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE   ||
            targetClass == VegaMatrix::IShift::TYPE                     ||
            targetClass == VegaMatrix::IRestorableShift::TYPE           ||
            targetClass == RootTimeVega::IShift::TYPE                   ||
            targetClass == RootTimeVega::IRestorableShift::TYPE         ||
            targetClass == VegaSkewParallel::IShift::TYPE               ||
            targetClass == VegaSkewParallel::IRestorableShift::TYPE     ||
            targetClass == VegaSkewPointwise::IShift::TYPE              ||
            targetClass == VegaSkewPointwise::IRestorableShift::TYPE    || 
            targetClass == CreditSpreadRhoParallel::Shift::TYPE         ||
            targetClass == CreditSpreadRhoPointwise::IShift::TYPE       ||
            targetClass == ITweakableWithRespectTo<Spot>::TYPE          ||
            targetClass == IRestorableWithRespectTo<Spot>::TYPE         ) {
            return false;
        }
    }

    return true;
}

void CorporateBond::getMarket(const IModel* model, const MarketData* market)
{
    static const string method = "CorporateBond::getMarket";

    bond->getMarket(model, market);

    if (ClosedFormCDSPSandFA::TYPE->isInstance(model) || 
        ClosedFormFA::TYPE->isInstance(model) ||
	FD1FE2C::TYPE->isInstance(model))  {
        createE2Csensitivities = true;
    } else {
        createE2Csensitivities = false;
    }

    CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());

    // if asset has been provided, make sure it's a firm asset
    if ( ncAsset ) {
        FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

        if (!firmAsset) {
            throw ModelException(method, "Underlying credit asset must be a FirmAsset");
        }

        firmAsset->calculateProcessedVol(bond->getMaturityDate());

    }

    // default to par curve DCC
    bondAccrualDCC = DayCountConventionSP(cdsParSpreads->dayCountConv());
        
#ifdef CDS_BACKWARD_COMPATIBILITY
    cdsParSpreads->setBadDayConvention(BDC);
    cdsParSpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
#endif
        
}

void CorporateBond::validatePop2Object()
{
    static const string method = "CorporateBond::validatePop2Object";
    try {
        BDC = BadDayConventionSP(BadDayConventionFactory::make(bdc));

        //DateTime initialAccrualDate = bond->getAccrualStartDate();
        //if (bond->getCoupons(initialAccrualDate)->size()==0) {
        //    throw ModelException(method, "The underlying bond has no cashflow");
        //}
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

void CorporateBond::Validate() 
{
    static const string method = "CorporateBond::Validate";

    try {
        // validate the asset specific stuff at the generic level
        validate();

        if ( koBeforeIssueDate == true ) {
            throw ModelException(method, "KO before issue date is not supported yet!");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// -- Pricing Methods -- //

bool CorporateBond::getPutLevel(const DateTime& putDate,
                                double&         putLevel) const
{
    static const string method = "CorporateBond::getPutLevel";
    
    bool isPuttable;
    if (!!putSchedule && putSchedule->length() > 0) {
        isPuttable = putSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                         bond->getRedemption(false),
                                                         bond->getMaturityDate(),
                                                         putDate, 
                                                         putLevel);

        if (isPuttable) {
            if (putAdjustForAccrued == true) {
                putLevel += bond->getAccruedAtDate(putDate);
            }
        }
    } else {
        isPuttable = false;
        putLevel   = 0.0;
    }

    return isPuttable;
}

bool CorporateBond::getCallLevel(const DateTime& callDate,
                                 double&         callLevel) const
{
    static const string method = "CorporateBond::getCallLevel";
    
    bool isCallable;
    if (!!callSchedule && callSchedule->length() > 0) {
        isCallable = callSchedule->interpolateCVBSchedule(Schedule::ADJUST_NONE,
                                                          bond->getRedemption(false),
                                                          bond->getMaturityDate(),
                                                          callDate, 
                                                          callLevel);

        if (isCallable) {
            if (callAdjustForAccrued == true) {
                callLevel += bond->getAccruedAtDate(callDate);
            }
        }
    } else {
        isCallable = false;
        callLevel   = 0.0;
    }

    return isCallable;	       
}

YieldCurveSP CorporateBond::getRiskyCurve() const
{
    // create a risky curve
    IYieldCurveSP riskyCurve = cdsParSpreads.get()->makeRiskyCurve(*discount.get());

    if (!riskyCurve.get())
    {
        throw ModelException("CorporateBond::getRiskyCurve", "failed to create risky curve");
    }

    return YieldCurveSP(dynamic_cast<YieldCurve*>(riskyCurve.get()));
}

double CorporateBond::parSpreadsValue(bool           priceE2C,
                                      bool           generateRiskyStream,
                                      ObjectArraySP& riskyStreamDetails)  const
{
    static const string method = "CorporateBond::parSpreadsValue";

    double   cashFlowRiskyPV    = 0.0;
    double   accruedPayment     = 0.0;
    double   defaultPayment     = 0.0;
    double   value              = 0.0;
    DateTime startDate          = valueDate;
    DateTime initialAccrualDate = bond->getAccrualStartDate();

    // create a 'key' for fast indexing into zero curve
    auto_ptr<YieldCurve::IKey> discFactorKey(discount.get()->logOfDiscFactorKey());

    // get cash flows from Bond
    CashFlowArraySP cashFlows               = bond->getCoupons(initialAccrualDate);
    CashFlowArraySP cashFlowsWithRedemption = bond->getCashFlows();
    CashFlowArraySP redemptionPayments      = bond->getRedemptionPayments();
    DateTime        maturity                = bond->getMaturityDate();

    DefaultRatesSP defRates;
    CDSHelper::CParSpreadDefaultRates psAnnDefRates = CDSHelper::CParSpreadDefaultRates(valueDate);
    DateTimeArraySP critDates;

    if (inDefault) {
        double  recovery = getBondRecovery();
        double  notional = bond->getNotional(defaultDate);

        value = recovery * notional;
        if ( payAccruedUponDefault) {
            value += recovery * bond->getAccruedAtDate(defaultDate);
        }
    } else {
        if (priceE2C) {
	    CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
            FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

// interpolate vol for firm asset
            firmAsset->calculateProcessedVol(maturity);
        
// calculate the clean liquidity spread curve
            DefaultRatesSP cleanLiquiditySpreads(
                new CDSHelper::CParSpreadDefaultRates(valueDate));
                
            try {
                if (!firmAsset->getLiquiditySpreadCurve().getName().empty())  {
                    ICDSParSpreadsSP liquiditySpreads =
                        firmAsset->getLiquiditySpreadCurve().getSP();
#ifdef CDS_BACKWARD_COMPATIBILITY
                    liquiditySpreads->setBadDayConvention(BDC);
                    liquiditySpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
#endif     
                    cleanLiquiditySpreads = 
                        liquiditySpreads->defaultRates();
                }
            } catch (exception &e) {
                throw ModelException(&e, "CorporateBond::basePrice: Failed to calculate clean liquidity spread curve");
            }
	    
            defRates = DefaultRatesSP(new CDSHelper::CFirmAssetDefaultRates(firmAsset, cleanLiquiditySpreads));
            critDates = DateTimeArraySP(new DateTimeArray(0));
            DateTime currentDate = valueDate;
            while (currentDate <= maturity) {
                critDates->push_back(currentDate);
                currentDate = currentDate.rollDate(1);
            }
        } else {
// set up default rates and critical dates
            defRates = cdsParSpreads->defaultRates();                      
            DateTimeArraySP defaults   = defRates->getDefaultDates();
            critDates  = CDSHelper::createTimeLine(
		valueDate, 
		startDate, 
		*cashFlowsWithRedemption, 
		*defaults, 
		discount.getSP());
        }
	
        if ( (!!callSchedule && callSchedule->length() > 0) || 
             (!!putSchedule  && putSchedule->length()  > 0) ) {
            // cannot price separate legs independently - need two-state discounting algorithm
	    
            // create a risky curve
            YieldCurveSP riskyCurve = getRiskyCurve();

            int         currentCFIdx = cashFlowsWithRedemption->size()-1;
            DateTime    currentDate  = (*cashFlowsWithRedemption)[currentCFIdx].date;
            double      currentValue = (*cashFlowsWithRedemption)[currentCFIdx].amount;
            DateTime    previousDate = currentDate;
            double      accruedIR    = 0.0;
            double      putLevel     = 0.0;
            double      callLevel    = 0.0;
            currentDate              = currentDate.rollDate(-1);
            --currentCFIdx;

            while (currentDate >= valueDate) {
                accruedIR    = bond->getAccruedAtDate(currentDate);
                currentValue = riskyCurve->riskyPV(
		    currentDate, 
		    previousDate,
		    currentValue,
		    bond->getNotional(currentDate) + accruedIR,
		    useBondRecovery, 
		    bondRecovery);

                if ( currentCFIdx>= 0 && 
                     (*cashFlowsWithRedemption)[currentCFIdx].date.getDate() == currentDate.getDate() ) {
                    currentValue += (*cashFlowsWithRedemption)[currentCFIdx].amount;
                    --currentCFIdx;
                }

                // max with put level if applicable
                if ( getPutLevel(currentDate, putLevel)) {
                    currentValue =Maths::max(putLevel, currentValue);
                }

                // min with call level if applicable
                if ( getCallLevel(currentDate, callLevel)) {
                    currentValue =Maths::min(callLevel, currentValue);
                }

                previousDate = currentDate;
                currentDate  = currentDate.rollDate(-1);
            }

            value = currentValue;
        } else {
            // a numerical integration of the different legs - these are independent as they only depend on the
            // default probabilities

            // 1.) calculate the risky pv of the coupons and redemption payments
	    
	    cashFlowRiskyPV =  CDSHelper::calculateRiskyCashFlowPV(
		valueDate,
		maturity,
		cashFlows,
		discount.getSP(),
		&(*defRates),
		true); // isE2C
	    
	    cashFlowRiskyPV +=  CDSHelper::calculateRiskyCashFlowPV(
		valueDate,
		maturity,
		redemptionPayments,
		discount.getSP(),
		&(*defRates),
		true); // isE2C
	    
            // 2.) calculate the pv's of the expected default payment and expected accrued settlement
	    //  calculateDefaultPayments(startDate, maturity, critDates, cashFlows, defRates, 
            //                        discFactorKey, defaultPayment, accruedPayment);

	    DoubleArraySP notionals = bond->getNotionalsAtDates(critDates);
 	    DoubleArraySP notionals2(new DoubleArray(notionals->size()));
 	    for (int i=0; i < (notionals2->size());i++) {
 		(*notionals2)[i] = (*notionals)[i] * getBondRecovery();
 	    }
	    
	    CDSHelper::calculateDefaultPayments(
		valueDate, 
		startDate, //effDate,
		cashFlows, //feePayments
		notionals2, //notional * (1. - getBondRecovery()), //bond->getNotionalsAtDates(critDates);
		payAccruedUponDefault, //accrueFee, 
		bond->getAccrualStartDate(), //swapEffectiveDate,
		maturity,//
		&(*defRates),  //
		&(*bondAccrualDCC), //dcc,
		false, //isE2C?false:startFromEffDate,
		true, //isE2C, 
		critDates,//
		discount.getSP(),
		accruedPayment,  //
		defaultPayment); //
	    
	    accruedPayment *= getBondRecovery();

            // add the different components and return price
            value = cashFlowRiskyPV + accruedPayment + defaultPayment;
        }

        // RISKY_STREAM_DETAILS
        try {
            if (generateRiskyStream) {
                CashFlowArraySP cashFlowsWithRedemption = bond->getCashFlows();
                CashFlowArraySP rates                   = bond->getCouponRates();
                CashFlowArraySP dayFracs                = bond->getDayCountFractions();

                if ( cashFlowsWithRedemption->size() != rates->size() ||
                     cashFlowsWithRedemption->size() != dayFracs->size()) {
                    throw ModelException("CorporateBond::parSpreadsValue",
                        "Internal inconsistency in Bond data - please contact Derivatives Research");
                }
                     

                riskyStreamDetails = ObjectArraySP(new ObjectArray(6));

                DateTimeArraySP dates(new DateTimeArray(cashFlows->size()));
                DoubleArraySP   payment(new DoubleArray(cashFlows->size()));
                DoubleArraySP   discFact(new DoubleArray(cashFlows->size()));
                DoubleArraySP   defaultProb(new DoubleArray(cashFlows->size()));
                DoubleArraySP   ratesUsed(new DoubleArray(cashFlows->size()));
                DoubleArraySP   dayCountFractions(new DoubleArray(cashFlows->size()));

                for (int i=0; i<cashFlows->size(); i++) {
                    (*dates)[i]             = (*cashFlowsWithRedemption)[i].date;
                    (*payment)[i]           = (*cashFlowsWithRedemption)[i].amount;
                    (*discFact)[i]          = exp((discFactorKey->calc(valueDate, (*cashFlowsWithRedemption)[i].date)));
                    (*defaultProb)[i]       = defRates->calcTotalDefaultPV(
                                                    valueDate < maturity ? valueDate: maturity, 
                                                    (*cashFlowsWithRedemption)[i].date);
                    (*ratesUsed)[i]         = (*rates)[i].amount;
                    (*dayCountFractions)[i] = (*dayFracs)[i].amount;
                }

                (*riskyStreamDetails)[0] = dates;
                (*riskyStreamDetails)[1] = payment;
                (*riskyStreamDetails)[2] = discFact;
                (*riskyStreamDetails)[3] = defaultProb;
                (*riskyStreamDetails)[4] = ratesUsed;
                (*riskyStreamDetails)[5] = dayCountFractions;
            }
        }
        catch (exception&) {
            riskyStreamDetails = ObjectArraySP(   );
        }
    }

    return value;
}

void CorporateBond::priceParSpreads(CResults* results, Control* control) const {
    static const string method = "CorporateBond::priceParSpreads";
    try {
        ObjectArraySP riskyStream;
        double        value        = parSpreadsValue(false, true, riskyStream);

        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing()) {

            // RISKY_STREAM_DETAILS
            OutputRequest* request = control->requestsOutput(OutputRequest::RISKY_STREAM_DETAILS);
            if (request && !!riskyStream) {
                results->storeRequestResult(request, riskyStream);
            }

            // CLEAN_PRICE
            if ( control->requestsOutput(OutputRequest::CLEAN_PRICE, request) ) {

                double accInt = 0.0;
                if (inDefault) {
                    if (payAccruedUponDefault) {
                        accInt = bond->getAccruedAtDate(defaultDate) * getBondRecovery();
                    }
                } else {
                    accInt = bond->getAccruedAtDate(valueDate);
                }

                double cleanPrice = value - accInt;
               results->storeRequestResult(request, cleanPrice);
            }


            DateTime maturity = bond->getMaturityDate();
            HolidaySP hols(Holiday::weekendsOnly());

            DefaultRatesSP psDefRates;
            psDefRates = cdsParSpreads->defaultRates();

            CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();
            IObjectSP currentSpread = cdsParSpreads->getCurrentSpreadOrUntweakable(valueDate, maturity);
            addRequests(control, results, cleanSpreadCurve, currentSpread);
        }

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void CorporateBond::addRequests(
    Control*                  control,
    Results*                  results,
    CashFlowArraySP           cleanSpreadCurve,
    IObjectSP                 currentSpread) const
{
    static const string method = "CorporateBond::addRequests";
    try
    {
        OutputNameConstSP cdspsOutputName(new OutputName(cdsParSpreads->getName()));

        // CLEAN_DEFAULT_SPREAD_CURVE
        OutputRequest* request = control->requestsOutput(OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE);
        if (request && cleanSpreadCurve.get()) {
            IObjectSP    cflows(new CashFlowList(cleanSpreadCurve.get()));
            results->storeRequestResult(request, cflows, cdspsOutputName);
        }

        // CURRENT_SPREAD
        request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
        if (request) {
            results->storeRequestResult(request, currentSpread);
        }

        // IND_CDS_PAR_SPREAD
        // This one is the same as CURRENT_SPREAD but goes in its own packet and is qualified by curve name
        request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
        if (request) {
            results->storeRequestResult(request, currentSpread, cdspsOutputName);
        }

        // PAYMENT_DATES
        request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request) {
            DateTimeArraySP dates(paymentDates());
            OutputRequestUtil::recordPaymentDates(control,results,dates.get()); 
        }

        // KNOWN_CASHFLOWS
        request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request)  {
            CashFlowArraySP cfl(knownCashFlows());
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    cfl.get()); 
        }

        // DEFAULT_PROBABILITY
        request = control->requestsOutput(OutputRequest::DEFAULT_PROBABILITY);
        if (request)  {
            if ( !!asset) {
                DateTime maturity = bond->getMaturityDate();
                CAsset*    ncAsset   = const_cast<CAsset*>(asset.get());
                FirmAsset* firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
                double defaultProb = 1. - firmAsset->CalcNoDefaultProb(valueDate, maturity);
                results->storeRequestResult(request, defaultProb);
            } else {
                results->storeNotApplicable(request);
            }
        }

        // IMPLIED_DEFAULT_PROBABILITY
        request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
        try {
            if (request)  {
                if (!!cdsParSpreads) {
                    DateTime maturity = bond->getMaturityDate();
                    // set the rolling effective date - for bootstrapping
                    
                    DateTime effDate = cdsParSpreads->spotDate(valueDate);
                    DefaultRatesSP psDefRates =
                        cdsParSpreads->defaultRates();
        
                    double defaultProb = 1.0 - psDefRates->calcDefaultPV(effDate, maturity);
                    results->storeRequestResult(request, defaultProb);
                } else {
                    results->storeNotApplicable(request);
                }
            }
        } catch (exception&) {
            results->storeNotApplicable(request);
        }
        
        // THETA_ACCRUED_INTEREST
        request = control->requestsOutput(OutputRequest::THETA_ACCRUED_INTEREST);
        if (request) {
            double thetaAI = calcThetaAccruedInterest();
            results->storeRequestResult(request, thetaAI);
        } 

        // ACCRUED_INTEREST
        request = control->requestsOutput(OutputRequest::ACCRUED_INTEREST);
        if (request) {
            double accInt = 0.;
            if (inDefault) {
                if (payAccruedUponDefault) {
                    accInt = bond->getAccruedAtDate(defaultDate) * getBondRecovery();
                }
            } else {
                // accrued should be to the settlement date
                DateTime aiDate = bond->settles(valueDate);
                accInt = bond->getAccruedAtDate(aiDate);
            }
            results->storeRequestResult(request, accInt);
        } 

        // RECOVERY_VALUE
        request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
        if (request) {
            double defaultValue;
            
            defaultValue = getBondRecovery() * notional;
            results->storeRequestResult(request, defaultValue);
        } 

        // COUPON_DUE 
        request = control->requestsOutput(OutputRequest::COUPON_DUE);
        if (request) {
            AccrualCalendarArraySP coupons = bond->couponDue(valueDate);
            if (coupons->empty()) {
                results->storeNotApplicable(request); 
            } else {
                results->storeRequestResult(request, coupons); 
            }
        }    

        // ACCRUAL_CALENDAR 
        request = control->requestsOutput(OutputRequest::ACCRUAL_CALENDAR);
        if (request) {
            AccrualCalendarArraySP accData = bond->accrualCalendar(valueDate);
            if (accData->empty()) {
                results->storeNotApplicable(request); 
            } else {
                results->storeRequestResult(request, accData); 
            }
        }    
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}


double CorporateBond::calcThetaAccruedInterest()const
{
    static const string method = "CorporateBond::calcThetaAccruedInterest";
    
    try {
        double  thetaAI = 0.;
        if (!inDefault) {
            HolidaySP hols(Holiday::weekendsOnly());
            BadDayConventionSP badDay(BadDayConventionFactory::make("Following"));
            DateTime tommorrow = valueDate.rollDate(1);
            tommorrow = badDay->adjust(tommorrow, hols.get());
        
            double          accruedInterestNow      = bond->getAccruedAtDate(valueDate);
            double          accruedInterestTomorrow = bond->getAccruedAtDate(tommorrow);
            DateTime        initialAccrualDate      = bond->getAccrualStartDate();
            CashFlowArraySP cashFlows               = bond->getCoupons(initialAccrualDate);
        
            double totalCoupons = 0.0;
            for (int i=0 ; i<cashFlows->size() ; ++i) {
                if ( ((*cashFlows)[i].date  > valueDate && (*cashFlows)[i].date <= tommorrow) &&
                     !(initialAccrualDate >= (*cashFlows)[i].date)                            &&
                     !(valueDate <= initialAccrualDate && initialAccrualDate <= (*cashFlows)[i].date)) {
                    totalCoupons += (*cashFlows)[i].amount;
                }
            }

            // Account for the case where today and tommorrow are in different fee periods
            thetaAI = (accruedInterestTomorrow + totalCoupons) - accruedInterestNow;
        }
        return thetaAI;
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}

void CorporateBond::priceFirmAsset(CResults* results, Control* control) const {
    static const string method = "CorporateBond::priceParSpreads";
    try {
        ObjectArraySP riskyStream;
        double        value        = parSpreadsValue(true, true, riskyStream);

        results->storePrice(value, discount->getCcy());

        if (control && control->isPricing()) {

            // RISKY_STREAM_DETAILS
            OutputRequest* request = control->requestsOutput(OutputRequest::RISKY_STREAM_DETAILS);
            if (request ) {
                results->storeRequestResult(request, riskyStream);
            }

            // CLEAN_PRICE
            if ( control->requestsOutput(OutputRequest::CLEAN_PRICE, request) ) {

                double accInt = 0.0;
                if (inDefault) {
                    if (payAccruedUponDefault) {
                        accInt = bond->getAccruedAtDate(defaultDate) * getBondRecovery();
                    }
                } else {
                    accInt = bond->getAccruedAtDate(valueDate);
                }

                double cleanPrice = value - accInt;
               results->storeRequestResult(request, cleanPrice);
            }

            DateTime maturity = bond->getMaturityDate();
            HolidaySP hols(Holiday::weekendsOnly());

            DefaultRatesSP psDefRates =
                cdsParSpreads->defaultRates();

            CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();
            IObjectSP currentSpread = cdsParSpreads->getCurrentSpreadOrUntweakable(valueDate, maturity);
            addRequests(control, results, cleanSpreadCurve, currentSpread);
        }

    } catch (exception& e) {
        throw ModelException(&e, method);
    }
}


double CorporateBond::getBondRecovery() const
{
    double recovery;
    if (!useBondRecovery) {
        recovery = cdsParSpreads->getRecovery();
    } else {
        recovery = bondRecovery;
    }
    return recovery;
}

// -- The model implementations and product definitions -- //

/** private class. */
class CorporateBondClosedFormCDSPS: public ClosedFormCDSPS::IProduct{
private:
    const CorporateBond* cf; // a reference

public:
    CorporateBondClosedFormCDSPS(const CorporateBond* cf): cf(cf){}

    void price(ClosedFormCDSPS* model,
               Control*         control, 
               CResults*        results) const{
        cf->priceParSpreads(results, control);
    }
};
    
/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPS::IProduct* CorporateBond::createProduct(
    ClosedFormCDSPS* model) const{
    return new CorporateBondClosedFormCDSPS(this);
}


/** private class */
class CorporateBondClosedFormFA: public ClosedFormFA::IProduct{
private:
    const CorporateBond* cf; // a reference

public:
    CorporateBondClosedFormFA(const CorporateBond* cf): cf(cf){}

    void price(ClosedFormFA* model,
               Control*         control, 
               CResults*        results) const{
        cf->priceFirmAsset(results, control);
    }
};
    
/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormFA::IProduct* CorporateBond::createProduct(
    ClosedFormFA* model) const{
    static const string method = "CorporateBond::createProduct";

    if (!!callSchedule && callSchedule->length() >0 ) {
        throw ModelException(method, "Can not price callable bonds in the closed for E2C model - use the finite difference model instead");
    }

    if (!!putSchedule && putSchedule->length() >0 ) {
        throw ModelException(method, "Can not price puttable bonds in the closed for E2C model - use the finite difference model instead");
    }

    return new CorporateBondClosedFormFA(this);
}


/** private class. */
class CorporateBondClosedFormCDSPSandFA: public ClosedFormCDSPSandFA::IProduct{
private:
    const CorporateBond* cf; // a reference

public:
    CorporateBondClosedFormCDSPSandFA(const CorporateBond* cf): cf(cf){}

    void price(ClosedFormCDSPSandFA* model,
               Control*              control, 
               CResults*             results) const{
        if (model->isPricePar()) {
            cf->priceParSpreads(results, control);
        } else {
            cf->priceFirmAsset(results, control);
        }
    }
};
    
/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPSandFA::IProduct* CorporateBond::createProduct(
    ClosedFormCDSPSandFA* model) const{
    static const string method = "CorporateBond::createProduct";

    if (!!callSchedule && callSchedule->length() >0 ) {
        throw ModelException(method, "Can not price callable bonds in the closed for E2C model - use the finite difference model instead");
    }

    if (!!putSchedule && putSchedule->length() >0 ) {
        throw ModelException(method, "Can not price puttable bonds in the closed for E2C model - use the finite difference model instead");
    }

    return new CorporateBondClosedFormCDSPSandFA(this);
}


// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime CorporateBond::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime CorporateBond::endDate(const Sensitivity* sensControl) const 
{
    return bond->getMaturityDate();
}

/** when do payments occur ? */
DateTimeArraySP CorporateBond::paymentDates() const 
{
    CashFlowArraySP cashFlows = bond->getCashFlows();
    DateTimeArraySP paydates(new DateTimeArray(cashFlows->size()));

    for (int i = 0; i < paydates->size(); i++)  {
        (*paydates)[i] = (*cashFlows)[i].date;
    }
    return paydates;
}

/** Returns all known cash flows */
CashFlowArraySP CorporateBond::knownCashFlows()const
{
    CashFlowArraySP cashFlows = bond->getCashFlows();
    return cashFlows;
}

/** IInstrumentAsAsset interface
    Returns a date after which the instrument can no longer be used as
    an asset  */
DateTime CorporateBond::maturityDate() const {
    return bond->getMaturityDate();
}
    
/** Returns the coupons that the instrument will make during
    its lifetime. This can include historic payments. */
// Note: Unlike ConvBond am returning ALL coupons
CashFlowArraySP CorporateBond::getCoupons() const {
    static const string method = "CorporateBond::getCoupons";
    try {
        CashFlowArraySP myCoupons = bond->getCashFlows();
        (*myCoupons)[myCoupons->size()-1].amount -= bond->getRedemption();
        return myCoupons;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns the dates on which the instrument has to be held in order to
    hold the right to the corresponding coupon as returned by 
    getCoupons() */
DateTimeArraySP CorporateBond::getExCouponDates() const {
    return bond->getExCouponDates();
}
    
/** Returns the yield curve used for discounting */
YieldCurveConstSP CorporateBond::getDiscount() const {
    return discount.getSP();
}
    
/** Returns the accured interest (if any) to date */
double CorporateBond::getAccrued() const {
    return bond->getAccruedAtDate(getValueDate());
}


/** Returns name identifying vol for vega parallel */
string CorporateBond::sensName(DeltaToCredit* shift) const 
{
    static const string method = "CorporateBond::sensName(DeltaToCredit)"; 

    if (createE2Csensitivities && !!asset) 
    {
        if (!FirmAsset::TYPE->isInstance(asset.get())) {
            throw  ModelException(method, "asset is not of type FirmAsset");
        }
        
        OutputNameArrayConstSP names(
            RiskProperty<Spot>().subjectNames(asset.getSP()));

        if (names->size() > 1) {
            throw ModelException(method,
                    "E2C currently not supported for baskets!");
        }

        return names->empty() ? "" : (*names)[0]->toString();
    }
    else {
        return "";
    }
}

/** Shifts the object using given shift */
bool CorporateBond::sensShift(DeltaToCredit* shift)
{
    static const string method = "CorporateBond::sensShift(DeltaToCredit)"; 

    try {
        if (createE2Csensitivities) {
            if (!FirmAsset::TYPE->isInstance(asset.get())) {
                throw  ModelException(method, "asset is not of type FirmAsset");
            }

            TweakOutcome outcome = PropertyTweakHypothesis<Spot>(
                shift->getShiftSize(), shift->getMarketDataName()).
                    applyTo_TweakOutcome(asset.getSP());

            shift->setInitialValue(outcome.oldValue());
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }

    return false; // none of our components has a DeltaToCredit type sensitivity
}

// for reflection
CorporateBond::CorporateBond(): Generic1FactorCredit(TYPE), callAdjustForAccrued(true), 
                                putAdjustForAccrued(true), bondRecovery(0), useBondRecovery(false),
                                payAccruedUponDefault(true), koBeforeIssueDate(false), 
                                inDefault(false), createE2Csensitivities(false), 
                                e2cBasePriceCalculated(false), e2cBasePrice(0.0)
{}

class CorporateBondHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorporateBond, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
        IMPLEMENTS(ClosedFormFA::IIntoProduct);
        IMPLEMENTS(ClosedFormCDSPSandFA::IIntoProduct);
        IMPLEMENTS(FD1FGeneric::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(DeltaToCredit::IShift);
        IMPLEMENTS(ObjectIteration::IOverride);
        IMPLEMENTS(IGetMarket);
        IMPLEMENTS(IInstrumentAsAsset);
        EMPTY_SHELL_METHOD(defaultCorporateBond);

        FIELD(bond,                                     "the bond");
        FIELD(callSchedule,                             "call schedule");
        FIELD_MAKE_OPTIONAL(callSchedule);
        FIELD(callAdjustForAccrued,              "true = settle accrued when bond is called");
        FIELD_MAKE_OPTIONAL(callAdjustForAccrued);
        FIELD(putSchedule,                              "put schedule");
        FIELD_MAKE_OPTIONAL(putSchedule);
        FIELD(putAdjustForAccrued,               "true = settle accrued when bond is put");
        FIELD_MAKE_OPTIONAL(putAdjustForAccrued);
        FIELD(useBondRecovery,                   "use passed swap recovery value");
        FIELD(bondRecovery,                      "defaults to par recovery rate unless useBondRecovery = true");
        FIELD_MAKE_OPTIONAL(bondRecovery);
        FIELD(payAccruedUponDefault,             "make accrued payment on default?");
        FIELD_MAKE_OPTIONAL(payAccruedUponDefault);
        FIELD(koBeforeIssueDate,                 "Knockout when bond defaults before issue date");
        FIELD(inDefault,                         "whether the bond has already defaulted");
        FIELD_MAKE_OPTIONAL(inDefault);
        FIELD(defaultDate,                       "the date on which the bond has defaulted");
        FIELD_MAKE_OPTIONAL(defaultDate);
        FIELD_MAKE_OPTIONAL(koBeforeIssueDate);
        FIELD(bdc,                               "bad day convention");
        FIELD(bondAccrualDCC,                           "swap accrual day count convention");
        FIELD(BDC, "bad day convention");
        FIELD(createE2Csensitivities,            "whether to create debt/equity sensitivities");
        FIELD(e2cBasePriceCalculated,            "internal field");
        FIELD(e2cBasePrice,                      "internal field");
        FIELD_MAKE_TRANSIENT(BDC);
        FIELD_MAKE_TRANSIENT(bondAccrualDCC);
        FIELD_MAKE_TRANSIENT(createE2Csensitivities);
        FIELD_MAKE_TRANSIENT(e2cBasePriceCalculated);
        FIELD_MAKE_TRANSIENT(e2cBasePrice);

        // look up field for use on recurse
        cdsParSpreadsField = clazz->getSuperClass()->getDeclaredField("cdsParSpreads");
        assetField         = clazz->getSuperClass()->getDeclaredField("asset");
    }

    static IObject* defaultCorporateBond(){
        return new CorporateBond();
    }
};

CClassConstSP const CorporateBond::TYPE = CClass::registerClassLoadMethod(
    "CorporateBond", typeid(CorporateBond), CorporateBondHelper::load);

bool  CorporateBondLoad() {
    return (CorporateBond::TYPE != 0);
}

DRLIB_END_NAMESPACE

