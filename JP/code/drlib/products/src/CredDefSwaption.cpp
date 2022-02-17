//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : CredDefSwaption.cpp
//
//   Description : Credit default swaption
//
//   Author      : André Segger
//
//   Date        : 06 June 2003
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CredDefSwap.hpp"
#include "edginc/Class.hpp"
#include "math.h"
#include "edginc/Maths.hpp"
#include "edginc/Nrfns.hpp"
#include "edginc/Format.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/ZeroCurve.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/Delta.hpp"
#include "edginc/DayCountConventionFactory.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/imsl.h"
#include "edginc/CDSHelper.hpp" 
#include "edginc/CredDefSwaption.hpp"
#include "edginc/Black.hpp"
#include "edginc/FlatVol.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/Spot.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/Spot.hpp"
#include <math.h>

DRLIB_BEGIN_NAMESPACE

static CFieldConstSP cdsParSpreadsField;
static CFieldConstSP assetField;

bool CredDefSwaption::recurse(const CFieldConstSP& field,
                              const CClassConstSP& targetClass) const
{
    // this gets called as part of the tweaking and allows us to specify 
    // whether the fields within this class should be tweaked or not. The
    // target class indicates what is being shifted
    if (field == cdsParSpreadsField) {
        if (targetClass == CreditSpreadRhoParallel::Shift::TYPE  ||
            targetClass == CreditSpreadRhoPointwise::IShift::TYPE) {
            return false;
        }
    }

    if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE                    ||
        targetClass == IRestorableWithRespectTo<VolParallel>::TYPE          ||
        targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE  ||
        targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE ||
        targetClass == VegaMatrix::IShift::TYPE                     ||
        targetClass == VegaMatrix::IRestorableShift::TYPE           ||
        targetClass == RootTimeVega::IShift::TYPE                   ||
        targetClass == RootTimeVega::IRestorableShift::TYPE         ||
        targetClass == VegaSkewParallel::IShift::TYPE               ||
        targetClass == VegaSkewParallel::IRestorableShift::TYPE     ||
        targetClass == VegaSkewPointwise::IShift::TYPE              ||
        targetClass == VegaSkewPointwise::IRestorableShift::TYPE    ) {
        return false;
    }

    if (field == assetField) {
        if (targetClass == ITweakableWithRespectTo<VolParallel>::TYPE                    ||
            targetClass == IRestorableWithRespectTo<VolParallel>::TYPE          ||
            targetClass == ITweakableWithRespectTo<VolPointwise>::TYPE  ||
            targetClass == IRestorableWithRespectTo<VolPointwise>::TYPE ||
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

void CredDefSwaption::getMarket(const IModel* model, const MarketData* market)
{

    market->GetReferenceDate(valueDate);

    if (ClosedFormCDSPSandFA::TYPE->isInstance(model) || 
        ClosedFormFA::TYPE->isInstance(model))  {
        createE2Csensitivities = true;
    } else {
        createE2Csensitivities = false;
    }

    // create a vol object
    string       volName = cdsParSpreads->getName();
    HolidaySP    hols(Holiday::noHolidays());
    TimeMetricSP timeMetric(new TimeMetric(1.0,hols.get()));
    FlatVolSP    flatVol(new FlatVol(volName, 
                                     valueDate,
                                     timeMetric.get(),
                                     fwdSwapSpreadVol));

    vol = CVolBaseWrapper(flatVol);
}

void CredDefSwaption::Validate() 
{
    static const string method = "CredDefSwaption::Validate";

    try 
    {
        // validate the asset specific stuff at the generic level
        validate();

        // May need to default swapRecovery to equal par spread recovery
        if (!useSwapRecovery)
        {
            swapRecovery = cdsParSpreads->getRecovery();
        }

        if (payAccruedFee)
        {
            swpAccrualDCC = DayCountConventionSP(DayCountConventionFactory::make(dcc));
        }
        else
        {
            // default to par curve DCC
            swpAccrualDCC = DayCountConventionSP(cdsParSpreads->dayCountConv());
        }

        BDC = BadDayConventionSP(BadDayConventionFactory::make(bdc));

#ifdef CDS_BACKWARD_COMPATIBILITY
        cdsParSpreads->setBadDayConvention(BDC);
        cdsParSpreads->setHolidays(HolidaySP(Holiday::weekendsOnly()));
#endif

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Returns name identifying vol for asset vega parallel */
string CredDefSwaption::sensName(CreditSpreadVegaParallel* shift) const
{
    string vegaName;

    OutputNameArrayConstSP names(
        RiskProperty<VolParallel>().subjectNames(vol.getSP()));

    if (names->size() == 1 ) {
        vegaName = (*names)[0]->toString();
    } else if (names->size() > 1 ) {
        throw ModelException("CredDefSwaption::sensName(CreditSpreadVegaParallel)",
                "CreditSpreadVegaParallel not supported for options on CDS baskets!");
    }
    return vegaName;
}

/** Shifts the object using given shift */
bool CredDefSwaption::sensShift(CreditSpreadVegaParallel* shift)
{
    try{
        PropertyTweakHypothesis<VolParallel>(shift->getShiftSize(),
                                             shift->getMarketDataName()).
            applyTo(vol.getSP());
    } catch (exception& e){
        throw ModelException(e, "CredDefSwaption::sensShift(CreditSpreadVegaParallel)");
    }
    return false; // none of our components has a CreditSpreadVegaParallel type sensitivity
}

// -- Pricing Methods -- //

void CredDefSwaption::priceParSpreads(CResults* results, Control* control) const{
    static const string method = "CredDefSwaption::priceParSpreads";
    try {
        double value    = 0.0;
        double variance = 0.0;
        double cdsValue = 0.0;
        double bpValue  = 0.0;
        DateTime maturity = swapMaturityDate;
        // set the rolling effective date - for bootstrapping
        DateTime effDate = cdsParSpreads->spotDate(valueDate);
        DateTime swapEffDate = (swapEffectiveDate>effDate)?swapEffectiveDate:effDate;

        DefaultRatesSP psDefRates;

        double fwdSpread = 0.0;
        //CashFlowArraySP feePayments;

        if ( valueDate <= optionMaturityDate ) {
            if (priceViaParCurves)
            {
                psDefRates = cdsParSpreads->defaultRates();
            }
            else // use credit spread methodology for calculating CREDIT_SPREAD_RHO's
            {
                buildCreditSpreadCurve(); // build if needed
                CDSHelper::CParSpreadDefaultRates psAnnDefRates(creditSpreadCurve.get(), discount.getSP(), 
                                                    valueDate, cdsParSpreads->getRecovery(), &maturity);
                psDefRates = psAnnDefRates.createFwdContDefCurve(valueDate);
            }

            // choose how to interpolate the vol - go for traditional route for now
            LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                                       swapStrike, 
                                                       valueDate, 
                                                       optionMaturityDate,
                                                       false));
            // interpolate the vol using our LN request
            CVolProcessedSP  procVol(vol->getProcessedVol(volRequest.get(), NULL));
            CVolProcessedBS* volBS = dynamic_cast<CVolProcessedBS*>(procVol.get());
            if ( !volBS ) {
                throw ModelException(method, "Processed vol object is of non Black-Scholes type");
            }
        
            // calculate the variance
            variance = volBS->CalcVar(valueDate, optionMaturityDate);

            fwdSpread = CDSHelper::CDSParSpread(valueDate,
                                                swapEffectiveDate,
                                                frequency,
                                                1.0,
                                                swapRecovery,
                                                payAccruedFee,
                                                swapEffectiveDate,
                                                maturity,
                                                discount.getSP(),
                                                false, /* isE2C */
                                                psDefRates.get(),
                                                swpAccrualDCC.get());

            // calculate the fair value in basis points
            value     = Black::price(!isCall, fwdSpread, swapStrike, 1.0, variance);

            // do Edgworth expansion if required
            if (doEdgworthExpansion) {
	        // calculate derivatives of lognormal distribution - Corrado/Tsu

                double d  = (log(fwdSpread/swapStrike) +.5 * variance) / sqrt(variance);
                double Q3 = 1/6.  * fwdSpread * sqrt(variance) * ((2.*sqrt(variance)-d)*N1Density(d) - variance * N1(d));
                double Q4 = 1/24. * fwdSpread * sqrt(variance) * 
                            ((d*d - 1.0 - 3*sqrt(variance)*(d-sqrt(variance)))*N1Density(d) + pow(sqrt(variance),3.0)*N1(d));

                value = value + skewness * Q3 + kurtosis * Q4;
            }

            bpValue   = value;
   
	        auto_ptr<YieldCurve::IKey> discFactorKey(discount.get()->logOfDiscFactorKey());
	        
	        /* create fee payments array */
	        CashFlowArraySP feePayments = CDSHelper::calculateFeePayments(
	            swapEffectiveDate,
	            maturity,
	            frequency,
	            swpAccrualDCC.get(), //dcc,
	            1.0,  /* notional */
	            1.0); /* spread */
	        
	        double accruedPayment, defaultPayment, cashFlowRiskyPV;
	        
	        CDSHelper::calculateDefaultPayments(
	            valueDate,
	            swapEffectiveDate,
	            feePayments,
	            1.0, //notional,
	            swapRecovery, //recovery,
	            payAccruedFee, //accrueFee,
	            swapEffectiveDate,
	            maturity,
	            discount.getSP(),
	            //discFactorKey,
	            psDefRates.get(), //defRates,
	            swpAccrualDCC.get(), //dcc,
	            false,/* isE2C */
	            true, /* startFromEffDate */ 
	            true, /*fullTimeLine */
	            accruedPayment,
	            defaultPayment,
	            cashFlowRiskyPV);
	    
	        value    *= notional * (cashFlowRiskyPV + accruedPayment);
	    
            if (!optionKOOnDefault && !isCall) {
                double defaultPV =  1.0 - psDefRates->calcDefaultPV(valueDate, swapEffectiveDate);
                value += (1-swapRecovery) * notional * defaultPV * discount->pv(swapEffectiveDate);
            } else if (!isCall && addProtectionForPut){
                // the quotes for single stock puts include short term protection until maturity ....
                // to get correct implieds for market quotes and ease booking of these instruments,
                // we include the price in the 

                double shortAccruedPayment, shortDefaultPayment, shortCashFlowPV;

	            /* create fee payments array */
	            CashFlowArraySP shortStubPayments = CDSHelper::calculateFeePayments(valueDate,
		                                                                            optionMaturityDate,
		                                                                            frequency,
		                                                                            swpAccrualDCC.get(), //dcc,
		                                                                            1.0,  // notional
		                                                                            0.0); // spread

	            CDSHelper::calculateDefaultPayments(valueDate,
		                                            valueDate,
		                                            shortStubPayments,
		                                            notional,
		                                            swapRecovery,
		                                            payAccruedFee,
		                                            valueDate,
		                                            optionMaturityDate,
		                                            discount.getSP(),
                                                        //discFactorKey,
		                                            psDefRates.get(),
		                                            swpAccrualDCC.get(),
		                                            false,
		                                            true,
		                                            true,
		                                            shortAccruedPayment,
		                                            shortDefaultPayment,
		                                            shortCashFlowPV);

                value += shortDefaultPayment;
            }

        } else {
            value = 0.0;
        }

        results->storePrice(value, discount->getCcy());
        if (control && control->isPricing() && valueDate <= optionMaturityDate ) {
	        CashFlowArraySP feePayments = CDSHelper::calculateFeePayments(
		    swapEffectiveDate,//effDate,
		    maturity,
		    frequency,
		    swpAccrualDCC.get(),//dcc,
		    notional,
		    swapStrike); //spread
	        
	        double accruedPayment, defaultPayment, cashFlowRiskyPV;
	    
	        /* create a 'key' for fast indexing into zero curve */
	        auto_ptr<YieldCurve::IKey> discFactorKey(discount->logOfDiscFactorKey());
	        
	        CDSHelper::calculateDefaultPayments(
		    valueDate,
		    swapEffectiveDate,//effDate,
		    feePayments,
		    notional,
		    swapRecovery,//recovery,
		    payAccruedFee,//accrueFee,
		    swapEffectiveDate,
		    maturity,
		    discount.getSP(),
		    //discFactorKey,
		    psDefRates.get(),//defRates,
		    swpAccrualDCC.get(),//dcc,
		    false, /* isE2C */
		    true,  /* startFromEffDate */ 
		    false, /* fullTimeLine */
		    accruedPayment,
		    defaultPayment,
		    cashFlowRiskyPV) ;
	    
	        defaultPayment *= exp(discFactorKey->calc(valueDate, swapEffectiveDate)); //effDate));
	        
	        cdsValue = cashFlowRiskyPV + accruedPayment - defaultPayment;
	    	        
            HolidaySP hols(Holiday::weekendsOnly());
            CashFlowArraySP cleanSpreadCurve = psDefRates->getCleanSpreadCurve();
            IObjectSP currentSpread = cdsParSpreads->getCurrentSpreadOrUntweakable(effDate, maturity);
            string ccyName = discount->getName();
            addRequests(control, results, cleanSpreadCurve, currentSpread, bpValue, cdsValue, fwdSpread);
         }
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

void CredDefSwaption::addRequests(Control* control,
                              Results* results,
                              CashFlowArraySP cleanSpreadCurve,
                              IObjectSP currentSpread,
                              const double bpValue,
                              const double cdsValue,
                              const double fwdSpread) const
{
    static const string method = "CredDefSwaption::addRequests";

    try
    {
        OutputNameConstSP cdpsOutputName(new OutputName(cdsParSpreads->getName()));

        // CLEAN_DEFAULT_SPREAD_CURVE
        OutputRequest* request = control->requestsOutput(OutputRequest::CLEAN_DEFAULT_SPREAD_CURVE);
        if (request && cleanSpreadCurve.get()) 
        {
            IObjectSP    cflows(new CashFlowList(cleanSpreadCurve.get()));
            results->storeRequestResult(request, cflows, cdpsOutputName);
        }

        // CURRENT_SPREAD
        request = control->requestsOutput(OutputRequest::CURRENT_SPREAD);
        if (request) 
        {
            results->storeRequestResult(request, currentSpread); 
        }

        // IND_CDS_PAR_SPREAD
        // This one is the same as CURRENT_SPREAD but goes in its own packet and is qualified by curve name
        request = control->requestsOutput(OutputRequest::IND_CDS_PAR_SPREAD);
        if (request) {
            results->storeRequestResult(request, currentSpread, cdpsOutputName);
        }

        // IND_VOL
        request = control->requestsOutput(OutputRequest::IND_VOL);
        if (request) 
        {
            results->storeRequestResult(request, fwdSwapSpreadVol); 
        }

        // IMPLIED_DEFAULT_PROBABILITY
        request = control->requestsOutput(OutputRequest::IMPLIED_DEFAULT_PROBABILITY);
        if (request)  {
            try {
                if (!!cdsParSpreads) {
                    // set the rolling effective date - for bootstrapping
                    DateTime effDate = cdsParSpreads->spotDate(valueDate);

                    DefaultRatesSP psDefRates = 
                        cdsParSpreads->defaultRates();
        
                    double defaultProb = 1.0 - psDefRates->calcDefaultPV(effDate, swapMaturityDate);
                    results->storeRequestResult(request, defaultProb);
                } else {
                    results->storeNotApplicable(request);
                }
            } catch (exception&) {
                results->storeNotApplicable(request);
            }
        }
        
        // ACCRUED_INTEREST
        request = control->requestsOutput(OutputRequest::ACCRUED_INTEREST);
        if (request)
        {
            double accInt = 0.0;
            results->storeRequestResult(request, accInt);
        } 

        // RECOVERY_VALUE
        request = control->requestsOutput(OutputRequest::RECOVERY_VALUE);
        if (request)
        {
            double recovery;
            double defaultValue;
            if (!useSwapRecovery) {
                recovery = cdsParSpreads->getRecovery();
            } else {
                recovery = swapRecovery;
            }

            defaultValue = recovery * notional;
            results->storeRequestResult(request, defaultValue);
        } 

        // BASIS POINT PRICE
        request = control->requestsOutput(OutputRequest::BASIS_POINT_VALUE);
        if (request) {
            results->storeRequestResult(request, bpValue);
        } 
        
        // CDS value
        request = control->requestsOutput(OutputRequest::CDS_VALUE);
        if (request) {
            results->storeRequestResult(request, cdsValue);
        } 

        // FORWARD_CDS_SPREAD
        request = control->requestsOutput(OutputRequest::FORWARD_CDS_SPREAD);
        if (request) {
            results->storeRequestResult(request, fwdSpread);
        } 
    }
    catch (exception& e) 
    {
        throw ModelException(&e, method);
    }
}

void CredDefSwaption::priceFirmAsset(CResults* results, Control* control) const{
    static const string method = "CredDefSwaption::priceFirmAsset";

    throw ModelException(method, "Not implemented yet");
}

// -- The model implementations and product definitions -- //

/** private class. */
class CredDefSwaptionClosedFormCDSPS: public ClosedFormCDSPS::IProduct{
private:
    const CredDefSwaption* cf; // a reference

public:
    CredDefSwaptionClosedFormCDSPS(const CredDefSwaption* cf): cf(cf){}

    void price(ClosedFormCDSPS* model,
               Control*         control, 
               CResults*        results) const{
        cf->priceParSpreads(results, control);
    }
};
    
/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormCDSPS::IProduct* CredDefSwaption::createProduct(
    ClosedFormCDSPS* model) const{
    return new CredDefSwaptionClosedFormCDSPS(this);
}


/** private class */
class CredDefSwaptionClosedFormFA: public ClosedFormFA::IProduct{
private:
    const CredDefSwaption* cf; // a reference

public:
    CredDefSwaptionClosedFormFA(const CredDefSwaption* cf): cf(cf){}

    void price(ClosedFormFA*    model,
               Control*         control, 
               CResults*        results) const{
        cf->priceFirmAsset(results, control);
    }
};
    
/** Implementation of ClosedFormCDSPS::IntoProduct interface */
ClosedFormFA::IProduct* CredDefSwaption::createProduct(
    ClosedFormFA* model) const{
    return new CredDefSwaptionClosedFormFA(this);
}


/** private class. */
class CredDefSwaptionClosedFormCDSPSandFA: public ClosedFormCDSPSandFA::IProduct{
private:
    const CredDefSwaption* cf; // a reference

public:
    CredDefSwaptionClosedFormCDSPSandFA(const CredDefSwaption* cf): cf(cf){}

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
ClosedFormCDSPSandFA::IProduct* CredDefSwaption::createProduct(
    ClosedFormCDSPSandFA* model) const{
    return new CredDefSwaptionClosedFormCDSPSandFA(this);
}


// -- Stuff every product needs to do -- //

/** what's today ? */
DateTime CredDefSwaption::getValueDate() const {
    return valueDate;
}

/** when to stop tweaking */
DateTime CredDefSwaption::endDate(const Sensitivity* sensControl) const 
{
    // return (*feePayments)[feePayments->size()-1].date;;
    return swapMaturityDate;
}

/** Returns name identifying yield curve for CREDIT_SPREAD_RHO_PARALLEL */
string CredDefSwaption::sensName(CreditSpreadRhoParallel* shift) const
{
    return cdsParSpreads->getName(); 
}

/** Shifts the object using given shift */
bool CredDefSwaption::sensShift(CreditSpreadRhoParallel* shift)
{
   static const string method = "CredDefSwaption::sensShift";
   try 
   {
       // Set to price via credit spreads
       priceViaParCurves = false;

       buildCreditSpreadCurve(); // build if needed
       // Do the shift
       creditSpreadCurve->sensShift(shift);

       return false; 
   }
   catch (exception &e) 
   {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void CredDefSwaption::sensRestore(CreditSpreadRhoParallel* shift)
{
    static const string method = "CredDefSwaption::sensRestore";
    try
    {
        creditSpreadCurve = CreditSpreadCurveSP(copy(csRestore.get()));
        priceViaParCurves = true;
    }
    catch (exception &e) 
    {
        throw ModelException(&e, method);
    }
}

/** Returns name identifying yield curve for CREDIT_SPREAD_RHO_POINTWISE */
string CredDefSwaption::sensName(CreditSpreadRhoPointwise* shift) const
{
    return cdsParSpreads->getName(); 
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryArrayConstSP CredDefSwaption::sensExpiries(CreditSpreadRhoPointwise* shift) const {
    buildCreditSpreadCurve(); // build if needed
    return creditSpreadCurve->getExpiries();
}

/** Shifts the object using given shift */
bool CredDefSwaption::sensShift(CreditSpreadRhoPointwise* shift)
{
   static const string method = "CredDefSwaption::sensShift";
   try 
   {
       // set to price via credit spreads
       priceViaParCurves = false;

       buildCreditSpreadCurve(); // build if needed
       // Do the shift
       creditSpreadCurve->sensShift(shift);
       return false; // dont want to go on a tweak the par curve rho which inherits from CreditSpreadCurve
   }
   catch (exception &e) 
   {
      throw ModelException(&e, method);
   }
}

/** Restores the object to its original form */
void CredDefSwaption::sensRestore(CreditSpreadRhoPointwise* shift)
{
    static const string method = "CredDefSwaption::sensRestore";
    try
    {
        creditSpreadCurve = CreditSpreadCurveSP(copy(csRestore.get()));
        priceViaParCurves = true;
    }
    catch (exception &e) 
    {
        throw ModelException(&e, method);
    }
}

// Derive a credit spread curve from the clean default curve derived from the par curve
void CredDefSwaption::buildCreditSpreadCurve()const
{
    if (creditSpreadCurve.get()){
        return; // already built
    }
    static const string method = "CredDefSwaption::buildCreditSpreadCurve";
    try
    {
        CreditSpreadCurveSP csCurve = cdsParSpreads->makeCreditSpreadCurve(valueDate, *discount.get());

        // copy it to the transient field
        creditSpreadCurve = csCurve;
        // make backup for restore too
        csRestore = CreditSpreadCurveSP(copy(csCurve.get()));
    
    }
    catch (exception &e) 
    {
        throw ModelException(&e, method);
    }
}

/** Returns name identifying vol for vega parallel */
string CredDefSwaption::sensName(DeltaToCredit* shift) const 
{
    static const string method = "CredDefSwaption::sensName(DeltaToCredit)"; 

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
bool CredDefSwaption::sensShift(DeltaToCredit* shift)
{
    static const string method = "CredDefSwaption::sensShift(DeltaToCredit)"; 

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
CredDefSwaption::CredDefSwaption(): Generic1FactorCredit(TYPE), doEdgworthExpansion(false), 
                                    skewness(0.0), kurtosis(0.0), useSwapRecovery(false),
                                    swapRecovery(0), payAccruedFee(true), addProtectionForPut(false),
                                    priceViaParCurves(true), createE2Csensitivities(false), e2cBasePriceCalculated(false), 
                                    e2cBasePrice(0.0)
{}

class CredDefSwaptionHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CredDefSwaption, clazz);
        SUPERCLASS(Generic1FactorCredit);
        IMPLEMENTS(ClosedFormCDSPS::IIntoProduct);
        IMPLEMENTS(ClosedFormFA::IIntoProduct);
        IMPLEMENTS(ClosedFormCDSPSandFA::IIntoProduct);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(CreditSpreadRhoParallel::RestorableShift);
        IMPLEMENTS(CreditSpreadRhoPointwise::IRestorableShift);
        IMPLEMENTS(CreditSpreadVegaParallel::IShift);
        IMPLEMENTS(DeltaToCredit::IShift);
        IMPLEMENTS(ObjectIteration::IOverride);
        IMPLEMENTS(IGetMarket);
        EMPTY_SHELL_METHOD(defaultCredDefSwaption);

        FIELD(optionIssueDate,           "option issue date");
        FIELD(optionMaturityDate,        "option maturity date");
        FIELD(optionKOOnDefault,         "whether option knocks out when underlying defaults");
        FIELD(optionType,                "(A)merican or (E)uropean");
        FIELD(fwdSwapSpreadVol,          "CDS Spread volatility");
        FIELD(swapStrike,                "CDS Spread volatility");
        FIELD(isCall,                    "true if option is a call");
        FIELD(doEdgworthExpansion,       "true=Edgworth expansion");
        FIELD_MAKE_OPTIONAL(doEdgworthExpansion);
        FIELD(skewness,                  "3rd central moment for Edgworth expansion");
        FIELD_MAKE_OPTIONAL(skewness);
        FIELD(kurtosis,                  "4th central moment for Edgworth expansion");
        FIELD_MAKE_OPTIONAL(kurtosis);
        FIELD(swapEffectiveDate,         "swap effective date");
        FIELD(swapMaturityDate,          "swap maturity date");
        FIELD(frequency,                 "swap frequency");
        FIELD(useSwapRecovery,           "use passed swap recovery value");
        FIELD(swapRecovery,              "defaults to par recovery rate unless useSwapRecovery = true");
        FIELD_MAKE_OPTIONAL(swapRecovery);
        FIELD(payAccruedFee,             "make accrued payment on default?");
        FIELD(addProtectionForPut,       "whether a put includes protection until maturity");
        FIELD_MAKE_OPTIONAL(addProtectionForPut);
        FIELD(dcc,                       "day count conv for accrual periods");
        FIELD_MAKE_OPTIONAL(dcc);
        FIELD(bdc,                       "bad day convention");
        FIELD(creditSpreadCurve,                "credit spread curve");
        FIELD(csRestore,                        "credit spread curve backup");
        FIELD(swpAccrualDCC,                    "swap accrual day count convention");
        FIELD(BDC, "bad day convention");
        FIELD(priceViaParCurves,         "TRUE = price via par curves");
        FIELD(createE2Csensitivities,    "whether to create debt/equity sensitivities");
        FIELD(e2cBasePriceCalculated,    "internal field");
        FIELD(e2cBasePrice,              "internal field");
        FIELD(vol,                       "internal field");

        FIELD_MAKE_TRANSIENT(creditSpreadCurve);
        FIELD_MAKE_TRANSIENT(csRestore);
        FIELD_MAKE_TRANSIENT(swpAccrualDCC);
        FIELD_MAKE_TRANSIENT(BDC);
        FIELD_MAKE_TRANSIENT(priceViaParCurves);
        FIELD_MAKE_TRANSIENT(createE2Csensitivities);
        FIELD_MAKE_TRANSIENT(e2cBasePriceCalculated);
        FIELD_MAKE_TRANSIENT(e2cBasePrice);
        FIELD_MAKE_TRANSIENT(vol);
        // FIELD_MAKE_TWEAKABLE(vol);

        // look up field for use on recurse
        cdsParSpreadsField = clazz->getSuperClass()->getDeclaredField("cdsParSpreads");
        assetField         = clazz->getSuperClass()->getDeclaredField("asset");
    }

    static IObject* defaultCredDefSwaption() {
        return new CredDefSwaption();
    }
};

CClassConstSP const CredDefSwaption::TYPE = CClass::registerClassLoadMethod(
    "CredDefSwaption", typeid(CredDefSwaption), CredDefSwaptionHelper::load);
bool  CredDefSwaptionLoad() {
    return (CredDefSwaption::TYPE != 0);
}




DRLIB_END_NAMESPACE

