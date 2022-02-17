//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GammaSwap.cpp
//
//   Description : Gamma Swap Instrument
//
//   Author      : Manos Venardos
//
//   Date        : 1 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/VegaMatrixLite.hpp"

DRLIB_BEGIN_NAMESPACE

/** Gamma Swap Instrument */
class GammaSwap: public Generic1Factor,
                 virtual public ClosedFormIntegrateLN::IIntoProduct,
                 virtual public LastSensDate,
                 virtual public ISupportVegaMatrixLite {
public:
    friend class GammaSwapProduct;

    static CClassConstSP const TYPE;

    /** Support VEGA_MATRIX_LITE */
    virtual void avoidVegaMatrixLite(const IModel* model);
    
    /** Instrument validation */
    virtual void Validate();

    /** Opportunity to grab samples */
    virtual void GetMarket(const IModel*          model,
                           const CMarketDataSP    market);

    /** Deals with instruments at or after last sample */
    virtual bool priceDeadInstrument(CControl* control,
                                     CResults* results) const;

    /** Indicates whether VEGA_MATRIX is sensible for this instrument */
    virtual bool avoidVegaMatrix(const IModel* model);

    /** returns all strikes on the vol surface to which this instrument is sensitive */
    virtual DoubleArraySP getSensitiveStrikes(OutputNameConstSP outputName,
                                              const IModel*     model);

    /** Returns rolls value date and sets initial spot for Theta,
        return true if sub objects need to be tweaked */
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Instrument into ClosedFormIntegrateLN product */
    ClosedFormIntegrateLN::IProduct* createProduct(ClosedFormIntegrateLN* model) const;

private:
    /** Validation + population of transient fields */
    virtual void validatePop2Object();

    //quadVar Term
    static double quadTerm(double          sumDiv,
                           double          fwd);
        
    /** For reflection */
    static void load(CClassSP& clazz);

    /** Empty shell method */
    static IObject* defaultGammaSwap();

    GammaSwap();

    // Input fields
    double                  strike;                 //!< Strike Vol
    int                     expectedN;              //!< Divisor for realized var
    int                     observationsPerYear;    //!< PPY
    bool                    dividendAdjusted;       //!< Whether to adjust returns for divs
    bool                    scaleByStrike;          //!< Whether to scale by 2*K
    string                  assetHistorySource;     //!< Designates asset history source
    IObservationBuilderSP   observationBuilder;     //!< Factory for observations
    // Transient fields
    DateTimeArraySP         obsDates;
    ObservationTypeArraySP  obsTypes;
    DoubleArraySP           obsSamples;
    ObservationSourceSP     assetHistorySourceObject;
};


void GammaSwap::avoidVegaMatrixLite(const IModel* model) {
    static const string method = "GammaSwap::avoidVegaMatrixLite";
    
    if (avoidVegaMatrix(model)) {
        // Check basic VEGA_MATRIX
        throw ModelException(method, "Instrument does not support VEGA_MATRIX and hence not VEGA_MATRIX_LITE");
    } else if (!SimpleEquity::TYPE->isInstance(asset.get())){
        // Allow LITE only for SimpleEquity
        throw ModelException(method, "Only SimpleEquity underlyings supported");
    }
}


void GammaSwap::Validate() {
    static const string routine = "GammaSwap::Validate";
    try {
        // Override fwdStarting flag on Generic1Factor
        fwdStarting = startDate > valueDate;

        // If forward-starting must start EOD today
        if (fwdStarting) {
            if (startDate.getDate() > valueDate.getDate()) {
                throw ModelException(routine,
                    "Forward starting of at least 1 day is not supported: startDate "
                    + startDate.toString() + " >> valueDate " + valueDate.toString());
            }
        }

        // check that settlement is cash
        if (instSettle->isPhysical() || instSettle->isMargin()) {
            throw ModelException(routine,
                                "Only cash settlement is allowed");
        }

        // Only simple equity is supported
        if (!SimpleEquity::TYPE->isInstance(asset.get())){
            throw ModelException(routine, "Only SimpleEquity (vanilla) is supported");
        }

        // Cross validation of market data
        AssetUtil::assetCrossValidate(asset.get(),
                                      false,
                                      valueDate,
                                      valueDate,
                                      discount,
                                      this);

    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void GammaSwap::GetMarket(const IModel*          model,
                          const CMarketDataSP    market) {
    static const string routine = "GammaSwap::GetMarket";
    try {
        // Call parent
        Generic1Factor::GetMarket(model, market);

        // Grab samples from asset history
        VarSwapUtilities::populateSamples(asset.get(),
            assetHistorySourceObject.get(), valueDate, *obsDates, *obsTypes, *obsSamples);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


bool GammaSwap::priceDeadInstrument(CControl* control,
                                    CResults* results) const {
    static const string routine = "GammaSwap::priceDeadInstrument";
    try {
        DateTime lastDate    = obsDates->back();
        DateTime settlement  = instSettle->settles(lastDate, asset.get());

        if (valueDate < lastDate) {
            // valueDate < lastDate

            // Alive instrument - nothing to do here
            return false;
        } else {
            // lastDate <= valueDate

            // Case 1: lastDate <= settlement <  valueDate
            // Case 2: lastDate <= valueDate  <= settlement

            // Compute output requests for information. Compute price as follows
            // Case 1: price = 0.0
            // Case 2: price = pv of known cashflows

            // Create dummy model on the fly (not really used as payoff is realized)
            ClosedFormIntegrateLNSP model(new ClosedFormIntegrateLN("VolPreferred"));
            // Create product and price it to record output requests and price
            auto_ptr<ClosedFormIntegrateLN::IProduct> product(createProduct(model.get()));
            product->price(model.get(), control, results);

            // Dead instrument
            return true;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


bool GammaSwap::avoidVegaMatrix(const IModel* model) {
    if (ClosedFormIntegrateLN::TYPE->isInstance(model)) {
        return false;
    } else {
        return true;
    }
}


DoubleArraySP GammaSwap::getSensitiveStrikes(OutputNameConstSP outputName,
                                             const IModel*     model) {
    static const string routine("GammaSwap::getSensitiveStrikes");
    try {
        if (avoidVegaMatrix(model)) {
            throw ModelException(routine, "VEGA_MATRIX is not valid for this instrument");
        }

        DoubleArraySP sensStrikes;
        if (ClosedFormIntegrateLN::TYPE->isInstance(model)) {
            ClosedFormIntegrateLN* cfi = dynamic_cast<ClosedFormIntegrateLN*>(const_cast<IModel*>(model));
            sensStrikes = cfi->sensitiveStrikes(outputName, asset.get(), valueDate, obsDates->back());
        } else {
            throw ModelException(routine, "VEGA_MATRIX is only supported for model ClosedFormIntegrateLN");
        }

        return sensStrikes;
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}


bool GammaSwap::sensShift(Theta* shift) {
    static const string routine = "GammaSwap::sensShift";
    try {
        // Fill in cashflows with spot or fwd
        VarSwapUtilities::thetaShiftCashFlows(shift, asset.get(), valueDate, *obsDates, *obsSamples);

        // Call parent
        return Generic1Factor::sensShift(shift);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


DateTime GammaSwap::endDate(const Sensitivity* sensControl) const {
    static const string routine = "GammaSwap::endDate";
    try {
        const DateTime& matDate = obsDates->back();
        const DateTime& instEnd = instSettle->settles(matDate, asset.get());
        const DateTime& assetEnd = asset->settleDate(matDate);
        const DateTime& end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
        return end;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

void GammaSwap::validatePop2Object() {
    static const string routine = "GammaSwap::validatePop2Object";
    try {
        if(Maths::isNegative(strike)) {
            throw ModelException(routine, "Strike is negative " + Format::toString(strike));
        }

        if(expectedN < 1) {
            throw ModelException(routine,
                "Expected N " + Format::toString(expectedN) + " must be positive ");
        }

        if(observationsPerYear < 1) {
            throw ModelException(routine,
                "Observations per year " + Format::toString(observationsPerYear) + " must be strictly positive");
        }

        if(scaleByStrike && Maths::isZero(strike)) {
            throw ModelException(routine,
                "Cannot scale by strike when strike is zero");
        }

        obsDates = observationBuilder->dateList();
        obsTypes = observationBuilder->obsTypes();
        obsSamples = DoubleArraySP(new DoubleArray(obsDates->size()));

        // Validate expectedN consistency with sample schedule
        int noReturns = obsDates->size() - 1;
        if (expectedN != noReturns) {
            throw ModelException(routine, "expectedN (" + Format::toString(expectedN) +
                ") does not match the number of returns in the sample schedule (" +
                Format::toString(noReturns) + ")");
        }

        // Override startDate on Generic1Factor
        startDate = obsDates->front();

        // build the observation source object
        assetHistorySourceObject
            = ObservationSourceSP(new ObservationSource(assetHistorySource));
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


void GammaSwap::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Gamma Swap Instrument");
    REGISTER(GammaSwap, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(ClosedFormIntegrateLN::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ISupportVegaMatrixLite);
    EMPTY_SHELL_METHOD(defaultGammaSwap);
    FIELD(strike, "Volatility strike");
    FIELD(expectedN, "Expected number of returns");
    FIELD(observationsPerYear, "Periods per year");
    FIELD_MAKE_OPTIONAL(observationsPerYear);
    FIELD(dividendAdjusted, "True: adjust returns for divs, False: use whole return");
    FIELD(scaleByStrike, "True: divide by 2*K, False: don't divide by 2*K");
    FIELD_MAKE_OPTIONAL(scaleByStrike);
    FIELD(observationBuilder, "Observation builder");
    FIELD(assetHistorySource, "Asset history source");
    FIELD_MAKE_OPTIONAL(assetHistorySource);
    // Transient fields
    FIELD(obsDates, "Observation dates");
    FIELD_MAKE_TRANSIENT(obsDates);
    FIELD(obsTypes, "Observation types");
    FIELD_MAKE_TRANSIENT(obsTypes);
    FIELD(obsSamples, "Observation samples");
    FIELD_MAKE_TRANSIENT(obsSamples);
    FIELD(assetHistorySourceObject, "Asset history source object");
    FIELD_MAKE_TRANSIENT(assetHistorySourceObject);
}


IObject* GammaSwap::defaultGammaSwap() {
    return new GammaSwap();
}


GammaSwap::GammaSwap(): Generic1Factor(TYPE), observationsPerYear(252),
scaleByStrike(true), assetHistorySource(IMarketObservable::DEFAULT_SOURCE) {}


CClassConstSP const GammaSwap::TYPE = CClass::registerClassLoadMethod(
    "GammaSwap", typeid(GammaSwap), GammaSwap::load);


////////////////////////////////////////////////////////////////////////////


/** ClosedFormIntegrateLN product */
class GammaSwapProduct: virtual public ClosedFormIntegrateLN::IProduct {
public:
    /** Constructor */
    GammaSwapProduct(const GammaSwap* inst);

    /** Price method: computes price and output requests.
        Can be called even for lastDate <= valueDate */
    virtual void price(ClosedFormIntegrateLN* model,
                       Control*               control,
                       CResults*              results);
private:
    const GammaSwap* inst;          //!< Instrument

    /** Prices a portfolio of vanillas with weights 2 / K */
    static double pricePortfolio(VanillaContractsRecorderSP     recorder,
                                 CAssetConstSP                  asset,
                                 const ClosedFormIntegrateLN*   model,
                                 const Control*                 control,
                                 const DateTime&                valueDAte,
                                 const DateTime&                maturity,
                                 double                         firstSample);

    /** Weight for call / put options portfolio */
    class GammaSwapIntegrandWeight: public Function1DDouble {
    public:
        /** Full constructor */
        GammaSwapIntegrandWeight(const Range& integrationDomain);

        /** Implements 1/k weights */
        virtual double operator()(double relativeStrike) const;
    };
    DECLARE_REF_COUNT(GammaSwapIntegrandWeight);
};


GammaSwapProduct::GammaSwapProduct(const GammaSwap* inst): inst(inst) {}


void GammaSwapProduct::price(ClosedFormIntegrateLN* model,
                             Control*               control,
                             CResults*              results) {
    static const string routine = "GammaSwapProduct::price";
    try {
        // Keep references for easier access
        const DateTime& valueDate = inst->valueDate;
        const DateTimeArray& obsDates = *inst->obsDates;
        DoubleArray& obsSamples = *inst->obsSamples;
        const DateTime& firstDate = obsDates.front();
        const DateTime& lastDate = obsDates.back();
        CAssetConstSP asset = inst->asset.getSP();
        double strike = inst->strike;
        int expectedN = inst->expectedN;
        int observationsPerYear = inst->observationsPerYear;
        InstrumentSettlementSP instSettle = inst->instSettle;
        const YieldCurveWrapper& discount = inst->discount;

        // Get option recorder
        VanillaContractsRecorderSP recorder = VanillaContractsRecorder::createVanillaOptionRecorder(control);

        // Define first sample
        double firstSample;
        if(firstDate < valueDate) {
            // Started so use historic fixing
            firstSample = obsSamples[0];
        } else if (firstDate == valueDate) {
            // Spot starting so use spot
            firstSample = asset->getSpot();
        } else {
            // Use forward at start date
            firstSample = asset->fwdValue(firstDate);
        }

        // PART 1: past and current contribution
        double pastFloatingVar = 0.0;
        double currentFloatingVar = 0.0;
        int iStep = 0;
        if(firstDate <= valueDate) {
            // Compute log-returns for past samples
            DoubleArraySP logReturns =  VarSwapUtilities::computeHistoricLogReturns(
                asset.get(), obsDates, obsSamples, valueDate, true, inst->dividendAdjusted, false);

            // Compute past realized variance
            for(iStep = 1; iStep < obsDates.size(); iStep++) {
                if(obsDates[iStep] < valueDate) {
                    double weight = obsSamples[iStep] / firstSample;
                    pastFloatingVar += weight * Maths::square((*logReturns)[iStep]);
                } else {
                    break;
                }
            }

            // Compute current contribution);
			if(valueDate <= lastDate){
				double weight = asset->getSpot() / firstSample;
				currentFloatingVar = weight * Maths::square((*logReturns)[iStep]);
			}
        }

        // PART 2: future contribution. Note may be fwd starting
        double futureFloatingVar = 0.0;
        if(valueDate < lastDate) {
            // PART A: Price terminal portfolios
            recorder->setScalingFactor(1.0);
            double terminalPortfolios = pricePortfolio(recorder, asset, model, control, valueDate, lastDate, firstSample);
            if(valueDate < firstDate) {
                // Forward starting case: put the timeline in reverse order
                recorder->setScalingFactor(-1.0);
                terminalPortfolios -= pricePortfolio(recorder, asset, model, control, valueDate, firstDate, firstSample);
            }

            // PART B: Price term structure of portfolios
            DateTimeArray datesTS;
            DoubleArray   weightsTS;
            
            //methodology used to compute prices with dividends
            string methodology(model->getDivMethodology());
            
            if( !CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX)
                && !CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_DEFAULT) ){
                    throw ModelException("divMethodology chosen for gamma swap is wrong");
            }

            //DIV_CONTINUOUS_APPROX - Manos implementation using logRatios of Fwd
            if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX)){    //None approach can't be handled with general methodology
                model->getPriceTSPortfolios(asset.get(), valueDate, obsDates, datesTS, weightsTS);
            }

            //other cases - Gad implementation using logRatios of ZC
            else{
                model->getPriceTSPortfolios(discount.get(), valueDate, obsDates, datesTS, weightsTS);
            }

            //termStructurePortfolios
            double termStructurePortfolios = 0.0;
            //one first term in 0
            if(! Maths::isZero(weightsTS.size()) ){
                //price for term in 0
                double this0Weight = weightsTS[0];
                recorder->setScalingFactor(-this0Weight);
                double this0Price = pricePortfolio(recorder, asset, model, control, valueDate, datesTS[0], firstSample);
                termStructurePortfolios -= (this0Weight * this0Price);
            }
            //then sum over rolling periods
            for(int iFutStep = 1; iFutStep < datesTS.size(); iFutStep++) {
                //aggregate dvidends
                double sumDiv = 0.0;

                //test for adding dividend in computation
                //nothing for Manos initial implementation 
                if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX) ){
                    //nothing
                }
                //add dividend for Gad implementation
                else{
                    //computation of dividend between steps taken in integral
                    DividendListSP divList = AssetUtil::getDiscreteDivs(asset.get(),
						       			    							valueDate,
                                                                        datesTS[iFutStep-1] <= valueDate ? valueDate : datesTS[iFutStep-1],
                                                                        datesTS[iFutStep] <= valueDate ? valueDate : datesTS[iFutStep],
                                                                        -1, //keep all dividends
										    			    			DividendCollector::DOLLAR_TO_YIELD); //transform dollar dividends into yield dividends
                    //discrete dividend amount
                    DoubleArrayConstSP divListAmount = (divList.get())->getDivAmounts();
                    
                    //number or dividends from valueDate to lastDate 
	    		    int NDiv = (divListAmount.get())->size();

                    for(int iDiv = 0; iDiv < NDiv; iDiv++){
                        sumDiv += (*divListAmount)[iDiv];
                    }
                }

                //modified weight to account for dividends
                double thisWeight = weightsTS[iFutStep] - sumDiv;
                //price
                recorder->setScalingFactor(-thisWeight);
                double thisPrice = pricePortfolio(recorder, asset, model, control, valueDate, datesTS[iFutStep], firstSample);
                
                //price
                termStructurePortfolios -= thisWeight * thisPrice;
            }

            // Total Floating Var before basis
            
            //total return
            futureFloatingVar = terminalPortfolios + termStructurePortfolios;
           
            //Quadratic variation term (flag inst->dividendAdjusted for the instrument)
            if( ! (inst->dividendAdjusted) ){
                //quadVar term
                double quadTermDiv = 0.0;
                //get discrete dividends from today till endDate rolling on obsDates
                for(int iObs = 0; iObs<obsDates.size()-1; iObs++){
                    //test on obsDates being used compare with 
                    DateTime Date1 = obsDates[iObs] <= valueDate ? valueDate : obsDates[iObs];
                    DateTime Date2 = obsDates[iObs+1] <= valueDate ? valueDate : obsDates[iObs+1];
                    
                    //dividend between relevant obsDates (after test)
                    DividendListSP divList = AssetUtil::getDiscreteDivs(asset.get(),
                                                                        valueDate,
                                                                        Date1, //obsDate[iObs]
						    					    			        Date2, //obsDate[iObs+1]
							    					    			    -1, //keep all dividends
								    					    		    DividendCollector::DOLLAR_TO_YIELD); //transform dollar dividends into yield dividends
                                        
                    //discrete dividend amount
                    DoubleArrayConstSP divListAmount = (divList.get())->getDivAmounts();
    		        //number or dividends collected above 
	    	        int NDiv = (divListAmount.get())->size();
    			    //aggregate dividends collected above
                    double div = 0.0; //individual dividend
                    double sumDiv = 0.0; //sum of dividends used in the quadVar function 
                    for(int iDiv = 0; iDiv < NDiv; iDiv++){
                        div = (*divListAmount)[iDiv];
                        sumDiv += log(1.0-div) * log(1.0-div);
                    }
			        
                    //forward value at samplingDate
			        //double fwd = asset->fwdValue(Date1); //t-1 sampling (not used for the moment but could be used if they decide to trade t-1 sampled gamma swaps)
                    double fwd = asset->fwdValue(Date2); //t sampling (default implementation for a gamma swap)
                    
                    //quadVar term
                    quadTermDiv =  GammaSwap::quadTerm(sumDiv, fwd);
        
                    //adding quadVar term (dividion by initial value for homogeneity)
                    futureFloatingVar += quadTermDiv / firstSample;        
                }
            }
        }

        // Compute value of legs
        double floatingLeg =
            ((double)observationsPerYear / (double)expectedN) *
            (pastFloatingVar + currentFloatingVar + futureFloatingVar);
        double fixedLeg = Maths::square(strike);

        // Compute price for non-settled instruments
        double fwdValuedPrice = 100.0 * inst->notional * (floatingLeg - fixedLeg);
        double optionScalingFactor = 100.0 * inst->notional * (double)observationsPerYear / (double)expectedN;

        if(inst->scaleByStrike) {
            fwdValuedPrice /=  2.0 * strike;
            optionScalingFactor /= 2.0 * strike;
        }

        double price = 0.0;
        double pv = 1.0;
        DateTime settlement  = instSettle->settles(lastDate, asset.get());
        if(valueDate <= settlement) {
            // PV and scale
            pv = instSettle->pv(valueDate, lastDate, discount.get(), asset.get());
            price = pv * fwdValuedPrice;
            optionScalingFactor *= pv;
		}
        recorder->scaleNbContracts(optionScalingFactor);

        // Store price
        results->storePrice(price, discount->getCcy());

        //Vega Matrix Lite
        if (control)
        {
            SensitivitySP sens(control->getCurrentSensitivity());
            VegaMatrixLiteSP vml(dynamic_cast<VegaMatrixLite *>(sens.get()));
            if (vml.get()) {
                VanillaInfo::storeVegaMatrix(
                    vml,
                    const_cast<GammaSwap*>(inst),
                    model,
                    recorder,
                    CAssetSP::constCast(asset),
                    valueDate,
                    lastDate,
                    results);

                return;
            }
        }

        // Output Requests now
        if (control && control->isPricing()) {
            // Compute past and future ExpectedN for vol requests
            VolRequestTime volReq;
            IVolProcessedSP vol(asset->getProcessedVol(&volReq));
            HolidayWrapper assetHols = vol->GetTimeMetric()->getHolidays();
            int pastExpectedN;
            int futureExpectedN;
            VarSwapUtilities::pastAndfutureExpN(obsDates, assetHols, valueDate, expectedN, pastExpectedN, futureExpectedN);

            // Compute all vol requests first
            double pastVar = pastFloatingVar + currentFloatingVar;
            double futureVar = futureFloatingVar;

            // TOTAL_VOL
            OutputRequest* request = control->requestsOutput(OutputRequest::TOTAL_VOL);
            if(request) {
                double ppyTime = (double)observationsPerYear / (double)expectedN;
                double totalVar = pastVar + futureVar;
                if(!Maths::isNegative(totalVar)) {
                    double totalVol = sqrt(totalVar * ppyTime);
                    results->storeRequestResult(request, totalVol);
                }
            }

            // Compute future vol
            double volFuture = 0.0;
            if(futureExpectedN && !Maths::isNegative(futureVar)) {
                double ppyFutureTime = (double)(observationsPerYear) / (double)(futureExpectedN);
                volFuture = sqrt(futureVar * ppyFutureTime);
            }

            // VOL_IN_FUTURE
            request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
            if (request) {
                results->storeRequestResult(request, volFuture);
            }

            // IND_VOL = VOL_IN_FUTURE
            request = control->requestsOutput(OutputRequest::IND_VOL);
            if (request) {
                results->storeRequestResult(request, volFuture);
            }

            // VOL_IN_PAST
            request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                double volPast = 0.0;
                if(pastExpectedN && !Maths::isNegative(pastVar)) {
                    double ppyPastTime = (double)(observationsPerYear) / (double)(pastExpectedN);
                    volPast = sqrt(pastVar * ppyPastTime);
                }
                results->storeRequestResult(request, volPast);
            }

            // PAST_WEIGHT
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                double pastWeight = (double)pastExpectedN / (double)(expectedN);
                results->storeRequestResult(request, pastWeight);
            }

            // STRIKE_VOL
            request = control->requestsOutput(OutputRequest::STRIKE_VOL);
            if (request) {
                results->storeRequestResult(request, strike);
            }

            // EXPECTED_N
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, expectedN);
            }

            // DISCOUNT_FACTOR
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request && valueDate <= settlement) {
                results->storeRequestResult(request, pv);
            }

            // PAYMENT_DATES
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray date(1, settlement);
                OutputRequestUtil::recordPaymentDates(control, results, &date);
            }

            // KNOWN_CASH_FLOWS
            request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
            if (request) {
                if (lastDate <= valueDate) {
                    CashFlow cf(settlement, fwdValuedPrice);
                    CashFlowArray cfl(1, cf);
                    OutputRequestUtil::recordKnownCashflows(control,
                        results,
                        discount->getCcy(),
                        &cfl);
                }
            }

            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                results,
                price,
                valueDate,
                discount.get(),
                asset.get(),
                inst->premiumSettle.get());

            // FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                lastDate,
                valueDate,
                asset.get());
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


double GammaSwapProduct::pricePortfolio(VanillaContractsRecorderSP      recorder,
                                        CAssetConstSP                   asset,
                                        const ClosedFormIntegrateLN*    model,
                                        const Control*                  control,
                                        const DateTime&                 valueDate,
                                        const DateTime&                 maturity,
                                        double                          firstSample) {
    static const string routine = "GammaSwapProduct::pricePortfolio";
    try {
        // I) Get integrator from model
        Range integrationDomain(OpenBoundary(0.0), Infinity(Infinity::Plus));
        Integrator1DSP integrator;
        model->getIntegrator(control, asset.get(), valueDate, maturity, recorder, integrator, integrationDomain);

        // II) Modify scaling factor
        double fwd = asset->fwdValue(maturity);
        double scalingFactor = recorder->getScalingFactor() * 2.0 * fwd / firstSample;
        recorder->setScalingFactor(scalingFactor);
        
        // III) Construct integrand and integrate
        GammaSwapIntegrandWeightConstSP optionWeights(new GammaSwapIntegrandWeight(integrationDomain));
        StaticReplicationIntegrandSP integrand(new StaticReplicationIntegrand(
            optionWeights, recorder, asset, valueDate, maturity, model->negativeFwdVarAllowed()));
        double integral = integrator->integrate(*integrand);
        double price = 2.0 * integral * fwd / firstSample;
        return price;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


GammaSwapProduct::GammaSwapIntegrandWeight::GammaSwapIntegrandWeight(const Range& integrationDomain):
Function1DDouble(integrationDomain) {}


double GammaSwapProduct::GammaSwapIntegrandWeight::operator()(double relativeStrike) const {
    double relativeStrikeWeight = 1.0 / relativeStrike;
    return relativeStrikeWeight;
}


////////////////////////////////////////////////////////////////////////////


ClosedFormIntegrateLN::IProduct* GammaSwap::createProduct(ClosedFormIntegrateLN* model) const {
	return new GammaSwapProduct(this);
}


bool GammaSwapLoad() {
    return (GammaSwap::TYPE != 0);
}


//quadVar Term
double GammaSwap::quadTerm(double          sumDiv,
                           double          fwd){
    static const string method("GammaSwap::quadTerm");
    try{
        //result
        double res;

        //forward component
        res = sumDiv * fwd;
        
        //result
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }                
}

DRLIB_END_NAMESPACE
