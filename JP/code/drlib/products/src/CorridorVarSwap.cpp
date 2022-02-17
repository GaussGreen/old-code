//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorridorVarSwap.cpp
//
//   Description : Corridor Variance Swap Instrument
//
//   Author      : Arnaud Jobert
//
//   Date        : 14 February 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/VolProcessed.hpp"
#include "edginc/Black.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/Array.hpp"
#include "edginc/VarSwapBasis.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/VegaMatrixLite.hpp"


DRLIB_BEGIN_NAMESPACE

/** Corridor Variance Swap Instrument */
class CorridorVarSwap: public Generic1Factor,
                       virtual public ISupportVegaMatrixLite,
                       virtual public ClosedFormIntegrateLN::IIntoProduct,
                       virtual public LegalTerms::Shift,
                       virtual public LastSensDate {
public:
    friend class CorridorVarSwapProduct;

    static CClassConstSP const TYPE;

    /** Support VEGA_MATRIX_LITE */
    virtual void avoidVegaMatrixLite(const IModel* model);
    
    /** Instrument validation */
    virtual void Validate();

    virtual bool sensShift(LegalTerms* shift);

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

    
    //checkRangeAtT used for sampling date
    static const string SAMPLING_T_MINUS_1;
    static const string SAMPLING_T;
	static const string SAMPLING_DEFAULT;
	
private:
    /** Validation + population of transient fields */
    virtual void validatePop2Object();

    //Newton Cotes method for numerical integration
    static double NewtonCotesCallIntegral(int                  n,
		                				  double               LMin,
							              double               LMax,
							              const DateTime&      divDate,
                                          const DateTime&      valueDate,
   							              double               fwd,
							              const CAsset*        asset,
                                          bool                 allowNegativeVar,
                                          VanillaContractsRecorderSP recorder);
 
    //brute force integration
    static double bruteForceIntegral(int                  n,
                                     double               LMin,
                                     double               LMax,
                                     const DateTime&      divDate,
                                     const DateTime&      valueDate,
                                     double               fwd,
                                     const CAsset*        asset,
                                     bool                 allowNegativeVar,
                                     VanillaContractsRecorderSP recorder);

    //quadVar Term
    static double quadTerm(double          L,
                           double          U,
                           const DateTime& divDate,
                           const DateTime& valueDate,
                           double          fwd,
                           const CAsset*   asset,
                           bool            allowNegativeVar,
                           VanillaContractsRecorderSP recorder);
    
    /** For reflection */
    static void load(CClassSP& clazz);

    /** Empty shell method */
    static IObject* defaultCorridorVarSwap();

    /** Default constructor */
    CorridorVarSwap();

    // Input fields

    double                  strike;                //!< Strike Vol
    int                     expectedN;             //!< Divisor for realized var
    int                     observationsPerYear;   //!< PPY
    bool                    dividendAdjusted;      //!< Whether to adjust returns for divs
    bool                    scaleByStrike;         //!< Whether to scale by 2*K
    double                  strikeRef;             //!< vega notional var swap scaling
    double                  LowerBarrier;          //!< Lower barrier
    double                  LowerEconBarrier;      //!< lower barrier for legal terms
    double                  UpperBarrier;          //!< Upper barrier
    double                  UpperEconBarrier;      //!< upper barrier for legal terms

    string                  assetHistorySource;     //!< Designates asset history source
    IObservationBuilderSP   observationBuilder;     //!< Factory for observations

    string	                checkRangeAtT;	        //!< Sampling used to weigh log-retruns
    bool                    percentageBarriers;     //!< True: barriers are % of first sample, False: absolute

    bool                    isCapped;
    double                  cap;
    double                  ExpRangeNForCap;
    
    // Transient fields
    DateTimeArraySP         obsDates;
    ObservationTypeArraySP  obsTypes;
    DoubleArraySP           obsSamples;
    ObservationSourceSP     assetHistorySourceObject; //!< object for asset history source
};


void CorridorVarSwap::avoidVegaMatrixLite(const IModel* model) {
    static const string method = "CorridorVarSwap::avoidVegaMatrixLite";
    
    if (avoidVegaMatrix(model)) {
        // Check basic VEGA_MATRIX
        throw ModelException(method, "Instrument does not support VEGA_MATRIX and hence not VEGA_MATRIX_LITE");
    } else if (!SimpleEquity::TYPE->isInstance(asset.get())) {
        // Allow LITE only for SimpleEquity
        throw ModelException(method, "Only SimpleEquity underlyings supported");
    }
}


void CorridorVarSwap::Validate() {
    static const string routine = "CorridorVarSwap::Validate";
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



// set barriers to be the economic (legal) ones
bool CorridorVarSwap::sensShift(LegalTerms* shift) {
        // just replace all 2 barriers with the corresponding economic one
        LowerBarrier = LowerEconBarrier;
        UpperBarrier = UpperEconBarrier;

        return true; // continue shifting
}


void CorridorVarSwap::GetMarket(const IModel*          model,
                          const CMarketDataSP    market) {
    static const string routine = "CorridorVarSwap::GetMarket";
    try {
        // Call parent
        Generic1Factor::GetMarket(model, market);

        // Grab samples from asset history
        VarSwapUtilities::populateSamples(asset.get(),
            assetHistorySourceObject.get(), valueDate,
            *obsDates, *obsTypes, *obsSamples);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


bool CorridorVarSwap::priceDeadInstrument(CControl* control,
                                          CResults* results) const {
    static const string routine = "CorridorVarSwap::priceDeadInstrument";
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


bool CorridorVarSwap::avoidVegaMatrix(const IModel* model) {
    if (ClosedFormIntegrateLN::TYPE->isInstance(model)) {
        return false;
    } else {
        return true;
    }
}


DoubleArraySP CorridorVarSwap::getSensitiveStrikes(OutputNameConstSP outputName,
                                                   const IModel*     model) {
    static const string routine("CorridorVarSwap::getSensitiveStrikes");
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


bool CorridorVarSwap::sensShift(Theta* shift) {
    static const string routine = "CorridorVarSwap::sensShift";
    try {
        // Fill in cashflows with spot or fwd
        VarSwapUtilities::thetaShiftCashFlows(shift, asset.get(), valueDate, *obsDates, *obsSamples);

        // Call parent
        return Generic1Factor::sensShift(shift);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


DateTime CorridorVarSwap::endDate(const Sensitivity* sensControl) const {
    static const string routine = "CorridorVarSwap::endDate";
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


void CorridorVarSwap::validatePop2Object() {
    static const string routine = "CorridorVarSwap::validatePop2Object";
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

        if(scaleByStrike && Maths::isZero(strikeRef)) {
            throw ModelException(routine,
                "Cannot scale by strike when strike Ref is zero");
        }

        // Validate against zero lower barrier
        if( Maths::isZero(LowerBarrier) ||
            Maths::isZero(LowerEconBarrier) ||
            Maths::isZero(UpperBarrier) ||
            Maths::isZero(UpperEconBarrier)) {
            throw ModelException(routine, "All input barriers must be greater than zero");
        }
        
        //Validate checkRangeAtT used to weigh log-returns
        if (   !CString::equalsIgnoreCase(checkRangeAtT, SAMPLING_DEFAULT)
            && !CString::equalsIgnoreCase(checkRangeAtT, SAMPLING_T_MINUS_1)
            && !CString::equalsIgnoreCase(checkRangeAtT, SAMPLING_T) ){
            throw ModelException("checkRangeAtT chosen is wrong");
        }

        // Lower different from Upper
        if(Maths::equals(LowerBarrier, UpperBarrier) ||
           Maths::equals(LowerEconBarrier, UpperEconBarrier)) {
            throw ModelException(routine, "Risk and Legal lower barriers must be differnet from upper barriers");
        }

        // Lower < Upper
        if(LowerBarrier > UpperBarrier || LowerEconBarrier > UpperEconBarrier) {
            throw ModelException(routine, "Risk and Legal lower barriers must be greater than upper barriers");
        }

        // Build dates & observations
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

        // Validate ExpRangeNForCap for the capped corridor: set to an arbitrary number between 0 and 1
        if (ExpRangeNForCap < 0.0 || ExpRangeNForCap > 1.0) {
            throw ModelException(routine, "ExpRangeNForCap must be set between 0 and 1");
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


void CorridorVarSwap::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Corridor Variance Swap Instrument");
    REGISTER(CorridorVarSwap, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(ISupportVegaMatrixLite);
    IMPLEMENTS(ClosedFormIntegrateLN::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(LegalTerms::Shift);
    EMPTY_SHELL_METHOD(defaultCorridorVarSwap);
    FIELD(strike, "Volatility strike");
    FIELD(expectedN, "Expected number of returns");
    FIELD(observationsPerYear, "Periods per year");
    FIELD_MAKE_OPTIONAL(observationsPerYear);
    FIELD(dividendAdjusted, "True: adjust returns for divs, False: use whole return");
    FIELD(scaleByStrike, "True: divide by 2*K, False: don't divide by 2*K");
    FIELD(strikeRef, "Reference strike");
    FIELD_MAKE_OPTIONAL(scaleByStrike);
    FIELD(LowerBarrier,"CorridorVarSwap lower barrier");
    FIELD(LowerEconBarrier,"CorridorVarSwap legal lower barrier");
    FIELD(UpperBarrier,"CorridorVarSwap upper barrier");
    FIELD(UpperEconBarrier,"CorridorVarSwap legal lower barrier");
    FIELD(observationBuilder, "Observation builder");
    FIELD(assetHistorySource, "Asset history source");
    FIELD_MAKE_OPTIONAL(assetHistorySource);
    FIELD(checkRangeAtT, "Sampling date for log-returns");
    FIELD_MAKE_OPTIONAL(checkRangeAtT);
    FIELD(percentageBarriers, "TRUE: barriers are a % of first sample. FALSE: absolute levels");
    FIELD_MAKE_OPTIONAL(percentageBarriers);
    FIELD(isCapped, "true: corridor payoff is capped, false: corridor payoff is not capped");
    FIELD_MAKE_OPTIONAL(isCapped);
    FIELD(cap, "Cap for corridor variance"); // now interpreted as a volatility number
    FIELD_MAKE_OPTIONAL(cap);
    FIELD(ExpRangeNForCap, "Expected number of returns for the cap");
    FIELD_MAKE_OPTIONAL(ExpRangeNForCap);

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


IObject* CorridorVarSwap::defaultCorridorVarSwap() {
    return new CorridorVarSwap();
}


CorridorVarSwap::CorridorVarSwap(): 
Generic1Factor(TYPE), observationsPerYear(252), scaleByStrike(true), 
assetHistorySource(IMarketObservable::DEFAULT_SOURCE), checkRangeAtT(CorridorVarSwap::SAMPLING_DEFAULT),
percentageBarriers(false),isCapped(false),cap(2.5),ExpRangeNForCap(1.0) {}


CClassConstSP const CorridorVarSwap::TYPE = CClass::registerClassLoadMethod(
    "CorridorVarSwap", typeid(CorridorVarSwap), CorridorVarSwap::load);


////////////////////////////////////////////////////////////////////////////


/** ClosedFormIntegrateLN product */
class CorridorVarSwapProduct: virtual public ClosedFormIntegrateLN::IProduct {
public:
    static const double PUT_SPREAD_WIDTH;       //!< Put spread width
    
    /** Constructor */
    CorridorVarSwapProduct(const CorridorVarSwap* inst);

    /** Price method: computes price and output requests.
        Can be called even for lastDate <= valueDate */
    virtual void price(ClosedFormIntegrateLN* model,
                       Control*               control,
                       CResults*              results);
private:
    const CorridorVarSwap* inst;                //!< Instrument

    /** Prices a portfolio of vanillas with weights 2 / K */
    static double pricePortfolio(VanillaContractsRecorderSP     recorder,
                                 CAssetConstSP                  asset,
                                 const ClosedFormIntegrateLN*   model,
                                 const Control*                 control,
                                 const DateTime&                valueDAte,
                                 const DateTime&                maturity,
                                 double                         lowerBarr,
                                 double                         upperBarr,
                                 bool                           allowNegativeVar);

    /** Weight for call / put options portfolio */
    class CorridorVarSwapIntegrandWeight: public Function1DDouble {
    public:
        /** Constructor */
        CorridorVarSwapIntegrandWeight(const Range& integrationDomain);

        /** Implements 1/k^2 weights */
        virtual double operator()(double relativeStrike) const;
    };
    DECLARE_REF_COUNT(CorridorVarSwapIntegrandWeight);
};

const double CorridorVarSwapProduct::PUT_SPREAD_WIDTH = 0.0001;


CorridorVarSwapProduct::CorridorVarSwapProduct(const CorridorVarSwap* inst): inst(inst) {}


void CorridorVarSwapProduct::price(ClosedFormIntegrateLN* model,
                                   Control*               control,
                                   CResults*              results) {
    static const string routine = "CorridorVarSwapProduct::price";
    try {
        // Keep references for easier access
        const DateTime& valueDate = inst->valueDate;
        const DateTimeArray& obsDates = *inst->obsDates;
        DoubleArray& obsSamples = *inst->obsSamples;
        const DateTime& firstDate = obsDates.front();
        const DateTime& lastDate = obsDates.back();
        CAssetConstSP asset = inst->asset.getSP();
        double strike = inst->strike;
        double strikeRef = inst->strikeRef;
        int expectedN = inst->expectedN;
        int observationsPerYear = inst->observationsPerYear;
        InstrumentSettlementSP instSettle = inst->instSettle;
        const YieldCurveWrapper& discount = inst->discount;
        double cap = inst->cap;
        double ExpRangeNForCap = inst->ExpRangeNForCap;
        

        // Predefine all barriers
        double instLowerBarrier = inst->LowerBarrier;
        double instUpperBarrier = inst->UpperBarrier;
        double instLowerEconBarrier = inst->LowerEconBarrier;
        double instUpperEconBarrier = inst->UpperEconBarrier;
        if(inst->percentageBarriers) {
            // Convert from % to absolute
            double scalingFactor = 0.0;
            if(firstDate <= valueDate) {
                // For started swaps scale by first sample
                scalingFactor = obsSamples[0];
            } else {
                // For forward staring scale by current spot (remember only allowed up to 1 day)
                scalingFactor = asset->getSpot();
            }
            instLowerBarrier *= scalingFactor;
            instUpperBarrier *= scalingFactor;
            instLowerEconBarrier *= scalingFactor;
            instUpperEconBarrier *= scalingFactor;
        }

        // Get option recorder
        VanillaContractsRecorderSP recorder = VanillaContractsRecorder::createVanillaOptionRecorder(control);
        double nbOptionsForRange = - Maths::square(strike) / (double)expectedN;
        double nbOptionsForVar = double(observationsPerYear) / (double)expectedN;

        
        //checkRangeAtT for sampling date
        string checkRangeAtT = inst->checkRangeAtT;

        // Compute past and future ExpectedN
        VolRequestTime volReq;
        IVolProcessedSP vol(asset->getProcessedVol(&volReq));
        HolidayWrapper assetHols = vol->GetTimeMetric()->getHolidays();
        int pastExpectedN;
        int futureExpectedN;
        VarSwapUtilities::pastAndfutureExpN(obsDates, assetHols, valueDate, expectedN, pastExpectedN, futureExpectedN);

        // PART 1: past contribution including partly past / partly future contribution
        double pastFloatingVar = 0.0;
        double futureFloatingVar = 0.0;
        double pastRange = 0.0;   // counting number of times historic spot was in range
        double futureRange = 0.0; // expected number of times spot will leave the range specified by barriers (using call/put spreads)
        int iStep = 1;
        if(firstDate <= valueDate) {
            // Compute log-returns for past samples
            DoubleArraySP logReturns = VarSwapUtilities::computeHistoricLogReturns(
                asset.get(), obsDates, obsSamples, valueDate, true, inst->dividendAdjusted, false);

            // Compute past realized variance
            double lowerBarr = instLowerEconBarrier; // cannot bend in the past
            double upperBarr = instUpperEconBarrier; // cannot bend in the past
            
            // SAMPLING DEFINITION: L < S(t-1) <= U
            if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_DEFAULT) 
                || CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T_MINUS_1) ){
                for(iStep = 0; iStep < obsDates.size() - 1; iStep++) {
                    const DateTime& thisStartDate = obsDates[iStep];
                    const DateTime& thisEndDate   = obsDates[iStep + 1];
                    if(thisStartDate < valueDate) {
                        // Started in the past
                        double sample = obsSamples[iStep];
                    bool inCorridor =
                            Maths::isPositive(sample - lowerBarr) && !Maths::isNegative(upperBarr - sample);
                        double weight = inCorridor? 1.0 : 0.0;
                        pastRange += weight;
                        double thisContribution = weight * Maths::square((*logReturns)[iStep + 1]);
                        pastFloatingVar += thisContribution;
                    }
                    else if(thisStartDate == valueDate) {
                        // Starting now and ending in the future
                        // Account as future range due to lag
                        double sample = asset->getSpot();
                        bool inCorridor = 
                            Maths::isPositive(sample - lowerBarr) && !Maths::isNegative(upperBarr - sample);
                        double weight = inCorridor? 1.0 : 0.0;
                        futureRange += weight;
                    }
                    else {
                        // Strictly future
                        break;
                    }
                }
            }
            
            // SAMPLING DEFINITION: L < S(t) <= U
            if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T) ){
                for(iStep = 0; iStep < obsDates.size() - 1; iStep++) {
                    const DateTime& thisStartDate = obsDates[iStep];
                    const DateTime& thisEndDate   = obsDates[iStep + 1];
                    if(thisEndDate < valueDate) {
                        // ending in the past
                        double sample = obsSamples[iStep+1];
                        bool inCorridor = 
                            Maths::isPositive(sample - lowerBarr) && !Maths::isNegative(upperBarr - sample);
                        double weight = inCorridor? 1.0 : 0.0;
                        pastRange += weight;
                        double thisContribution = weight * Maths::square((*logReturns)[iStep + 1]);
                        pastFloatingVar += thisContribution;
                    }
                    else if(thisEndDate == valueDate) {
                        // Starting in the past and ending today
                        // Account as past range
                        double sample = asset->getSpot();
                        bool inCorridor = 
                            Maths::isPositive(sample - lowerBarr) && !Maths::isNegative(upperBarr - sample);
                        double weight = inCorridor? 1.0 : 0.0;
                        pastRange += weight;
                        double thisContribution = weight * Maths::square((*logReturns)[iStep + 1]);
                        pastFloatingVar += thisContribution;
                    }
                    else if(thisEndDate > valueDate && thisStartDate <= valueDate ) {
                        // Starting now or in the past and ending strictly in the future
                        // Accounts as future, future range is acomputed using probability / contribution is computed using a log return
                        double lowerBarr = instLowerBarrier;
                        double upperBarr = instUpperBarrier;
                        double delta = CorridorVarSwapProduct::PUT_SPREAD_WIDTH;
                        DoubleArray strikes(4);
                        DoubleArray weights(4);
                        strikes[0] = lowerBarr * (1 - delta);
                        strikes[1] = lowerBarr * (1 + delta);
                        strikes[2] = upperBarr * (1 - delta);
                        strikes[3] = upperBarr * (1 + delta);
                        weights[0] = 1.0 / (2.0 * delta * lowerBarr);
                        weights[1] = 1.0 / (-2.0 * delta * lowerBarr);
                        weights[2] = 1.0 / (-2.0 * delta * upperBarr);
                        weights[3] = 1.0 / (2.0 * delta * upperBarr);
                        //vol request
                        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, lastDate, false));
                        volRequest->allowNegativeFwdVar(model->negativeFwdVarAllowed());
                        //fwd
                        double fwd = asset->fwdValue(thisEndDate);
                        // Compute probabilities
                        // Price difference of put spreads to approximate probability of being
                        // in corridor for a given date
                        const DateTime& thisDate = thisEndDate;
                        CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                        double thisYearFrac = volBS->calcTradingTime(valueDate, thisDate);
                        recorder->setScalingFactor(nbOptionsForRange);
                        double thisProb = 0.0;
                        for(int iStrike = 0; iStrike < 4; iStrike++) {
                            volRequest->setStrike(strikes[iStrike]);
                            volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
                            double vol = volBS->CalcVol(valueDate, thisDate);
                            double var = Maths::square(vol) * thisYearFrac;
                            thisProb +=  BlackPrice(false, fwd, strikes[iStrike], 1.0, vol, thisYearFrac, 
                                weights[iStrike], thisDate, recorder);
                        }
                        futureRange += thisProb;
                        
                        //contribution
                        futureFloatingVar += thisProb * Maths::square( log(obsSamples[iStep] / asset->getSpot()) );
                    }
                    else{
                        // Strictly future
                        break;
                    }
                }
            }
        }

        // PART 2: future contribution. Note may be fwd starting

        // PART 2A: Future Digitals
        // preliminaries needed to compute additional fixed leg component due to call/put spreads, i.e. futureRange:
        double lowerBarr = instLowerBarrier;
        double upperBarr = instUpperBarrier;
        double delta = CorridorVarSwapProduct::PUT_SPREAD_WIDTH;
        DoubleArray strikes(4);
        DoubleArray weights(4);
        strikes[0] = lowerBarr * (1 - delta);
        strikes[1] = lowerBarr * (1 + delta);
        strikes[2] = upperBarr * (1 - delta);
        strikes[3] = upperBarr * (1 + delta);
        weights[0] = 1.0 / (2.0 * delta * lowerBarr);
        weights[1] = 1.0 / (-2.0 * delta * lowerBarr);
        weights[2] = 1.0 / (-2.0 * delta * upperBarr);
        weights[3] = 1.0 / (2.0 * delta * upperBarr);

        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, lastDate, false));
        volRequest->allowNegativeFwdVar(model->negativeFwdVarAllowed());

        // build up daily maturities and forwards
        DateTimeArray remainingDates;
        // SAMPLING DEFINITION: L < S(t-1) <= U
        if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_DEFAULT) 
            || CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T_MINUS_1) ){
            int k = iStep;
            for (; k < obsDates.size(); k++) {
                // Only count if the return is strictly in the future
                // If it's partly past partly future it has been counted already in the
                // pastRange
                if(obsDates[k - 1] > valueDate) {
                    remainingDates.push_back(obsDates[k]);
                }
            }
            //fwds
            DoubleArray remainingFwds(remainingDates.size());
            asset->fwdValue(remainingDates, remainingFwds);
            // Compute probabilities
            for (k = 0; k < remainingDates.size(); k++) {
                // Price difference of put spreads to approximate probability of being
                // in corridor for a given date
                const DateTime& thisDate = remainingDates[k];
                CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                double thisYearFrac = volBS->calcTradingTime(valueDate, thisDate);
                double fwd = remainingFwds[k];
                recorder->setScalingFactor(nbOptionsForRange);
                for(int iStrike = 0; iStrike < 4; iStrike++) {
                    volRequest->setStrike(strikes[iStrike]);
                    volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
                    double vol = volBS->CalcVol(valueDate, thisDate);
                    double var = Maths::square(vol) * thisYearFrac;
                    double thisPrice = BlackPrice(false, fwd, strikes[iStrike], 1.0, vol, thisYearFrac, 
                         weights[iStrike], thisDate, recorder);
                    futureRange += thisPrice;
                }
            }
        }
        // SAMPLING DEFINITION: L < S(t) <= U
        if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T) ){
            int k = iStep;
            for (; k < obsDates.size(); k++) {
                // Only count if the return is strictly in the future
                // If it's partly past partly future it has been counted already in the
                // pastRange
                if(obsDates[k-1] > valueDate) {
                    remainingDates.push_back(obsDates[k]);
                }
            }
            //fwds
            DoubleArray remainingFwds(remainingDates.size());
            asset->fwdValue(remainingDates, remainingFwds);
            // Compute probabilities
            for (k = 0; k < remainingDates.size(); k++) {
                // Price difference of put spreads to approximate probability of being
                // in corridor for a given date
                const DateTime& thisDate = remainingDates[k];
                CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                double thisYearFrac = volBS->calcTradingTime(valueDate, thisDate);
                double fwd = remainingFwds[k];
                recorder->setScalingFactor(nbOptionsForRange);
                for(int iStrike = 0; iStrike < 4; iStrike++) {
                    volRequest->setStrike(strikes[iStrike]);
                    volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
                    double vol = volBS->CalcVol(valueDate, thisDate);
                    double var = Maths::square(vol) * thisYearFrac;
                    double thisPrice = BlackPrice(false, fwd, strikes[iStrike], 1.0, vol, thisYearFrac, 
                         weights[iStrike], thisDate, recorder);
                    futureRange += thisPrice;
                }
            }
        }

        // PART 2B: Future Corridor variance
        VarSwapBasis::VarSwapBasisProcSP basisProc;
        if(valueDate < lastDate) {
            //negativeVariance allowed or not for volRequest
            bool allowNegativeVariance = model->negativeFwdVarAllowed();

            // PART I: Price terminal portfolios
            recorder->setScalingFactor(nbOptionsForVar);
            double terminalPortfolios = pricePortfolio(recorder, asset, model, control, 
                valueDate, lastDate, lowerBarr, upperBarr, allowNegativeVariance);
            if(valueDate < firstDate) {
                // Forward starting case: put the timeline in reverse order
                recorder->setScalingFactor(-nbOptionsForVar);
                terminalPortfolios -= pricePortfolio(recorder, asset, model, control, 
                    valueDate, firstDate, lowerBarr, upperBarr, allowNegativeVariance);
            }

            // II: Price term structure of portfolios
            DateTimeArray datesTS;
            DoubleArray   weightsTS;

            //methodology used to compute prices with dividends
            string methodology(model->getDivMethodology());
            
            //DIV_CONTINUOUS_APPROX - Manos initial implementation using logRatios of Fwd
            if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX) ){
                model->getPriceTSPortfolios(asset.get(), valueDate, obsDates, datesTS, weightsTS);
            }
            //other cases - Gad implementation using logRatios of ZC
            else{
                model->getPriceTSPortfolios(discount.get(), valueDate, obsDates, datesTS, weightsTS);
            }

            double termStructurePortfolios = 0.0;
            // SAMPLING DEFINITION SAMPLING_T_MINUS_1: L < S(t-1) <= U
            if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_DEFAULT) 
                || CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T_MINUS_1) ){
                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, lastDate, false));
                volRequest->allowNegativeFwdVar(allowNegativeVariance);
                DoubleArray fwds(datesTS.size());
                asset->fwdValue(datesTS, fwds);
                recorder->setScalingFactor(nbOptionsForVar);
                for(int iFutStep = 0; iFutStep < datesTS.size(); iFutStep++) {
                    double thisWeight = weightsTS[iFutStep];

                    
                    const DateTime& thisDate = datesTS[iFutStep];
                    double fwd = fwds[iFutStep];
                    
                    // Price option at Lower Barrier
                    volRequest->setStrike(lowerBarr);
                    CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                    double thisYearFrac = volBS->calcTradingTime(valueDate, thisDate);
                    double volLow = volBS->CalcVol(valueDate, thisDate);
                    double varLow = Maths::square(volLow) * thisYearFrac;
                    bool isCall = lowerBarr < fwd ? false : true;
                    
                    double lowContribution = BlackPrice(isCall, fwd, lowerBarr, 1.0, volLow, thisYearFrac,
                        -2.0 * thisWeight / lowerBarr, thisDate, recorder);

                    // Price option at Upper Barrier
                    volRequest->setStrike(upperBarr);
                    volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
                    double volHigh = volBS->CalcVol(valueDate, thisDate);
                    double varHigh = Maths::square(volHigh) * thisYearFrac;
                    isCall = upperBarr < fwd ? false : true;
                    
                    double highContribution = BlackPrice(isCall, fwd, upperBarr, 1.0, volHigh, thisYearFrac,
                        2.0 * thisWeight / upperBarr, thisDate, recorder);
                    

                    // Add discontinuity contributions
                    termStructurePortfolios += lowContribution + highContribution;
                }
            }
            // SAMPLING DEFINITION: L < S(t) <= U
            if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T) ){
                LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, lastDate, false));
                volRequest->allowNegativeFwdVar(allowNegativeVariance);
                DoubleArray fwds(datesTS.size());
                asset->fwdValue(datesTS, fwds);
                recorder->setScalingFactor(nbOptionsForVar);
                for(int iFutStep = 0; iFutStep < datesTS.size(); iFutStep++) {
                    double thisWeight = weightsTS[iFutStep];
                    
                    const DateTime& thisDate = datesTS[iFutStep];
                    double fwd = fwds[iFutStep];
                    
                    // Price option at Lower Barrier
                    volRequest->setStrike(lowerBarr);
                    CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
                    double thisYearFrac = volBS->calcTradingTime(valueDate, thisDate);
                    double volLow = volBS->CalcVol(valueDate, thisDate);
                    double varLow = Maths::square(volLow) * thisYearFrac;
                    bool isCall = lowerBarr < fwd ? false : true;

                    double lowContribution = BlackPrice(isCall, fwd, lowerBarr, 1.0, volLow, thisYearFrac,
                        -2.0 * thisWeight / lowerBarr, thisDate, recorder);

                    // Price option at Upper Barrier
                    volRequest->setStrike(upperBarr);
                    volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
                    double volHigh = volBS->CalcVol(valueDate, thisDate);
                    double varHigh = Maths::square(volHigh) * thisYearFrac;
                    isCall = upperBarr < fwd ? false : true;
                    
                    double highContribution = BlackPrice(isCall, fwd, upperBarr, 1.0, volHigh, thisYearFrac,
                        2.0 * thisWeight / upperBarr, thisDate, recorder);
                    
                    // Add discontinuity contributions
                    termStructurePortfolios += lowContribution + highContribution;
                }
            }
            
			// II': Price term structure of portfolios in the presence of dividends
            
			//contribution of options maturing at dividend dates using shift in barriers
			double divStripesOptionValues = 0.0;
            //DIV_CONTINUOUS_APPROX: we just add the quad Var term corresponding to dividends
            if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX) ){
                //Nothing to do, the div effect has been accounted in using fwds instead of pv in portfolioTS
                //(Manos implementation vs Gad implementation above)
            }
            //other cases: we need to add the dividend effect due to the shift in barriers in the repplication formula
            else{
                //get discrete dividends from today till endDate
			    DividendListSP divList = AssetUtil::getDiscreteDivs(asset.get(),
									    							valueDate,                        
                                                                    firstDate <= valueDate ? valueDate : firstDate,
											    					lastDate,
												    				-1, //keep all dividends
													    			DividendCollector::DOLLAR_TO_YIELD); //transform dollar dividends into yield dividends
                //discrete dividend dates
                DateTimeArrayConstSP divListDates = (divList.get())->getExDivDates(); 
			    //discrete dividend amount
                DoubleArrayConstSP divListAmount = (divList.get())->getDivAmounts();
			
    			//number or dividends from valueDate to lastDate 
	    		int NDiv = (divListDates.get())->size();
			
		    	//transform divListDates in DateTimeArray
                DateTimeArray divListDatesArray(NDiv);//ant div
                for(int j = 0; j< NDiv; j++){
						divListDatesArray[j] = DateTime((*divListDates)[j].getDate(), DateTime::BEFORE_EX_DIV_TIME);
                }
            
			    //forward value at dividend dates (ant Div / post div)
			    DoubleArray fwdsAtDiv(NDiv);
                asset->fwdValue(divListDatesArray, fwdsAtDiv);
                
			    //definition of points used for the integrals			
			    double lowerBarr = instLowerBarrier;
			    double upperBarr = instUpperBarrier;

                //order of integral approximation
			    int intOrder;

                //DIV_BRUTE_FORCE
                if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_BRUTE_FORCE) ){
                    //order
                    intOrder = 1000;

                    //rolling over dividends
			        for(int iDiv = 0; iDiv < NDiv; iDiv++){
				        double d = (*divListAmount)[iDiv];
				        const DateTime& dDate = divListDatesArray[iDiv];
                        double fwd = fwdsAtDiv[iDiv];
                        //stripes for upperBarrier using bruteForce integration
                        recorder->setScalingFactor(-2.0 * nbOptionsForVar);
                        divStripesOptionValues -= 2.0 * CorridorVarSwap::bruteForceIntegral(intOrder,
                                                                                            upperBarr,
                                                                                            upperBarr / (1.0 - d),
                                                                                            dDate,
                                                                                            valueDate,
                                                                                            fwd,
                                                                                            asset.get(),
                                                                                            allowNegativeVariance,
                                                                                            recorder);

                        //stripes for lowerBarrier using bruteForce integration
                        recorder->setScalingFactor(2.0 * nbOptionsForVar);
                        divStripesOptionValues += 2.0 * CorridorVarSwap::bruteForceIntegral(intOrder,
                                                                                            lowerBarr,
                                                                                            lowerBarr / (1.0 - d),
                                                                                            dDate,
                                                                                            valueDate,
                                                                                            fwd,
                                                                                            asset.get(),
                                                                                            allowNegativeVariance,
                                                                                            recorder);
                    }
                }
                //other cases
                else{
                    //DIV_DEFAULT
                    if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_DEFAULT) ){
                        intOrder = 5;
                    }
                    //different orders 2..7
                    else{
                        char* endPtr;
                        intOrder = (int)strtol(methodology.c_str(), &endPtr, 10 /* base 10 */);
                    }

                    //rolling over dividends
			        for(int iDiv = 0; iDiv < NDiv; iDiv++){
				        double d = (*divListAmount)[iDiv];
				        const DateTime& dDate = divListDatesArray[iDiv];
                        double fwd = fwdsAtDiv[iDiv];
                        
        				//stripes for uppperBarrier using NC
		        		recorder->setScalingFactor(-2.0 * nbOptionsForVar);
                        divStripesOptionValues -= 2.0 * CorridorVarSwap::NewtonCotesCallIntegral(intOrder,
				        												                         upperBarr,
						        										                         upperBarr / (1.0 - d),
                                                                                                 dDate,
                                                                                                 valueDate,
												        				                         fwd,
                                                                                                 asset.get(),
                                                                                                 allowNegativeVariance,
                                                                                                 recorder);

				        //stripes for lowerBarrier using NC
				        recorder->setScalingFactor(2.0 * nbOptionsForVar);
                        divStripesOptionValues += 2.0 * CorridorVarSwap::NewtonCotesCallIntegral(intOrder,
									        							                         lowerBarr,
											        					                         lowerBarr / (1.0 - d),
													        			                         dDate,
                                                                                                 valueDate,
                                                                                                 fwd,
                                                                                                 asset.get(),
                                                                                                 allowNegativeVariance,
                                                                                                 recorder);
                    }
                }
            }
			
            // Total Floating Var before basis
            //total return
            //Manos initial implementation
            if( CString::equalsIgnoreCase(methodology, ClosedFormIntegrateLN::DIV_CONTINUOUS_APPROX) ){
                futureFloatingVar = terminalPortfolios + termStructurePortfolios;
            }
            //other cases
            else{
                //must add divStripeOptionValues corresponding to dividend component
                futureFloatingVar = terminalPortfolios + termStructurePortfolios + divStripesOptionValues;
            }

            //Quadratic variation term (flag inst->dividendAdjusted for the instrument)
            double quadTermDiv = 0.0;
            if( ! (inst->dividendAdjusted) ){
                for(int iObs = 0; iObs < obsDates.size()-1; iObs++){
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
			        double fwd1 = asset->fwdValue(Date1); //t-1 sampling (default implementation for a corridor)
                    double fwd2 = asset->fwdValue(Date2); //t sampling
                    
                    //barrier Levels
                    double lowerBarr = instLowerBarrier;
			        double upperBarr = instUpperBarrier;
                    
                    recorder->setScalingFactor(sumDiv * nbOptionsForVar);
                    //if T-1 sampling, obsDate[iObs] is counting
                    if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_DEFAULT) || CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T_MINUS_1)  ){
                        //quadVar term
                        quadTermDiv +=  sumDiv * CorridorVarSwap::quadTerm(lowerBarr,
                                                                           upperBarr,
                                                                           Date1,
                                                                           valueDate,
                                                                           fwd1,
                                                                           asset.get(),
                                                                           allowNegativeVariance,
                                                                           recorder);
                    }
                    if( CString::equalsIgnoreCase(checkRangeAtT, CorridorVarSwap::SAMPLING_T)  ){
                        //quadVar term
                        quadTermDiv +=  sumDiv * CorridorVarSwap::quadTerm(lowerBarr,
                                                                           upperBarr,
                                                                           Date2,
                                                                           valueDate,
                                                                           fwd2,
                                                                           asset.get(),
                                                                           allowNegativeVariance,
                                                                           recorder);
                    }
                }
                //Add quadTermDiv to futureFloatingVar
                futureFloatingVar += quadTermDiv;
            }

            /// III) Corridor Swap Basis
            if (model->getUseBasis()) {
                // Obtain Variance Swap basis
                DateTimeArray dates(1, lastDate);
                if(valueDate < firstDate) {
                // Forward starting case: put the timeline in reverse order
                    dates.push_back(firstDate);
                }
                VarSwapBasis::VarSwapBasisRequestSP basisRequest(new VarSwapBasis::VarSwapBasisRequest(dates, discount.getSP()));
                IVolProcessedSP procVol;
                try {
                    IVolProcessedSP tmp(asset->getProcessedVol(basisRequest.get()));
                    procVol = tmp;
                } catch(exception&) {
                    //  Ignore all errors and proceed without basis
                    //  Could be that VolSpline threw an error like: don't know what VarSwapBasisRequest is
                    basisProc = VarSwapBasis::VarSwapBasisProcSP(   );
                }

                if(procVol.get()) {
                    VarSwapBasis::VarSwapBasisProcError* error = dynamic_cast<VarSwapBasis::VarSwapBasisProcError*>(procVol.get());
                    if(error) {
                        throw ModelException(const_cast<ModelException&>(error->getException()), routine);
                    }
                    basisProc = VarSwapBasis::VarSwapBasisProcSP::dynamicCast(procVol);
                }

                // We now have a processed basis
                if(basisProc.get()) {
                    //intermediate quantities needed in basis computation, i.e. futureExpectedN;
                    ATMVolRequest volReq;
                    CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
                    double tradYear = vol->GetTimeMetric()->yearFrac(valueDate, lastDate);
                    double varSwapVolBasis = basisProc->interpVolBasis(tradYear);

                    // Disallow cutoff methodology
                    double priceCutoff = basisProc->interpPriceCutoff(tradYear);
                    if(!Maths::isZero(priceCutoff)) {
                        throw ModelException("Price cutoff not supported. Either remark basis or set useBasis to false for this instrument.");
                    }
                    double skewCutoff  = basisProc->interpSkewCutoff(tradYear);
                    if(!Maths::isZero(skewCutoff)) {
                        throw ModelException("Skew cutoff not supported. Either remark basis or set useBasis to false for this instrument.");
                    }

                    // Compute corridor basis adjustment
                    double corridorVolBasis = 0.0;
                    if(!Maths::isZero(futureRange * futureFloatingVar)) {
                        double tmpVolFuture = sqrt((double)observationsPerYear / futureRange * futureFloatingVar);
                        double corridorVolBasis = varSwapVolBasis * futureRange / futureExpectedN;
	                    double corridorBasisVar = (double)futureExpectedN / (double)observationsPerYear * (2.0 * tmpVolFuture * corridorVolBasis + Maths::square(corridorVolBasis));
	                    futureFloatingVar += corridorBasisVar;
                    }
                }
            }
        }

        // Compute value of legs
        double floatingLeg =
            ((double)observationsPerYear / (double)expectedN) *
            (pastFloatingVar + futureFloatingVar);
        if (inst->isCapped) {
            if(Maths::isZero(ExpRangeNForCap)) {
                    throw ModelException("ExpRangeNForCap should be set to a non-zero probability number; otherwise price will be zero");
            }
            floatingLeg = Maths::min(floatingLeg, ExpRangeNForCap * Maths::square(cap));
            // cap is now interpreted as a volatility number instead of a vol multiplier (i.e. new cap = old cap * strike)
        }

        double totalRange = pastRange + futureRange;
        double fixedLeg = Maths::square(strike) * totalRange / (double)expectedN;

        // Compute price for non-settled instruments
        double fwdValuedFloatingLeg = 100.0 * inst->notional * floatingLeg;
        double fwdValuedFixedLeg = 100.0 * inst->notional * fixedLeg;
        double optionScalingFactor = 100.0 * inst->notional;
        if(inst->scaleByStrike) {
            fwdValuedFloatingLeg /=  2.0 * strikeRef;
            fwdValuedFixedLeg /=  2.0 * strikeRef;
            optionScalingFactor /= 2.0 * strikeRef;
        }

        double fwdValuedPrice = fwdValuedFloatingLeg - fwdValuedFixedLeg;

        double floatingLegPrice = 0.0;
        double fixedLegPrice = 0.0;
        double price = 0.0;
        double pv = 1.0;
        DateTime settlement  = instSettle->settles(lastDate, asset.get());
        if((valueDate <= settlement)) {
            // PV and scale
            pv = instSettle->pv(valueDate, lastDate, discount.get(), asset.get());
            floatingLegPrice = pv * fwdValuedFloatingLeg;
            fixedLegPrice = pv * fwdValuedFixedLeg;
            price = floatingLegPrice - fixedLegPrice;
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
                    const_cast<CorridorVarSwap*>(inst),
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
            OutputRequest* request;
            
            // Compute all vol requests first
            // TOTAL_VOL
            double totalVar = pastFloatingVar + futureFloatingVar;
            if (!Maths::isZero(totalRange) && !Maths::isNegative(totalVar)) {
                double totalvol = sqrt((double)observationsPerYear / totalRange * totalVar);
                request = control->requestsOutput(OutputRequest::TOTAL_VOL);
                if(request) {
                    results->storeRequestResult(request, totalvol);
                }
            }

            // Compute future vol
            double volFuture = 0.0;
            bool volFutureDefined;
            if(futureExpectedN) {
                // There exists future
                if(!Maths::isZero(futureRange * futureFloatingVar) && !Maths::isNegative(futureFloatingVar)) {
                    volFuture = sqrt((double)(observationsPerYear) / futureRange * futureFloatingVar);
                    volFutureDefined = true;
                } else {
                    // Not defined
                    volFutureDefined = false;
                }
            } else {
                // No future
                volFutureDefined = true;
            }

            if(volFutureDefined) {
                // VOL_IN_FUTURE
                request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
                if(request) {
                    results->storeRequestResult(request, volFuture);
                }

                // IND_VOL
                request = control->requestsOutput(OutputRequest::IND_VOL);
                if (request) {
                    results->storeRequestResult(request, volFuture);
                }
            }

            // Compute past vol
            double volPast = 0.0;
            // If there are no past samples in the corridor, VOL_IN_PAST is now
            // reported as 0
            if(pastExpectedN) {
                if(!Maths::isZero(pastRange) && !Maths::isNegative(pastFloatingVar)) {
                    // Well defined
                    volPast = sqrt((double)(observationsPerYear) / pastRange * pastFloatingVar);
                }
            } 

            // VOL_IN_PAST
            request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                results->storeRequestResult(request, volPast);
            }

            // PAST_WEIGHT
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                if(!Maths::isZero(totalRange)) {
                    double pastWeight = pastRange / totalRange;
                    results->storeRequestResult(request, pastWeight);
                }
            }

            //EXPECTED_PCT_TIME_IN_RANGE
            request = control->requestsOutput(OutputRequest::EXPECTED_PCT_TIME_IN_RANGE);
            if (request && !Maths::isZero(futureExpectedN)) {
                double expectedPctTimeInRange = futureRange / futureExpectedN;
                results->storeRequestResult(request, expectedPctTimeInRange);
            }

            // CORRIDOR_VARIANCE_VALUE
            request = control->requestsOutput(OutputRequest::CORRIDOR_VARIANCE_VALUE);
            if (request) {
                results->storeRequestResult(request, floatingLegPrice);
            }

            // CORRIDOR_ACCRUAL_VALUE
            request = control->requestsOutput(OutputRequest::CORRIDOR_ACCRUAL_VALUE);
            if (request) {
                results->storeRequestResult(request, -fixedLegPrice);
            }

            // Record basis information at maturity i.e. ignore forward starting VarSwaps for here
            if(basisProc.get()) {
                ATMVolRequest volReq;
                CVolProcessedBSSP vol(asset->getProcessedVol(&volReq));
                double tradYear = vol->GetTimeMetric()->yearFrac(valueDate, lastDate);
                basisProc->recordRequests(tradYear, control, results);
            }

            //STRIKE_VOL
            request = control->requestsOutput(OutputRequest::STRIKE_VOL);
            if (request) {
               results->storeRequestResult(request, strike);
            }

            //STRIKE_REF
            request = control->requestsOutput(OutputRequest::STRIKE_REF);
            if (request) {
               results->storeRequestResult(request, strikeRef);
            }

            //EXPECTED_N
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, expectedN);
            }

            //DISCOUNT_FACTOR
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request && valueDate <= settlement) {
                results->storeRequestResult(request, pv);
            }

            // BARRIER LEVEL
            request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
            if (request) {
                // Report barriers within a window
                DateTime upperBarrierDate = BarrierLevel::barrierWindow(valueDate);
                BarrierLevelArraySP reportLevels(new BarrierLevelArray(0));
                for(int iObsDate = 0; iObsDate < obsDates.size(); iObsDate++) {
                    const DateTime& thisObsDate = obsDates[iObsDate];
                    if(valueDate <= thisObsDate && thisObsDate <= upperBarrierDate) {
                        // Down barrier
                        BarrierLevel bl(false, thisObsDate, instLowerEconBarrier, false);
                        reportLevels->push_back(bl);
                        
                        // Up barrier
                        BarrierLevel bh(true,  thisObsDate, instUpperEconBarrier, false);
                        reportLevels->push_back(bh);
                    }
                }
                if (!reportLevels->empty()) {
                    OutputRequestUtil::recordBarrierLevels(control,
                                                           results,
                                                           asset->getTrueName(),
                                                           reportLevels.get());
                }
            }

            //PAYMENT_DATES
            request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
            if (request) {
                DateTimeArray date(1, settlement);
                OutputRequestUtil::recordPaymentDates(control, results, &date);
            }

            //KNOWN_CASH_FLOWS
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

            //DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                results,
                price,
                valueDate,
                discount.get(),
                asset.get(),
                inst->premiumSettle.get());

            //FWD_AT_MAT
            InstrumentUtil::recordFwdAtMat(control,
                results,
                lastDate,
                valueDate,
                asset.get());
        }
    }  catch(exception& e) {
      throw ModelException(e, routine);
    }
}


double CorridorVarSwapProduct::pricePortfolio(VanillaContractsRecorderSP    recorder,
                                              CAssetConstSP                 asset,
                                              const ClosedFormIntegrateLN*  model,
                                              const Control*                control,
                                              const DateTime&               valueDate,
                                              const DateTime&               maturity,
                                              double                        lowerBarr,
                                              double                        upperBarr,
                                              bool                          allowNegativeVar) {
    static const string routine = "CorridorVarSwapProduct::pricePortfolio";
    try {
        // I) Get integrator from model
        double fwd = asset->fwdValue(maturity);
        ClosedBoundary upper(upperBarr / fwd);
        ClosedBoundary lower(lowerBarr / fwd);
        Range integrationDomain(lower, upper);
        Integrator1DSP integrator;
        model->getIntegrator(control, asset.get(), valueDate, maturity, recorder,
            integrator, integrationDomain);


        // II) Modify scaling factor
        double scalingFactor = recorder->getScalingFactor() * 2.0;
        recorder->setScalingFactor(scalingFactor);

        // III) Construct integrand and integrate
        CorridorVarSwapIntegrandWeightConstSP optionWeights(new CorridorVarSwapIntegrandWeight(integrationDomain));
        StaticReplicationIntegrandSP integrand(new StaticReplicationIntegrand(
            optionWeights, recorder, asset, valueDate, maturity, allowNegativeVar));
        double integral = integrator->integrate(*integrand);
        double price = 2.0 * integral;
        return price;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}


CorridorVarSwapProduct::CorridorVarSwapIntegrandWeight::CorridorVarSwapIntegrandWeight(const Range& integrationDomain):
Function1DDouble(integrationDomain) {}


double CorridorVarSwapProduct::CorridorVarSwapIntegrandWeight::operator()(double relativeStrike) const {
    double relativeStrikeWeight = 1.0 / (Maths::square(relativeStrike));
    return relativeStrikeWeight;
}

//sampling used to weigh log-retruns
const string CorridorVarSwap::SAMPLING_DEFAULT = string("Default");
const string CorridorVarSwap::SAMPLING_T_MINUS_1 = string("SamplingT-1");
const string CorridorVarSwap::SAMPLING_T = string("SamplingT");


////////////////////////////////////////////////////////////////////////////


ClosedFormIntegrateLN::IProduct* CorridorVarSwap::createProduct(ClosedFormIntegrateLN* model) const {
	return new CorridorVarSwapProduct(this);
}


bool CorridorVarSwapLoad() {
    return (CorridorVarSwap::TYPE != 0);
}

//utility functions for integration

//quadVar Term
double CorridorVarSwap::quadTerm(double          L,
                                 double          U,
                                 const DateTime& divDate,
                                 const DateTime& valueDate,
                                 double          fwd,
                                 const CAsset*   asset,
                                 bool            allowNegativeVar,
                                 VanillaContractsRecorderSP recorder){
    static const string method("CorridorVarSwap::quadTerm");
    try{
        //price of digitals
        double delta = CorridorVarSwapProduct::PUT_SPREAD_WIDTH;
        DoubleArray strikes(4);
        DoubleArray weights(4);
        strikes[0] = L * (1 - delta);
        strikes[1] = L * (1 + delta);
        strikes[2] = U * (1 - delta);
        strikes[3] = U * (1 + delta);
        weights[0] = 1.0 / (2.0 * delta * L);
        weights[1] = 1.0 / (-2.0 * delta * L);
        weights[2] = 1.0 / (-2.0 * delta * U);
        weights[3] = 1.0 / (2.0 * delta * U);

        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, divDate, false));
        volRequest->allowNegativeFwdVar(allowNegativeVar);

        CVolProcessedBSSP volBS(asset->getProcessedVol(volRequest.get()));
        double yearFrac = volBS->calcTradingTime(valueDate, divDate);
       
        // Compute probabilities
        // Price difference of put spreads to approximate probability of being in corridor for a given date
        double res = 0.0;
        for(int iStrike = 0; iStrike < 4; iStrike++) {
            volRequest->setStrike(strikes[iStrike]);
            volBS = CVolProcessedBSSP(asset->getProcessedVol(volRequest.get()));
            double vol = volBS->CalcVol(valueDate, divDate);
            double var = Maths::square(vol) * yearFrac;
            bool isCall = (strikes[iStrike]) < fwd ? false : true;
            double thisPrice = BlackPrice(isCall, fwd, strikes[iStrike], 1.0, vol, yearFrac,
                weights[iStrike], divDate, recorder);
            res += thisPrice;
        }

        //forward component
        bool boolFwd = Maths::isPositive(fwd - L) && Maths::isPositive(U - fwd);
        if(boolFwd){
            res += 1.0;
        } 
        
        //result
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }                
}

//Newton Cotes method for numerical integration
double CorridorVarSwap::NewtonCotesCallIntegral(int                  n,
							                    double               LMin,
							                    double               LMax,
							                    const DateTime&      divDate,
                                                const DateTime&      valueDate,
   							                    double               fwd,
							                    const CAsset*        asset,
                                                bool                 allowNegativeVar,
                                                VanillaContractsRecorderSP recorder){
    static const string method("CorridorVarSwap::NewtonCotesCallIntegral");
    try{
        //definition of strikes
        DoubleArray strikesNC(n);

        //definition of weights
        DoubleArray coeffNC(n);

        //main coefficient in Newton Cotes
        double hNC = (LMax - LMin) / ((double)n - 1.0);
        double resCoeff;
        
        //strikes
        for(int iNC = 0; iNC < n; iNC++){
            strikesNC[iNC] = LMin + (double)iNC * hNC;
        }

       //order
       switch(n){
       //Trapezium, n = 2
       case 2:
           //resCoeff
           resCoeff = hNC / 2.0;

	       //definitions of coeffNC
	       coeffNC[0] = 1.0;
	       coeffNC[1] = 1.0;

	       break;
           
       //Simpson, n = 3
       case 3:
           //resCoeff
           resCoeff = hNC / 3.0;
           
           //definitions of coeffNC
	       coeffNC[0] = 1.0;
	       coeffNC[1] = 4.0;
	       coeffNC[2] = 1.0;
           
           break;
           
           //Simpson 3/8, n = 4
       case 4:
           //resCoeff
           resCoeff = 3.0 * hNC / 8.0;
           
           //definitions of coeffNC
	       coeffNC[0] = 1.0;
	       coeffNC[1] = 3.0;
	       coeffNC[2] = 3.0;
	       coeffNC[3] = 1.0;
           
           break;

           //Boole, n = 5
        case 5:
            //resCoeff
	        resCoeff = 2.0 * hNC / 45.0;

	        //definitions of coeffNC
	       coeffNC[0] = 7.0;
	       coeffNC[1] = 32.0;
	       coeffNC[2] = 12.0;
	       coeffNC[3] = 32.0;
	       coeffNC[4] = 7.0;
	       
	       break;

           //6 points rule, n = 6
        case 6:
            //resCoeff
            resCoeff = 5.0 * hNC / 288.0;

            //definition of coeffNC
            coeffNC[0] = 19.0;
            coeffNC[1] = 75.0;
            coeffNC[2] = 50.0;
            coeffNC[3] = 50.0;
            coeffNC[4] = 75.0;
            coeffNC[5] = 19.0;
            
            break;

           //7 points rule, n = 7
        case 7:
            //resCoeff
	        resCoeff = 1.0 * hNC / 140.0;

	        //definitions of coeffNC
	       coeffNC[0] = 41.0;
	       coeffNC[1] = 216.0;
	       coeffNC[2] = 27.0;
	       coeffNC[3] = 272.0;
	       coeffNC[4] = 27.0;
	       coeffNC[5] = 216.0;
	       coeffNC[6] = 41.0;
	       
           break;

           //other cases
        default:
            throw ModelException("if even not implemented yet, if odd error", method);
        }

        //sum of option prices
        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, divDate, false));
        volRequest->allowNegativeFwdVar(allowNegativeVar);
        CVolProcessedBSSP volBS( asset->getProcessedVol(volRequest.get()) );
        double yearFrac = volBS->calcTradingTime(valueDate, divDate);
        double res = 0.0;
        for(int iOption = 0; iOption < n; iOption++){
	        volRequest->setStrike(strikesNC[iOption]);
	        volBS = CVolProcessedBSSP ( asset->getProcessedVol(volRequest.get()) );
	        double vol = volBS->CalcVol(valueDate, divDate);
            double var = Maths::square(vol) * yearFrac;
	        bool isCall = (strikesNC[iOption]) < fwd ? false : true;
            
            double thisPrice = BlackPrice(isCall, fwd, strikesNC[iOption], 1.0, vol, yearFrac,
                resCoeff * coeffNC[iOption] / (strikesNC[iOption] * strikesNC[iOption]), divDate, recorder);
            res += thisPrice;
        }

        //result
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

//bruteForce integration
double CorridorVarSwap::bruteForceIntegral(int                  n,
                                           double               LMin,
                                           double               LMax,
                                           const DateTime&      divDate,
                                           const DateTime&      valueDate,
                                           double               fwd,
                                           const CAsset*        asset,
                                           bool                 allowNegativeVar,
                                           VanillaContractsRecorderSP recorder){
    static const string method("CorridorVarSwap::bruteForceIntegral");
    try{
        //definition of strikes
        double hNC = (LMax - LMin) / ((double)n - 1.0);

        //strikes
        DoubleArray strikesNC(n);
        for(int iNC = 0; iNC < n; iNC++){
            strikesNC[iNC] = LMin + (double)iNC * hNC;
        }

        //sum of option prices
        LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(1.0, valueDate, divDate, false));
        volRequest->allowNegativeFwdVar(allowNegativeVar);
        CVolProcessedBSSP volBS( asset->getProcessedVol(volRequest.get()) );
        double yearFrac = volBS->calcTradingTime(valueDate, divDate);
        double res = 0.0;
        for(int iOption = 0; iOption < n; iOption++){
	        volRequest->setStrike(strikesNC[iOption]);
	        volBS = CVolProcessedBSSP( asset->getProcessedVol(volRequest.get()) );
	        double vol = volBS->CalcVol(valueDate, divDate);
            double var = Maths::square(vol) * yearFrac;
	        bool isCall = (strikesNC[iOption]) < fwd ? false : true;
	        
            double thisPrice = BlackPrice(isCall, fwd, strikesNC[iOption], 1.0, vol, yearFrac,
                hNC / (strikesNC[iOption] * strikesNC[iOption]), divDate, recorder);
            res += thisPrice;
        }

        //result
        return res;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}



DRLIB_END_NAMESPACE
