//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarianceOption.cpp
//
//   Description : Variance Option contract
//
//   Author      : Andrew McCleery
//
//   Date        : 21 Apr 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Black.hpp"
#include "edginc/Generic1Factor.hpp"
#include "edginc/FourierEngine.hpp"
#include "edginc/FourierProcessSVJ.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/VarSwapUtilities.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/ObservationBuilder.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/VolRequestTime.hpp"

DRLIB_BEGIN_NAMESPACE

// Variance Option
class VarianceOption: public Generic1Factor,
                 virtual public FourierEngine::IIntoProduct,
                 virtual public Theta::Shift,
                 virtual public LastSensDate
{
public:
    static CClassConstSP const TYPE;

    /** Instrument validation */
    virtual void Validate();

    /** Retrieve market data */
    void GetMarket(const IModel*       model,
                   const CMarketDataSP market);

    virtual bool priceDeadInstrument(CControl* control, CResults* results) const;

    // Rolls value date and set initial spot for Theta
    virtual bool sensShift(Theta* shift);

    virtual DateTime endDate(const Sensitivity* sensControl) const;

    /** Implementation of FourierEngine::IntoProduct interface */
    virtual FourierProduct* createProduct(const FourierEngine* model) const;

protected:
    friend class VarianceOptionFP;

private:
    // Validation and population of transient fields
    void validatePop2Object();

    /** For reflection */
    static void load(CClassSP& clazz);

    /** Empty shell method */
    static IObject* defaultVarianceOption();

    // constructor
    VarianceOption();

    // not implemented
    VarianceOption(const VarianceOption& rhs);
    VarianceOption& operator=(const VarianceOption& rhs);

    // Input fields
    double                  strike;                 //!< Strike Vol
    double                  refVol;                 //!< Used in Payoff divisor to scale notional
    int                     expectedN;              //!< Divisor for realized var
    int                     observationsPerYear;    //!< PPY
    bool                    dividendAdjusted;       //!< Whether to adjust returns for divs
    bool                    scaleByRefVol;          //!< Whether to scale by 2*refVol
    string                  assetHistorySource;     //!< Designates asset history source
    IObservationBuilderSP   observationBuilder;     //!< Factory for observations
    bool                    isCall;                 //!< Call/Put flag

    // Transient fields
    DateTimeArraySP         obsDates;
    ObservationTypeArraySP  obsTypes;
    DoubleArraySP           obsSamples;
    ObservationSourceSP     assetHistorySourceObject; //!< object for asset history source
};

/** Instrument validation */
void VarianceOption::Validate() {
    static const string routine = "VarianceOption::Validate";
    try {
        // Override fwdStarting flag on Generic1Factor
        fwdStarting = startDate > valueDate;

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

/** Retrieve market data */
void VarianceOption::GetMarket(const IModel*          model,
                          const CMarketDataSP    market) {
    static const string routine = "VarianceOption::GetMarket";
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

bool VarianceOption::priceDeadInstrument(CControl* control,
                                    CResults* results) const {
    static const string routine = "VarianceOption::priceDeadInstrument";
    try {
        DateTime lastDate    = obsDates->back();
        DateTime settlement  = instSettle->settles(lastDate, asset.get());

        if (valueDate < lastDate) {
            // Case 1: valueDate < lastDate <= settlement
            // Alive instrument - nothing to do here
            return false;
        }
        else if (valueDate.isGreaterOrEqual(settlement)) {
            // Case 2: lastDate <= settlement <  valueDate
            // price = 0.0
            results->storePrice(0.0, discount->getCcy());
            return true;
        } else {
            // Case 3: lastDate <= valueDate  <= settlement
            // price = pv of known cashflows
            double histVar = VarSwapUtilities::realisedVar(
                  asset.get(), *obsDates.get(), *obsSamples.get(), valueDate, true, dividendAdjusted, false);

            double totalYears = (double)expectedN/(double)observationsPerYear;
            double volPast = sqrt(histVar / totalYears);

            double price = Maths::max(volPast*volPast - strike*strike, 0.0);
            if (scaleByRefVol) {
                price /= (2.0*strike);
            }
            double pv = instSettle->pv(valueDate,
                                       obsDates->back(),
                                       discount.get(),
                                       asset.get());
            price *= pv*notional*100.0;
            results->storePrice(price, discount->getCcy());                        

            if (control && control->isPricing()) {
                OutputRequest* request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
                if (request) {
                    results->storeRequestResult(request, volPast);
                }
                request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
                if (request) {
                    results->storeRequestResult(request, 1.0);
                }
                request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
                if (request) {
                    results->storeRequestResult(request, 0.0);
                }
                request = control->requestsOutput(OutputRequest::TOTAL_VOL);
                if (request) {
                    results->storeRequestResult(request, volPast);
                }
                request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
                if (request) {
                    results->storeRequestResult(request, pv);
                }

                // KNOWN_CASHFLOWS
                request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
                if (request) {
                    DateTime paymentDate = instSettle->settles(obsDates->back(), asset.get());
                    CashFlow cf(paymentDate, price / pv);
                    CashFlowArray cfl(1, cf);
                    OutputRequestUtil::recordKnownCashflows(control,
                        results,
                        discount->getCcy(),
                        &cfl);
                }

                // PAYMENT_DATES
                request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
                if (request) {
                    DateTime paymentDate = instSettle->settles(obsDates->back(), asset.get());
                    DateTimeArray date(1, paymentDate);
                    OutputRequestUtil::recordPaymentDates(control,results,&date);
                }

                // Expected N
                request = control->requestsOutput(OutputRequest::EXPECTED_N);
                if (request) {
                    results->storeRequestResult(request, expectedN);
                }

                // Option Implied Vol
                request = control->requestsOutput(OutputRequest::VAR_OPTION_IND_VOL);
                if (request) {
                    results->storeRequestResult(request, 0.0);
                }
            }
            // Dead instrument
            return true;
        }
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

// Rolls value date and set initial spot for Theta
bool VarianceOption::sensShift(Theta* shift) {
    static const string routine = "VarianceOption::sensShift";
    try {
        // Fill in cashflows with spot or fwd
        VarSwapUtilities::thetaShiftCashFlows(shift, asset.get(), valueDate, *obsDates, *obsSamples);

        // Call parent
        return Generic1Factor::sensShift(shift);
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

DateTime VarianceOption::endDate(const Sensitivity* sensControl) const {
    static const string routine = "VarianceOption::endDate";
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

void VarianceOption::validatePop2Object() {
    static const string routine = "VarianceOption::validatePop2Object";
    try {
        if(Maths::isNegative(strike)) {
            throw ModelException(routine, "Strike cannot be negative (" + Format::toString(strike) + ")");
        }

        if(observationsPerYear < 1) {
            throw ModelException(routine,
                "Observations per year " + Format::toString(observationsPerYear) + " must be strictly positive");
        }
 
        if(scaleByRefVol && Maths::isZero(refVol)) {
            throw ModelException(routine,
                "Cannot scale by refVol when refVol is zero or not supplied");
        }

        obsDates = observationBuilder->dateList();
        obsTypes = observationBuilder->obsTypes();
        obsSamples = DoubleArraySP(new DoubleArray(obsDates->size()));

        // expectedN defaulting and validation
        int noReturns = obsDates->size() - 1;
        if (expectedN < 1) {
            if (0 == expectedN) {
                // Not initialised
                expectedN = noReturns;
            }
            else {
                throw ModelException(routine,
                    "Expected N " + Format::toString(expectedN) + " must be positive ");
            }
        }
        else {
            // Validate expectedN consistency with sample schedule
            if (expectedN != noReturns) {
                throw ModelException(routine, "expectedN (" + Format::toString(expectedN) +
                    ") does not match the number of returns in the sample schedule (" +
                    Format::toString(noReturns) + ")");
            }
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

void VarianceOption::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Variance Option Instrument");
    REGISTER(VarianceOption, clazz);
    SUPERCLASS(Generic1Factor);
    IMPLEMENTS(FourierEngine::IIntoProduct);
    IMPLEMENTS(LastSensDate);
    EMPTY_SHELL_METHOD(defaultVarianceOption);
    FIELD(strike, "Volatility strike");
    FIELD(refVol, "Reference vol for notional scaling");
    FIELD_MAKE_OPTIONAL(refVol);
    FIELD(expectedN, "Expected number of returns");
    FIELD_MAKE_OPTIONAL(expectedN);
    FIELD(observationsPerYear, "Periods per year");
    FIELD_MAKE_OPTIONAL(observationsPerYear);
    FIELD(dividendAdjusted, "Whether to adjust returns for divs");
    FIELD(scaleByRefVol, "True: divide by 2*refVol, False: don't divide by 2*refVol");
    FIELD_MAKE_OPTIONAL(scaleByRefVol);
    FIELD(assetHistorySource, "Asset history source");
    FIELD_MAKE_OPTIONAL(assetHistorySource);
    FIELD(observationBuilder, "Observation builder");
    FIELD(isCall, "Call/Put flag");
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

IObject* VarianceOption::defaultVarianceOption() {
    return new VarianceOption();
}

VarianceOption::VarianceOption(): Generic1Factor(TYPE), refVol(0.0),
expectedN(0), observationsPerYear(252), scaleByRefVol(true),
assetHistorySource(IMarketObservable::DEFAULT_SOURCE) {}


CClassConstSP const VarianceOption::TYPE = CClass::registerClassLoadMethod(
    "VarianceOption", typeid(VarianceOption), VarianceOption::load);

////////////////////////////////////////////////////////////////////////////

// ** using Fourier method for SVJ **
template <class Process, class Product>
class VarianceOptionIntegrand:public Function1DDouble {
public:
    VarianceOptionIntegrand(const  Process& process,
                    const  Product& product,
                    const  DateTime& maturity,
					bool   isVolSwap,
                    double strike):
        Function1DDouble(Range(OpenBoundary(0.0), Infinity(Infinity::Plus))),
        process(process), product(product), matdate(maturity), isVolSwap(isVolSwap), strike(strike) {}

    double operator()(double  x) const {
		if(isVolSwap) {
			// based on Matytsin 2000 Merrill Lynch presentation for VolatilitySwap
			return (1.0 - exp(process.cumulant(product, -x, matdate) + x * strike).real()) / (x * sqrt(x));
		}
		else {
#ifdef RG
			static double omega = +0.05;    // XXX make this an input
			// based on similar derivation to that in 'Single Integrand Approach' (Section 2.1.1.)
			// Venardos (2002) 'The Theory Behind the Fourier Engine'
			Complex z(omega, x);
			return (exp(process.cumulant(product, z, matdate) + z * strike) / (z * z)).real();
#else
			// based on Matytsin 2000 Merrill Lynch presentation
			Complex z(0.0, x);
			return ( (1.0 - exp(process.cumulant(product, z, matdate)-z*strike))/(x * x) ).real();
#endif
		}
    }

private:

    const  Process& process;
    const  Product& product;
    DateTime        matdate;
	bool            isVolSwap;
    double          strike;
};


class VarianceOptionFP: public FourierProduct,
                public FourierProductIntegrator1D,
                public StFourierProductQuadVar,
                public FwdStFourierProductQuadVar {
private:
    const VarianceOption*    inst;

    DateTime			maturity;
    //double				mult;
    double				years;
    double				totalYears;
    double				strike;
    double				volPast;
    double				pastWeight;

public:
    // equivalent to InstIntoFourierProduct
    VarianceOptionFP(const VarianceOption* inst);

    virtual void validate(FourierProcess* process) const;

    virtual void price(const FourierEngine* model,
                       Control*             control,
                       Results*             results);

    /** Constructs integrands for var option */
    Function1DDoubleArrayConstSP Integrand(const FourierEngine * model,
                                           const Integrator1D*  integrator);

    const DateTime& getStartDate() const {
        return inst->startDate;
    }

    /** Post process method for VolSwapFP */
    virtual void postResults(const FourierEngine* model,
                             const Integrator1D*  integrator,
                             const FourierProductIntegrator1D::IntegralArray& integrals,
                             CControl*            control,
                             CResults*            results);
};


// equivalent to InstIntoFourierProduct
VarianceOptionFP::VarianceOptionFP(const VarianceOption* inst):
    FourierProduct(inst->asset.get(),
                   inst->valueDate,
                   inst->discount.get(),
                   inst->instSettle.get()),
                   inst(inst),
                   maturity(inst->obsDates->back()){
}

void VarianceOptionFP::validate(FourierProcess* process) const {
    const static string method("VarianceOptionFP::validate");
    int i = 0;

    for(; i < mAsset->NbAssets(); i++) {
        string CcyTtreatment(mAsset->assetGetCcyTreatment(i));
        if(CcyTtreatment == "S") {
            throw ModelException(method,
                                 "Currency Struck assets not supported.");
        }
    }
    // give the process the chance to do some extra validation
    // and to initialize its transient fields
    process->validate(this);
}

void VarianceOptionFP::price(const FourierEngine* model,
                             Control*             control,
                             Results*             results){

    // valueDate >= matDate is taken care of here
    if(inst->priceDeadInstrument(control, results)){
        return; // dead instrument priced
    }

    FourierProduct::price(model,
                          control,
                          results);
}

/** Constructs integrands for var option */
Function1DDoubleArrayConstSP VarianceOptionFP::Integrand(const FourierEngine * model,
                                                    const Integrator1D*  integrator) {
    static const string method = "VarianceOptionFP::Integrand";

    try{
        // Keep references for easier access
        const DateTime& valueDate = inst->valueDate;
        CAssetConstSP asset = inst->asset.getSP();
        const DateTimeArray& obsDates = *inst->obsDates;
        DoubleArray& obsSamples = *inst->obsSamples;
        int expectedN = inst->expectedN;
        int observationsPerYear = inst->observationsPerYear;

        const FourierProcess& process = model->getProcess();

        // STRIKE = (cap*strike)^2 * ExpN / PPY - PASTVAR - FUTURE_DIV_VAR

        // Compute past realized variance
        double histVar = VarSwapUtilities::realisedVar(
                asset.get(), obsDates, obsSamples, valueDate, true, inst->dividendAdjusted, false);

        // Compute past and future ExpectedN
        int pastExpectedN;
        int futureExpectedN;
        VolRequestTime volReq;
        IVolProcessedSP vol(asset->getProcessedVol(&volReq));
        HolidayWrapper assetHols = vol->GetTimeMetric()->getHolidays();
        VarSwapUtilities::pastAndfutureExpN(obsDates, assetHols, valueDate, expectedN, pastExpectedN, futureExpectedN);

        //double annPastVar = varPast * pastWeight * totalYears;
        // pastWeight, totalYears, strike, volPast are transient
        pastWeight = (double)pastExpectedN/(double)expectedN;
        totalYears = (double)expectedN/(double)observationsPerYear;
        strike = (inst->strike * inst->strike) * totalYears - histVar;
        if (obsDates.front() >= valueDate) {
            volPast = 0.0;
        } else {
            volPast = sqrt(histVar / (pastWeight*totalYears));
        }

        if (!(inst->dividendAdjusted) && Maths::isPositive(1.0 - pastWeight)) {

            // now add effect due to dividends.  Note: discreteDivsAdjustment computes the effect of the future dividend
            // on the total variance, to avoid confusion with both past / future weights, and switching from a continuous
            // integration to a discrete summand.  Thus we have to normalize to determine modification to
            // varFuture
            double effectFromDivs = VarSwapUtilities::futureDiscreteDivsAdjustment(asset.get(), valueDate,
                obsDates.front(), maturity, observationsPerYear, expectedN + 1);

            strike -= effectFromDivs * totalYears;
        }

        // STRIKE = (cap*strike)^2 * ExpN / PPY - PASTVAR - FUTURE_DIV_VAR

        Function1DDoubleArraySP functions(new Function1DDoubleArray(1));

        if (inst->fwdStarting){
            // populate transient field
            years = process.getTimeMetric().yearFrac(obsDates.front(), maturity);

            const FwdStFourierProductQuadVar& thisProd = *this;
            const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);

            if(!thisProc) {
                throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface." );
            }

            (*functions)[0] = Function1DDoubleSP(new VarianceOptionIntegrand<FwdStFourierProcessQuadVar, FwdStFourierProductQuadVar>
                                                      (*thisProc,
                                                       thisProd,
                                                       maturity,
												       false,
                                                       strike));
        } else {
            // populate transient field
            years = process.getTimeMetric().yearFrac(valueDate, maturity);

            const StFourierProductQuadVar& thisProd = *this;
            const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);

            if(!thisProc) {
                throw ModelException(method, "Process does not support StFourierProcessQuadVar interface." );
            }

            (*functions)[0] = Function1DDoubleSP(new VarianceOptionIntegrand<StFourierProcessQuadVar, StFourierProductQuadVar>
                                                      (*thisProc,
                                                       thisProd,
                                                       maturity,
												       false,
                                                       strike));
        }

        return functions;
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Post process method for VolSwapFP */
void VarianceOptionFP::postResults(const FourierEngine* model,
                              const Integrator1D*  integrator,
                              const FourierProductIntegrator1D::IntegralArray& integrals,
                              CControl*            control,
                              CResults*            results) {
    static const string method = "VarianceOptionFP::postResults";
    double price = 0.0;

    try {
        const FourierProcess& process = model->getProcess();
        double varHM;
        if (inst->fwdStarting){
            const FwdStFourierProductQuadVar& thisProd = *this;
            const FwdStFourierProcessQuadVar* thisProc = dynamic_cast<const FwdStFourierProcessQuadVar*>(&process);
            if(!thisProc) {
                throw ModelException(method, "Process does not support FwdStFourierProcessQuadVar interface" );
            }
            varHM = thisProc->expectation(thisProd, maturity);
        } else {
            const StFourierProductQuadVar& thisProd = *this;
            const StFourierProcessQuadVar* thisProc = dynamic_cast<const StFourierProcessQuadVar*>(&process);
            if(!thisProc) {
                throw ModelException(method, "Process does not support StFourierProcessQuadVar interface" );
            }
            varHM = thisProc->expectation(thisProd, maturity);
        }
        //future realized variance or quadratic variation (NOT NORMALIZED BY TIME, ie a sum of logReturnsSQ)
        double futVar = varHM * years;
		price = integrals[0] / Maths::PI;
		//price = price +.5*(varHM * years - strike);
		price = price +.5*(futVar - strike);
        price /= totalYears;

        // See if option price is positive
        if (price <= 0.0) {
            if(price > -0.01) {
                price = 0.0;
            } else {
                throw ModelException("(Scaleless) Price of variance option is negative " +
                                     Format::toString(price) +
                                     ". Check integrator parameters e.g. tolerance, nbSteps.");
            }
        }

        // E (AvgTotalVar)

        // Do dividend adjustments
        double divVar = 0.0;
        if (!inst->dividendAdjusted) {
              // No dividend adjustment = use full return = add the effect of divs to vanillas
              divVar = VarSwapUtilities::futureDivsVariance(
                  inst->asset.get(), inst->valueDate, inst->valueDate, maturity);
        }
        
        //pastVar is realized Var accumulated in the past (volPast is normalized, hence multiply by time)
        //pastWeight * totalYears = discrete fraction of time from first obsDate to valueDate
        double pastVar = volPast * volPast * pastWeight * totalYears;
        
        //futVar is future integrated variance or future quadratic variation
        //obtained using expectation from the relevant stochVol model (expectation is normalized, hence multiply by time)  
        //years = continuous fraction of time from first obsDate to maturity
        /*
        double futVar = varHM * years;
        */

        //divVar is future variance due to future dividends (specific to quadVar)
        //(no need to multiply by time, sum log(Retun)^2 has the dimension of a variance)
        
        //total realized + future var (or quadVar) normalized by fraction of years
        //double unscaledFwd   = ((varHM + divVar) * years + volPast * volPast * pastWeight * totalYears) / totalYears;
        double unscaledFwd   = (pastVar + futVar + divVar) / totalYears;
        
        double pv = inst->instSettle->pv(inst->valueDate,
                                         maturity,
                                         inst->discount.get(),
                                         inst->asset.get());

        if(!inst->isCall) {
            price -= unscaledFwd - inst->strike*inst->strike;
        }

        // Keep some info aside before we apply notions, scaling by 2K and pv
        // Call: E ( max(AvgTotalVar - Strike^2,0) )
        // Put: E ( max(Strike^2 - AvgTotalVar,0) )
        double unscaledPrice = price;

        price *= pv*inst->notional*100.0;

        if (inst->scaleByRefVol) {
            price /= (2.0*inst->refVol);
        }

        results->storePrice(price, discount->getCcy());

        // take care of additional outputs
        if (control && control->isPricing()) {
            // Past Vol
            OutputRequest* request = control->requestsOutput(OutputRequest::VOL_IN_PAST);
            if (request) {
                results->storeRequestResult(request, volPast);
            }

            // Past Weight
            request = control->requestsOutput(OutputRequest::PAST_WEIGHT);
            if (request) {
                results->storeRequestResult(request, pastWeight);
            }

            /* Compute Past and Future Vols */ //see above for definitions of var terms
            
            //contributionPast is pastVar normalized by ExpNPast / PPY = volPast * volPast * pastWeight
            double contributionPast = volPast*volPast*pastWeight;
            //contributionFuture
            double contributionFuture = 0.0;
            if (pastWeight < 1.0) {
                // Note if pricing on maturity date we have numPastReturns = expectedN (pastWeight = 1)
                // The way time is apportioned means that VOL_IN_FUTURE = 0 and TOTAL_VOL = VOL_IN_PAST in this case
                
                //future realized variance or quadratic variation is futVar = varHM * years
                //future contribution due to dividends is divVar

                //contributionFuture = (varHM * years + divVar) normalized by ExpN / PPY
                contributionFuture = (futVar + divVar) * (double)inst->observationsPerYear/(double)(inst->expectedN);
            }

            // Future Vol
            request = control->requestsOutput(OutputRequest::VOL_IN_FUTURE);
            if (request) {

                double volFuture = 0.0;
                if (pastWeight < 1.0) {
                    // expN / futureN
                    double adj = 1.0 / (1.0 - pastWeight);
                    volFuture = sqrt(contributionFuture*adj);
                }

                results->storeRequestResult(request, volFuture);
            }

			// IND_VOL:
            double totalVolSquared = contributionPast + contributionFuture;

            // Total Vol
            request = control->requestsOutput(OutputRequest::TOTAL_VOL);
            if (request) {
                results->storeRequestResult(request, sqrt(totalVolSquared));
            }

            // Discount
            request = control->requestsOutput(OutputRequest::DISCOUNT_FACTOR);
            if (request) {
                results->storeRequestResult(request, pv);
            }

            // Expected N
            request = control->requestsOutput(OutputRequest::EXPECTED_N);
            if (request) {
                results->storeRequestResult(request, inst->expectedN);
            }
/*
            OLD CODE, kept for the moment
            
            // Option Implied Vol
            request = control->requestsOutput(OutputRequest::VAR_OPTION_IND_VOL);
            if (request) {
                if(Maths::isPositive(years)) {
                    double variance = 0.0;
                    try {
                        bool success = Black::impliedVariance(inst->isCall,                                          // option call or put
                                                              unscaledFwd,                                     // forward price
                                                              inst->strike*inst->strike,                                       // strike
                                                              1.0,                                             // pv
                                                              0.3 * 0.3 * years,                               // initial var guess
                                                              unscaledPrice,                                   // option price
                                                              2.0 * 0.3 * 1.0e-5 * years,                      // var accuracy
                                                              variance);

                        if(success) {
                            double impliedVol = sqrt(variance / years);
                            results->storeRequestResult(request, impliedVol);
                        }
                    } catch(exception&) {
                        // Don't bother
                    }
                }
            }
*/
            //GAD, correction for IND_VOL
            // Option Consistent Implied Vol (based on future realized variance / past realized variance accounted for by modifying strike)
            request = control->requestsOutput(OutputRequest::VAR_OPTION_IND_VOL);
            if (request) {
                if(Maths::isPositive(years)) {
                    double variance = 0.0;
                    try{
                        //consistent values
                        //consistent fwd is based on future realized variance
                        double consFwdPrice = (unscaledFwd - pastVar / totalYears) / (1.0 - pastWeight) ; //past var has been taken out and put in the strike
                        //consistent strike accounts for past realized variance
                        double consStrike = (inst->strike*inst->strike - pastVar / totalYears) / (1.0 - pastWeight) ;   
                        //consistent price is scaled by pastN
                        double consUnscaledPrice = unscaledPrice / (1.0 - pastWeight);

                        //test on strike, if negative option is worthless and we give zero implied vol
                        if(Maths::isPositive(consStrike)){
                            bool success = Black::impliedVariance(inst->isCall,                                  // option call or put
                                                                consFwdPrice,//changed                           // forward price
                                                                consStrike, //changed                           // strike
                                                                1.0,                                             // pv
                                                                0.3 * 0.3 * years,                               // initial var guess
                                                                consUnscaledPrice, //changed                     // option price
                                                                2.0 * 0.3 * 1.0e-5 * years,                      // var accuracy
                                                                variance);
                            if(success) {
                                double impliedVol = sqrt(variance / years);
                                results->storeRequestResult(request, impliedVol);
                            }
                        }
                    }
                    catch(exception&) {
                        // Don't bother
                    }
                }
            }
            //end of GAD

			// payment_dates
			request = control->requestsOutput(OutputRequest::PAYMENT_DATES);
			if (request) {
				DateTime paymentDate = inst->instSettle->settles(maturity, inst->asset.get());
				DateTimeArray date(1, paymentDate);
				OutputRequestUtil::recordPaymentDates(control,results,&date);
			}

            // DELAY_PRICE
            InstrumentUtil::delayPriceHelper(control,
                                             results,
                                             price,
                                             inst->valueDate,
                                             inst->discount.get(),
                                             inst->asset.get(),
                                             inst->premiumSettle.get());

            try {
                // FWD_AT_MAT
                InstrumentUtil::recordFwdAtMat(control,
                                               results,
                                               maturity,
                                               inst->valueDate,
                                               inst->asset.get());
            } catch(exception& e) {
                // Don't rethrow
                ModelException ee = ModelException::addTextToException(e, method + ": Computation of forward price failed");

                IObjectSP fwdObj(new Untweakable(ee));

                results->storeRequestResult(control->requestsOutput(OutputRequest::FWD_AT_MAT),
                                            fwdObj,
                                            OutputNameSP(new OutputName(inst->asset.get()->getName())));
            }
        }
    }

    catch (exception& e){
        throw ModelException(e, method);
    }
}

FourierProduct* VarianceOption::createProduct(const FourierEngine* model) const {
    static const string routine("VarianceOption::createProduct");
    try {
        return new VarianceOptionFP(this);

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

bool VarianceOptionLoad()
{
    return (VarianceOption::TYPE != 0);
}

DRLIB_END_NAMESPACE
