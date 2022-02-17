//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VarSwapHedgingSupport.cpp
//
//   Author      : Romain Garagnon
//
//   Description   Computes the P&L distribution for a trading strategy based on:
//                 1. a long or short position on a variance swap
//                 2. the replicating portfolio of call and put options
//                 This includes:
//                 1. transaction costs for options
//                 2. transaction costs for stocks and futures (delta-hedge)
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/SampleList.hpp"
#include "edginc/VolVarSwap.hpp"
#include "edginc/Black.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE

//////// Instrument class ////////
class VarSwapHedgingSupport: public Generic1Factor, 
                             virtual public LastSensDate,
                             virtual public IMCIntoProduct {
protected:
    
    bool                 positionLS; // short or long position on variance swap

	// variance swap
	SampleListSP         samples; // sample dates and past samples if any
    double               strikeVol; // contract strike in volatility units
    int                  observationsPerYear; // number of sample points per year
    bool                 subtractMeanVol; // if true, the mean volatility is subtracted
    bool                 dontScaleByStrike; // if true, the notional is not divided by 2 * strikeVol
    bool                 noDivAdj; // if true, the additional effect E[sum(log(1-D)^2)] due to discrete dividends to variance is added
    bool                 divAdjOnExDate; // add div on exDate sample or subtract from exDate - 1 sample
    ImpliedIntegrationSP varSwapModel; // model used to price the variance swap
    VarianceSwapSP       varSwap; // variance swap instrument

    // replicating portfolio
    DoubleArray          strikes; // strikes of the options used to replicate the variance swap
    DoubleArraySP        impliedVols; // implied volatilities at start date of these options
    DoubleArraySP        bidVols; // bid implied volatility - mid implied volatility
    DoubleArraySP        askVols; // ask implied volatility - mid implied volatility
    BoolArray            callputFlags; // flags indicating whether the option is a call or a put
    DoubleArraySP        hedges; // number of options (calls and puts) in the replicating portfolio
    double               lowStrike; // lower strike used for stripping the log-contract
    double               highStrike; // upper strike used for stripping the log-contract

    double               stockTransCost; // transaction cost for stocks (per 1mm notional)
    double               futureTransCost; // transaction cost for futures (per 1mm notional)

    // transient (to get implied vols only)
    MarketObjectArray    assetObjects;
            
public:
    static CClassConstSP const TYPE;
    friend class VarSwapHedgingSupportProd;

    /** validation */
    virtual void Validate() {
        static const string routine("VarSwapHedgingSupport::Validate");

        // call Generic1Factor validate()
        validate();

        AssetUtil::assetCrossValidate(asset.get(),
                                      fwdStarting,
                                      startDate,
                                      valueDate,
                                      discount,
                                      this);

        // calculate implied volatilities at start date for the different options
        impliedVols = impliedVolatilities();
        
        // create variance swap instrument
        string payoffType = "FORWARD"; // payoff type = FORWARD (variance swap)
        double cap = 1.0; // shouldn't be used
        int numPastReturns = 0; // number of past returns
        int numTotalReturns = 0;  // total number of returns
        numReturns(numPastReturns, numTotalReturns);

        CashFlowArraySP mySamples = CashFlow::createCashFlowArray(samples.get()->getDates(),
                                                                  samples.get()->getValues());

        CAssetSP myAsset = CAssetSP::dynamicCast(assetObjects[0]);
        CAssetWrapper myAssetVS(copy(myAsset.get()));

        varSwap = VarianceSwapSP(
                      new VarianceSwap(valueDate,
                                       startDate, 
                                       fwdStarting, 
                                       1.0, // initial spot 
                                       ccyTreatment,
                                       instSettle, 
                                       premiumSettle, 
                                       myAssetVS, 
                                       discount,
                                       (*mySamples.get()), 
                                       strikeVol, 
                                       observationsPerYear,
                                       subtractMeanVol, 
                                       payoffType, 
                                       dontScaleByStrike,
                                       cap, 
                                       noDivAdj,
                                       divAdjOnExDate,
                                       false, // is vanilla = false
                                       numPastReturns,
                                       numTotalReturns));
        
        // compute the number of calls and puts
        hedges = computeHedges();
    }

    /** validation */
    void validatePop2Object() {
        static const string routine("VarSwapHedgingSupport::validatePop2Object");
        if (fwdStarting)
		    throw ModelException(routine, "fwdStarting flag has to be set to false.");
        
        if (oneContract)
		    throw ModelException(routine, "oneContract flag has to be set to false.");

        // checks that start date = value date = first sample date
        if (!startDate.equals(valueDate))
		    throw ModelException(routine, "startDate and valueDate must coincide.");

        if (!startDate.equals(samples.get()->getFirstDate()))
		    throw ModelException(routine, "startDate and first sample date must coincide.");
        
		// checks that the asset is not currency struck
		if (ccyTreatment == CAsset::CCY_TREATMENT_STRUCK)
			throw ModelException(routine, "Ccy struck type not supported.");
		
		// checks that the asset is not currency protected
		if (ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED)
			throw ModelException(routine, "Ccy protected type not supported.");
		
		// checks there is at least one past sample
        if (samples.get()->getDates().size() == 0)
		    throw ModelException(routine, "There must be at least one sample date.");

        // checks the number of strikes in the replicating portfolio
        int nberOptions = strikes.size();
        if (nberOptions < 2) {
            throw ModelException(routine, "There must be at least two strikes in the replication portfolio.");
        }

        // checks if the number of strikes is the same than the number of call / put flags
        if (nberOptions != callputFlags.size())
		    throw ModelException(routine, "The number of strikes and call put flags "
                                          "should be the same.");

        // checks if the number of strikes is the same than the number bid / ask vols
        if (nberOptions != bidVols.get()->size())
		    throw ModelException(routine, "The number of strikes and bid volatilities "
                                          "should be the same.");

        if (nberOptions != askVols.get()->size())
		    throw ModelException(routine, "The number of strikes and ask volatilities "
                                          "should be the same.");

        // checks there are at least one call and one put
        int iStrk;
        int indexFirstCall; // index giving the position of the first call of the first call
        bool foundCall = false; // flag indicating if there is a call
        for (iStrk = 0 ; iStrk < nberOptions ; iStrk++) {
            if (callputFlags[iStrk])
            {
                indexFirstCall = iStrk;
                foundCall = true;
            }
        }
        
        if (!foundCall) // no call has been found
            throw ModelException(routine, "At least one call is needed.");
        if (indexFirstCall == 0) // no put has been found
            throw ModelException(routine, "At least one put is needed.");
        
        // check that there are not put inserted in the calls
        // check that there are no call inserted in the puts
        for (iStrk = indexFirstCall ; iStrk < nberOptions ; iStrk++) {
            if (!callputFlags[iStrk])
                throw ModelException(routine, "The " + Format::toString(iStrk + 1) + "-th strike "
                                              "corresponding to a put is inserted between two calls .");
        }
        
        // checks that the strikes are positive
        for (iStrk = 0 ; iStrk < nberOptions ; iStrk++) {
            if (!Maths::isPositive(strikes[iStrk]))
                throw ModelException(routine, "The " + Format::toString(iStrk + 1) + "-th strike is negative.");
        }

        // checks that the strikes are in increasing order
        for (iStrk = 1 ; iStrk < nberOptions ; iStrk++) {
            if (Maths::isPositive(strikes[iStrk - 1] - strikes[iStrk]) &&
                (callputFlags[iStrk - 1] || !callputFlags[iStrk]))
                throw ModelException(routine, "The " + Format::toString(iStrk) + "-th strike is "
                                              " greater than the " + Format::toString(iStrk + 1) + "-th strike.");
        }
    }

    // initiate GetMarket 
    void GetMarket(const IModel* model, const CMarketDataSP market)
    {
        market->GetReferenceDate(valueDate);

        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, discount, asset);
        discount.getData(model, market);
        instSettle->getMarket(model, market.get());
        if (premiumSettle.get())
        {
            premiumSettle->getMarket(model, market.get());
        }

        // copy of the asset with BS compatible volatility type
        MarketDataFetcherSP mdf(new MarketDataFetcherLN("VolSpline"));
        NonPricingModel dummyModel(mdf);
        MarketObjectSP myAsset(market->GetData(&dummyModel, asset.getName(), CAsset::TYPE));
        assetObjects.push_back(myAsset);
    }

    /** calculate implied volatilities at start date for the different options */
    DoubleArraySP impliedVolatilities() const {
        static const string method = "VarSwapHedgingSupport::impliedVolatilities";

        int nberStrikes = strikes.size(); // nber of strikes
        DoubleArraySP impliedVols = DoubleArraySP(new DoubleArray(nberStrikes));

        DateTime matDate = samples.get()->getDates().back(); // maturity date
        
        int iStrk;
        for (iStrk = 0 ; iStrk < strikes.size() ; iStrk++)
        {
            // forward price at start date
            double fwdAtStart = fwdStarting ? asset->fwdValue(startDate) : initialSpot;
            
            // choose how to interpolate the vol
            LinearStrikeVolRequestSP volRequest(
                new LinearStrikeVolRequest(fwdAtStart * strikes[iStrk], 
                                           startDate, 
                                           matDate,
                                           fwdStarting));

            // interpolate the vol for each strike
            CAssetSP myAsset = CAssetSP::dynamicCast(assetObjects[0]);
            CVolProcessedSP vol(myAsset->getProcessedVol(volRequest.get()));

            // cast to the type of vol we are expecting
            CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);

            // this should never happen if our get market data has worked properly
            if (!vol)
                throw ModelException(method, "No Black-Scholes volatility");
            
            // calculate the volatility
            double volatility = volBS->CalcVol(startDate, matDate);
                if (!Maths::isPositive(volatility))
                    throw ModelException(method, "The volatility calculated at date " + startDate.toString() +
                                                 " for strike " + Format::toString(strikes[iStrk]) + " is negative.");
                
            (*impliedVols)[iStrk] = volatility;
        }

        return impliedVols;
    }

    /** computes the number of past returns and the total number of returns */
    void numReturns(int& numPastReturns, int& numTotalReturns) const {
        static const string method = "VarSwapHedgingSupport::numReturns";
        try {
            DateTimeArray sampleDates = samples.get()->getDates(); // sample dates
            
            // total number of returns
            numTotalReturns = samples.get()->getDates().size() - 1;
            
            // all the sample dates are in the past
            if (valueDate >= sampleDates.back()) {
                numPastReturns = numTotalReturns; // all the returns are in the past
            }
            else if (valueDate > sampleDates.front()) {

                HolidayConstSP marketHols = AssetUtil::getHoliday(asset.get()); // market holidays

                // find the next sample date    
                int iDate = 0;
                while (iDate < sampleDates.size() && sampleDates[iDate] < valueDate) {
                    iDate++;
                }

                // expected number of returns from next sample date
                int numFutureReturns = marketHols->businessDaysDiff(sampleDates[iDate], sampleDates.back());

                // number of returns to next sample is defined as:
                // (total num returns expected on instrument creation) - (expected num returns from next sample)
                numPastReturns = numTotalReturns - numFutureReturns;
                numPastReturns = Maths::max(numPastReturns, 1);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    /** computes the number of options (calls and puts) in the replicating portfolio */
    DoubleArraySP computeHedges() const {
        static const string method = "VarSwapHedgingSupport::computeHedges";

        int nberOptions = strikes.size(); // number of different options in the replicating portfolio
        
        // number of options (calls and puts) in the replicating portfolio
        DoubleArraySP hedges = DoubleArraySP(new DoubleArray(nberOptions));
        
        DateTime maturity = samples.get()->getDates().back(); // maturity of the variance swap
        double initSpot = asset->getSpot(); // initial spot price
        HolidayConstSP marketHols = AssetUtil::getHoliday(asset.get()); // market holidays
        
        VolRequestTimeSP volTimeReq(new VolRequestTime()); // vol request for time
		CVolProcessedSP volTimeProcessed(asset->getProcessedVol(volTimeReq.get())); 
		TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric(); // time metric
        double timeToMaturity =  timeMetric->yearFrac(samples->getFirstDate(),
                                                      samples->getLastDate());// time to maturity in years
        
        int iStrk;
        int nberCalls = 0; // number of call options
        int nberPuts = 0; // number of put options
        for (iStrk = 0 ; iStrk < nberOptions ; iStrk++) {
            if (callputFlags[iStrk]) {
                nberCalls++;
            }
            else {
                nberPuts++;
            }
        }
        
        // computation of the hedges for the calls
        DoubleArray strikeCalls(nberCalls + 1); // strikes for calls including highStrike)
        DoubleArray payoffCalls(nberCalls + 1); // corresponding payoffs
        DoubleArray callDelta(nberCalls); // call delta for the different strikes
        
        // strikes and payoffs for calls
        for (iStrk = 0 ; iStrk < strikeCalls.size() - 1 ; iStrk++) {
            strikeCalls[iStrk] = strikes[nberPuts + iStrk];
        }
        strikeCalls[strikeCalls.size() - 1] = highStrike;
        
        for (iStrk = 0 ; iStrk < strikeCalls.size() ; iStrk++) {
            // the payoff function is given by
            // - (2 / T) * ((log(ST / S0) - ((ST - S0) / S0))
            payoffCalls[iStrk] = -2.0 * (log(strikeCalls[iStrk]) - (strikeCalls[iStrk] - 1.0)) / timeToMaturity;
        }

        // call deltas
        for (iStrk = 0 ; iStrk < strikeCalls.size() - 1 ; iStrk++) {
            callDelta[iStrk] = (payoffCalls[iStrk+1] - payoffCalls[iStrk]) /
                (strikeCalls[iStrk+1] - strikeCalls[iStrk]) / initSpot;
        }
        
        // call hedges
        (*hedges)[nberPuts] = callDelta[0];
        for (iStrk = 1 ; iStrk < nberCalls ; iStrk++) {
            (*hedges)[nberPuts + iStrk] = callDelta[iStrk] - callDelta[iStrk-1];
        }

        // computation of the hedges for the puts
        DoubleArray strikePuts(nberPuts + 1); // strikes for puts including lowStrike)
        DoubleArray payoffPuts(nberPuts + 1); // corresponding payoffs
        DoubleArray putDelta(nberPuts); // put delta for the different strikes

        // strikes and payoffs for puts
        strikePuts[0] = lowStrike;
        for (iStrk = 1 ; iStrk < strikePuts.size() ; iStrk++) {
            strikePuts[iStrk] = strikes[iStrk-1];
        }
        
        for (iStrk = 0 ; iStrk < strikePuts.size() ; iStrk++) {
            payoffPuts[iStrk] = -2.0 * (log(strikePuts[iStrk]) - (strikePuts[iStrk] - 1)) / timeToMaturity;
        }

        // put deltas
        for (iStrk = 1 ; iStrk < strikePuts.size() ; iStrk++) {
            putDelta[iStrk-1]  = (payoffPuts[iStrk-1] - payoffPuts[iStrk]) /
                (strikePuts[iStrk-1] - strikePuts[iStrk]) / initSpot;
        }

        // put hedges
        (*hedges)[nberPuts-1] = -putDelta[nberPuts-1];
        for (iStrk = nberPuts-1 ; iStrk > 0 ; iStrk--) {
            (*hedges)[iStrk-1] = putDelta[iStrk] - putDelta[iStrk-1];
        }

        return hedges;
    }

    /** implement MonteCarlo::IntoProduct interface */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const;

    /** return the date at which tweaking stops (default implementation assumes product can be priced with a MC) */
    DateTime endDate(const Sensitivity*) const {
        DateTime end = samples.get()->getDates().back();
        return end;
    }

private:
    VarSwapHedgingSupport(): Generic1Factor(TYPE), 
                                 dontScaleByStrike(false),
                                 noDivAdj(false),
                                 divAdjOnExDate(false) {}

    // for reflection
    VarSwapHedgingSupport(const VarSwapHedgingSupport& rhs); // not implemented
    VarSwapHedgingSupport& operator=(const VarSwapHedgingSupport& rhs); // not implemented

    static IObject* defaultVarSwapHedgingSupport() {
        return new VarSwapHedgingSupport();
    }

    //// roll through time (setting historic values)
    bool sensShift(Theta* theta) {
        return false;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){

	    REGISTER(VarSwapHedgingSupport, clazz);
        SUPERCLASS(Generic1Factor);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(LastSensDate);   
        IMPLEMENTS(Theta::Shift);
        EMPTY_SHELL_METHOD(defaultVarSwapHedgingSupport);
        FIELD(positionLS, "short or long position on variance swap");
        FIELD(samples, "sample dates and past samples if any");
        FIELD(strikeVol, "contract strike in volatility units");
        FIELD(observationsPerYear, "number of sample points per year");
        FIELD(subtractMeanVol, "if true, the mean volatility is subtracted");
        FIELD(dontScaleByStrike, "if true, the notional is not divided by 2 * strikeVol");
        FIELD_MAKE_OPTIONAL(dontScaleByStrike);
        FIELD(noDivAdj, "if true, the additional effect E[sum(log(1-D)^2)] due to discrete dividends to variance is added");
        FIELD_MAKE_OPTIONAL(noDivAdj);
        FIELD(divAdjOnExDate, "true - add to ex-date sample; false - subtract from previous sample");
        FIELD_MAKE_OPTIONAL(divAdjOnExDate);   
        FIELD(varSwapModel, "model to price the variance swap");
        FIELD(varSwap, "variance swap used in payoff function")
        FIELD_MAKE_TRANSIENT(varSwap);
        FIELD(strikes, "strikes of the options used to replicate the variance swap");
        FIELD(impliedVols, "implied volatilities at start date of these options");
        FIELD_MAKE_TRANSIENT(impliedVols);
        FIELD(bidVols, "bid implied volatility - mid implied volatility");
        FIELD(askVols, "ask implied volatility - mid implied volatility");
        FIELD(callputFlags, "flags indicating whether the option is a call or a put");
        FIELD(hedges, "number of options (calls and puts) in the replicating portfolio");
        FIELD_MAKE_TRANSIENT(hedges);
        FIELD(lowStrike, "lower strike used for stripping the log-contract");
        FIELD(highStrike, "upper strike used for stripping the log-contract");
        FIELD(stockTransCost, "transaction cost for stocks (per 1mm notional)");
        FIELD(futureTransCost, "transaction cost for futures (per 1mm notional)");
        FIELD(assetObjects, "");
        FIELD_MAKE_TRANSIENT(assetObjects);

        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
};
// end of class VarSwapHedgingSupport

//////// product class //////////
class VarSwapHedgingSupportProd : public IMCProduct, virtual public IMCProductLN {
private:

    const VarSwapHedgingSupport* inst; // reference to original instrument
    DateTimeArray                simulationDates; // simulation dates
	DoubleArray                  histSpots; // all historical spots
    DateTimeArray                histDates; // all historical dates
    DoubleArray                  pvFactor; // pv factor to last date
    DoubleArray                  growthFactor; // growth factor to settlement date of last date
    DoubleArray                  yearFrac; // time in years to last date
    
    VarianceSwapSP               varSwapMC; // variance swap instrument
    CControlSP                   control; // control used to compute price and delta of variance swap
    
    double                       pathValue; // store the path value
    double                       plDeltaHedge; // store the P&L related to delta-hedge
    double                       transactionCost; // store the transaction cost
    // for preservation of the past
    double                       pathValueSoFar; // store the path value
    double                       plDeltaHedgeSoFar; // store the P&L related to delta-hedge
    double                       transactionCostSoFar; // store the transaction cost
    
	class VarSwapHedgingSupportPrices;
    typedef refCountPtr<VarSwapHedgingSupportPrices> VarSwapHedgingSupportPricesSP;
    
    class VarSwapHedgingSupportPrices: public IMCPricesSimple {
    public:
        
        enum PricesIndices {
            PL_AVERAGE = 0,
            NB_PRICES
        };

		/** Adds supplied P&L to compute the P&L average of the hedging strategy */
        void add(double price, int index) {
            simplePrices[index]->add(price);
        }

        /** Adds supplied P&L to the array representing the P&L distribution  */
        void addPathValue(double price, int indexPath) {
            (*pathValues)[indexPath] = price;
        }

        /** Adds supplied delta-hedge P&L to the array representing the P&L distribution  */
        void addPathDeltaHedgeValue(double price, int indexPath) {
            (*pathDeltaHedgeValues)[indexPath] = price;
        }

        /** Adds supplied spots path to the array of spots paths  */
        void addSpotsPath(const double* spotPath, int nberSpots, int indexPath) {
            int iSpot;
            for (iSpot = 0 ; iSpot < nberSpots ; iSpot++)
            {
                (*spotPaths)[indexPath][iSpot] = spotPath[iSpot];
            }
        }

        /** Returns the averaged P&L together with its standard error */
        virtual void getResult(double* result, double* resultStdErr) const {
            simplePrices[PL_AVERAGE]->getResult(result, resultStdErr);
        }

        /** Returns the P&L distribution */
        DoubleArraySP getPathValues() const {
            return pathValues;
        }

        /** Returns the delta-hedge P&L distribution */
        DoubleArraySP getPathDeltaHedgeValues() const {
            return pathDeltaHedgeValues;
        }

        /** Returns the spots paths */
        DoubleArrayArraySP getSpotsPaths() const {
            return spotPaths;
        }

        double getResults(int index) const {
            double pathValue, dummy;
            simplePrices[index]->getResult(&pathValue, &dummy);
			return pathValue;
        }

        /** Returns true if the path, identified by pathIdx, should be
            repriced. If false is returned then there will be no add()
            method called for this path and the IMCPrices object must
            take any appropriate action */
        virtual bool repriceForGreek(int pathIdx) {
            return true;
        }

        /** Support product-level caching across tweaks */
        virtual int storagePerPath(IMCProduct* product) const {
            return 0;
        }

        virtual void configureCache(const IntArray& changedAssets) {
		}

        /** Returns a deep copy of this object */
        IMCPrices* clone() const {
            VarSwapHedgingSupportPricesSP copy(new VarSwapHedgingSupportPrices());
            copy->simplePrices.resize(simplePrices.size());
            for (unsigned int iPrice = 0 ; iPrice < copy->simplePrices.size() ; iPrice++) {
                IMCPricesSP temp(this->simplePrices[iPrice]->clone());
                copy->simplePrices[iPrice] = DYNAMIC_POINTER_CAST<IMCPricesSimple>(temp);
            }
            return copy.get();
        }

        VarSwapHedgingSupportPrices(int NbIter, int NbSubSamples, int nberSamples):
        simplePrices(NB_PRICES) {
            pathValues = DoubleArraySP(new DoubleArray(NbIter)); // resize pathValues
            pathDeltaHedgeValues = DoubleArraySP(new DoubleArray(NbIter)); // resize pathDeltaHedgeValues
            // resize spotPaths
            spotPaths = DoubleArrayArraySP(new DoubleArrayArray(NbIter));
            int iIter;
            for (iIter = 0 ; iIter < NbIter ; iIter++) {
                (*spotPaths)[iIter].resize(nberSamples);
            }
            for (unsigned int iPrice = 0 ; iPrice < simplePrices.size() ; iPrice++) {
                simplePrices[iPrice] = IMCPricesSimpleSP(new MCPricesSimple(NbIter, NbSubSamples));
            }
        }

    private:
        VarSwapHedgingSupportPrices() {
		}

		/** adds supplied price to this set of IMCPrices */
        virtual void add(double price) {
            throw ModelException("IMCPrices::add(double)",
                                 "internal error");
        }

        /** Returns the last price stored. Undefined behaviour if no
            prices yet stored */
        virtual double lastPrice() const {
            throw ModelException("IMCPrices::lastPrice",
                                 "internal error");
        }

        /** On pricing run returns MAX(x, 0.0). It should be used for applying
            the 'final point of optionality' within the product. This allows
            QuickGreeks type of IMCPrices not to apply the max when doing first
            order greeks (this may sound strange but see the docs for why) */
        virtual double maxWithZero(double x) const {
            throw ModelException("IMCPrices::maxWithZero",
                                 "internal error");
        }

        /** Reset this object so that it can be used for the same operation
            again. Normally, a new IMCPrices object is created for each
            pricing run. However, for quick x gamma, it is important to use
            the same one for each pair of assets */
        virtual void reset() {
            throw ModelException("IMCPrices::reset",
                                 "internal error");
        }

        /** Ease cloning */
        virtual IMCPrices* emptyConstructor() const {
            throw ModelException("IMCPrices::emptyConstructor",
                                 "internal error");
        }

        vector<IMCPricesSimpleSP> simplePrices;
        DoubleArraySP pathValues;
        DoubleArraySP pathDeltaHedgeValues;
        DoubleArrayArraySP spotPaths;
    };

public:
    
	virtual IMCPrices* createOrigPrices(int nbIter,
                                     int nbSubSamples,
                                     int mode) {
        return new VarSwapHedgingSupportPrices(nbIter, nbSubSamples, inst->samples.get()->getDates().size());
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen) const {
        // empty
    }

    /** equivalent to InstIntoMCProduct. */
    VarSwapHedgingSupportProd(const VarSwapHedgingSupport* inst,
                              IRefLevelSP refLevel,
                              const SimSeriesSP& simSeries) :
        IMCProduct(inst->asset.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  refLevel,
                  simSeries,
                  inst->samples,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
                  inst(inst) {
                    
		simulationDates = simSeries->getAllDates();

        histDates = inst->samples.get()->getDates();
        histSpots = inst->samples.get()->getValues();

        // compute pv factors from value date to settlement date of current date
		pvFactor.resize(histDates.size());
        
        int iDate;
		for (iDate = 0 ; iDate < histDates.size() ; iDate++)
		{
            DateTime settleDate = inst->instSettle->settles(histDates[iDate], NULL);
			if (inst->valueDate <= histDates[iDate])
                pvFactor[iDate] = inst->instSettle->pv(inst->valueDate,
                                                       settleDate,
                                                       inst->discount.get(),
                                                       inst->asset.get());
            else
                pvFactor[iDate] = 999999.9; // this number should never be used
		}

        // compute growth factors from payment dates to last date
		growthFactor.resize(histDates.size());
        DateTime lastSettleDate = inst->instSettle->settles(histDates.back(), NULL);

        for (iDate = 0 ; iDate < histDates.size() ; iDate++)
		{
            if (inst->valueDate <= histDates[iDate])
                growthFactor[iDate] = 1.0 / inst->instSettle->pv(histDates[iDate],
                                                                 lastSettleDate,
                                                                 inst->discount.get(),
                                                                 inst->asset.get());
            else
                growthFactor[iDate] = 999999.9; // this number should never be used
		}

        // compute year fraction from payment dates to last date
		yearFrac.resize(histDates.size());

        VolRequestTimeSP volTimeReq(new VolRequestTime()); // vol request for time
		CVolProcessedSP volTimeProcessed(getMultiFactors()->factorGetProcessedVol(0, volTimeReq.get())); 
		TimeMetricConstSP timeMetric = volTimeProcessed->GetTimeMetric(); // time metric
        
        for (iDate = 0 ; iDate < histDates.size() ; iDate++)
		{
			if (inst->valueDate <= histDates[iDate])
                yearFrac[iDate] = timeMetric->yearFrac(histDates[iDate], histDates.back());
            else
                yearFrac[iDate] = 999999.9; // this number should never be used
		}

        // create control used to compute price and delta of variance swap
        control = CControlSP(Control::makeFromFlags("", 0.0));
        
        // initialize path value
        pathValue = 0.0;
        pathValueSoFar = 0.0;

        // initialize P&L related to delta-hedge
        plDeltaHedge = 0.0;
        plDeltaHedgeSoFar = 0.0;

        // initialize transaction cost
        transactionCost = 0.0;
        transactionCostSoFar = 0.0;
        
        // create a copy of the original variance swap
        varSwapMC = VarianceSwapSP(copy(inst->varSwap.get()));
    }

    /** compute the P&L of the hedging strategy for a given path */
	void payoff(const IPathGenerator*  pathGen, IMCPrices&  prices);

    /** compute the payoff of the variance swap at maturity */
    double payoffVarSwapAtMat(const double* path, int nberSamples) const;

	/** output statistics of P&L distribution and P&L distribution itself */
    virtual void recordExtraOutput(Control* control, Results* results, const IMCPrices&) const;
    
	/** interpolate volatility for the log-normal path generator */
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen, int iAsset) const;

};
// end of class VarSwapHedgingSupportProd

/** compute the payoff of the variance swap at maturity */
double VarSwapHedgingSupportProd::payoffVarSwapAtMat(const double* pastSamples, int nberSamples) const {

    static const string routine("VarSwapHedgingSupportProd::payoffVarSwapAtMat");
    try {
        double sumLogSqrReturn = 0.0; // sum of the square log-returns

        int iSample = 0;
        for (iSample = 1 ; iSample < nberSamples ; iSample++)
        {
            sumLogSqrReturn += ::pow(log(pastSamples[iSample]/pastSamples[iSample-1]), 2);
        }
        sumLogSqrReturn /= (nberSamples - 1); // 1-day realized volatility
        sumLogSqrReturn *= inst->observationsPerYear; // annualized realized variance

        double scaleFactor = 100.0 / (inst->dontScaleByStrike ? 1.0 : (2.0 * inst->strikeVol));
        
        return scaleFactor * (sumLogSqrReturn - inst->strikeVol * inst->strikeVol);
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** compute the P&L of the hedging strategy for a given path */
void VarSwapHedgingSupportProd::payoff(const IPathGenerator*  pathGen, IMCPrices& prices) {

	static const string routine("VarSwapHedgingSupportProd::payoff");
	try {
		VarSwapHedgingSupportPrices& myprices = static_cast<VarSwapHedgingSupportPrices&>(prices);

        // path details
        const double* path = pathGen->Path(0,0); // access path
		int pathLength = histDates.size();
        int beginIdx = pathGen->begin(0);
		int endIdx = pathGen->end(0);
        
        // set a flag isDoingPast
        bool doingPast = pathGen->doingPast();
        
        // preservation of the past
		if (!doingPast)
		{
			pathValue = pathValueSoFar; // path value
            plDeltaHedge = plDeltaHedgeSoFar; // P&L related to delta-hedge
            transactionCost = transactionCostSoFar; // transaction cost
        }
        
        // scaling as done for variance swap
        double scaleFactor = 100.0 / (inst->dontScaleByStrike ? 1.0 : (2.0 * inst->strikeVol));
        
        // main loop
        int Idx;
        for (Idx = beginIdx ; Idx < endIdx ; Idx++)
		{
            // at this t = 0 for a short position on variance swap
            if (doingPast && (Idx == 0))
            {
                // compute the initial price and delta of the variance swap
                Results* results = inst->varSwapModel->Run(varSwapMC.get(), control.get());
                double priceVarSwap = results->retrievePrice(); // price of the variance swap
                
                // compute the price of the replicating portfolio
                // double fwd = varSwapMC.get()->getAssetFwdValue(histDates.back());
                double fwd = inst->asset->fwdValue(histDates.back());
                double initSpot = inst->initialSpot;
                double variance, impliedVol;

                DoubleArraySP pricesCallPut = DoubleArraySP(new DoubleArray(inst->strikes.size()));
                double sumPriceCallPut = 0.0;
                int iStrk;
                for (iStrk = 0 ; iStrk < inst->strikes.size() ; iStrk++)
                {
                    // implied volatility incorporating bid-ask
                    impliedVol = (*inst->impliedVols)[iStrk] +
                        (inst->positionLS ? (*inst->bidVols)[iStrk] : (*inst->askVols)[iStrk]);
                    // total variance to maturity
                    variance = yearFrac[0] * pow(impliedVol, 2.0);
                
                    // price of the call or put
                    (*pricesCallPut)[iStrk] = Black::price(inst->callputFlags[iStrk],
                                                           fwd, 
                                                           inst->strikes[iStrk] * initSpot,
                                                           pvFactor[pvFactor.size()-1],
                                                           variance);

                    sumPriceCallPut += ((*inst->hedges)[iStrk]) * ((*pricesCallPut)[iStrk]);
                }
                sumPriceCallPut *= scaleFactor;

                // compute the price of 1/S0 shares
                double priceShare = scaleFactor * (2.0 / yearFrac[0]) *
                    (1.0 + ((inst->positionLS ? -1.0 : 1.0) * (1 / path[0]) * (inst->stockTransCost / 100.0)));

                // compute the price of 1/S0 futures struck at S0
                double notionalFuture = scaleFactor * (2.0 / yearFrac[0]);
                
                double priceFuture = scaleFactor * pvFactor[pvFactor.size()-1] *
                    (2.0 / yearFrac[0]) * ((fwd - initSpot) / initSpot) +
                    (inst->positionLS ? 1.0 : -1.0) *
                    (notionalFuture / 1000000.0) * inst->futureTransCost;
                
                // cost of the strategy at time t = 0;
                pathValue = (inst->positionLS ? -1.0 : 1.0) * growthFactor[Idx] *
                    (priceVarSwap + priceFuture - sumPriceCallPut - priceShare);
            }

            if (!doingPast && (Idx >= 1))
            {
                double transactionCost = 0.0;
                // long position on variance swap
                if (inst->positionLS)
                {
                    transactionCost = ((path[Idx] > path[Idx-1]) ? 1.0 : -1.0) * (inst->stockTransCost / 100.0);
                }
                // short position on variance swap
                else
                {
                    transactionCost = ((path[Idx] > path[Idx-1]) ? -1.0 : 1.0) * (inst->stockTransCost / 100.0);
                }
                
                // buy (2 / (T * St) - 2 / (T * St-1)) shares
                plDeltaHedge += growthFactor[Idx] * scaleFactor *
                    (2.0 / yearFrac[0]) * ((1.0 / path[Idx-1]) - (1.0 / path[Idx])) * (path[Idx] + transactionCost);
            }
            
            // at t = T for a short position on variance swap
            if (!doingPast && (Idx == endIdx - 1))
            {
                // compute the intrasic value of the variance swap
                double priceVarSwap = payoffVarSwapAtMat(path, histDates.size()); // price of the variance swap

                // compute the intrasic of the replicating portfolio
                double initSpot = inst->initialSpot;

                DoubleArraySP pricesCallPut = DoubleArraySP(new DoubleArray(inst->strikes.size()));
                double sumPriceCallPut = 0.0;
                int iStrk;
                for (iStrk = 0 ; iStrk < inst->strikes.size() ; iStrk++)
                {
                    (*pricesCallPut)[iStrk] = (inst->callputFlags[iStrk] ?
                        max(path[Idx] - (inst->strikes[iStrk] * initSpot), 0.0) :
                        max((inst->strikes[iStrk] * initSpot) - path[Idx], 0.0));

                    sumPriceCallPut += ((*inst->hedges)[iStrk]) * ((*pricesCallPut)[iStrk]);
                }
                sumPriceCallPut *= scaleFactor;
                
                // compute the price of 1/ST shares
                double priceShare = scaleFactor * (2.0 / yearFrac[0]) *
                    (1.0 + ((inst->positionLS ? 1.0 : -1.0) * (1 / path[Idx]) * (inst->stockTransCost / 100.0)));

                // compute the price of 1/S0 futures struck at S0
                double priceFuture = scaleFactor * (2.0 / yearFrac[0]) * ((path[Idx] - initSpot) / initSpot);

                // update path value
                pathValue += (inst->positionLS ? -1.0 : 1.0) * growthFactor[Idx] *
                    (-priceVarSwap - priceFuture + sumPriceCallPut + priceShare + plDeltaHedge);
            }
        }

        // preservation of the past 
		if (doingPast)
		{
			pathValueSoFar = pathValue; // path value
            plDeltaHedgeSoFar = plDeltaHedge; // P&L related to delta-hedge
            transactionCostSoFar = transactionCost; // transaction cost
        }

        // preservation of the past
		if (!doingPast)
		{
            myprices.add(pathValue, 0);
            myprices.addPathValue(pathValue, pathGen->getPathIndex());
            myprices.addPathDeltaHedgeValue(plDeltaHedge, pathGen->getPathIndex());
            myprices.addSpotsPath(path, pathLength, pathGen->getPathIndex());
        }
    }
	catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/* interpolate volatility for the log-normal path generator */
CVolRequestLNArray VarSwapHedgingSupportProd::getVolInterp(const IMCPathGenerator* pathGen,
                                                           int                     iAsset) const {
    DateTime imntStartDate = inst->fwdStarting ? inst->startDate : inst->valueDate;
    CVolRequestLNArray reqarr(1); // one interpolation level per asset here

    reqarr[0] = CVolRequestLNSP(
		new LinearStrikeVolRequest(inst->asset->getSpot(), //ATM
								   imntStartDate, 
								   inst->endDate(0),
								   inst->fwdStarting));
    return reqarr;
}

/** output statistics of P&L distribution and P&L distribution itself */
void VarSwapHedgingSupportProd::recordExtraOutput(Control* control, Results* results, const IMCPrices& prices) const {

	const VarSwapHedgingSupportPrices& myprices = static_cast<const VarSwapHedgingSupportPrices&>(prices);

    if (control && control->isPricing()) {
        OutputRequest* request =
			control->requestsOutput(OutputRequest::PL_AVERAGE);
        if (request) {
			double plAverage = myprices.getResults(0); // P&L average
			results->storeRequestResult(request, plAverage * pvFromPaymentDate());
        }
        request = control->requestsOutput(OutputRequest::PL_DISTRIBUTION);
        if (request) {
            DoubleArraySP plDistribution = myprices.getPathValues(); // P&L distribution
            // sort this distribution by increasing order
            // Algorithm::shellSort(*plDistribution);
            // sort this distribution by decreasing order
            // int sizeDist = plDistribution->size();
            // DoubleArraySP plDistributionSorted = DoubleArraySP(new DoubleArray(sizeDist));
            // int iPath = 0;
            // for (iPath = 0 ; iPath < sizeDist ; iPath++)
            // {
            //     (*plDistributionSorted)[iPath] = (*plDistribution)[sizeDist-1-iPath] * pvFromPaymentDate();
            // }
            results->storeRequestResult(request, plDistribution);
        }
        request = control->requestsOutput(OutputRequest::PATHS);
        if (request) {
            DoubleArrayArraySP spotsPaths = myprices.getSpotsPaths(); // spots paths
            results->storeRequestResult(request, spotsPaths);
        }
        request = control->requestsOutput(OutputRequest::PL_DELTA_HEDGE);
        if (request) {
            DoubleArraySP plDeltaHedgeDistribution = myprices.getPathDeltaHedgeValues(); // P&L Delta-hedge distribution
            results->storeRequestResult(request, plDeltaHedgeDistribution);
        }
        
        request = control->requestsOutput(OutputRequest::HEDGES);
        if (request) {
            results->storeRequestResult(request, inst->hedges);
        }
    }
}

/** create a MC payoff product */
IMCProduct* VarSwapHedgingSupport::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says
    // which assets need which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(1)); /* create empty one */
    simSeries->addDates(samples.get()->getDates());

    IRefLevelSP refLevel(IRefLevel::Util::makeTrivialAverage(samples.get()->getFirstDate()));
    // VarSwapHedgingSupport* ptr = const_cast<VarSwapHedgingSupport*>(this);
    return new VarSwapHedgingSupportProd(this, refLevel, simSeries);
}

CClassConstSP const VarSwapHedgingSupport::TYPE = CClass::registerClassLoadMethod(
    "VarSwapHedgingSupport", typeid(VarSwapHedgingSupport), VarSwapHedgingSupport::load);

/** force linker to include this file (avoid having header file) */
extern bool VarSwapHedgingSupportLoad()
{
    return true && VarSwapHedgingSupport::TYPE;
}

DRLIB_END_NAMESPACE
