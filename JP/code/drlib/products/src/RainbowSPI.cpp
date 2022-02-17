//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RainbowSPI.cpp
//
//   Description : Rainbow SPI
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIDynamicBasket.hpp"
#include "edginc/SPILockIn.hpp"
#include "edginc/SyntheticPortfolioInsurance.hpp"
#include "edginc/Algorithm.hpp"

DRLIB_BEGIN_NAMESPACE

class RainbowSPI : public SyntheticPortfolioInsurance {
public:
    static CClassConstSP const TYPE;
    friend class RainbowSPIMC;
    friend class RainbowSPIMCSV;

    RainbowSPI(): SyntheticPortfolioInsurance(TYPE), numSPIs(0) {
        isRainbowSPI = true;
    } // for reflection

    // validate the basic instrument
    void validatePop2Object(){
        static const string routine("RainbowSPI::validatePop2Object");
        if (assetType.size() < 1) {
            throw ModelException(routine, 
                                 "Must supply array of asset types for Rainbow SPI");
        }
        if (rainbowWeights.size() < 1) {
            throw ModelException(routine, 
                                 "Must supply array of rainbow weights for Rainbow SPI");
        }
        if (participations.size() < 1) {
            throw ModelException(routine, 
                                 "Must supply array of asset participations for Rainbow SPI");
        }
        if (rainbowWeights.size() != assetType.size() ||
            rainbowWeights.size() != participations.size()) {
            throw ModelException(routine, 
                                 "Array of asset types must be same length as array of rainbow weights and "
                                 "array of participations for Rainbow SPI");
        }
        if (rainbowWeights.size() == 1) {
            throw ModelException(routine, 
                                 "Must have at least 2 assets for Rainbow SPI");
        }

        // check the weights matrix is consistently sized. There should be a column 
        // for each underlying containing the weight that underlying has for 
        // each asset (SPI or otherwise)
        if (weights.numRows() != assetType.size()) {
            throw ModelException(routine, 
                                 "Matrix of weights must have as many rows as there are rainbow assets");
        }
        if (weights.numCols() != assets.get()->numFactors()) {
            throw ModelException(routine, 
                                 "Matrix of weights must have as many columns as there are underlyings");
        }

        int i;
        // let's have a look at what assets we have and how they group
        for (i = 0; i < assetType.size(); i++) {
            if (assetType[i] != SPI_ASSET_TYPE_SPI &&
 //               assetType[i] != SPI_ASSET_TYPE_ALPHA &&
                assetType[i] != SPI_ASSET_TYPE_NON_SPI) {
                throw ModelException(routine, 
                                     "Asset type[" + Format::toString(i) +"] (" + assetType[i] +
                                     ") is unrecognised. Must be either " +
                                     SPI_ASSET_TYPE_SPI + /*", " + SPI_ASSET_TYPE_ALPHA +*/
                                     " or " + SPI_ASSET_TYPE_NON_SPI);
            }
            if (assetType[i] != SPI_ASSET_TYPE_NON_SPI) {
                numSPIs++;
            }
        }
        if (numSPIs == 0) {
            throw ModelException(routine, "At least one SPI Asset needed for"
                                    " a Rainbow SPI");
        }
        // we're going to insist that the SPI assets come first just to 
        // make our lives easier
        for (i = 0; i < numSPIs; i++) {
            if (assetType[i] == SPI_ASSET_TYPE_NON_SPI) {
                throw ModelException(routine, "All SPI assets must come before"
                        " non SPI assets in the list");
            }
        }

        // for now we don't allow individual SPI assets to interact
        // for this reason we want deductions from baskets to be independent
        // thus we disallow paid fees unless on notional and coupons unless they're fixed
        // Given payoff amounts are rainbowed we also insist rainbow weights sum to 100%
        // ALSO WE HAVEN'T EXPOSED THE ALPHA ASSET YET
        double sumRbowWghts = 0.0;
        for (i = 0; i < rainbowWeights.size(); i++) {
            sumRbowWghts += rainbowWeights[i];
        }
        if (!Maths::equals(sumRbowWghts, 1.0)) {
            throw ModelException(routine, 
                                 "Rainbow weights must add up to 100%");
        }

        if (numSPIs > 1) {
            if (basket->feesSPI->feeType() == SPI_FEES_TYPE_KNP) {
                throw ModelException(routine, 
                                    "Cannot have paid equity fees for Rainbow SPI");
            }

            if (!basket->couponsSPI->getCouponsSPI()->onlyFixedCoupons()) {
                throw ModelException(routine, 
                                    "Rainbow SPI can only have fixed coupons (i.e. minCpn = maxCpn)");
            }
        }

        SyntheticPortfolioInsurance::validatePop2Object();
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below

private:
    RainbowSPI(const RainbowSPI& rhs); // not implemented
    RainbowSPI& operator=(const RainbowSPI& rhs); // not implemented

    static IObject* defaultRainbowSPI(){
        return new RainbowSPI();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RainbowSPI, clazz);
        SUPERCLASS(SyntheticPortfolioInsurance);
        EMPTY_SHELL_METHOD(defaultRainbowSPI);
        FIELD(assetType,           "Asset types SPI or non-SPI");
        FIELD(weights,              "Asset weights");
        FIELD(rainbowWeights,          "Rainbow weights");
        FIELD(participations,        "Participations");
        FIELD(numSPIs,        "How many SPIs are we running");
        FIELD_MAKE_TRANSIENT(numSPIs);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    StringArray             assetType;
    DoubleMatrix            weights;
    DoubleArray             rainbowWeights;
    DoubleArray             participations;
    int                     numSPIs;
};

class RainbowSPIReport : public CObject {
public:
     static CClassConstSP const TYPE;
     int                 historySize; // how many points (helps limit mem alloc)
     DateTimeArray       historySampleDate;
     DoubleArray         historyPayoff;
     DoubleArray         assetLevels;
     SPIReportArraySP    spiReports;

     RainbowSPIReport(int numRbwAssets, SPIReportArraySP reports) :
          CObject(TYPE), spiReports(reports) {
          historySize = (*spiReports)[0]->historySampleDate.size();
          historySampleDate = DateTimeArray(historySize);
          historyPayoff = DoubleArray(historySize, 0.0);
          assetLevels = DoubleArray(numRbwAssets, 0.0);
     } 

     void downsize() {
          int i;
          for (i = 0; i < spiReports->size(); i++) {
               (*spiReports)[i]->downsize();
          }
          historySize = (*spiReports)[0]->historySize;
          historySampleDate.resize(historySize);
          historyPayoff.resize(historySize);
          // note, for last step if we are beyond maturity we have already added
          // the final payoff so we don't want to overwrite it hence the +=
          for (i = 0; i < historySize; i++) {
               historySampleDate[i] = (*spiReports)[0]->historySampleDate[i];
               historyPayoff[i] += (*spiReports)[0]->historyPayoff[i];
          }
     }

     void addMaturityData(double finalPayoff, DoubleArray assetLvls) {
          for (int i = 0; i < assetLevels.size(); i++) {
               assetLevels[i] = assetLvls[i];
          }
          historyPayoff[historySize-1] = finalPayoff;
     }

private:        
     RainbowSPIReport() : CObject(TYPE), historySize(0) {} // for reflection

     RainbowSPIReport(const RainbowSPIReport& rhs); // not implemented
     RainbowSPIReport& operator=(const RainbowSPIReport& rhs); // not implemented

     static IObject* defaultSPIReport() {
          return new RainbowSPIReport();
     }
        
     static void load(CClassSP& clazz) {
          REGISTER(RainbowSPIReport, clazz);
          SUPERCLASS(CObject);
          EMPTY_SHELL_METHOD(defaultSPIReport);
          FIELD(historySampleDate, "historySampleDate");
          FIELD(historyPayoff, "Payoff");
          FIELD(assetLevels, "Final levels of all rainbow assets");
          FIELD(spiReports, "Reports for each SPI asset");;
          FIELD(historySize, "historySize");
          FIELD_MAKE_TRANSIENT(historySize);
          clazz->setPublic(); // make visible to EAS/spreadsheet
     }
};
typedef smartPtr<RainbowSPIReport> RainbowSPIReportSP;

class RainbowSPIMC : public SyntheticPortfolioInsuranceMC {
public:
    // copies of data from instrument
    StringArray             assetType;
    DoubleMatrix            weights;
    DoubleArray             rainbowWeights;
    DoubleArray             participations;
    RainbowSPIReportSP      rainbowReport;

private:
    // a run-time dynamic basket/report for each SPI algorithm
    int                     numSPIs;
    int                     numRbwAssets;
    vector<SPIRunTimeSP>    dynBasksRT;
    SPIReportArraySP        spiReports;
    int                     nbAssets; // convenient
    
    int gapRiskContributor() {
        // note gap risk is never negative
        double maxTotal = -1.0;
        int maxIdx = -1;
        for (int spiIdx = 0; spiIdx < numSPIs; spiIdx++) {
            double total = 0.0;
            for(int i=0; i<dynBasksRT[spiIdx]->gapRiskByBucket.size(); i++) {
                total += dynBasksRT[spiIdx]->gapRiskByBucket[i];
            }
            if (total > maxTotal) {
                maxIdx = spiIdx;
                maxTotal = total;
            }
        }
        return maxIdx;
    }

public:
    RainbowSPIMC(const RainbowSPI* inst, const SimSeriesSP& simSeries) :
        SyntheticPortfolioInsuranceMC(inst, simSeries),
        assetType(inst->assetType),
        weights(inst->weights),
        rainbowWeights(inst->rainbowWeights),
        participations(inst->participations),
        numSPIs(inst->numSPIs),
        numRbwAssets(inst->rainbowWeights.size()),
        dynBasksRT(0),
        nbAssets(getNumAssets()) {
        int numDates = inst->basket->getRebalanceDates().size();

        // build run time object and report for first SPI
        dynBasksRT.push_back(dynBaskRT);
        spiReports = SPIReportArraySP(new SPIReportArray(numSPIs));
        (*spiReports)[0] = SPIReportSP(new SPIReport(numDates));

        // build run time objects and reports for other SPIs
        for (int i = 1; i < numSPIs; i++) {
           dynBasksRT.push_back(SPIRunTimeSP(new SPIRunTime(inst, lockIn, settlement,
                                        &paymentDate, feeNotifDates)));
            // take copies?
            (*spiReports)[i] = SPIReportSP(new SPIReport(numDates));
        }
        rainbowReport = RainbowSPIReportSP(new RainbowSPIReport(numRbwAssets, 
                                                                spiReports));
    }

    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static         const string routine("RainbowSPIMC::payoff");
        int            beginIdx = pathGen->begin(0); // same for all assets
        int            endIdx   = pathGen->end(0);
        PricesSPI&     myPrices = static_cast<PricesSPI&>(prices);
        bool           doingPast = pathGen->doingPast();
        double         strike = Maths::isZero(inst->strike)? 
                                dynBaskRT->lockIn->getInitialLockIn() : inst->strike;
        double         futurePayoff = 0.0; // running total of payoff in future for rainbow
        double         thisFuturePayoff; // running total of payoff in future for single SPI

        int j, spiIdx;
        int startIdx = Maths::max(dynBaskRT->iStepFirstRebal, beginIdx);
        bool finishedEarly = false;

        for (spiIdx = 0; spiIdx < numSPIs && !finishedEarly; spiIdx++) {
            thisFuturePayoff = 0;

            dynBasksRT[spiIdx]->init(startIdx, doingPast);

            for(int iStep=startIdx; iStep<endIdx && !finishedEarly; iStep++) {
                // currently don't do out performance
                // we also ignore publication lag
                dynBasksRT[spiIdx]->E[0] = 0.0;
                for (int j = 0; j < nbAssets; j++) {
                    dynBasksRT[spiIdx]->E[0] += weights[j][spiIdx] *
                            pathGen->Path(j,0)[iStep] / pathGen->refLevel(j, 0);
                }
                
                finishedEarly = dynBasksRT[spiIdx]->singleSPIStep(doingPast, iStep, 
                                                            endIdx, thisFuturePayoff);
                if (inst->debugOn || doingPast) {
                    (*spiReports)[spiIdx]->dumpReportStep(iStep, finishedEarly, 
                                                        dynBasksRT[spiIdx].get());
                }
            }
            // set the futurePayoff. Note each spi must generate the same flows
            // before maturity so no paid fees or unequal coupons
            // should validate against this at construction so
            // if we get no match something has gone badly wrong
            if (spiIdx == 0) {
                futurePayoff = thisFuturePayoff;
            } else {
                if (!Maths::equals(futurePayoff, thisFuturePayoff)) {
                    throw ModelException("RainbowSPIMC::payoff", 
                                "Internal error - component SPIs are generating "
                                "different future fees or coupons");
                }
            }
            // preserve values for past 
            if (doingPast){ 
                // all known things are known by this point
                // note we treat unsettled cash as per future payoff comment above
                if (spiIdx == 0) {
                    unsettledCash = 
                        dynBasksRT[spiIdx]->instKnownCashFlows->getUnpaidValue();
                } else {
                    double thisUnsettledCash = 
                        dynBasksRT[spiIdx]->instKnownCashFlows->getUnpaidValue();
                    if (!Maths::equals(unsettledCash, thisUnsettledCash)) {
                        throw ModelException("RainbowSPIMC::payoff", 
                                    "Internal error - component SPIs are generating "
                                    "different past fees or coupons");
                    }
                }
                dynBasksRT[spiIdx]->setSoFar();
                // for vol interp; only needed if we have started
                if (endIdx>dynBasksRT[spiIdx]->iStepFirstRebal) {
                    dynBasksRT[spiIdx]->Einit[0] = 0.0;
                    for (int j = 0; j < nbAssets; j++) {
                        dynBasksRT[spiIdx]->Einit[0] += weights[j][spiIdx] *
                                pathGen->Path(j,0)[0] / pathGen->refLevel(j, 0);
                    }
                }
            }
        }

        if (!doingPast) {
            // gap risk is a statistic too. Only called for future.
            // we add the contributions from the spi giving the biggest total
            int gapRiskIdx = gapRiskContributor();
            myPrices.addGapRisk(dynBasksRT[gapRiskIdx]->gapRiskByBucket);
        }

        if (!doingPast || !hasFuture()) {
            double finalPayoff = inst->matCashFlow;
            // Price the call option.
            // BL here is the value at mat date, B has been updated for final lag adj
            // With the advent of independent strike we no longer know that max(B,BL) > strike
            // so the extra max(,0) is needed.

            DoubleArray rbowPerfs(numRbwAssets, 0.0);
            DoubleArray assetLevels(numRbwAssets, 0.0);
            for (j = 0; j < numSPIs; j++) {
                double finalBasket = dynBasksRT[j]->sumB / inst->averageOutDates.size();
                assetLevels[j] = Maths::max(finalBasket,dynBasksRT[j]->BL);
                rbowPerfs[j] = participations[j] * (assetLevels[j] - strike);
            }
            for (j = numSPIs; j < numRbwAssets; j++) {
                assetLevels[j] = 0.0; 
                for (int k = 0; k < nbAssets; k++) {
                    assetLevels[j] += weights[k][j] *
                            pathGen->Path(k,0)[endIdx-1] / pathGen->refLevel(k, 0);
                }
                rbowPerfs[j] = participations[j] * (assetLevels[j] - strike);
            }
            Algorithm::shellSort(rbowPerfs);
            for(j = 0; j < numRbwAssets; j++) {
                finalPayoff += rainbowWeights[j] * rbowPerfs[j];
            }
            if(inst->debugOn || doingPast) {
                rainbowReport->addMaturityData(finalPayoff, assetLevels); 
            }
            if (doingPast) {
                // another known cashflow!
                DateTime payDate = settlement->settles(inst->basket->getRebalanceDates()[endIdx-1],0);
                // Nice that SPIKnownFlows internally copes with avoiding a double-count
                // note also we use the object for the first SPI
                // see the comment above 
                dynBasksRT[0]->instKnownCashFlows->addFlow(payDate, finalPayoff * inst->notional);
                unsettledCash = dynBasksRT[0]->instKnownCashFlows->getUnpaidValue();
            } else {
                futurePayoff += finalPayoff;
            }
            futurePayoff *= inst->notional;
            // Add any fee cash flows/coupons which are notified but not yet paid.
            // Note that we provide values at final settlement date.
            myPrices.add(futurePayoff + unsettledCash);
        }
        for (spiIdx = 0; spiIdx < numSPIs; spiIdx++) {
            if (dynBasksRT[spiIdx]->dynBask->hasBFHistory) {
                dynBasksRT[spiIdx]->dynBask->bondFloorToday = 
                    dynBasksRT[spiIdx]->bondFloor->getLevelToday();
            }
        }
    }

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const {
        // Only for pricing run
        if (control->isPricing()) {
            const PricesSPI& myPrices = static_cast<const PricesSPI&>(prices);
            OutputRequest* request;
            request = control->requestsOutput(OutputRequest::SPI_GAP_RISK);
            if (request) {
                results->storeRequestResult(request, myPrices.getGapRisk());
            }
            request = control->requestsOutput(OutputRequest::SPI_GAP_RISK_PROFILE);
            if (request) {
                const CashFlowListSP getGapRiskProfile(myPrices.getGapRiskBucketted());
                results->storeRequestResult(request, getGapRiskProfile);
            }
            request = control->requestsOutput(OutputRequest::SPI_BOND);
            if (request) {
                double bondClosing;
                bool haveValue = dynBaskRT->getBondClosing(bondClosing);
                results->storeRequestResult(request, haveValue?
                                            IObjectSP(CDouble::create(bondClosing)):
                                            IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR_LEVEL);
            if (request) {
                double bondFloorLevel = dynBaskRT->dynBask->bondFloorToday;
                results->storeRequestResult(request, bondFloorLevel > 0.0?
                                            IObjectSP(CDouble::create(bondFloorLevel)):
                                            IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_LOAN_COST_RATE);
            if (request) {
                // This is for processing : the value will be calculated at EOD and written 
                // back to the input in time for the O/N batch.
                CashFlowSP lcr = dynBaskRT->loanCost->getLoanCostRateForPyramid();
                results->storeRequestResult(request, 
                                            !!lcr ? IObjectSP(lcr.get()) : IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_REPORT);
            if (request) {
                // Restrict results to those which are meaningful - perhaps just the past
                // or only those up to any early redemption, for example
                rainbowReport->downsize();
                results->storeRequestResult(request, rainbowReport);
            }
        }
    }
};

class RainbowSPIMCSV : public SyntheticPortfolioInsuranceMCSV {
public:
    // copies of data from instrument
    StringArray             assetType;
    DoubleMatrix            weights;
    DoubleArray             rainbowWeights;
    DoubleArray             participations;
    RainbowSPIReportSP      rainbowReport;

private:
    // a run-time dynamic basket/report for each SPI algorithm
    int                     numSPIs;
    int                     numRbwAssets;
    vector<SPIRunTimeSP>    dynBasksRT;
    SPIReportArraySP        spiReports;
    int                     nbAssets; // convenient
    
    

    int gapRiskContributor() {
        // note gap risk is never negative
        double maxTotal = -1.0;
        int maxIdx = -1;
        for (int spiIdx = 0; spiIdx < numSPIs; spiIdx++) {
            double total = 0.0;
            for(int i=0; i<dynBasksRT[spiIdx]->gapRiskByBucket.size(); i++) {
                total += dynBasksRT[spiIdx]->gapRiskByBucket[i];
            }
            if (total > maxTotal) {
                maxIdx = spiIdx;
                maxTotal = total;
            }
        }
        return maxIdx;
    }

public:
    RainbowSPIMCSV(const RainbowSPI* inst, const SimSeriesSP& simSeries) :
        SyntheticPortfolioInsuranceMCSV(inst, simSeries),
        assetType(inst->assetType),
        weights(inst->weights),
        rainbowWeights(inst->rainbowWeights),
        participations(inst->participations),
        numSPIs(inst->numSPIs),
        numRbwAssets(inst->rainbowWeights.size()),
        dynBasksRT(0),
        nbAssets(getNumAssets()) {
        int numDates = inst->basket->getRebalanceDates().size();

        // build run time object and report for first SPI
        dynBasksRT.push_back(dynBaskRT);
        spiReports = SPIReportArraySP(new SPIReportArray(numSPIs));
        (*spiReports)[0] = SPIReportSP(new SPIReport(numDates));

        // build run time objects and reports for other SPIs
        for (int i = 1; i < numSPIs; i++) {
           dynBasksRT.push_back(SPIRunTimeSP(new SPIRunTimeSV(inst, lockIn, settlement,
                                                              feeNotifDates)));
            // take copies?
            (*spiReports)[i] = SPIReportSP(new SPIReport(numDates));
        }
        rainbowReport = RainbowSPIReportSP(new RainbowSPIReport(numRbwAssets, 
                                                                spiReports));
    }

    void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        static const string routine = "SyntheticPortfolioInsuranceMCSV::pathGenUpdated";
        
        try {
            // parent
            SyntheticPortfolioInsuranceMCSV::pathGenUpdated(newPathGen);
            // then our own
            for(int i=0; i<numSPIs; i++) {
                dynBasksRT[i]->pathGenUpdated(newPathGen);
            }
        } catch (exception& e) {
            throw ModelException(e, routine);
        }
    }

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    void collectStateVars(IStateVariableCollectorSP svCollector) const{
        // parent
        SyntheticPortfolioInsuranceMCSV::collectStateVars(svCollector);
        // then our own
        for(int i=0; i<numSPIs; i++) {
            dynBasksRT[i]->collectStateVars(svCollector);
        }
    }

    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices) {
        static         const string routine("RainbowSPIMCSV::payoff");
        int            beginIdx = spotSV->path(0/*iAsset*/).begin(); // same for all assets
        int            endIdx   = spotSV->path(0/*iAsset*/).end();
        PricesSPI&     myPrices = static_cast<PricesSPI&>(prices);
        bool           doingPast = pathGen->doingPast();
        double         strike = Maths::isZero(inst->strike)? 
                                dynBaskRT->lockIn->getInitialLockIn() : inst->strike;
        double         futurePayoff = 0.0; // running total of payoff in future for rainbow
        double         thisFuturePayoff; // running total of payoff in future for single SPI

        int j, spiIdx;
        int startIdx = Maths::max(dynBaskRT->iStepFirstRebal, beginIdx);
        bool finishedEarly = false;

        for (spiIdx = 0; spiIdx < numSPIs && !finishedEarly; spiIdx++) {
            thisFuturePayoff = 0;

            dynBasksRT[spiIdx]->init(startIdx, doingPast);

            for(int iStep=startIdx; iStep<endIdx && !finishedEarly; iStep++) {
                // currently don't do out performance
                // we also ignore publication lag
                dynBasksRT[spiIdx]->E[0] = 0.0;
                for (int j = 0; j < nbAssets; j++) {
                    dynBasksRT[spiIdx]->E[0] += weights[j][spiIdx] *
                        spotSV->path(j)[iStep] / refLevelSV->refLevel(j);
                }
                
                finishedEarly = dynBasksRT[spiIdx]->singleSPIStep(doingPast, iStep, 
                                                            endIdx, thisFuturePayoff);
                if (inst->debugOn || doingPast) {
                    (*spiReports)[spiIdx]->dumpReportStep(iStep, finishedEarly, 
                                                        dynBasksRT[spiIdx].get());
                }
            }
            // set the futurePayoff. Note each spi must generate the same flows
            // before maturity so no paid fees or unequal coupons
            // should validate against this at construction so
            // if we get no match something has gone badly wrong
            if (spiIdx == 0) {
                futurePayoff = thisFuturePayoff;
            } else {
                if (!Maths::equals(futurePayoff, thisFuturePayoff)) {
                    throw ModelException("RainbowSPIMCSV::payoff", 
                                "Internal error - component SPIs are generating "
                                "different future fees or coupons");
                }
            }
            // preserve values for past 
            if (doingPast){ 
                // all known things are known by this point
                // note we treat unsettled cash as per future payoff comment above
                if (spiIdx == 0) {
                    unsettledCash = 
                        dynBasksRT[spiIdx]->instKnownCashFlows->getUnpaidValue();
                } else {
                    double thisUnsettledCash = 
                        dynBasksRT[spiIdx]->instKnownCashFlows->getUnpaidValue();
                    if (!Maths::equals(unsettledCash, thisUnsettledCash)) {
                        throw ModelException("RainbowSPIMCSV::payoff", 
                                    "Internal error - component SPIs are generating "
                                    "different past fees or coupons");
                    }
                }
                dynBasksRT[spiIdx]->setSoFar();
                // for vol interp; only needed if we have started
                if (endIdx>dynBasksRT[spiIdx]->iStepFirstRebal) {
                    dynBasksRT[spiIdx]->Einit[0] = 0.0;
                    for (int j = 0; j < nbAssets; j++) {
                        dynBasksRT[spiIdx]->Einit[0] += weights[j][spiIdx] *
                            spotSV->path(j)[0] / refLevelSV->refLevel(j);
                    }
                }
            }
        }

        if (!doingPast) {
            // gap risk is a statistic too. Only called for future.
            // we add the contributions from the spi giving the biggest total
            int gapRiskIdx = gapRiskContributor();
            myPrices.addGapRisk(dynBasksRT[gapRiskIdx]->gapRiskByBucket);
        }

        if (!doingPast || !hasFuture()) {
            double finalPayoff = inst->matCashFlow;
            // Price the call option.
            // BL here is the value at mat date, B has been updated for final lag adj
            // With the advent of independent strike we no longer know that max(B,BL) > strike
            // so the extra max(,0) is needed.

            DoubleArray rbowPerfs(numRbwAssets, 0.0);
            DoubleArray assetLevels(numRbwAssets, 0.0);
            for (j = 0; j < numSPIs; j++) {
                double finalBasket = dynBasksRT[j]->sumB / inst->averageOutDates.size();
                assetLevels[j] = Maths::max(finalBasket,dynBasksRT[j]->BL);
                rbowPerfs[j] = participations[j] * (assetLevels[j] - strike);
            }
            for (j = numSPIs; j < numRbwAssets; j++) {
                assetLevels[j] = 0.0; 
                for (int k = 0; k < nbAssets; k++) {
                    assetLevels[j] += weights[k][j] *
                            spotSV->path(k)[endIdx-1] / refLevelSV->refLevel(k);
                }
                rbowPerfs[j] = participations[j] * (assetLevels[j] - strike);
            }
            Algorithm::shellSort(rbowPerfs);
            for(j = 0; j < numRbwAssets; j++) {
                finalPayoff += rainbowWeights[j] * rbowPerfs[j];
            }
            if(inst->debugOn || doingPast) {
                rainbowReport->addMaturityData(finalPayoff, assetLevels); 
            }
            if (doingPast) {
                // another known cashflow!
                DateTime payDate = settlement->settles(inst->basket->getRebalanceDates()[endIdx-1],0);
                // Nice that SPIKnownFlows internally copes with avoiding a double-count
                // note also we use the object for the first SPI
                // see the comment above 
                dynBasksRT[0]->instKnownCashFlows->addFlow(payDate, finalPayoff * inst->notional);
                unsettledCash = dynBasksRT[0]->instKnownCashFlows->getUnpaidValue();
            } else {
                futurePayoff += finalPayoff * dfSV->firstDF();
            }
            futurePayoff *= inst->notional;
            // Add any fee cash flows/coupons which are notified but not yet paid.
            // Note that we provide values at final settlement date.
            myPrices.add(futurePayoff + unsettledCash);
        }
        for (spiIdx = 0; spiIdx < numSPIs; spiIdx++) {
            if (dynBasksRT[spiIdx]->dynBask->hasBFHistory) {
                dynBasksRT[spiIdx]->dynBask->bondFloorToday = 
                    dynBasksRT[spiIdx]->bondFloor->getLevelToday();
            }
        }
    }

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const {
        // Only for pricing run
        if (control->isPricing()) {
            const PricesSPI& myPrices = static_cast<const PricesSPI&>(prices);
            OutputRequest* request;
            request = control->requestsOutput(OutputRequest::SPI_GAP_RISK);
            if (request) {
                results->storeRequestResult(request, myPrices.getGapRisk());
            }
            request = control->requestsOutput(OutputRequest::SPI_GAP_RISK_PROFILE);
            if (request) {
                const CashFlowListSP getGapRiskProfile(myPrices.getGapRiskBucketted());
                results->storeRequestResult(request, getGapRiskProfile);
            }
            request = control->requestsOutput(OutputRequest::SPI_BOND);
            if (request) {
                double bondClosing;
                bool haveValue = dynBaskRT->getBondClosing(bondClosing);
                results->storeRequestResult(request, haveValue?
                                            IObjectSP(CDouble::create(bondClosing)):
                                            IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR_LEVEL);
            if (request) {
                double bondFloorLevel = dynBaskRT->dynBask->bondFloorToday;
                results->storeRequestResult(request, bondFloorLevel > 0.0?
                                            IObjectSP(CDouble::create(bondFloorLevel)):
                                            IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_LOAN_COST_RATE);
            if (request) {
                // This is for processing : the value will be calculated at EOD and written 
                // back to the input in time for the O/N batch.
                CashFlowSP lcr = dynBaskRT->loanCost->getLoanCostRateForPyramid();
                results->storeRequestResult(request, 
                                            !!lcr ? IObjectSP(lcr.get()) : IObjectSP(new NotApplicable()));
            }
            request = control->requestsOutput(OutputRequest::SPI_REPORT);
            if (request) {
                // Restrict results to those which are meaningful - perhaps just the past
                // or only those up to any early redemption, for example
                rainbowReport->downsize();
                results->storeRequestResult(request, rainbowReport);
            }
        }
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* RainbowSPI::createProduct(const MonteCarlo* model) const {
    static const string method("RainbowSPI::createProduct");
    try{
        // we need to create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                     one */
        DateTimeArray instDates = DateTime::merge(lockIn->getLockInSPINoInit()->getEssentialDates(),
                                                  averageOutDates);
        simSeries->addDates(basket->prepareCoarseRebalDates(instDates,
                                                            valueDate));
        
        if (model->stateVarUsed()) {
            return new RainbowSPIMCSV(this, simSeries);
        }
        return new RainbowSPIMC(this, simSeries);
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const RainbowSPI::TYPE = CClass::registerClassLoadMethod(
    "RainbowSPI", typeid(RainbowSPI), RainbowSPI::load);

// * for class loading (avoid having header file) */
bool RainbowSPILoad() {
	return (RainbowSPI::TYPE !=0 &&
                SyntheticPortfolioInsurance::TYPE != 0 &&
                SPIDynamicBasket::TYPE != 0 &&
                SPILockInWrapper::TYPE != 0 &&
                SPIBondWrapper::TYPE != 0 &&
                SPIFeesWrapper::TYPE != 0 &&
                SPIAlgorithmWrapper::TYPE != 0 &&
                SPILoanCostWrapper::TYPE != 0);
}

CClassConstSP const RainbowSPIReport::TYPE = CClass::registerClassLoadMethod(
    "RainbowSPIReport", typeid(RainbowSPIReport), RainbowSPIReport::load);

DRLIB_END_NAMESPACE
