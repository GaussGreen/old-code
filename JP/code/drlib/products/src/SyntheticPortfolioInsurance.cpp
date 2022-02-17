//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SyntheticPortfolioInsurance.cpp
//
//   Description : Synthetic Portfolio Insurance (aka SPI)
//
//   Date        : Jan 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Delta.hpp" // to override gamma
#include "edginc/CrossGamma.hpp" // to override gamma
#include "edginc/SensMgr.hpp"
#include "edginc/SPIDynamicBasket.hpp"
#include "edginc/SPILockIn.hpp"
#include "edginc/SyntheticPortfolioInsurance.hpp"

DRLIB_BEGIN_NAMESPACE

/** Validate instrument having aquired market data */
void SyntheticPortfolioInsurance::Validate(){
    // check the dates first. Some schedules may need to 
    // be built from date builders which require market 
    // data (hols/assets)
    validateDates();

    // need to apply settlement to fees (it may have its own)
    basket->feesSPI->getFeesSPI()->Validate(instSettle.get());
    basket->couponsSPI->getCouponsSPI()->Validate(instSettle.get());

    // validate bond floor - basket needs to know valueDate to check
    // hence we check it here rather than in validatePop2Object
    basket->validateBondFloor(valueDate);

    // now defer to the parent for the rest
    GenericNFBase::Validate();
}

/** Validate dates in instrument having aquired market data 
    in case dates need to be built from date builders whicb
    need market data*/
void SyntheticPortfolioInsurance::validateDates() {
    static const string routine("SyntheticPortfolioInsurance::validateDates");
    try {
        // we know this is going to be a MultiAsset surely?
        const IAsMultiFactors* asMulti = dynamic_cast<const IAsMultiFactors*>(assets.get());
        if (!asMulti){
            throw ModelException("IMultiFactors view not supported by"
                                    " object of type "+
                                    assets.get()->getClass()->getName());
        }
        IMultiFactorsSP mAsset(asMulti->asMultiFactors());
        ObservationSourceArraySP obs(pastValues->getSources(mAsset.get()));
        basket->validateDates(ccyHols.get(), mAsset.get(), *obs, excludeAssetHols);

        // at this stage we can guarantee rebalance dates will be there now
        // Another optional parameter - averageOutDates. If none provided
        // use a single date at the final basket definition date.
        if (averageOutDates.empty()) {
            averageOutDates.push_back(basket->getFinalDate());
        }
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
}

// validate the basic instrument
void SyntheticPortfolioInsurance::validatePop2Object(){
    static const string routine("SyntheticPortfolioInsurance::validatePop2Object");
    GenericNFBase::validatePop2Object();

    try {
        // Using inst settlement for intra-instrument payments means we must forbid
        // a fixed cash settlement date
        if (!CashSettlePeriod::TYPE->isInstance(instSettle.get())) {
            throw ModelException(routine,"SPI only supports instrument settlement of CashSettlePeriod!");
        }
        basket->algorithm->getAlgorithmSPI()->crossValidate(cutoff->getCutoffSPI(), 
                                                            basket->feesSPI->getFeesSPI(),
                                                            isRainbowSPI);
    }
    catch (exception& e) {
        // this catch is to try to prevent vc6 from crashing owing to
        // misoptimisation
        throw ModelException(e, routine);
    }

    if (!Maths::equals(dayCountBasis, 360.0) &&
        !Maths::equals(dayCountBasis, 365.0)) {
        throw ModelException(routine, 
                                "dayCountBasis must be 360 or 365");
    }

    if (Maths::isNegative(strike)) {
        throw ModelException(routine,
                                "Strike (" + Format::toString(strike) +
                                ") cannot be negative");
    } 
    // lock ins and coupons together create a potential circularity so we disallow them
    if (!lockIn->getLockInSPINoInit()->doesNothing() && 
            !basket->couponsSPI->doesNothing() ) {
        throw ModelException(routine,
                                "Cannot specify a lock in and a coupon together");
    }

    // This is since we're introducing an optional param and wish to have
    // decent backwards compatible behaviour
    if (Maths::isZero(strike)) {
        // unfortunately rather verbose and a duplication of code in
        // SyntheticPortfolioInsuranceMC constructor ... XXX
        strike = lockIn->getLockInSPINoInit()->getInitialLockIn();
    }
}
   
void SyntheticPortfolioInsurance::GetMarket(const IModel*          model, 
                const CMarketDataSP    market) {
    try {
        GenericNFBase::GetMarket(model, market); // parent
        IBondSPISP bondSPI(basket->bond->getBondSPI());
        bondSPI->getYieldCurveData(model, market);
        if (!ccyHols.isEmpty()) {
            ccyHols.getData(model, market);
        }
    } catch (exception &e) {
        throw ModelException(e, "SyntheticPortfolioInsurance::GetMarket", "failed");
    }
}

// returns the array of all sampling dates the instrument will ever need
// excluding (possibly) the ref level dates (handled in GenericNFBase)
const DateTimeArray SyntheticPortfolioInsurance::samplingDates() const {
    return basket->getAllRebalanceDates();
}

// for DeltaTPlusN we need to move forward lag many rebalance periods
int SyntheticPortfolioInsurance::deltaOffset(const Holiday* hols) const  {
    DateTime lagDate;

    if (valueDate >= basket->getLastRebalDate()) {
        // before final date we move to final date
        if (valueDate < basket->getFinalDate()) {
            lagDate = basket->getFinalDate();
        }
        else { // just default to Delta T+1
            return 1;  
        }
    }
    else {
        int i;
        const DateTimeArray& rebalDates = basket->getRebalanceDates();

        for(i=1; i<rebalDates.size(); i++) {
            if (rebalDates[i] > valueDate) {
                break;
            }
        }

        // rebalDate[i-1] will be the most recent rebal date (or if we've
        // not started it will be the first valuation date)
        // now shift so we've moved lag many rebal dates
        // shouldn't be able to under/overrun this array given checks above
        int lag = (i == 1 && basket->overrideInitialLag) ? 0 : basket->maxExecLagDays;
        lagDate = rebalDates[i + lag - 1];
    }

    int offset = hols->businessDaysDiff(valueDate, lagDate);
    
    // lagDate is the rebalance date at which the most recent potential
    // rebalance will play out. If, for example, we're valuing SOD and the rebal is EOD 
    // we need to go forward another day otherwise we'll end up before the rebal

    if (valueDate.getTime() < lagDate.getTime()) {
        offset++;
    }

    // can't go backwards - only possible if no lag or override and not daily rebalancing
    return offset > 0 ? offset : 0;
}

//// roll through time (setting historic values)
bool SyntheticPortfolioInsurance::sensShift(Theta* theta){
    const Theta::Util thetaUtil = theta->getUtil(valueDate);
    // use valueDate before it changes
    basket->roll(thetaUtil,
                    discount.get(), 
                    dayCountBasis);

    // store some cached data for rolling bond floor if necessary
    const DateTime& origDate = thetaUtil.getOriginalValueDate();
    const DateTime& newDate = thetaUtil.getNewValueDate();

    if (newDate > origDate && !newDate.equals(origDate, false)) {
        // this is a genuine roll so we need to stop using the Euribor curve
        basket->bond->getBondSPI()->setYC(discount);

        // and cache the yield curve for the bond floor if necessary
        if (basket->hasBFHistory) {
            basket->cachedDFQuot = basket->getDiscountFactorsQuotients(false, discount.get());
            // also store today's bond floor as historic if it's a rebal date as it uses a different curve
            const DateTime& nextDailyRebalDate = basket->getNextRebalanceDate(origDate);
            if (nextDailyRebalDate.getDate()==origDate.getDate()) {
                CashFlowArray *BFHistory = basket->bondFloorHistory.get();
                for (int i = 0; i < BFHistory->size(); i++) {
                    if ((*BFHistory)[i].date.equals(origDate, false)) {
                        (*BFHistory)[i].amount = basket->bondFloorToday;
                        basket->lastHistoricBFIndex++;
                    }
                }
            }
        }
    }

    GenericNFBase::sensShift(theta); // and then call parent's method
    return true; // continue to tweak components which implement Theta
}

// used by the EDR_GET_BOND_FLOOR_HISTORY addin function to retrospectively
// fill in the bond floor history given a set of correct yield curves
// for historic dates
CashFlowArraySP SyntheticPortfolioInsurance::getBondFloorHistory(MonteCarlo* mc,
                                    DateTimeArray baseDates,
                                    YieldCurveArray ycArray) {
    static const string routine = "SyntheticPortfolioInsurance::getBondFloorHistory";

    if (!basket->hasBFHistory) {
        throw ModelException(routine, "Cannot calculate bond floor history for an SPI\n"
                                    "which doesn't have one");                                            
    }

    // check curves fall on succesive rebalance dates in the past
    const DateTimeArray& rebalDates = basket->getAllRebalanceDates();
    int size = ycArray.size();
    if (baseDates[0].getDate() > valueDate.getDate()) {
        throw ModelException(routine, "A curve has been supplied with a base date (" +
                                        baseDates[0].toString() +
                                        ") which is in the future");                                            
    }

    int i;
    int firstRebalIndex = -1;
    for (i = 0; i < rebalDates.size(); i++) {
        if (rebalDates[i].equals(baseDates[0], false)) {
            break;
        }
    }
    if (i >= rebalDates.size()) {
        throw ModelException(routine, "The first yield curve has a base date (" +
                                        baseDates[0].toString() +
                                        ") which is not amongst the rebalance dates");                                            
    }
    firstRebalIndex = i;

    for (i = 1; i < size; i++) {
        if (baseDates[i].getDate() > valueDate.getDate()) {
            throw ModelException(routine, "A curve has been supplied with a base date (" +
                                            baseDates[i].toString() +
                                            ") which is in the future");                                            
        }
        if (!(baseDates[i].equals(rebalDates[i + firstRebalIndex], false))) {
            throw ModelException(routine, "The curve with base date (" +
                                    baseDates[i].toString() +
                                    ") must be on the next rebalance date after the previous yield curve");                                            
        }
    }

    // now we're ready to calculate the bond floors at successive dates
    CControlSP ctrl(new Control(SensitivityArrayConstSP(   ), OutputRequestArrayConstSP(   ),0,""));
    OutputRequestSP bfLevelRequest(new OutputRequest(OutputRequest::SPI_BOND_FLOOR_LEVEL));
    OutputNameSP requestOutName(new OutputName(bfLevelRequest->getRequestName()));
    ctrl->addRequest(bfLevelRequest);
    ResultsSP results(new Results());
    HolidaySP holidays(Holiday::noHolidays());

    int lastDate = valueDate.getDate();
    for (i = 0; i < size; i++) {
        // set the date to the historical date everywhere
        ThetaSP thetaShift(new Theta(baseDates[i].getDate() - lastDate, holidays));
        thetaShift->applyScenario(IObjectSP::attachToRef(this));
        lastDate = baseDates[i].getDate();

        // Clear out BFHistory for given date and validate
        // this resets the lastHistoricIndex flag and stops a memory error
        (*(basket->bondFloorHistory))[i+firstRebalIndex].amount = 0.0;
        basket->validateBondFloor(baseDates[i]);

        // overwrite curve in instrument
        basket->bond->getBondSPI()->setYC(ycArray[i]);

        // get bond floor by calling price function
        mc->Price(this, ctrl.get(), results.get());

        // stick result in the bond floor history before moving on
        double bfLevel = results->retrieveScalarGreek(bfLevelRequest->getPacketName(),
                                                                        requestOutName);
        (*(basket->bondFloorHistory))[i+firstRebalIndex].amount = bfLevel;
    }

    CashFlowArraySP result(basket->bondFloorHistory);

    return result;
}

// we validate the SPI before we call getBondFloorHistory it will fail if it detects zeros in historic BF
// we should make sure any dates we're calculating have something non-zero
void SyntheticPortfolioInsurance::setRequiredBFDatesToValid(DateTimeArray baseDates) {
        
    if (basket->bondFloorHistory.get() && !basket->bondFloorHistory->empty()) {
        const DateTimeArray& rebalDates = basket->getAllRebalanceDates();
        int BFSize = basket->bondFloorHistory->size();
        int j;
        for (int i = 0; i < baseDates.size(); i++) {
            for (j = 0; j < rebalDates.size(); j++) {
                if (rebalDates[j].equals(baseDates[i], false)) {
                    break;
                }
            }
            if (j >= rebalDates.size()) {
                throw ModelException("SyntheticPortfolioInsurance::setRequiredBFDatesToValid", 
                                        "Curve has a base date  (" +
                                                baseDates[i].toString() +
                                                ") which is not amongst the rebalance dates");                                            
            }
            if (j < BFSize) {
                (*(basket->bondFloorHistory))[j].amount = 1.0;
            }
        }
    }
}

// implementation of TargetRedemption::IEventHandler interface
void SyntheticPortfolioInsurance::getEvents(const TargetRedemption* tarn, IModel* model, 
                const DateTime& eventDate, EventResults* events) const {
    static const string method = "SyntheticPortfolioInsurance::getEvents";

    // check whether we have TARN feature at all
    if (basket->couponsSPI->getCouponsSPI()->hasEarlyRedemption()) {
        try {
            MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
            if (mc) {
                auto_ptr<IMCProduct> prod(createProduct(mc));
                MCPathGeneratorSP past = prod->runPast(mc);
                prod->retrieveEvents(events);
            } else {
                throw ModelException(method, 
                        "Internal error - expected Monte Carlo model for SPI pricing");
            }
        } catch (exception& e) {
            throw ModelException(e, method);
        }
    }
}

// implementation of SPIFixing::IEventHandler interface
void SyntheticPortfolioInsurance::getEvents(const SPIFixing* fixing, IModel* model, 
                const DateTime& eventDate, EventResults* events) const {
    static const string method = "SyntheticPortfolioInsurance::getEvents";

    try {
        int i = -1;
        const DateTimeArray& rebalDates = basket->getAllRebalanceDates();
        try {
            i = eventDate.find(rebalDates);
        } catch (exception&) {
            // only do fixings if this is a rebalance date
            return;
        }

        int numDates = rebalDates.size();

        // add ZC bond closing event
        events->addEvent(new SPIFixing(eventDate, SPIFixing::ZC_BOND_FIXING, 
                                basket->bond->getBondSPI()->getBondPriceToday(eventDate)));

        // add MM rate event - make it contingent on NEEDING it?
        // note in call to getAccrueFactor we ignore discount and use correct curve
        if (i < numDates - 1) {
            ILoanCostSPISP loanCost = basket->loanCost->getLoanCostSPI();
            loanCost->init(eventDate, basket->bond->getBondSPI()->getYC(),
                           dayCountBasis, rebalDates);
            loanCost->getAccrueFactor(discount.get(), eventDate,
                                        rebalDates[i+1], dayCountBasis);
            CashFlowSP lcr = loanCost->getLoanCostRateForPyramid();
            if (lcr.get()) {
                events->addEvent(new SPIFixing(eventDate, 
                                               SPIFixing::MM_RATE_FIXING, 
                                               lcr->amount));
            }
        }

        if (basket->hasBFHistory) {
            // add BF fixing event if needed
            // note BF is path dependent so must run past first
            MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
            if (mc) {
                auto_ptr<IMCProduct> prod(createProduct(mc));
                MCPathGeneratorSP past = prod->runPast(mc);
                events->addEvent(new SPIFixing(eventDate, 
                                            SPIFixing::BOND_FLOOR_FIXING, 
                                            basket->bondFloorToday));
            } else {
                throw ModelException(method, 
                        "Internal error - expected Monte Carlo model for SPI pricing");
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

SyntheticPortfolioInsurance::SyntheticPortfolioInsurance(CClassConstSP clazz): GenericNFBase(clazz), averageOutDates(0),
        strike(0.0), matCashFlow(0.0), ccyHols(0), excludeAssetHols(false),  
        isRainbowSPI(false), debugOn(false)/*, weightedBasket(false)*/ {} // for reflection

SyntheticPortfolioInsurance::SyntheticPortfolioInsurance(): GenericNFBase(TYPE), averageOutDates(0),
        strike(0.0), matCashFlow(0.0), ccyHols(0), excludeAssetHols(false), 
        isRainbowSPI(false), debugOn(false)/*, weightedBasket(false)*/ {} // for reflection

IObject* SyntheticPortfolioInsurance::defaultSyntheticPortfolioInsurance(){
    return new SyntheticPortfolioInsurance();
}

/** Invoked when Class is 'loaded' */
void SyntheticPortfolioInsurance::load(CClassSP& clazz){
    REGISTER(SyntheticPortfolioInsurance, clazz);
    SUPERCLASS(GenericNFBase);
    IMPLEMENTS(IMCIntoProduct);
    //IMPLEMENTS(IDeltaTPlusN);
    IMPLEMENTS(DeltaTPlusN::IDeltaTPlusNImnt);
    IMPLEMENTS(TargetRedemption::IEventHandler);
    IMPLEMENTS(SPIFixing::IEventHandler);
    EMPTY_SHELL_METHOD(defaultSyntheticPortfolioInsurance);
    FIELD(basket,                      "Handle to Dynamic Basket");
    FIELD(averageOutDates,      "Basket averaging out dates");
    FIELD_MAKE_OPTIONAL(averageOutDates);
    FIELD(lockIn,                      "Handle to LockIn");
    FIELD(dayCountBasis,        "Day Count Basis");
    FIELD(cutoff,                      "Handle to Cutoff");
    FIELD(isCall,               "True=>Call, else Put");
    FIELD(strike,               "Strike");
    FIELD_MAKE_OPTIONAL(strike);
    FIELD(matCashFlow,          "Payoff=Call or Put + matCashFlow");
    FIELD_MAKE_OPTIONAL(matCashFlow);
    FIELD(ccyHols, "Currency (bond floor) holidays");
    FIELD_MAKE_OPTIONAL(ccyHols);
    FIELD(excludeAssetHols, "Should we exclude asset holidays from the rebalance dates?");
    FIELD_MAKE_OPTIONAL(excludeAssetHols);
    FIELD(isRainbowSPI,         "Is this a rainbow SPI");
    FIELD_MAKE_TRANSIENT(isRainbowSPI);
    FIELD(gapRiskBucketOffsets, "Offsets for bucketting of Gap Risk Profile, e.g. 1M");
    FIELD_MAKE_OPTIONAL(gapRiskBucketOffsets);
    FIELD(debugOn,              "True=>enables Debug sheet with variable values from final iteration");
    FIELD_MAKE_OPTIONAL(debugOn);
//    FIELD(tempWeights,              "weights");
//    FIELD_MAKE_OPTIONAL(tempWeights);
//    FIELD(weightedBasket,              "True=>enables us to have weighted basket as the asset");
//    FIELD_MAKE_OPTIONAL(weightedBasket);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

/*****************************************************************************/
SPIGapEventStatistics::SPIGapEventStatistics() : CObject(TYPE),
        numGapEvents(0), gapEventProbability(0.0), expectedConditionalValue(0.0) {}

void SPIGapEventStatistics::evaluateFinalLevel(bool isCall, double finalBasket, 
                                            double finalLockIn, double strike){
    if (isCall) {
        double effStrike = Maths::max(finalLockIn, strike);
        if (Maths::isNegative(finalBasket - effStrike)){
            numGapEvents++;
            expectedConditionalValue += (finalBasket - effStrike);
        }
    } else {
        if (Maths::isNegative(strike - finalBasket)) {
            numGapEvents++;
            expectedConditionalValue += strike - finalBasket;
        }
    }
}

void SPIGapEventStatistics::reset() {
    gapEventProbability = 0.0;
    expectedConditionalValue = 0.0;
    numGapEvents = 0;
}

void SPIGapEventStatistics::finalise(int numIters, double discFactor) {
    gapEventProbability = (double) numGapEvents/numIters;
    expectedConditionalValue *= discFactor;
    expectedConditionalValue /= ((numGapEvents > 0) ? numGapEvents : 1.0);    
}

// The CombinableResult interface
/** scale by factor x */
void SPIGapEventStatistics::scale(double x) {
    gapEventProbability *= x;
}

/** add an object to this result */
void SPIGapEventStatistics::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const SPIGapEventStatistics& resultToAdd = 
            dynamic_cast<const SPIGapEventStatistics&>(static_cast<const IObject&>(x));

    // gapEventProbablity works in the usual STATISTICAL way
    gapEventProbability += scaleFactor * resultToAdd.gapEventProbability;           

    // a little bit of cheating for the conditional expectation
    expectedConditionalValue *= numGapEvents;
    expectedConditionalValue += (resultToAdd.expectedConditionalValue *
                                 resultToAdd.numGapEvents);
    numGapEvents += resultToAdd.numGapEvents;
    expectedConditionalValue /= ((numGapEvents > 0) ? numGapEvents : 1.0);    
}

IObject* SPIGapEventStatistics::defaultStatistics(){
    return new SPIGapEventStatistics();
}

/** Invoked when Class is 'loaded' */
void SPIGapEventStatistics::load(CClassSP& clazz){
    REGISTER(SPIGapEventStatistics, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultStatistics);
    IMPLEMENTS(CombinableResult);
    FIELD(numGapEvents, "Number of times basket ends 'out of the money'");
    FIELD_MAKE_TRANSIENT(numGapEvents);
    FIELD(gapEventProbability, "Probability of ending 'out of the money'");
    FIELD(expectedConditionalValue, "Expected final basket value given it ends out of the money");
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const SPIGapEventStatistics::TYPE = CClass::registerClassLoadMethod(
    "SPIGapEventStatistics", typeid(SPIGapEventStatistics), 
    SPIGapEventStatistics::load);

/***************************************************************************/

// Provides gap risk too
// and expect cutoff info
PricesSPI::PricesSPI(int                    NbIter,
              int                    NbSubSamples,
              const DateTimeArray&   gapRiskBucketDates):
        MCPricesSimple(NbIter, NbSubSamples),
        sumCutoffDate(0), numCutoff(0),
        gapRiskBuckets(gapRiskBucketDates.size(), 0.0),
        gapRiskBucketDates(gapRiskBucketDates) {
    gapEventStats = SPIGapEventStatisticsSP(new SPIGapEventStatistics());
}

/** Clears out SumSubPrices and resets iSubSample */
void PricesSPI::reset() {
    MCPricesSimple::reset();
    for(int iBucket=0; iBucket<gapRiskBuckets.size(); iBucket++) {
        gapRiskBuckets[iBucket] = 0.0;
    }
    sumCutoffDate = 0;
    numCutoff = 0;
    gapEventStats->reset();
}

int PricesSPI::storagePerPath(IMCProduct* product) const {
    return 0; // no caching
}
void PricesSPI::configureCache(const IntArray& changedAssets) {
    ; // empty - no caching
}

void PricesSPI::addGapRisk(const DoubleArray& gapRiskByBucket) {
    // bucketting handled in payoff - more efficient there
    for(int i=0; i<gapRiskBuckets.size(); i++) {
        gapRiskBuckets[i] += gapRiskByBucket[i];
    }
}

void PricesSPI::addCutoffData(bool isCutoff,
                    const DateTime& when) {
    if (isCutoff) {
        sumCutoffDate += when.getDate();
        numCutoff++;
    }
}

void PricesSPI::addGapEventStats(bool isCall, double finalBasket,
                                 double lockIn, double strike) {
    gapEventStats->evaluateFinalLevel(isCall, finalBasket, lockIn, strike);
}

double PricesSPI::getGapRisk() const {
    double gapRisk = 0.0;
    for (int i=0; i<gapRiskBuckets.size(); i++) {
        gapRisk += gapRiskBuckets[i];
    }
    return gapRisk/NbIter;
}

CashFlowListSP PricesSPI::getGapRiskBucketted() const {
    CashFlowArraySP gapRiskBucketResult(new CashFlowArray(gapRiskBucketDates.size()));
    for(int i=0; i<gapRiskBucketDates.size(); i++) {
        (*gapRiskBucketResult)[i].date = gapRiskBucketDates[i];
        (*gapRiskBucketResult)[i].amount = gapRiskBuckets[i]/NbIter;
    }
    CashFlowListSP gapRiskBuckets(new CashFlowList(gapRiskBucketResult.get()));
    return gapRiskBuckets;
}

double PricesSPI::getExpectedCutoff() const {
    return (double)numCutoff / (double)NbIter;
}

DateTimeSP PricesSPI::getExpectedCutoffDate() const {
    int expCutoffDate = (numCutoff>0 ? sumCutoffDate/numCutoff : 0);
    return DateTimeSP(new DateTime(expCutoffDate, DateTime::START_OF_DAY_TIME));
}

SPIGapEventStatisticsSP PricesSPI::getGapEventStats(double discFactor) const {
    gapEventStats->finalise(NbIter, discFactor);
    //take copy o/w it might get trampled on by the Greeks
    return SPIGapEventStatisticsSP(copy(gapEventStats.get()));
}

PricesSPI::~PricesSPI() {}

IMCPrices* PricesSPI::emptyConstructor() const{
    DateTimeArray   gapRiskBucketDates(0);
    return new PricesSPI(1, 1, gapRiskBucketDates);
}

/*****************************************************************************/

SPIReport::SPIReport(int numDates): 
     CObject(TYPE), historySize(0) {
    historySampleDate = DateTimeArray(numDates);
    historyDynBasket = DoubleArray(numDates, 0.0);
    historyBond = DoubleArray(numDates, 0.0);
    historyBF = DoubleArray(numDates, 0.0);
    historyIsRebal = BoolArray(numDates);
    historyIsCutoff = BoolArray(numDates);
    historyUE = DoubleArray(numDates, 0.0);
    historyUC = DoubleArray(numDates, 0.0);
    historySE = DoubleArray(numDates, 0.0);
    historyTE = DoubleArray(numDates, 0.0);
    historySC = DoubleArray(numDates, 0.0);
    historyBL = DoubleArray(numDates, 0.0);
    historynZ = DoubleArray(numDates, 0.0);
    historynE0 = DoubleArray(numDates, 0.0);
    historynE1 = DoubleArray(numDates, 0.0);
    historyE0 = DoubleArray(numDates, 0.0);
    historyE1 = DoubleArray(numDates, 0.0);
    historyA0 = DoubleArray(numDates, 0.0);
    historyA1 = DoubleArray(numDates, 0.0);
    historyLB = DoubleArray(numDates, 0.0);
    historyAccumulatedFee = DoubleArray(numDates, 0.0);
    historyFeeAmount = DoubleArray(numDates, 0.0);
    historyCouponAmt = DoubleArray(numDates, 0.0);
    historyCumDynBasket = DoubleArray(numDates, 0.0);
    historyCumSE = DoubleArray(numDates, 0.0);
    historyCumTE = DoubleArray(numDates, 0.0);
    historyPayoff = DoubleArray(numDates, 0.0);
    historyOutPerfAsset0 = DoubleArray(numDates, 0.0);
    historyOutPerfAsset1 = DoubleArray(numDates, 0.0);
    historyRC0 = DoubleArray(numDates, 0.0);
    historyRC1 = DoubleArray(numDates, 0.0);
}

void SPIReport::downsize() {
    historySampleDate.resize(historySize);
    historyDynBasket.resize(historySize);
    historyBond.resize(historySize);
    historyBF.resize(historySize);
    historyIsRebal.resize(historySize);
    historyIsCutoff.resize(historySize);
    historyUE.resize(historySize);
    historyUC.resize(historySize);
    historySE.resize(historySize);
    historyTE.resize(historySize);
    historySC.resize(historySize);
    historyBL.resize(historySize);
    historynZ.resize(historySize);
    historynE0.resize(historySize);
    historynE1.resize(historySize);
    historyE0.resize(historySize);
    historyE1.resize(historySize);
    historyA0.resize(historySize);
    historyA1.resize(historySize);
    historyLB.resize(historySize);
    historyAccumulatedFee.resize(historySize);
    historyCouponAmt.resize(historySize);
    historyFeeAmount.resize(historySize);
    historyCumDynBasket.resize(historySize);
    historyCumSE.resize(historySize);
    historyCumTE.resize(historySize);
    historyPayoff.resize(historySize);
    historyOutPerfAsset0.resize(historySize);
    historyOutPerfAsset1.resize(historySize);
    historyRC0.resize(historySize);
    historyRC1.resize(historySize);
}

void SPIReport::dumpReportStep(int iStep, 
                               bool redeemedEarly,
                               SPIRunTime* rtSPI) {
    int flag = rtSPI->sampleFlags[iStep];
    if (redeemedEarly) {
        // arbitrary ... not yet assigned this step
        historyIsRebal[iStep] = false;
        historyIsCutoff[iStep] = false;
        historyUE[iStep] = 0.;
        historyUC[iStep] = 0.;
        historySE[iStep] = 0.;
        historyTE[iStep] = 0.;
        historySC[iStep] = 0.;
        historynZ[iStep] = 0.;
        historynE0[iStep] = 0.;
        if (rtSPI->numAlgAssets>1) {
            historynE1[iStep] = 0.;
        }
    } else {
        historyIsRebal[iStep] = rtSPI->isRebal;
        historyIsCutoff[iStep] = rtSPI->isCutoff;
        historyUE[iStep] = rtSPI->UE;
        historyUC[iStep] = rtSPI->UC;
        historySE[iStep] = rtSPI->SE;
        historyTE[iStep] = rtSPI->TE;
        historySC[iStep] = rtSPI->SC;
        historynZ[iStep] = rtSPI->nZ;
        historynE0[iStep] = rtSPI->nE[0];
        if (rtSPI->numAlgAssets>1) {
            historynE1[iStep] = rtSPI->nE[1];
        }
    }

    historySize = iStep+1; // will then end up being num steps in last run
    historySampleDate[iStep] = rtSPI->dynBask->getRebalanceDates()[iStep];
    historyBF[iStep] = rtSPI->BF;
    historyBL[iStep] = rtSPI->BL;
    historyDynBasket[iStep] = rtSPI->B;
    historyBond[iStep] = rtSPI->Z(iStep);
    historyE0[iStep] = rtSPI->E[0];
    if (flag & FEE_NOTIFICATION || redeemedEarly) {
        historyAccumulatedFee[iStep] = rtSPI->paidFeeToPay;
        historyPayoff[iStep] = rtSPI->paidFeeToPay;
    } else {
        historyAccumulatedFee[iStep] =  rtSPI->accumulatedFeeToPay;
        historyPayoff[iStep] = 0.0;
    }
    historyFeeAmount[iStep] = rtSPI->currentFeeAmount;
    historyA0[iStep] = rtSPI->A[0];
    historyLB[iStep] = rtSPI->LB;
    historyRC0[iStep] = rtSPI->RC[0];
    if (rtSPI->numAlgAssets>1) {
        historyE1[iStep] = rtSPI->E[1];
        historyA1[iStep] = rtSPI->A[1];
        historyRC1[iStep] = rtSPI->RC[1];
    }
    historyCumSE[iStep] = rtSPI->couponsRT->getCumSE();
    historyCumTE[iStep] = rtSPI->couponsRT->getCumTE();
    historyCouponAmt[iStep] = rtSPI->couponAmt;
    historyPayoff[iStep] += rtSPI->couponAmt;
    historyPayoff[iStep] += rtSPI->payoff;
    // crap behaviour but it's what it used to do and it would break tests to change it
    double cumB = rtSPI->couponsRT->getCumB();
    if (Maths::isZero(cumB)) {
        historyCumDynBasket[iStep] = rtSPI->B;
    } else {
        historyCumDynBasket[iStep] = cumB;
    }
}

SPIReport::SPIReport() : CObject(TYPE), historySize(0) {}

IObject* SPIReport::defaultReport(){
    return new SPIReport();
}

/** Invoked when Class is 'loaded' */
void SPIReport::load(CClassSP& clazz){
    REGISTER(SPIReport, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultReport);
    FIELD(historySampleDate,                 "historySampleDate");
    FIELD(historyDynBasket, "Basket");
    FIELD(historyBond, "Bond");
    FIELD(historyBF, "Bond Floor");
    FIELD(historyIsRebal, "Rebalance?");
    FIELD(historyIsCutoff, "Cutoff?");
    FIELD(historyUE, "Unbalanced Exposure");
    FIELD(historyUC, "Unbalanced Crash");
    FIELD(historySE, "Sustainable Exposure");
    FIELD(historyTE, "Target Exposure");
    FIELD(historySC, "Sustainable Crash");
    FIELD(historyBL, "Locked-In ");
    FIELD(historynZ, "Bond Alloc");
    FIELD(historynE0, "Asset#1 Alloc");
    FIELD(historynE1, "Asset#2 Alloc");
    FIELD(historyE0, "Asset#1 Level");
    FIELD(historyE1, "Asset#2 Level");
    FIELD(historyA0, "Asset#1 Lag Adj");
    FIELD(historyA1, "Asset#2 Lag Adj");
    FIELD(historyLB, "Loan Balance");
    FIELD(historyAccumulatedFee, "Accumulated Paid Fee");
    FIELD(historyFeeAmount, "Fee Amount");
    FIELD(historyCouponAmt, "Coupon");
    FIELD(historyCumDynBasket, "Cum Basket");
    FIELD(historyCumSE, "Cum Sustainable Exposure");
    FIELD(historyCumTE, "Cum Target Exposure");
    FIELD(historyOutPerfAsset0, "Outperformance Asset#1 Level");
    FIELD(historyOutPerfAsset1, "Outperformance Asset#2 Level");
    FIELD(historyPayoff, "Payoff");
    FIELD(historyRC0, "Asset#1 Rebalance Cost");
    FIELD(historyRC1, "Asset#2 Rebalance Cost");
    FIELD(historySize, "historySize");
    FIELD_MAKE_TRANSIENT(historySize);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

/*****************************************************************************/

/** equivalent to InstIntoMCProduct. Need to call parent's constructor */
SyntheticPortfolioInsuranceMC::SyntheticPortfolioInsuranceMC(const SyntheticPortfolioInsurance*   inst,
                                const SimSeriesSP&                   simSeries):
    IMCProduct(inst->assets.get(),
                inst->valueDate,
                inst->discount.get(),
                inst->refLevel,
                simSeries,
                inst->pastValues,
                inst->instSettle.get(),
                simSeries->getLastDate()), // sets payment date
    unsettledCash(0.),
    today_haveValues(false), today_dynBasket(0.0),
    today_Z(0.0), today_nZ(0.0), today_E(0),
    today_nE(0), today_UE(0.0), today_SE(0.0),
    today_TE(0.0), today_BF(0.0), today_BL(0.0), 
    haveDebugITE(false), debugITE(0.0),
    inst(inst) {

    static const string routine("SyntheticPortfolioInsuranceMC::SyntheticPortfolioInsuranceMC");
    try{
        int             i;
        int             numDates = inst->basket->getRebalanceDates().size();
        const DateTime& lastRebalDate = inst->basket->getLastRebalDate();
        
        lockIn = inst->lockIn->getLockInSPI(inst->basket->getRebalanceDates(),
                                            lastRebalDate);
        
        // locate fee payment dates. Helps to have the actual rebalance dates (at least as far as the back is 
        // concerned).
        DateTimeArray trueRebalDates(inst->basket->getRebalanceDates().size()-inst->basket->maxPubLagDays-inst->basket->maxExecLagDays);
        for(i=0; i<trueRebalDates.size(); i++) {
            trueRebalDates[i] = inst->basket->getRebalanceDates()[i+inst->basket->maxPubLagDays];
        }
        feeNotifDates = inst->basket->feesSPI->getFeesSPI()->getNotificationDates(trueRebalDates);
        if (feeNotifDates.size()>0 &&
            !feeNotifDates.back().equals(lastRebalDate)) {
            throw ModelException(routine,
                                    "Last fee payment date (" + feeNotifDates.back().toString() + 
                                    ") must match final rebalance date (" + lastRebalDate.toString() +
                                    ")");
        }
        
        // For extra outputs (PAYMENT_DATES & KNOWN_CASHFLOWS) which 
        DateTimeArray feePaymentDates = inst->basket->feesSPI->getFeesSPI()->getPaymentDates();

        // plus the maturity flow - "merge" it in (don't duplicate)
        const DateTime& matSettleDate = settlement->settles(inst->basket->getRebalanceDates().back(), 
                                                            0); // asset optional
        if (feePaymentDates.size()<1 ||
            matSettleDate>feePaymentDates.back()) {
            feePaymentDates.push_back(matSettleDate);
        }
        // And any coupons
        instPaymentDates = DateTime::merge(feePaymentDates,
                                            inst->basket->couponsSPI->getCouponsSPI()->getPaymentDates());
        
        // average dates : none until rebalancing has started and final date must be last average date
        if (!DateTime::isSubset(inst->basket->getRebalanceDates(),inst->averageOutDates)) {
            throw ModelException(routine,
                                    "Average out dates must be a subset of the rebalance dates");
        }
        if (inst->averageOutDates.back() != inst->basket->getRebalanceDates().back()) {
            throw ModelException(routine,
                                    "Final average out date must be the same as the final entered 'rebalance' date");
        }
        if (inst->averageOutDates[0] < trueRebalDates[0]) {
            throw ModelException(routine,
                                    "Cannot average out before rebalancing has started!");
        }
        dynBaskRT = SPIRunTimeSP(new SPIRunTime(inst, lockIn, settlement,
                                        &paymentDate, feeNotifDates));

        numSimAssets = getNumAssets();
        if (!inst->isRainbowSPI) {
            // We may simulate more assets than we run through the algorithm. Currently
            // cope with 1 or 2 risky assets in the SPI algorithm, but can allow any number
            // of simulated assets, so long as we know how to aggregate them into 1 or 2
            // "algo risky assets"
            // For now we support only a single risky asset and an outperformance 
            // May want a nicer way to identify this ...
            if (numSimAssets != dynBaskRT->numAlgAssets) {
                if (/*!inst->weightedBasket && */numSimAssets != 2) {
                    throw ModelException(routine,
                                            "SPI currently supports at most 2 simulated assets, but " + 
                                            Format::toString(numSimAssets) + " supplied");
                }
                if (dynBaskRT->numAlgAssets > 1) {
                    throw ModelException(routine,
                                            "Simulating 2 assets is only supported for a 1-factor algorithm");
                }
            }
        }

        // Used to be "if (inst->debugOn) {"
        // but now being used for client vals so always allocate, though in payoff
        // don't always write (only for past and if debugOn)
        report = SPIReportSP(new SPIReport(numDates));

        double IB = inst->basket->initialBasketLevel;
        double IBF = dynBaskRT->bondFloor->getLevel(lockIn->getInitialLockIn(),
                                            dynBaskRT->iStepFirstRebal);
        double ISC = dynBaskRT->algo->sustCrash(1.-IBF/IB, dynBaskRT->iStepFirstRebal);
        double ISE = Maths::max((1.-IBF/IB)/ISC, 0.);
        debugITE = dynBaskRT->algo->targetExp(ISE);
        haveDebugITE = true;

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

// override hasFuture so that early termination can be captured in the past
bool SyntheticPortfolioInsuranceMC::hasFuture() const {
    return IMCProduct::hasFuture() && !dynBaskRT->terminatedEarly;
}

// Satisfy IHandlePaymentEvents interface
// for SPI specific record of events (pay dates, cashflows etc)
void SyntheticPortfolioInsuranceMC::recordEvents(Control* control,
                    Results* results) {
    static const string method("SyntheticPortfolioInsuranceMC::recordEvents");
    try {
        // PAYMENT_DATES is a list of all dates on which payments may occur
        // including past and potential future dates.
        // For SPI this means all paid fee payment dates and the final 
        // maturity flow
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            OutputRequestUtil::recordPaymentDates(control,results,&instPaymentDates);
        }
        
        // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
        // and be supplied for all past cash flows, and any future ones
        // that are determined.
        // For SPI this means any past paid fees, and possibly a current 
        // fee due to be paid today, plus possibly the final option flow
        // XXX Need to check which ccy this is in.
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished() && 
            dynBaskRT->instKnownCashFlows->getFlows()) {
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    dynBaskRT->instKnownCashFlows->getFlows()); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

IMCPrices* SyntheticPortfolioInsuranceMC::createOrigPrices(int  nbIter,
                                    int  nbSubSamples,
                                    int  mode) {
    // XXX don't support anything other than very simple "mode" so
    // XXX should I validate against mode!=0?
    return new PricesSPI(nbIter, nbSubSamples, dynBaskRT->gapRiskBucketDates);
}

/** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
void SyntheticPortfolioInsuranceMC::initialiseLN(const  MCPathGenerator*  pathGen)const{
    // empty
}

/** Called within the simulation loop */
void SyntheticPortfolioInsuranceMC::payoff(const MCPathGenerator*  pathGen,
            IMCPrices&                prices) {
    static         const string routine("SyntheticPortfolioInsuranceMC::payoff");
    int            beginIdx = pathGen->begin(0); // same for all assets
    int            endIdx   = pathGen->end(0);
    PricesSPI&     myPrices = static_cast<PricesSPI&>(prices);
    bool           doingPast = pathGen->doingPast();
    double         strike = Maths::isZero(inst->strike)? 
                                dynBaskRT->lockIn->getInitialLockIn() : inst->strike;
    bool           outPerf = false;
    double         perf1 = 0.0, perf2 = 0.0;
    double         futurePayoff = 0.0; // running total of payoff in future

    // First few steps are required just to collect equity levels for subsequent 
    // publication lag offset use.
    // This should work even with the call for the past, because the beginIdx/endIdx are
    // based on an absolute scale.
    int startIdx = Maths::max(dynBaskRT->iStepFirstRebal, beginIdx);
    
    dynBaskRT->init(startIdx, doingPast);

    // May want a nicer way to identify this ...
    if (numSimAssets != dynBaskRT->numAlgAssets) {
        outPerf = true;
    }

    bool finishedEarly = false;
    for(int iStep=startIdx; iStep<endIdx && !finishedEarly; 
                                            iStep++) {
        for(int iAsset=0; iAsset<dynBaskRT->numAlgAssets; iAsset++) {
            // Account for any publication lag here, so E[] is implicitly lagged for
            // the pub lag before anyone else has a look-in.
            // Note the lag can be different for different assets
            int lag = inst->basket->getPubLag(iAsset, iStep);

            if (!outPerf) {
                dynBaskRT->E[iAsset] = pathGen->Path(iAsset,0)[iStep-lag] / pathGen->refLevel(iAsset, 0); 
            } else {
                // Need to collapse several simulated assets into 1 or 2 risky assets...
                // For now we support a single risky asset and an outperformance - validated earlier
//                if (inst->weightedBasket) {
//                    dynBaskRT->E[iAsset] = 0.0;
//                    for (int k = 0; k < inst->tempWeights.size(); k++) {
//                        dynBaskRT->E[iAsset] += inst->tempWeights[k] *
//                            pathGen->Path(/*iAsset*/k,0)[iStep-lag] / pathGen->refLevel(/*iAsset*/k, 0);
//                    }
//                } else {
                    perf1 = pathGen->Path(/*iAsset*/0,0)[iStep-lag] / pathGen->refLevel(/*iAsset*/0, 0);
                    perf2 = pathGen->Path(/*iAsset*/1,0)[iStep-lag] / pathGen->refLevel(/*iAsset*/1, 0);
                    dynBaskRT->E[iAsset] = 1.0 + perf1 - perf2;
//                }
            }
        }
        finishedEarly = dynBaskRT->singleSPIStep(doingPast, iStep, endIdx, futurePayoff);
        if (inst->debugOn || doingPast) {
            report->dumpReportStep(iStep, finishedEarly, dynBaskRT.get());
            if (outPerf) {
                // report both underlyings behind the outperformance
                report->historyOutPerfAsset0[iStep] = perf1;
                report->historyOutPerfAsset1[iStep] = 1.0 + perf1 - dynBaskRT->E[0];
            }
        }
    }

    // preserve values for past 
    if (doingPast){ 
        // all known things are known by this point
        unsettledCash = dynBaskRT->instKnownCashFlows->getUnpaidValue();
        dynBaskRT->setSoFar();
        // for vol interp; only needed if we have started
        if (endIdx>dynBaskRT->iStepFirstRebal) {
            int iAsset;
            for(iAsset=0; iAsset<dynBaskRT->numAlgAssets; iAsset++) {
                dynBaskRT->Einit[iAsset] = pathGen->Path(iAsset,0)[0] / pathGen->refLevel(iAsset, 0);
            }
        }
        // some extra outputs : values "today" which means at the latest past date.
        // Only if we have any past, of course!
        if (!today_haveValues && endIdx>0) {
            today_date = inst->basket->getRebalanceDates()[endIdx-1];
            today_dynBasket = dynBaskRT->B;
            today_Z = dynBaskRT->Z(endIdx-1); // used for Client Vals
            today_nZ = dynBaskRT->nZ;
            today_E = dynBaskRT->E; // these will be any pub-lagged values
            today_nE = dynBaskRT->nE;
            today_UE = dynBaskRT->UE;
            today_SE = dynBaskRT->SE;
            today_TE = dynBaskRT->TE;
            today_BF = dynBaskRT->BF;
            today_BL = dynBaskRT->BL;
            today_haveValues = true;
        }
    } else {
        // gap risk is a statistic too. Only called for future.
        myPrices.addGapRisk(dynBaskRT->gapRiskByBucket);
        if (inst->debugOn) {
            // if there's a cutoff use that date, otherwise no cutoff 
            myPrices.addCutoffData(dynBaskRT->cutoffStep>0,
                                    inst->basket->getRebalanceDates()[dynBaskRT->cutoffStep]);
        }
    }
    if (!doingPast || !hasFuture()) {
        if (!finishedEarly) {
            double finalPayoff = inst->matCashFlow;
            double finalBasket = dynBaskRT->sumB / inst->averageOutDates.size();
            // Price the call option.
            // BL here is the value at mat date, B has been updated for final lag adj
            // With the advent of independent strike we no longer know that max(B,BL) > strike
            // so the extra max(,0) is needed.
            if (inst->isCall) {
                finalPayoff += Maths::max(Maths::max(finalBasket,dynBaskRT->BL) - strike, 0.0);
            } else {
                finalPayoff += Maths::max(dynBaskRT->BL - finalBasket, 0.);
            }
            // work out the 'gap risk' stats
            myPrices.addGapEventStats(inst->isCall, finalBasket, dynBaskRT->BL, strike);
            if(inst->debugOn || doingPast) {
                // record payoff - not the fees and coupons which are reported separately
                // note additon in case there were any fees/coupons at this date
                report->historyPayoff[endIdx-1] += finalPayoff;
            }
            if (doingPast) {
                // another known cashflow!
                DateTime payDate = settlement->settles(inst->basket->getRebalanceDates()[endIdx-1],0);
                // Nice that SPIKnownFlows internally copes with avoiding a double-count
                dynBaskRT->instKnownCashFlows->addFlow(payDate, finalPayoff * inst->notional);
                unsettledCash = dynBaskRT->instKnownCashFlows->getUnpaidValue();
            } else {
                futurePayoff += finalPayoff;
            }
        }
        futurePayoff *= inst->notional;
        // Add any fee cash flows/coupons which are notified but not yet paid.
        // Note that we provide values at final settlement date.
        myPrices.add(futurePayoff + unsettledCash);
    }
    if (dynBaskRT->dynBask->hasBFHistory) {
        dynBaskRT->dynBask->bondFloorToday = dynBaskRT->bondFloor->getLevelToday();
    }
}

// for the LogNormal path generator
CVolRequestLNArray SyntheticPortfolioInsuranceMC::getVolInterp(const MCPathGenerator* pathGen,
                                int                     iAsset) const {
    const DateTime&    startDate = getRefLevel()->getAllDates().front();
    const DateTime&    today = getToday();
    bool               fwdStarting = startDate.isGreater(today);
    CVolRequestLNArray reqarr(1);
    int                endIdx = pathGen->end(0);
    double             interp;

    // Separate fwdStarting from endIdx since if we roll from
    // fwdStarting true->false, endIdx may still not be active (i.e.
    // first rebal date after ref level set) and we still do not know
    // what the initial exposure will be. For most cases the ref level 
    // date and first rebal date should coincide, so this is no functional change.
    if (endIdx==0) {
#if FWD_START_VOL_INTERP
        // prepare values for nE etc as if at the first date
        PricesSPI myPrices(1, 1);
        payoff(pathGen, myPrices);
        interp = dynBaskRT->volInterp(0/*endIdx*/, dynBaskRT->BLSoFar, inst->strike);
#else
        interp = 1.0; // ATM for now
#endif
    } else {
        interp = dynBaskRT->volInterp(endIdx, 
                                        dynBaskRT->bondFloor->getLevel(dynBaskRT->BLSoFar, endIdx-1),
                                        inst->strike);
    }
    if (!fwdStarting) {
        interp *= pathGen->refLevel(iAsset, 0);
    }
    reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interp,
                                                                startDate,
                                                                getSimSeries()->getLastDate(),
                                                                fwdStarting));

    return reqarr;
}

const string SyntheticPortfolioInsuranceMC::DEBUG_ITE = "DEBUG_ITE";

/** invoked after final simulated path is run. */
void SyntheticPortfolioInsuranceMC::recordExtraOutput(CControl*     control,
                                Results*      results,
                                const IMCPrices& prices) const{
    // Only for pricing run
    if (control->isPricing()) {

        // if gamma or x-gamma requested simply set them to 0
        int i;
        DeltaSP delta = DeltaSP(dynamic_cast<Delta*>(control->sensitivityRequested(Delta::TYPE).get()));
        if (delta.get()) {
            SensMgr sensMgr(inst->assets);
            OutputNameArrayConstSP names = delta->hasOverrideNames()? 
                delta->overrideNames() : sensMgr.allNames(delta.get());
            names = OutputName::trim(names); // remove duplicates/empties
            for(int i=0; i<names->size(); i++) {
                results->storeScalarGreek(0.0, Delta::SECOND_ORDER_NAME, (*names)[i]);

                // store divisor used (needed for cross derivatives)
                results->storeScalarGreek(delta->getShiftSize(),
                                            Delta::NAME+
                                            Results::SHIFT_SIZE_POSTFIX,
                                            (*names)[i]);
            }
        }
        CrossGammaSP xgamma = CrossGammaSP(dynamic_cast<CrossGamma*>(control->sensitivityRequested(CrossGamma::TYPE).get()));
        if (xgamma.get()) {
            SensMgr sensMgr(inst->assets);
            Delta deltaForNames(Delta::DEFAULT_SHIFT);
            OutputNameArrayConstSP names = xgamma->hasOverrideNames()? 
                xgamma->overrideNames() : sensMgr.allNames(&deltaForNames);
            names = OutputName::trim(names); // remove duplicates/empties
            // Note need to do both ways round
            for (i = 0; i < names->size(); i++){
                for(int j=0; j<names->size(); j++) {
                    if (i!=j) {
                        OutputNameSP name(new OutputName((*names)[i].get(), 
                                                            (*names)[j].get()));
                        results->storeScalarGreek(0.0, CrossGamma::NAME, name);
                    }
                }
            }
        }

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
        request = control->requestsOutput(OutputRequest::SPI_DYN_BASKET);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_LEVELS);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DoubleArray(today_E)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_ALLOCS);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DoubleArray(today_nE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_TODAY);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_Z)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND);
        if (request) {
            double bondClosing;
            bool haveValue = dynBaskRT->getBondClosing(bondClosing);
            results->storeRequestResult(request, haveValue?
                                        IObjectSP(CDouble::create(bondClosing)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_ALLOC);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_nZ)):
                                        IObjectSP(new NotApplicable()));
        }
        double equityNAV = 0.0;
        for(i=0; i<today_nE.size(); i++) {
            equityNAV += today_E[i] * today_nE[i];
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_NAV);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(equityNAV)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_NAV_PCT);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(equityNAV/today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        double bondNAV = today_Z * today_nZ;
        request = control->requestsOutput(OutputRequest::SPI_BOND_NAV);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(bondNAV)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_NAV_PCT);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(bondNAV/today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_UNBAL_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_UE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_SUST_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_SE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_TARGET_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_TE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_BF)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR_LEVEL);
        if (request) {
            double bondFloorLevel = dynBaskRT->dynBask->bondFloorToday;
            results->storeRequestResult(request, bondFloorLevel > 0.0?
                                        IObjectSP(CDouble::create(bondFloorLevel)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_LOCKED_IN_VALUE);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_BL)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_REPORT_DATE);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DateTime(today_date)):
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
        if (request && !inst->isRainbowSPI) {
            // Restrict results to those which are meaningful - perhaps just the past
            // or only those up to any early redemption, for example
            report->downsize();
            results->storeRequestResult(request, report);
        }

        request = control->requestsOutput(OutputRequest::SPI_GAP_EVENT_STATS);
        if (request && !inst->isRainbowSPI) {
            results->storeRequestResult(request, 
                            myPrices.getGapEventStats(pvFromPaymentDate()));
        }

        // Used to be "if (inst->debugOn) {" but now always do this
        {
            const PricesSPI& myPrices = static_cast<const PricesSPI&>(prices);
            DateTimeSP expCutoffDate = myPrices.getExpectedCutoffDate();
            results->storeGreek(expCutoffDate, 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName("EXPECTED_CUTOFF_DATE")));
            results->storeGreek(IObjectSP(CDouble::create(myPrices.getExpectedCutoff())), 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName("EXPECTED_CUTOFF")));

            results->storeGreek(haveDebugITE?
                                IObjectSP(CDouble::create(debugITE)):
                                IObjectSP(new NotApplicable()), 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName(DEBUG_ITE)));
        }

#if 0
        // dump out filtered sim dates
        DateTimeArraySP allSimDates = DateTimeArraySP(new DateTimeArray(getSimSeries()->getAllDates()));
        results->storeGreek(allSimDates, 
                            Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("SIM_DATES")));
#endif

    }
}

// collects events which should be sitting there
// called after past is run
void SyntheticPortfolioInsuranceMC::retrieveEvents(EventResults* events) const {
    if (dynBaskRT->terminalStep > 0) {
        const DateTime redemptionDate =
                    inst->basket->getRebalanceDates()[dynBaskRT->terminalStep];
        events->addEvent(new TargetRedemption(redemptionDate, 
                                              dynBaskRT->terminalCoupon, 
                                              dynBaskRT->terminalTotalCoupon,
                                              dynBaskRT->targetCoupon, 
                                              dynBaskRT->bonusCoupon,
                                              dynBaskRT->terminalTotalRedemption,
                                              TargetRedemption::NOT_APPLICABLE));
    }
}

/*****************************************************************************/
/** start of SV */
/*************************************************************************/

/** equivalent to InstIntoMCProduct. Need to call parent's constructor */
SyntheticPortfolioInsuranceMCSV::SyntheticPortfolioInsuranceMCSV(const SyntheticPortfolioInsurance*   inst,
                                                                 const SimSeriesSP&                   simSeries):
    MCProductClient(inst->assets.get(),
                    inst->valueDate,
                    inst->discount.get(),
                    inst->refLevel,
                    simSeries,
                    inst->pastValues,
                    inst->instSettle.get(),
                    simSeries->getLastDate()), // sets payment date
    unsettledCash(0.),
    today_haveValues(false), today_dynBasket(0.0),
    today_Z(0.0), today_nZ(0.0), today_E(0),
    today_nE(0), today_UE(0.0), today_SE(0.0),
    today_TE(0.0), today_BF(0.0), today_BL(0.0), 
    haveDebugITE(false), debugITE(0.0),
    inst(inst),
    spotGen(new SVGenSpot(simSeries)),
    refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                  getToday())),
    dfGen(new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                           inst->instSettle, 
                           inst->basket->getRebalanceDates().back())) 
{

    static const string routine("SyntheticPortfolioInsuranceMCSV::SyntheticPortfolioInsuranceMCSV");
    try{
        int             i;
        int             numDates = inst->basket->getRebalanceDates().size();
        const DateTime& lastRebalDate = inst->basket->getLastRebalDate();
        
        lockIn = inst->lockIn->getLockInSPI(inst->basket->getRebalanceDates(),
                                            lastRebalDate);
        
        // locate fee payment dates. Helps to have the actual rebalance dates (at least as far as the back is 
        // concerned).
        DateTimeArray trueRebalDates(inst->basket->getRebalanceDates().size()-inst->basket->maxPubLagDays-inst->basket->maxExecLagDays);
        for(i=0; i<trueRebalDates.size(); i++) {
            trueRebalDates[i] = inst->basket->getRebalanceDates()[i+inst->basket->maxPubLagDays];
        }
        feeNotifDates = inst->basket->feesSPI->getFeesSPI()->getNotificationDates(trueRebalDates);
        if (feeNotifDates.size()>0 &&
            !feeNotifDates.back().equals(lastRebalDate)) {
            throw ModelException(routine,
                                    "Last fee payment date (" + feeNotifDates.back().toString() + 
                                    ") must match final rebalance date (" + lastRebalDate.toString() +
                                    ")");
        }
        
        // For extra outputs (PAYMENT_DATES & KNOWN_CASHFLOWS) which 
        DateTimeArray feePaymentDates = inst->basket->feesSPI->getFeesSPI()->getPaymentDates();

        // plus the maturity flow - "merge" it in (don't duplicate)
        const DateTime& matSettleDate = settlement->settles(inst->basket->getRebalanceDates().back(), 
                                                            0); // asset optional
        if (feePaymentDates.size()<1 ||
            matSettleDate>feePaymentDates.back()) {
            feePaymentDates.push_back(matSettleDate);
        }
        // And any coupons
        DateTimeArray couponPayDates = inst->basket->couponsSPI->getCouponsSPI()->getPaymentDates();
        instPaymentDates = DateTime::merge(feePaymentDates, couponPayDates);
        
        // average dates : none until rebalancing has started and final date must be last average date
        if (!DateTime::isSubset(inst->basket->getRebalanceDates(),inst->averageOutDates)) {
            throw ModelException(routine,
                                    "Average out dates must be a subset of the rebalance dates");
        }
        if (inst->averageOutDates.back() != inst->basket->getRebalanceDates().back()) {
            throw ModelException(routine,
                                    "Final average out date must be the same as the final entered 'rebalance' date");
        }
        if (inst->averageOutDates[0] < trueRebalDates[0]) {
            throw ModelException(routine,
                                    "Cannot average out before rebalancing has started!");
        }
        dynBaskRT = SPIRunTimeSP(new SPIRunTimeSV(inst, lockIn, settlement, 
                                                  feeNotifDates));

        numSimAssets = getNumAssets();
        if (!inst->isRainbowSPI) {
            // We may simulate more assets than we run through the algorithm. Currently
            // cope with 1 or 2 risky assets in the SPI algorithm, but can allow any number
            // of simulated assets, so long as we know how to aggregate them into 1 or 2
            // "algo risky assets"
            // For now we support only a single risky asset and an outperformance 
            // May want a nicer way to identify this ...
            if (numSimAssets != dynBaskRT->numAlgAssets) {
                if (numSimAssets != 2) {
                    throw ModelException(routine,
                                            "SPI currently supports at most 2 simulated assets, but " + 
                                            Format::toString(numSimAssets) + " supplied");
                }
                if (dynBaskRT->numAlgAssets > 1) {
                    throw ModelException(routine,
                                            "Simulating 2 assets is only supported for a 1-factor algorithm");
                }
            }
        }

        // Used to be "if (inst->debugOn) {"
        // but now being used for client vals so always allocate, though in payoff
        // don't always write (only for past and if debugOn)
        report = SPIReportSP(new SPIReport(numDates));

        double IB = inst->basket->initialBasketLevel;
        double IBF = dynBaskRT->bondFloor->getLevel(lockIn->getInitialLockIn(),
                                            dynBaskRT->iStepFirstRebal);
        double ISC = dynBaskRT->algo->sustCrash(1.-IBF/IB, dynBaskRT->iStepFirstRebal);
        double ISE = Maths::max((1.-IBF/IB)/ISC, 0.);
        debugITE = dynBaskRT->algo->targetExp(ISE);
        haveDebugITE = true;

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


void SyntheticPortfolioInsuranceMCSV::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    static const string routine = "SyntheticPortfolioInsuranceMCSV::pathGenUpdated";

    try {
        spotSV = spotGen->getSpotSV(newPathGen);
        refLevelSV = refLevelGen->getRefLevelSV(refLevelSV, newPathGen);
        dfSV = dfGen->getSVDiscFactor(newPathGen);
        dynBaskRT->pathGenUpdated(newPathGen);
    } catch (exception& e) {
        throw ModelException(e, routine);
    }
};

/** Appends 'true' (ie non derived) state variable generators
    required to the supplied collector.*/
void SyntheticPortfolioInsuranceMCSV::collectStateVars(IStateVariableCollectorSP svCollector) const{
    // ask for a reference level State Variable
    svCollector->append(refLevelGen.get());
    svCollector->append(spotGen.get());
    svCollector->append(dfGen.get());
    dynBaskRT->collectStateVars(svCollector);
};

// override hasFuture so that early termination can be captured in the past
bool SyntheticPortfolioInsuranceMCSV::hasFuture() const {
    return IMCProduct::hasFuture() && !dynBaskRT->terminatedEarly;
}

// Satisfy IHandlePaymentEvents interface
// for SPI specific record of events (pay dates, cashflows etc)
void SyntheticPortfolioInsuranceMCSV::recordEvents(Control* control,
                    Results* results) {
    static const string method("SyntheticPortfolioInsuranceMCSV::recordEvents");
    try {
        // PAYMENT_DATES is a list of all dates on which payments may occur
        // including past and potential future dates.
        // For SPI this means all paid fee payment dates and the final 
        // maturity flow
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            OutputRequestUtil::recordPaymentDates(control,results,&instPaymentDates);
        }
        
        // KNOWN_CASHFLOWS should have dates a subset of PAYMENT_DATES
        // and be supplied for all past cash flows, and any future ones
        // that are determined.
        // For SPI this means any past paid fees, and possibly a current 
        // fee due to be paid today, plus possibly the final option flow
        // XXX Need to check which ccy this is in.
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished() && 
            dynBaskRT->instKnownCashFlows->getFlows()) {
            OutputRequestUtil::recordKnownCashflows(control,
                                                    results,
                                                    discount->getCcy(),
                                                    dynBaskRT->instKnownCashFlows->getFlows()); 
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

IMCPrices* SyntheticPortfolioInsuranceMCSV::createOrigPrices(int  nbIter,
                                    int  nbSubSamples,
                                    int  mode) {
    // XXX don't support anything other than very simple "mode" so
    // XXX should I validate against mode!=0?
    return new PricesSPI(nbIter, nbSubSamples, dynBaskRT->gapRiskBucketDates);
}

/** Use this opportunity to do any LogNormal driven initialisation
    of the instrument before the main MC loop. e.g closed form barrier adjustment */
void SyntheticPortfolioInsuranceMCSV::initialiseLN(const  MCPathGenerator*  pathGen)const{
    // empty
}

/** Called within the simulation loop */
void SyntheticPortfolioInsuranceMCSV::payoff(const MCPathGenerator*  pathGen,
                                             IMCPrices&              prices) {
    static         const string routine("SyntheticPortfolioInsuranceMCSV::payoff");
    int            beginIdx = spotSV->path(0/*iAsset*/).begin(); // same for all assets
    int            endIdx   = spotSV->path(0/*iAsset*/).end();
    PricesSPI&     myPrices = static_cast<PricesSPI&>(prices);
    double         strike = Maths::isZero(inst->strike)? 
                                dynBaskRT->lockIn->getInitialLockIn() : inst->strike;
    bool           outPerf = false;
    double         perf1 = 0.0, perf2 = 0.0;
    double         futurePayoff = 0.0; // running total of payoff in future

    // First few steps are required just to collect equity levels for subsequent 
    // publication lag offset use.
    // This should work even with the call for the past, because the beginIdx/endIdx are
    // based on an absolute scale.
    int startIdx = Maths::max(dynBaskRT->iStepFirstRebal, beginIdx);
    
    dynBaskRT->init(startIdx, doingPast());

    // May want a nicer way to identify this ...
    if (numSimAssets != dynBaskRT->numAlgAssets) {
        outPerf = true;
    }

    bool finishedEarly = false;
    for(int iStep=startIdx; iStep<endIdx && !finishedEarly; 
                                            iStep++) {
        for(int iAsset=0; iAsset<dynBaskRT->numAlgAssets; iAsset++) {
            // Account for any publication lag here, so E[] is implicitly lagged for
            // the pub lag before anyone else has a look-in.
            // Note the lag can be different for different assets
            int lag = inst->basket->getPubLag(iAsset, iStep);

            if (!outPerf) {
                const SVPath& path = spotSV->path(iAsset);
                dynBaskRT->E[iAsset] = path[iStep-lag] / refLevelSV->refLevel(iAsset); 
            } else {
                // Need to collapse several simulated assets into 1 or 2 risky assets...
                // For now we support a single risky asset and an outperformance - validated earlier
                perf1 = spotSV->path(/*iAsset*/0)[iStep-lag] / refLevelSV->refLevel(/*iAsset*/0);
                perf2 = spotSV->path(/*iAsset*/1)[iStep-lag] / refLevelSV->refLevel(/*iAsset*/1);
                dynBaskRT->E[iAsset] = 1.0 + perf1 - perf2;
            }
        }
        finishedEarly = dynBaskRT->singleSPIStep(doingPast(), iStep, endIdx, futurePayoff);
        if (inst->debugOn || doingPast()) {
            report->dumpReportStep(iStep, finishedEarly, dynBaskRT.get());
            if (outPerf) {
                // report both underlyings behind the outperformance
                report->historyOutPerfAsset0[iStep] = perf1;
                report->historyOutPerfAsset1[iStep] = 1.0 + perf1 - dynBaskRT->E[0];
            }
        }
    }

    // preserve values for past 
    if (doingPast()){ 
        // all known things are known by this point
        unsettledCash = dynBaskRT->instKnownCashFlows->getUnpaidValue();
        dynBaskRT->setSoFar();
        // for vol interp; only needed if we have started
        if (endIdx>dynBaskRT->iStepFirstRebal) {
            int iAsset;
            for(iAsset=0; iAsset<dynBaskRT->numAlgAssets; iAsset++) {
                dynBaskRT->Einit[iAsset] = spotSV->path(iAsset)[0] / refLevelSV->refLevel(iAsset);
            }
        }
        // some extra outputs : values "today" which means at the latest past date.
        // Only if we have any past, of course!
        if (!today_haveValues && endIdx>0) {
            today_date = inst->basket->getRebalanceDates()[endIdx-1];
            today_dynBasket = dynBaskRT->B;
            today_Z = dynBaskRT->Z(endIdx-1); // used for Client Vals
            today_nZ = dynBaskRT->nZ;
            today_E = dynBaskRT->E; // these will be any pub-lagged values
            today_nE = dynBaskRT->nE;
            today_UE = dynBaskRT->UE;
            today_SE = dynBaskRT->SE;
            today_TE = dynBaskRT->TE;
            today_BF = dynBaskRT->BF;
            today_BL = dynBaskRT->BL;
            today_haveValues = true;
        }
    } else {
        // gap risk is a statistic too. Only called for future.
        myPrices.addGapRisk(dynBaskRT->gapRiskByBucket);
        if (inst->debugOn) {
            // if there's a cutoff use that date, otherwise no cutoff 
            myPrices.addCutoffData(dynBaskRT->cutoffStep>0,
                                    inst->basket->getRebalanceDates()[dynBaskRT->cutoffStep]);
        }
    }
    if (!doingPast() || !hasFuture()) {
        if (!finishedEarly) {
            double finalPayoff = inst->matCashFlow;
            double finalBasket = dynBaskRT->sumB / inst->averageOutDates.size();
            // Price the call option.
            // BL here is the value at mat date, B has been updated for final lag adj
            // With the advent of independent strike we no longer know that max(B,BL) > strike
            // so the extra max(,0) is needed.
            if (inst->isCall) {
                finalPayoff += Maths::max(Maths::max(finalBasket,dynBaskRT->BL) - strike, 0.0);
            } else {
                finalPayoff += Maths::max(dynBaskRT->BL - finalBasket, 0.);
            }
            // work out the 'gap risk' stats
            myPrices.addGapEventStats(inst->isCall, finalBasket, dynBaskRT->BL, strike);
            if(inst->debugOn) {
                // record payoff - not the fees and coupons which are reported separately
                // note additon in case there were any fees/coupons at this date
                report->historyPayoff[endIdx-1] += finalPayoff;
            }
            if (doingPast()) {
                // another known cashflow!
                DateTime payDate = settlement->settles(inst->basket->getRebalanceDates()[endIdx-1],0);
                // Nice that SPIKnownFlows internally copes with avoiding a double-count
                dynBaskRT->instKnownCashFlows->addFlow(payDate, finalPayoff * inst->notional);
                unsettledCash = dynBaskRT->instKnownCashFlows->getUnpaidValue();
            } else {
                futurePayoff += finalPayoff * dfSV->firstDF();
            }
        }
        futurePayoff *= inst->notional;
        // Add any fee cash flows/coupons which are notified but not yet paid.
        double myPayoff = futurePayoff + unsettledCash;
        myPrices.add(myPayoff);
    }
    if (dynBaskRT->dynBask->hasBFHistory) {
        dynBaskRT->dynBask->bondFloorToday = dynBaskRT->bondFloor->getLevelToday();
    }
}

// for the LogNormal path generator
CVolRequestLNArray SyntheticPortfolioInsuranceMCSV::getVolInterp(const MCPathGenerator* pathGen,
                                int                     iAsset) const {
    const DateTime&    startDate = getRefLevel()->getAllDates().front();
    const DateTime&    today = getToday();
    bool               fwdStarting = startDate.isGreater(today);
    CVolRequestLNArray reqarr(1);
    int                endIdx = spotSV->path(0/*iAsset*/).end();
    double             interp;

    // Separate fwdStarting from endIdx since if we roll from
    // fwdStarting true->false, endIdx may still not be active (i.e.
    // first rebal date after ref level set) and we still do not know
    // what the initial exposure will be. For most cases the ref level 
    // date and first rebal date should coincide, so this is no functional change.
    if (endIdx==0) {
#if FWD_START_VOL_INTERP
        // prepare values for nE etc as if at the first date
        PricesSPI myPrices(1, 1);
        payoff(pathGen, myPrices);
        interp = dynBaskRT->volInterp(0/*endIdx*/, dynBaskRT->BLSoFar, inst->strike);
#else
        interp = 1.0; // ATM for now
#endif
    } else {
        interp = dynBaskRT->volInterp(endIdx, 
                                        dynBaskRT->bondFloor->getLevel(dynBaskRT->BLSoFar, endIdx-1),
                                        inst->strike);
    }
    if (!fwdStarting) {
      interp *= refLevelSV->refLevel(iAsset);
    }
    reqarr[0] = CVolRequestLNSP(new LinearStrikeTSVolRequest(interp,
							     startDate,
							     getSimSeries()->getLastDate(),
							     fwdStarting));

    return reqarr;
}

const string SyntheticPortfolioInsuranceMCSV::DEBUG_ITE = "DEBUG_ITE";

/** invoked after final simulated path is run. */
void SyntheticPortfolioInsuranceMCSV::recordExtraOutput(CControl*     control,
                                Results*      results,
                                const IMCPrices& prices) const{
    // Only for pricing run
    if (control->isPricing()) {

        // if gamma or x-gamma requested simply set them to 0
        int i;
        DeltaSP delta = DeltaSP(dynamic_cast<Delta*>(control->sensitivityRequested(Delta::TYPE).get()));
        if (delta.get()) {
            SensMgr sensMgr(inst->assets);
            OutputNameArrayConstSP names = delta->hasOverrideNames()? 
                delta->overrideNames() : sensMgr.allNames(delta.get());
            names = OutputName::trim(names); // remove duplicates/empties
            for(int i=0; i<names->size(); i++) {
                results->storeScalarGreek(0.0, Delta::SECOND_ORDER_NAME, (*names)[i]);

                // store divisor used (needed for cross derivatives)
                results->storeScalarGreek(delta->getShiftSize(),
                                            Delta::NAME+
                                            Results::SHIFT_SIZE_POSTFIX,
                                            (*names)[i]);
            }
        }
        CrossGammaSP xgamma = CrossGammaSP(dynamic_cast<CrossGamma*>(control->sensitivityRequested(CrossGamma::TYPE).get()));
        if (xgamma.get()) {
            SensMgr sensMgr(inst->assets);
            Delta deltaForNames(Delta::DEFAULT_SHIFT);
            OutputNameArrayConstSP names = xgamma->hasOverrideNames()? 
                xgamma->overrideNames() : sensMgr.allNames(&deltaForNames);
            names = OutputName::trim(names); // remove duplicates/empties
            // Note need to do both ways round
            for (i = 0; i < names->size(); i++){
                for(int j=0; j<names->size(); j++) {
                    if (i!=j) {
                        OutputNameSP name(new OutputName((*names)[i].get(), 
                                                            (*names)[j].get()));
                        results->storeScalarGreek(0.0, CrossGamma::NAME, name);
                    }
                }
            }
        }

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
        request = control->requestsOutput(OutputRequest::SPI_DYN_BASKET);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_LEVELS);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DoubleArray(today_E)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_ALLOCS);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DoubleArray(today_nE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_TODAY);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_Z)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND);
        if (request) {
            double bondClosing;
            bool haveValue = dynBaskRT->getBondClosing(bondClosing);
            results->storeRequestResult(request, haveValue?
                                        IObjectSP(CDouble::create(bondClosing)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_ALLOC);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_nZ)):
                                        IObjectSP(new NotApplicable()));
        }
        double equityNAV = 0.0;
        for(i=0; i<today_nE.size(); i++) {
            equityNAV += today_E[i] * today_nE[i];
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_NAV);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(equityNAV)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_EQUITY_NAV_PCT);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(equityNAV/today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        double bondNAV = today_Z * today_nZ;
        request = control->requestsOutput(OutputRequest::SPI_BOND_NAV);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(bondNAV)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_NAV_PCT);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(bondNAV/today_dynBasket)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_UNBAL_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_UE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_SUST_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_SE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_TARGET_EXPO);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_TE)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_BF)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_BOND_FLOOR_LEVEL);
        if (request) {
            double bondFloorLevel = dynBaskRT->dynBask->bondFloorToday;
            results->storeRequestResult(request, bondFloorLevel > 0.0?
                                        IObjectSP(CDouble::create(bondFloorLevel)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_LOCKED_IN_VALUE);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(CDouble::create(today_BL)):
                                        IObjectSP(new NotApplicable()));
        }
        request = control->requestsOutput(OutputRequest::SPI_REPORT_DATE);
        if (request) {
            results->storeRequestResult(request, today_haveValues?
                                        IObjectSP(new DateTime(today_date)):
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
        if (request && !inst->isRainbowSPI) {
            // Restrict results to those which are meaningful - perhaps just the past
            // or only those up to any early redemption, for example
            report->downsize();
            results->storeRequestResult(request, report);
        }

        request = control->requestsOutput(OutputRequest::SPI_GAP_EVENT_STATS);
        if (request && !inst->isRainbowSPI) {
            results->storeRequestResult(request, 
                            myPrices.getGapEventStats(pvFromPaymentDate()));
        }

        // Used to be "if (inst->debugOn) {" but now always do this
        {
            const PricesSPI& myPrices = static_cast<const PricesSPI&>(prices);
            DateTimeSP expCutoffDate = myPrices.getExpectedCutoffDate();
            results->storeGreek(expCutoffDate, 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName("EXPECTED_CUTOFF_DATE")));
            results->storeGreek(IObjectSP(CDouble::create(myPrices.getExpectedCutoff())), 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName("EXPECTED_CUTOFF")));

            results->storeGreek(haveDebugITE?
                                IObjectSP(CDouble::create(debugITE)):
                                IObjectSP(new NotApplicable()), 
                                Results::DEBUG_PACKET, 
                                OutputNameSP(new OutputName(DEBUG_ITE)));
        }

#if 0
        // dump out filtered sim dates
        DateTimeArraySP allSimDates = DateTimeArraySP(new DateTimeArray(getSimSeries()->getAllDates()));
        results->storeGreek(allSimDates, 
                            Results::DEBUG_PACKET, 
                            OutputNameSP(new OutputName("SIM_DATES")));
#endif

    }
}

// collects events which should be sitting there
// called after past is run
void SyntheticPortfolioInsuranceMCSV::retrieveEvents(EventResults* events) const {
    if (dynBaskRT->terminalStep > 0) {
        const DateTime redemptionDate =
                    inst->basket->getRebalanceDates()[dynBaskRT->terminalStep];
        events->addEvent(new TargetRedemption(redemptionDate, 
                                              dynBaskRT->terminalCoupon, 
                                              dynBaskRT->terminalTotalCoupon,
                                              dynBaskRT->targetCoupon, 
                                              dynBaskRT->bonusCoupon,
                                              dynBaskRT->terminalTotalRedemption,
                                              TargetRedemption::NOT_APPLICABLE));
    }
}

/*************************************************************************/
/** end of SV */
/*************************************************************************/

CClassConstSP const SPIReport::TYPE = CClass::registerClassLoadMethod(
    "SPIReport", typeid(SPIReport), 
    SPIReport::load);

// work around for msvc 7 bug
//typedef SyntheticPortfolioInsuranceMC::ReportArray ReportArray;
DEFINE_TEMPLATE_TYPE(SPIReportArray);

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* SyntheticPortfolioInsurance::createProduct(const MonteCarlo* model) const {
    static const string method("SyntheticPortfolioInsurance::createProduct");
    try{
        // we need to create a SimSeries object which says which assets need
        // which dates to be simulated
        SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                     one */
        DateTimeArray instDates = DateTime::merge(lockIn->getLockInSPINoInit()->getEssentialDates(),
                                                  averageOutDates);
        simSeries->addDates(basket->prepareCoarseRebalDates(instDates,
                                                            valueDate));
        
        if (/*weightedBasket || */assets->NbAssets()==1 || assets->NbAssets()==2) {
             if (model->stateVarUsed()) {
                  return new SyntheticPortfolioInsuranceMCSV(this, simSeries);
             }
             return new SyntheticPortfolioInsuranceMC(this, simSeries);
        } else {
            throw ModelException(method, 
                                 "Only 1 or 2 assets supported.");
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

CClassConstSP const SyntheticPortfolioInsurance::TYPE = CClass::registerClassLoadMethod(
    "SyntheticPortfolioInsurance", typeid(SyntheticPortfolioInsurance), SyntheticPortfolioInsurance::load);

// * for class loading (avoid having header file) */
bool SyntheticPortfolioInsuranceLoad() {
    return (SyntheticPortfolioInsurance::TYPE != 0 &&
            SPIDynamicBasket::TYPE != 0 &&
            SPILockInWrapper::TYPE != 0 &&
            SPIBondWrapper::TYPE != 0 &&
            SPIFeesWrapper::TYPE != 0 &&
            SPIAlgorithmWrapper::TYPE != 0 &&
            SPILoanCostWrapper::TYPE != 0);
}

IObjectSP SPIAddin::getBondFloorHistory(SPIAddin* params){
    static const string routine = "SPIAddin::getBondFloorHistory";
    try {
        // Take a copy of the instrument as we're going to mess with it
        SyntheticPortfolioInsuranceSP spiCopy(copy(params->spi.get()));
        params->model->getInstrumentAndModelMarket(params->market.get(), spiCopy.get());

        // before we validate we should make sure that the dates for which we require
        // the bond floor history have a value present o/w it will fail
        // first validate dates in case they are parametric.
        spiCopy->validateDates();
        spiCopy->setRequiredBFDatesToValid(params->baseDates);
        spiCopy->Validate();

        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(params->model.get());

        // make sure all the yield curves are ok before we proceed
        int numCurves = params->ycArray.size();
        if (numCurves == 0) {
            throw ModelException("No yield curves have been supplied for calculating past bond floors");
        }
        if (numCurves != params->baseDates.size()) {
            throw ModelException("Different number of base dates and yield curves supplied");
        }
        for (int i = 0; i < numCurves; i++) {
            params->ycArray[i]->getMarket(params->model.get(), params->market.get());
        }

        return spiCopy->getBondFloorHistory(mc, params->baseDates, params->ycArray);
    } 
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** for reflection */
SPIAddin::SPIAddin():  CObject(TYPE){}

/** Invoked when Class is 'loaded' */
void SPIAddin::load(CClassSP& clazz){
    REGISTER(SPIAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSPIAddin);
    FIELD(spi, "SPI object");
    FIELD(ycArray, "array of historic yield curves");
    FIELD(baseDates, "base dates for the historic yield curves");
    FIELD(model, "Model object");
    FIELD(market, "Market cache");
    Addin::registerInstanceObjectMethod("GET_BOND_FLOOR_HISTORY",
                                        Addin::UTILITIES,
                                        "Calculates and returns the bond floor history for an SPI retrospectively",
                                        TYPE,
                                        false,
                                        Addin::expandMulti,
                                        (Addin::ObjMethod*)getBondFloorHistory);
}

IObject* SPIAddin::defaultSPIAddin(){
    return new SPIAddin();
}

CClassConstSP const SPIAddin::TYPE = CClass::registerClassLoadMethod(
    "SPIAddin", typeid(SPIAddin), load);

DRLIB_END_NAMESPACE
