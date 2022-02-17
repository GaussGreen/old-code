//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIDynamicBasket.cpp
//
//   Description : Dynamic basket interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/SPIDynamicBasket.hpp"

DRLIB_BEGIN_NAMESPACE

SPIDynamicBasket::SPIDynamicBasket(CClassConstSP clazz): CObject(clazz), maxNumSkipDates(0), skipping(false),
    iStepEarlyLagAdj(0), iStepLateLagAdj(0),
    couponsSPI(0), overrideInitialLag(false), hasRebalCost(false),
    pastRebalCostRate(0.0), futureRebalCostRate(0.0), 
    overrideInitialRebalCost(false), overrideFinalRebalCost(false), hasBFHistory(false), 
    bondFloorHistory(0), lastHistoricBFIndex(-1), cachedDFQuot(0), bondFloorToday(0.0) {} // for reflection

void SPIDynamicBasket::validatePop2Object(){
    static const string routine = "SPIDynamicBasket::validatePop2Object";

    if (numPubLagDays.size() != numExecLagDays.size()) {
        throw ModelException(routine,
                             "Must have equal length (=Num Assets) arrays for numPubLagDays (#=" +
                             Format::toString(numPubLagDays.size()) + ") and numExecLagDays (" +
                             Format::toString(numExecLagDays.size()) + ")");
    }
    maxExecLagDays = numExecLagDays[0];
    double minExecLagDays = maxExecLagDays;
    maxPubLagDays = numPubLagDays[0];
    for(int iAsset=1; iAsset<numPubLagDays.size(); iAsset++) {
        maxExecLagDays = Maths::max(maxExecLagDays, numExecLagDays[iAsset]);
        minExecLagDays = Maths::min(maxExecLagDays, numExecLagDays[iAsset]);
        maxPubLagDays = Maths::max(maxPubLagDays, numPubLagDays[iAsset]);
    }
    if (!Maths::isPositive(initialBasketLevel)) {
        throw ModelException(routine,
                             "initialBasketLevel (" + Format::toString(initialBasketLevel) +
                             ") must be positive");
    }

    if (overrideInitialLag && hasRebalCost) {
        throw ModelException(routine,
                             "cannot have rebalancing costs and use overrideInitialLag flag");
    }
    if (minExecLagDays == 0 && hasRebalCost) {
        throw ModelException(routine,
                             "cannot have rebalancing costs with zero execution lag");
    }

    // Tests suggest safe only up to biweekly, else differences can be large
    if (maxNumSkipDates<0 ||
        maxNumSkipDates>10) {
        throw ModelException(routine,
                             "Max number of skip dates (" + Format::toString(maxNumSkipDates) +
                             ") must be between 0 and 10 (= #business days)");
    }

    // New parameter which may be empty, so offer a safe placeholder
    if (!couponsSPI) {
        couponsSPI = SPICouponsWrapperSP(new SPICouponsWrapper());
    }
}

// has to be done after market data is fetched because it may have to build
// the dates from a DateBuilder using holiday info
// Called from the Validate method on the SPI but also called in the
// bond floor tool at the start because we need to get the dates before we
// Validate properly (in case of zero value in bond floor)
void SPIDynamicBasket::validateDates(const Holiday*       ccyHols, 
                                     const IMultiFactors* assets,
                                     const ObservationSourceArray& sources,
                                     bool           excludeAssetHols) {
    static const string routine = "SPIDynamicBasket::validateDates";

    if (!rebalanceDates.get()) {
        throw ModelException(routine,
                             "No rebalanceDates!");
    }

    // first let's pull the rebalance dates from the (possibly) parametrized form
    rebalanceDates->constructDates(assets, sources, ccyHols, excludeAssetHols);
    theRebalanceDates = rebalanceDates->dates();

    if (theRebalanceDates->empty()) {
        throw ModelException(routine,
                             "No rebalanceDates!");
    }
    if (theRebalanceDates->size() <= maxExecLagDays + maxPubLagDays) {
        throw ModelException(routine,
                             "Given lags require more than " + 
                             Format::toString(maxPubLagDays) + " + "  +
                             Format::toString(maxExecLagDays) + " = " +
                             Format::toString(maxExecLagDays + maxPubLagDays) +
                             " rebalance dates, but only " +
                             Format::toString(theRebalanceDates->size()) + 
                             " given");
    }
    for(int i=1; i<theRebalanceDates->size(); i++) {
        if ((*theRebalanceDates)[i-1]>=(*theRebalanceDates)[i]) {
            throw ModelException(routine,
                                 "Potential rebalance dates must be strictly increasing but [" +
                                 Format::toString(i) + "]=" + (*theRebalanceDates)[i-1].toString() +
                                 " is not before [" + Format::toString(i+1) + "]=" +
                                 (*theRebalanceDates)[i].toString());
        }
    }        
}

void SPIDynamicBasket::validateBondFloor(const DateTime today) {
    static const string routine = "SPIDynamicBasket::validateBondFloor";

    // Check whether we need a bond floor history (don't allow one if we don't need it)
    // we don't need one if we have a linear bond floor. Otherwise
    // will need one if we have past dates and we have guaranteed coupons or fees on notional
    if ( !bond->getBondSPI()->getLinearBondFloor() &&
        (couponsSPI->getCouponsSPI()->guaranteedCoupons()
          || feesSPI->feeType() == SPI_FEES_TYPE_PNL )) {
        int numPastDates;
        int i;
        for (numPastDates = 0; numPastDates < theRebalanceDates->size()
                        && (*theRebalanceDates)[numPastDates] < today; 
                        numPastDates++) {/*empty loop*/}
        
        if (!bondFloorHistory || bondFloorHistory->empty()) {
            throw ModelException(routine,
                 "A bond floor history must be supplied in the presence of guaranteed coupons or fees on notional");
        }
        if (bondFloorHistory->size() > theRebalanceDates->size()) {
            throw ModelException(routine,
                 "Bond floor history contains more dates than there are rebalance dates");
        }
        for (i = 0; i < bondFloorHistory->size(); i++) {
            if ((*theRebalanceDates)[i] != (*bondFloorHistory)[i].date) {
                throw ModelException(routine, "Date in bond floor history (" 
                    + (*bondFloorHistory)[i].date.toString() +
                    ") does not match corresponding rebalance date (" +
                    (*theRebalanceDates)[i].toString() + ")");
            }
        }
        if (bondFloorHistory->size() < numPastDates) {
            throw ModelException(routine,
                 "Bond floor history does not cover all the past rebalance dates");
        }
        for (i = 0; i < numPastDates; i++) {
            if (!Maths::isPositive((*bondFloorHistory)[i].amount)) {
                throw ModelException(routine, "Invalid value in bond floor history for date (" 
                    + (*bondFloorHistory)[i].date.toString() + ")");
            }
        }

        lastHistoricBFIndex = numPastDates - 1;
        hasBFHistory = true;
    }
    else {
        // stop people putting in a bond floor history if we don't need one
        if (!!bondFloorHistory && !bondFloorHistory->empty()) {
            throw ModelException(routine,
                "A bond floor history must NOT be supplied if\n"
                "there are no guaranteed coupons or fees on notional\n" 
                "or if a linear bond floor is supplied");
        }
    }
}

// For performance try to remove potential rebalance dates so
// long as any gap left is <= maxNumSkipDates. 
// Constrain to have some unskipped at the start and at the end.
// Any past dates are of course unskipped.
const DateTimeArray& SPIDynamicBasket::prepareCoarseRebalDates(const DateTimeArray& someEssentialDates,
                                             const DateTime&      today) {
    sparseRebalanceDates = DateTimeArray(0);
    HolidayConstSP hols = HolidayConstSP(Holiday::weekendsOnly());
    if (maxNumSkipDates>0) {
        DateTimeArray moreEssentialDates = DateTime::merge(someEssentialDates,
                                                           bond->getBondSPI()->getEssentialDates());
        DateTimeArray someMoreEssentialDates = DateTime::merge(moreEssentialDates,
                                                               feesSPI->getFeesSPI()->getEssentialDates());
        DateTimeArray essentialDates = DateTime::merge(someMoreEssentialDates,
                                                       couponsSPI->getCouponsSPI()->getEssentialDates());

        // Somewhat arbitrary but designed to give time for lag to resolve at the beginning
        // and final days to be accurate leading up to maturity. No lag adjustment is performed
        // between these dates. In order to turn off pub lag during skipping it is best to
        // have daily dates after today at least as far out as the pub lag and then we don't
        // ignore true levels - just simulated ones.
        // First, ...
        int iFirstFut;
        for(iFirstFut = 0; iFirstFut < theRebalanceDates->size() && (*theRebalanceDates)[iFirstFut]<=today; iFirstFut++) {
            // ...find rebalance date near today.
        }
        // then go maxExecLagDays after with buffer of 4 periods
        int iStepEarliestSkip = Maths::min(iFirstFut+maxExecLagDays+4, theRebalanceDates->size()-1);
        const DateTime& earliestSkipDate = (*theRebalanceDates)[iStepEarliestSkip];

        int numDailyBack = Maths::max(20, maxExecLagDays+1); // at least cover final lag period
        int iStepLatestSkip = Maths::max(0, theRebalanceDates->size()-numDailyBack);
        const DateTime& latestSkipDate = (*theRebalanceDates)[iStepLatestSkip];

        int iEssential;
        // start 'e' at first essential date before latestSkipDate
        for (iEssential=essentialDates.size()-1; iEssential>=0 && 
                 essentialDates[iEssential]>=latestSkipDate; iEssential--) {/* finding iEssential*/}
        
        // Not sure how to get hold of holiday info so account for weekends
        // only for now. Start from the back so as today moves the sim dates
        // remain the same.
        for(int iRebal=theRebalanceDates->size()-1; iRebal>=0; iRebal--) {
            const DateTime& thisDate = (*theRebalanceDates)[iRebal];
            // check essential dates first since we're relying on ordering and that they are subset
            if (iEssential>=0 && thisDate == essentialDates[iEssential]) { 
                sparseRebalanceDates.insert(sparseRebalanceDates.begin(), thisDate); // all essential dates
                iEssential--;
            } else if (thisDate >= latestSkipDate) {
                sparseRebalanceDates.insert(sparseRebalanceDates.begin(), thisDate); // include all dates at end
            } else if (thisDate <= earliestSkipDate) {
                sparseRebalanceDates.insert(sparseRebalanceDates.begin(), thisDate); // include all dates past and at start 
            } else {
                int diff = hols->businessDaysDiff(thisDate, sparseRebalanceDates.front());
                if (diff>maxNumSkipDates) {
                    sparseRebalanceDates.insert(sparseRebalanceDates.begin(), thisDate); // no gaps too long
                }
            }
        }
        // sanity check - also avoids somewhat cryptic error messages later on e.g. if lock-in has
        // a Sunday late on, but rebalanceDates are entered without the Sunday, then the algorithm
        // here means that many earlier essential dates are also missed. Consequently the error message 
        // from lock-in may be "early lock-in date is missing" when in fact the full rebalance dates
        // contain it, though the sparse ones don't. I'm not going to change the algorithm now just to
        // get easier error msgs.
        if (iEssential >= 0) { 
            throw ModelException("SPIDynamicBasket::prepareCoarseRebalDates",
                                 "Essential date " + 
                                 essentialDates[iEssential].toString() + 
                                 " is not in the rebalance dates!\n" +
                                 "Try turning off SkipDates and rerunning \n" +
                                 "for a more helpful error message, or\n" +
                                 "check lock-ins, bond, fee payment dates.");
        }

        // only set the skipping flag if we are actually skipping.
        // someone might have 10 day skipping but monthly rebalance dates. 
        // If code thinks it's skipping in this case then GAP_RISK will be incorrectly scaled 
        // and lag adjustment will be incorrectly ignored
        if (sparseRebalanceDates.size() < theRebalanceDates->size()) {
            skipping = true;
            // reset these to be in terms of sparseRebalanceDates
            for(iStepEarlyLagAdj=0; iStepEarlyLagAdj<sparseRebalanceDates.size()-1 &&
                    sparseRebalanceDates[iStepEarlyLagAdj]<=earliestSkipDate; iStepEarlyLagAdj++) {
                // finding iStepEarlyLagAdj. Note < size()-1 so final result always in range
            }
            iStepLateLagAdj = Maths::max(0, sparseRebalanceDates.size()-numDailyBack+maxExecLagDays);

            // Note that skipping dates requires adjustment to gap risk measure,
            // so we cache some useful info here (= num bus days between sim dates)
            daysBetweenSimDates = IntArray(sparseRebalanceDates.size(), 1);
            // first value is 1 so that we don't lose first date's gap risk
            for(int i=1; i<sparseRebalanceDates.size(); i++) {
                // XXX what if we have something other than daily rebal and some are skipped 
                // we may only want to scale by the number of periods not days
                // this might be tricky to do - currently we only have daily or monthly rebal (the latter can't skip)
                // so it's not an issue now. Needs some thought for the future maybe
                daysBetweenSimDates[i] = hols->businessDaysDiff(sparseRebalanceDates[i-1], 
                                                                sparseRebalanceDates[i]);
            }
        }
    } else {
        sparseRebalanceDates = *theRebalanceDates;
    }

    return sparseRebalanceDates;
}

const DateTimeArray& SPIDynamicBasket::getRebalanceDates() const {
    if (sparseRebalanceDates.empty()) {
        throw ModelException("SPIDynamicBasket::getRebalanceDates",
                             "Internal error!");
    }
    return sparseRebalanceDates;
}

const DateTimeArray& SPIDynamicBasket::getAllRebalanceDates() const {
    return *(theRebalanceDates.get());
}

// Dealing with actual - not sparse - rebalance dates 
const DateTime& SPIDynamicBasket::getNextRebalanceDate(const DateTime& from) const {
    for(int i=0; i<theRebalanceDates->size(); i++) {
        if ((*theRebalanceDates)[i]>=from) {
            return (*theRebalanceDates)[i];
        }
    }
    throw ModelException("SPIDynamicBasket::getNextRebalanceDate",
                         "No rebalance dates after " + from.toString());
}

const DateTime& SPIDynamicBasket::getLastRebalDate() const {
    return sparseRebalanceDates[sparseRebalanceDates.size()-maxExecLagDays-1];
}

// Date on which basket is finally known
const DateTime& SPIDynamicBasket::getFinalDate() const {
    return theRebalanceDates->back();
}

bool SPIDynamicBasket::doLagAdj(int iStep) const {
    return maxExecLagDays>0 &&
        (!skipping ||
         iStep < iStepEarlyLagAdj ||
         iStep > iStepLateLagAdj);
}

bool SPIDynamicBasket::doRebalCost(int iStep) const {
    return hasRebalCost &&
        (!skipping ||
         iStep < iStepEarlyLagAdj ||
         iStep > iStepLateLagAdj);
}

int SPIDynamicBasket::getPubLag(int iAsset,
              int iStep) {
    if (doLagAdj(iStep)) {
        return numPubLagDays[iAsset];
    }
    return 0;     // No lag if skipping
}

void SPIDynamicBasket::scaleGapRisk(int     iStep,
                  double& gaprisk) const {
    if (skipping) { // else no need
        gaprisk *= daysBetweenSimDates[iStep];
    }
}

    // if we've rolled forward we've stored DF quotients between all rebalance dates
// now recast these as a set of DFQuotients between the current sparse rebal dates
const DoubleArray SPIDynamicBasket::unpickDiscountFactorQuotients() const{

    DoubleArray DFQuot;
    int currentRebalIdx = 0;

    // note first sparse date will always be the first rebal date but just to be sure
    if ((*theRebalanceDates)[0] != sparseRebalanceDates[0]) {
        throw ModelException("SPIDynamicBasket::unpickDiscountFactorQuotients", 
                            "Internal error");
    }
    for (int i = 0; i < sparseRebalanceDates.size() - 1 ; i++) {
        double DFQ = 1.0;
        while (currentRebalIdx < theRebalanceDates->size() - 1 && 
                    (*theRebalanceDates)[currentRebalIdx+1] <= sparseRebalanceDates[i+1]) {
            DFQ *= cachedDFQuot[currentRebalIdx++];
        }
        DFQuot.push_back(DFQ);
    }

    return DFQuot;
}

const DoubleArray SPIDynamicBasket::getDiscountFactorsQuotients(bool sparse, const YieldCurve* disc) const{
    const DateTimeArray dfDates = sparse ? sparseRebalanceDates : *(theRebalanceDates.get());
    DoubleArray DFQuot;

    for (int i = 0; i < dfDates.size() - 1 ; i++) {
        DFQuot.push_back(disc->pv(dfDates[i], dfDates[i+1]));
    }

    return DFQuot;
}

// Update contents for a time shift
void SPIDynamicBasket::roll(const Theta::Util&   thetaUtil,
          const YieldCurve*    disc,
          double               dayCountBasis) {
    bond->getBondSPI()->roll(thetaUtil, disc);
    loanCost->getLoanCostSPI()->roll(thetaUtil, disc, dayCountBasis);
}

// note there's also a protected default constructor for use with the derived class
SPIDynamicBasket::SPIDynamicBasket(): CObject(TYPE), maxNumSkipDates(0), skipping(false),
    iStepEarlyLagAdj(0), iStepLateLagAdj(0),
    couponsSPI(0), overrideInitialLag(false), hasRebalCost(false),
    pastRebalCostRate(0.0), futureRebalCostRate(0.0), 
    overrideInitialRebalCost(false), overrideFinalRebalCost(false), hasBFHistory(false), 
    bondFloorHistory(0), lastHistoricBFIndex(-1), cachedDFQuot(0), bondFloorToday(0.0) {} // for reflection

IObject* SPIDynamicBasket::defaultSPIDynamicBasket(){
    return new SPIDynamicBasket();
}

/** Invoked when Class is 'loaded' */
void SPIDynamicBasket::load(CClassSP& clazz){
    REGISTER(SPIDynamicBasket, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSPIDynamicBasket);
    FIELD(bond,             "SPI Bond Wrapper ");
    FIELD(feesSPI,          "SPI fees Wrapper");
    FIELD(couponsSPI,       "SPI Coupons Wrapper");
    FIELD_MAKE_OPTIONAL(couponsSPI);
    FIELD(algorithm,        "algorithm");
    FIELD(loanCost,         "loanCost");
    FIELD(numPubLagDays, "numPubLagDays");
    FIELD(numExecLagDays, "numExecLagDays");
    FIELD(overrideInitialLag, "overrideInitialLag");
    FIELD_MAKE_OPTIONAL(overrideInitialLag);
    FIELD(initialBasketLevel, "initialBasketLevel");
    FIELD(rebalanceDates, "rebalanceDates");
    FIELD(maxNumSkipDates,      "Will try to skip dates leaving no gaps bigger than this # business days");
    FIELD_MAKE_OPTIONAL(maxNumSkipDates);
    FIELD(bondFloorHistory,"bond floor history for past closings");
    FIELD_MAKE_OPTIONAL(bondFloorHistory); // we may not need it
    FIELD(hasRebalCost, "is there a cost associated with rebalancing?");
    FIELD_MAKE_OPTIONAL(hasRebalCost);
    FIELD(pastRebalCostRate, "the percentage cost associated with rebalancing in the past");
    FIELD_MAKE_OPTIONAL(pastRebalCostRate);
    FIELD(futureRebalCostRate, "the percentage cost associated with rebalancing in the future");
    FIELD_MAKE_OPTIONAL(futureRebalCostRate);
    FIELD(overrideInitialRebalCost, "prevent rebalance cost being deducted on first rebalance date");
    FIELD_MAKE_OPTIONAL(overrideInitialRebalCost);
    FIELD(overrideFinalRebalCost, "prevent rebalance cost being deducted on last rebalance date");
    FIELD_MAKE_OPTIONAL(overrideFinalRebalCost);
    //transient fields
    FIELD(daysBetweenSimDates, "daysBetweenSimDates");
    FIELD(skipping, "are we actualy skipping dates");
    FIELD(theRebalanceDates, "rebalance dates after unpcaking parametric form");
    FIELD(maxExecLagDays, "maxExecLagDays");
    FIELD(maxPubLagDays, "maxPubLagDays");
    FIELD(iStepEarlyLagAdj, "iStepEarlyLagAdj");
    FIELD(iStepLateLagAdj, "iStepLateLagAdj");
    FIELD(hasBFHistory, "is there a bond floor history");
    FIELD(lastHistoricBFIndex, "index of last rebalance date before today (before rolling)");
    FIELD(cachedDFQuot,"cached DF quotients for rolling bond floor");
    FIELD(bondFloorToday,"cached today value for rolling bond floor");
    FIELD_MAKE_TRANSIENT(theRebalanceDates); // built from the date builder
    FIELD_MAKE_TRANSIENT(daysBetweenSimDates);
    FIELD_MAKE_TRANSIENT(maxExecLagDays);
    FIELD_MAKE_TRANSIENT(maxPubLagDays);
    FIELD_MAKE_TRANSIENT(iStepEarlyLagAdj);
    FIELD_MAKE_TRANSIENT(iStepLateLagAdj);
    FIELD_MAKE_TRANSIENT(skipping);
    FIELD_MAKE_TRANSIENT(hasBFHistory);
    FIELD_MAKE_TRANSIENT(lastHistoricBFIndex);
    FIELD_MAKE_TRANSIENT(cachedDFQuot);
    FIELD_MAKE_TRANSIENT(bondFloorToday);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIDynamicBasket::TYPE = CClass::registerClassLoadMethod(
    "SPIDynamicBasket", typeid(SPIDynamicBasket), SPIDynamicBasket::load);

/****************************************************************************/

SPIDynamicBasketSuper::SPIDynamicBasketSuper(): SPIDynamicBasket(TYPE) {} // for reflection

IObject* SPIDynamicBasketSuper::defaultSPIDynamicBasketSuper(){
    return new SPIDynamicBasketSuper();
}

/** Invoked when Class is 'loaded' */
void SPIDynamicBasketSuper::load(CClassSP& clazz){
    REGISTER(SPIDynamicBasketSuper, clazz);
    SUPERCLASS(SPIDynamicBasket);
    EMPTY_SHELL_METHOD(defaultSPIDynamicBasketSuper);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
CClassConstSP const SPIDynamicBasketSuper::TYPE = CClass::registerClassLoadMethod(
    "SPIDynamicBasketSuper", typeid(SPIDynamicBasketSuper), SPIDynamicBasketSuper::load);

DRLIB_END_NAMESPACE
