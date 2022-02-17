//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : QuantoCDSParSpreads.cpp
//
//   Description : Holds the par spreads in a different currency to the name's
//                 native currency (eg IBM par spreads in euros)
//
//   Author      : Mark A Robson
//
//   Date        : November 26, 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QuantoCDSParSpreads.hpp"
#include "edginc/CDSParSpreadsBase.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/RiskyLogOfDiscFactorKey.hpp"
#include "edginc/Actual360.hpp"
#include "edginc/Actual365F.hpp"

DRLIB_BEGIN_NAMESPACE
QuantoCDSParSpreads::IAlgorithm::IAlgorithm() {}
QuantoCDSParSpreads::IAlgorithm::~IAlgorithm(){}

QuantoCDSParSpreads::IAlgorithmBuilder::~IAlgorithmBuilder(){}
void QuantoCDSParSpreads::IAlgorithmBuilder::load(CClassSP& clazz){
    REGISTER_INTERFACE(IAlgorithmBuilder, clazz);
    EXTENDS(IObject);
}

CClassConstSP const QuantoCDSParSpreads::IAlgorithmBuilder::TYPE =
CClass::registerInterfaceLoadMethod(
    "QuantoCDSParSpreads::IAlgorithmBuilder", typeid(IAlgorithmBuilder), load);


QuantoCDSParSpreads::~QuantoCDSParSpreads()
{}

/** Override CObject default to preserve the IAlgorithm and cached data */
IObject* QuantoCDSParSpreads::clone() const{
    QuantoCDSParSpreads* myCopy = STATIC_CAST(QuantoCDSParSpreads,
                                              CObject::clone());
    myCopy->algorithm = algorithm;
    myCopy->quantodSpreadCurve = quantodSpreadCurve;
    myCopy->theDefaultRates = theDefaultRates;
    myCopy->theHashCode = theHashCode;  // POD default copy constructor
    return myCopy;
}

/** Override CObject::fieldsUpdated() to invalidate cached data */
void QuantoCDSParSpreads::fieldsUpdated(const CFieldArray& fields){
    theDefaultRates.reset();
    quantodSpreadCurve.reset();
    theHashCode.valid = false;
}

/** Returns the name of this object.*/
string QuantoCDSParSpreads::getName() const{
    return name;
}

string QuantoCDSParSpreads::getParSpreadsName() const {
    return parSpreads.getName();
}


/** little utility helper - should have getCorrelationName on MarketObject
    and put this on Correlation and just take mo's rather than strings etc */
static CorrelationSP getCorrelation(const IModel*        model,
                                    const string&        mo1,
                                    CClassConstSP        clazz1,
                                    const string&        mo2,
                                    CClassConstSP        clazz2,
                                    const MarketData*    market){
    const string& corrName = market->getCorrelationName(mo1, mo2);
    CorrelationBaseSP corr(model->getCorrelation(corrName,
                                                 clazz1,
                                                 clazz2,
                                                 Correlation::TYPE,
                                                 market));
    return CorrelationSP::dynamicCast(corr);
}

void QuantoCDSParSpreads::getMarket(const IModel*     model,
                                    const MarketData* market){
    parSpreads.getData(model, market);
    try {
        fx.getData(model, market);
        // Now get the 'growth' YC and correlations
        FXAssetSP fx = this->fx.getSP(); // for ease
        // duplicate these yield curves and call setProjectionCurve() on them
        MarketObjectSP domMO(market->GetData(model, fx->getYCName(),
                                             IYieldCurve::TYPE));
        domYC = IYieldCurveSP::dynamicCast(domMO);
        domYC->setProjectionCurve();
        MarketObjectSP forMO(market->GetData(model, fx->getRiskCcy()->getName(),
                                             IYieldCurve::TYPE));
        forYC = IYieldCurveSP::dynamicCast(forMO);
        forYC->setProjectionCurve();

        // now for the correlations
        ICDSParSpreadsConstSP cds(parSpreads.getSP()); // for ease
        corrDomFor = getCorrelation(model, domYC->getCcy(), domYC->getClass(),
                                    forYC->getCcy(), forYC->getClass(), market);
        corrFxDom = getCorrelation(model, fx->getName(), fx->getClass(),
                                   domYC->getCcy(), domYC->getClass(), market);
        corrFxFor = getCorrelation(model, fx->getName(), fx->getClass(),
                                   forYC->getCcy(), forYC->getClass(), market);
        corrFxCDS = getCorrelation(model,
                                   cds->getCorrelationName(), cds->getClass(),
                                   fx->getName(), fx->getClass(),
                                   market);
        corrCDSDom = getCorrelation(model,
                                    cds->getCorrelationName(), cds->getClass(),
                                    domYC->getCcy(), domYC->getClass(),
                                    market);
        corrCDSFor = getCorrelation(model,
                                    cds->getCorrelationName(), cds->getClass(),
                                    forYC->getCcy(), forYC->getClass(),
                                    market);
    }
    catch (ModelException& e) {
        ModelException f;
        try {
            f = ModelException(e, "QuantoCDSParSpreads::getMarket()", string() +
                "Couldn't retrieve all necessary correlations "
                "(did you mean to book name '" +
                parSpreads->getCorrelationName() + "' with " +
                fx->getBaseCcyIsoCode() + " and " + fx->getRiskCcyIsoCode() +
                " yield curves?)");
        }
        catch (...) {
            throw e;
        }
        throw f;
    }

    // now look up how to do the quanto adjustment
    if (!IAlgorithmBuilder::TYPE->isInstance(model)){
        throw ModelException("QuantoCDSParSpreads::getMarket",
                             model->getClass()->getName()+ " model does not "
                             "support quantod CDS Par Spread Curves");
    }
    const IAlgorithmBuilder& builder =
        dynamic_cast<const IAlgorithmBuilder&>(*model);
    algorithm = builder.cdsQuantoAlgorithm(this);
}

ParSpreadCurveConstSP QuantoCDSParSpreads::getParSpreadCurve() const {
    try {
        if (!quantodSpreadCurve) {
            ExpiryArrayConstSP expiries = getParSpreadsExpiries();
            DoubleArray quantodSpreads(expiries->size());

            DefaultRatesSP rates = defaultRates();
            DateTime effectiveDate =
                    parSpreads->spotDate(rates->getValueDate());

            for (int e = 0; e < expiries->size(); ++e) {
                quantodSpreads[e] = rates->cdsParSpread(
                    effectiveDate,
                    parSpreads->getSwapFrequency(),
                    1, // notional
                    parSpreads->getRecovery(),
                    parSpreads->isFeeAccrued(),
                    effectiveDate, // accruedEffectiveDate
                    (*expiries)[e]->toDate(effectiveDate),
                    fx->getBaseCcy(), // discount in our visible (new) ccy
                    false, // isE2C
                    parSpreads->dayCountConv());
            }

            quantodSpreadCurve.reset(new ParSpreadCurve(
                getName() + "_QUANTOD_SPREADS",
                ExpiryArraySP::constCast(expiries),
                quantodSpreads));
        }
        return quantodSpreadCurve;
    }
    catch (exception& e) {
        throw ModelException(e, "QuantoCDSParSpreads::getParSpreadCurve()");
    }
}

/** Quanto'd par spreads, computed from defaultRates() */
DoubleArrayConstSP QuantoCDSParSpreads::getParSpreads() const {
    try {
        return getParSpreadCurve()->getParSpreads();
    }
    catch (exception& e) {
        throw ModelException(e, "QuantoCDSParSpreads::getParSpreads()");
    }
}

/** Allows a shift to multiple points on the curve simultaneously */
// tenorShifts is expected to be the same size as the expiries/spreads
// shifts are additive
// modifies the object
// This method is outside the usual sensitivity framework
void QuantoCDSParSpreads::bucketShift(const DoubleArray& tenorShifts)
{
    parSpreads->bucketShift(tenorShifts);
}

/** Expiries (relative: eg. 1M) */
ExpiryArrayConstSP QuantoCDSParSpreads::getParSpreadsExpiries() const {
    return parSpreads->getParSpreadsExpiries();
}

bool QuantoCDSParSpreads::isFeeAccrued() const {
    return parSpreads->isFeeAccrued();
}

int QuantoCDSParSpreads::getSwapFrequency() const {
    return parSpreads->getSwapFrequency();
}

//// Returns the clean spreads in the new currency
DefaultRatesSP QuantoCDSParSpreads::defaultRates() const{
    if (!theDefaultRates) {
        // Obtain the clean spread associated to this CDSParSpread
        // from the cache
        theDefaultRates = getCachedDefaultRates(this, QUANTO_CDS_PAR_SPREADS);

        // cache IR grid points where we interpolate
        algorithm->cacheGridPoints(domYC, forYC);
    }

    return theDefaultRates;
}

/** return the bootstrapped dates */
DateTimeArray QuantoCDSParSpreads::zeroDates() const
{
    DefaultRatesSP dR = defaultRates();
    DateTimeArraySP dates = dR->getDefaultDates();
    return *(dates.get());
}

// accessor methods for logOfDiscFactorKey
IDiscountCurve::IKey* QuantoCDSParSpreads::getDiscountKey() const
{
    return parSpreads->getDiscountKey();
}

DefaultRates::IKey* QuantoCDSParSpreads::getRiskyKey() const
{
    DefaultRatesSP dR = defaultRates();
    return dR->logOfDefaultPVKey();
}

/** Returns a key used to optimise repeated calculations of discount
    factors/forward rate. The calc method for this key returns the 
    natural logarithm of the discount factor (or equivalently the
    product of the forward rate (continuous, Act/365F) and the negative
    year fraction (Act/365F) betweeen the two dates.
    The default implementation has no performance improvements. */
IDiscountCurve::IKey* QuantoCDSParSpreads::logOfDiscFactorKey() const
{
    return new RiskyLogOfDiscFactorKey(this);
}

/** Returns the recovery rate */
double QuantoCDSParSpreads::getRecovery() const{
    return parSpreads->getRecovery();
}

/** Returns true if the name has defaulted */
bool QuantoCDSParSpreads::defaulted() const{
    return parSpreads->defaulted();
}

/** Returns the date that this name defaulted. If a default has not
    happened an 'empty' DateTime is returned */
const DateTime& QuantoCDSParSpreads::getDefaultDate() const{
    return parSpreads->getDefaultDate();
}

/** return CDS Par Spreads currency [isocode] */
string QuantoCDSParSpreads::getCcy() const{
    return fx->getBaseCcyIsoCode();
}

/** return CDS Par Spreads discount yield curve name (not the isocode) */
string QuantoCDSParSpreads::getYCName() const{
    return fx->getYCName();
}

/** This currently fails as it's not clear what the right behaviour is.
    One possibility is to just use the vol of the curve that was quanto'd */
IVolProcessed* QuantoCDSParSpreads::getProcessedVol(
    const CVolRequest* volRequest) const{
    throw ModelException("QuantoCDSParSpreads::getProcessedVol",
                         "Not currently supported");
}

/** Returns the name of the curve that was quanto'd - not ideal but then
    a) we currently haven't a way of specifying such a correlation and
    b) you probably want the model to consider the separate factors
    anyway */
string QuantoCDSParSpreads::getCorrelationName() const{
    return parSpreads->getCorrelationName();
}

/** Returns the date when to stop tweaking after lastDate */
const DateTime QuantoCDSParSpreads::stopTweaking(
    const DateTime& lastDate) const{
    // as a first order approximation lets just try the underlying curve
    return parSpreads->stopTweaking(lastDate);
}

/** Part of ICreditCurve interface - currently just passes call down to
    original (ie unquantoed) curve */
double QuantoCDSParSpreads::getCurrentSpread(
    const DateTime& valueDate,
    const DateTime& maturityDate) const{
    return getCurrentSpread(valueDate, maturityDate,
                            getBadDayConvention().get(),
                            getHolidays().get());
}

BadDayConventionConstSP QuantoCDSParSpreads::getBadDayConvention() const {
    return parSpreads->getBadDayConvention();
}

HolidayConstSP QuantoCDSParSpreads::getHolidays() const {
    return parSpreads->getHolidays();
}

/** Part of ICreditCurve interface - currently just passes call down to
    original (ie unquantoed) curve */
double QuantoCDSParSpreads::getCurrentSpread(const DateTime& valueDate,
                                             const DateTime& maturityDate,
                                             const BadDayConvention* bdc,
                                             const Holiday* hols) const
{
    return getParSpreadCurve()->getCurrentSpread(valueDate, maturityDate,
                                                 bdc, hols);
}

/** ParSpreadMaxTweakSize support */
string QuantoCDSParSpreads::sensName(ParSpreadMaxTweakSize* shift) const {
    return name;
};

DoubleArrayConstSP QuantoCDSParSpreads::sensMaxTweakSize(ParSpreadMaxTweakSize* shift) const {
    return calculateTweakSizes();
}

/** Calculates the maximum tweak sizes for this ICDSParSpreads.
Currently just passes this down to original (ie unquantoed) curve */
DoubleArrayConstSP QuantoCDSParSpreads::calculateTweakSizes() const{
    return parSpreads->calculateTweakSizes();
}

/** Returns the day count convention used for par CDS */
DayCountConvention* QuantoCDSParSpreads::dayCountConv() const{
    return parSpreads->dayCountConv();
}

/** Returns the effective date */
DateTime QuantoCDSParSpreads::spotDate(const DateTime& valueDate) const{
    return parSpreads->spotDate(valueDate);
}

#ifdef CDS_BACKWARD_COMPATIBILITY
/** Set the bad day convention */
void QuantoCDSParSpreads::setBadDayConvention(BadDayConventionSP bdc){
    parSpreads->setBadDayConvention(bdc);
}

/** Set the holidays */
void QuantoCDSParSpreads::setHolidays(HolidaySP holidays){
    parSpreads->setHolidays(holidays);
}
#endif

/** builds a QuantoCDSParSpreads object with the names of the relevant
    market data filled in. Need to call getMarket after this in order
    for the market data to be populated */
QuantoCDSParSpreads::QuantoCDSParSpreads(const string& parSpreadsName,
                                         const string& fxName):
    MarketObject(TYPE), name(parSpreadsName),
    parSpreads(parSpreadsName), fx(fxName)
{
    // Need to initialise the hash code buffer (cache) to invalid
    theHashCode.valid = false;
}

//// default constructor
QuantoCDSParSpreads::QuantoCDSParSpreads(): MarketObject(TYPE) {
    // Need to initialise the hash code buffer (cache) to invalid
    theHashCode.valid = false;
}

IObject* QuantoCDSParSpreads::defaultConstructor(){
    return new QuantoCDSParSpreads();
}

void QuantoCDSParSpreads::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(QuantoCDSParSpreads, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(ICDSParSpreads);
    IMPLEMENTS(ParSpreadMaxTweakSize::IShift);
    IMPLEMENTS(Duration::IParHandlerWithClosedForm);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Name for this par spread curve");
    FIELD_NO_DESC(domYC);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(domYC);
    FIELD_NO_DESC(forYC);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(forYC);
    FIELD(parSpreads, "Par Spreads in native currency");
    FIELD(fx, "FX to go from native to new currency");
    FIELD_NO_DESC(corrDomFor);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrDomFor);
    FIELD_NO_DESC(corrFxDom);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrFxDom);
    FIELD_NO_DESC(corrFxFor);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrFxFor);
    FIELD_NO_DESC(corrFxCDS);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrFxCDS);
    FIELD_NO_DESC(corrCDSDom);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrCDSDom);
    FIELD_NO_DESC(corrCDSFor);
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrCDSFor);
}


//------------------------------
// Cache-related methods
//------------------------------

/** Returns a DefaultRates object from the cache.
 * This method collects from the cache the default rates associated to the
 * ICDSParSpreads passed as parameter (calculating them if not there in
 * the first place).
 * It is typically invoked from the "defaultRates" method if required (there
 * can be another layer of caching there) and "external classes" are
 * strongly DISCOURAGED from using this method directly: the "defaultRates"
 * method should be used instead.
 * Here we are in an QuantoCDSParSpreads so the cache will be kept in the
 * parSpreads - redirect the call to it */
const DefaultRatesSP QuantoCDSParSpreads::getCachedDefaultRates(
    const ICDSParSpreads* cds,
    const TypeOfEntry     entryType) const
{
    return parSpreads->getCachedDefaultRates(cds, entryType);
}
/** Provides a facility for previously constructed DefaultRates
 *  objects to be added to the cache
 * Here we are in an QuantoCDSParSpreads so the cache will be kept in the 
 * parSpreads - redirect the call to it */
const bool QuantoCDSParSpreads::cacheDefaultRates(const ICDSParSpreads* cds,
                                                  const TypeOfEntry     entryType,
                                                  const DefaultRatesSP  entry)
{
    return parSpreads->cacheDefaultRates(cds, entryType, entry);
}


/** Returns a DefaultRates object.
 * This method does NOT use the (potentially availabe) caching mechanism
 * which avoids the slow process of calculating default rates everytime.
 * It is invoked from inside the cache to calculate the default rates
 * for the first time and its direct use from any external classes is
 * strongly DISCOURAGED: the "defaultRates" method should be used instead.
 * Here we are in QuantoCDSParSpreads so the way to compute the clean spreads
 * is bootstrapping... so do that */
const DefaultRatesSP QuantoCDSParSpreads::computeDefaultRatesForCache() const {
    return CONST_POINTER_CAST<DefaultRates>(
            algorithm->currencyAdjust(domYC, forYC, parSpreads.getSP(),
                                      fx.getSP(),
                                      corrFxCDS, corrFxDom, corrFxFor,
                                      corrCDSDom, corrCDSFor, corrDomFor));
}


/** Resets the pointer to the cache - The cache will not be used
 * in this object from this point onwards */
void QuantoCDSParSpreads::resetCache() const {
    parSpreads->resetCache();
}

/** Hash code function used for caching - the cache needs improved
 * performance compared with the default "hashCode" function in CObjet: only
 * the required components are hashed (see comment in equalToOpt below) */
int QuantoCDSParSpreads::hashCodeOpt() const {
    if (!theHashCode.valid) {
        theHashCode.hc = hashCode();
        // The algorithm is not registered, so needs to be hashed by hand.
        theHashCode.hc ^= algorithm->hashCodeOpt();
        theHashCode.valid = true;
    }
    return theHashCode.hc;
}


/** Comparison function used for caching - the cache needs improved
 * performance compared with the default "equalTo" function in CObjet: only
 * the required components are compared. This means that the objects being
 * compared are not required to be identical: only equal as far as the
 * default rates calculation is concerned (e.g, if the volType is not used
 * in the calculation, the equalToOpt method should not fail if the volTypes
 * of both objects are different) */
bool QuantoCDSParSpreads::equalToOpt(const ICDSParSpreads* cds) const{
    bool equal = equalTo(cds);
    if (equal) {
        // The algorithm is not registered, so needs to be compared by hand.
        // However, only bother doing it if everything else is equal so far

        // If here, the static cast will succeed
        const QuantoCDSParSpreads* quantoCds = STATIC_CAST(QuantoCDSParSpreads, cds);
        equal = algorithm->equalToOpt(quantoCds->algorithm.get());
    }
    return equal;
}


//------------------------------
// Duration::IParHandler methods
//------------------------------
//return benchmarks for par instruments
//these will determine the durations to calculate unless specified
ExpiryArrayConstSP QuantoCDSParSpreads::getParBenchmarks() const {
    return getParSpreadsExpiries();
}

//return closed form solution
ExpiryResultArraySP QuantoCDSParSpreads::getDuration(const Duration* duration) const {
    return CDSParSpreadsBase::CDSGetDuration(duration,
                                             getParBenchmarks(),
                                             getYCName(),
                                             getName());
}

//------------------------------------------
//  IBadDayAdjuster method
//------------------------------------------
/** Returns "date" bad day adjusted using the bad day convention
 * and holidays in this object */
DateTime QuantoCDSParSpreads::badDayAdjust(const DateTime& date) const {
    return parSpreads->badDayAdjust(date);
}

/** Add a number of business days to a date */
DateTime QuantoCDSParSpreads::addBusinessDays(const DateTime& from, int busDays) const {
    return parSpreads->addBusinessDays(from, busDays);
}

double QuantoCDSParSpreads::survivalProb(
                            const DateTime& d1, 
                            const DateTime& d2) const
{
    return this->defaultRates()->calcTotalDefaultPV(d1, d2);
}



/**Probability of no default in [this->baseDate(),dt], given no default at base date.*/
double QuantoCDSParSpreads::survivalProb(const DateTime& dt) const
{
    return survivalProb(getValueDate(), dt);
}

/*TODO: this currently return PV with no delay - need to implement fixed date recovery.*/
double QuantoCDSParSpreads::protectionPV(const DateTime& paymentDate, 
                                         const DateTime& startDt, 
                                         const DateTime& endDt,
                                         RecoveryType    recTyp,
                                         const DateTime& recoveryDate) const
{
    if (recTyp == IDiscountCurveRisky::RECOVER_0)
    {
        return .0;
    }

    double  protectionValue = (1 - survivalProb(startDt, endDt)) * this->getYieldCurve()->pv(paymentDate, recoveryDate);

    if (recTyp == IDiscountCurveRisky::RECOVER_R)
    {
        protectionValue *= getRecovery();
    }
    else if (recTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
    {
        protectionValue *= (1-getRecovery());
    }

    return protectionValue;
}


/**Returns the value at paymentDate (and conditional on no default before then) 
 * of a contingent claim that pays recTyp (1, R, 1-R or, trivially, 0)
 * in case of default between startDate and endDate, and zero otherwise.
 * The protection pay-out is paid at date recoveryTiming.toDate(defaultDate).
 * This allows a recovery delay if the Expiry is a relative date, or an
 * absolute recovery if it's a fixed date (e.g. for use in recovery-at-maturity
 * calculations for unconditional settlement. If the Expiry pointer is null,
 * then no delay.
 * Whether integration is done continuously or discretely and what, 
 * if any approximations,
 * are used is up to the curve implementation. */
double QuantoCDSParSpreads::protectionPV(const DateTime& paymentDate, 
                                         const DateTime& startDt, 
                                         const DateTime& endDt,
                                         RecoveryType    recTyp,
                                         double          recoveryDelay) const
{   
    if (recTyp == IDiscountCurveRisky::RECOVER_0)
    {
        return .0;
    }

    DateTimeArray               noDates = DateTimeArray();
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual360());
    double                      protectionValue = .0;
    double                      accrualValue = .0;

    CDSPricer   cdsPricer 
        = CDSPricer(.0,                     // double couponRate,
                    noDates,                // paymentDates,
                    this->defaultRates(),   // defaultRates,
                    1.0,                    // notional,
                    0.0,                    // recovery,
                    false,                  // payAccruedFeed,
                    dcc,                    // swpAccrualDCC
                    paymentDate,            // valueDate,
                    startDt,                // protectionStartDate,
                    endDt,                  // protectionEndDate,
                    startDt,                // accruedStartDate,
                    getYieldCurve(),
                    getPrepayCurve());
    
    if (recoveryDelay == .0)
    {
        cdsPricer.calculateDefaultPayment(  false, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::NO_DELAY_DISCOUNTING);
    }
    else
    {
        cdsPricer.calculateDefaultPayment(  false, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::DELAYED_DISCOUNTING,
                                            recoveryDelay);
    }

    if (recTyp == IDiscountCurveRisky::RECOVER_R)
    {
        protectionValue *= getRecovery();
    }
    else if (recTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
    {
        protectionValue *= (1-getRecovery());
    }

    return protectionValue;
}

/*TODO: this currently return PV with no delay - need to implement fixed date recovery.*/
double QuantoCDSParSpreads::annuityPV(
                        const CashFlowArray&    payments,
                        const DateTime&         paymentDate,
                        RecoveryType            accruedRecTyp,
                        const DateTime&         recoveryDate,
                        DateTime                accrueStartDate) const 
{
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual365F());
    CashFlowArray*              paymentsCopy = new CashFlowArray(payments);
    CashFlowArraySP             paymentsSP = CashFlowArraySP(paymentsCopy);
    double                      protectionValue = .0;
    double                      accrualValue = .0;
    double                      couponValue = .0;

    if (payments.empty())
    {
        return .0;
    }
    if (accrueStartDate.empty())
      {
	accrueStartDate = payments[0].date;
      }
    if (accruedRecTyp != IDiscountCurveRisky::RECOVER_0)
    {
        CDSPricer   cdsPricer 
            = CDSPricer(paymentsSP,                // paymentDates,
                        this->defaultRates(),   // defaultRates,
                        1.0,                    // notional,
                        0.0,                    // recovery,
                        true,                   // payAccruedFeed,
                        dcc,                    // swpAccrualDCC
                        paymentDate,            // valueDate,
                        accrueStartDate,            // protectionStartDate,
                        payments[payments.size()-1].date,            // protectionEndDate,
                        accrueStartDate,            // accruedStartDate,
                        getYieldCurve());
        
        cdsPricer.calculateDefaultPayment(  true, 
                                            protectionValue, 
                                            accrualValue,
                                            CDSPricer::NO_TIME_DISCOUNTING);
        
        
        if (accruedRecTyp == IDiscountCurveRisky::RECOVER_R)
        {
            accrualValue *= getRecovery();
        }
        else if (accruedRecTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
        {
            accrualValue *= (1-getRecovery());
        }
    }
    
    couponValue = pv(payments, paymentDate);

    return couponValue + accrualValue * this->getYieldCurve()->pv(paymentDate, recoveryDate);
}

/**Returns the value at paymentDate (and conditional on no default before then) 
 * of a sequence of payments, with simple linear accrued-interest
 * recovery (so the claim on a coupon, C, payable at the end of accrual period (S,T) given
 * default at t is C(t-S)/(T_S)) in case of default defined by accruedRecType as for 
 * protectionPV. This allows the curve to compute default-accrual PVs itself
 * in a way which reflects the type of curve (e.g. flat forwards, etc.)
 * The accrual periods are the intervals between consecutive payment dates. To have a first
 * accrual period, you should set the first payment to zero, and then the first date is
 * the beginning of the first accrual period.
 * Payments before paymentDate are ignored, except to the extent that they affect accrued
 * interest due.
 * Unless recType==RECOVER_0, the cash-flow dates MUST BE IN INCREASING ORDER, or the
 * default accrual calculations will not be very meaningful! if recType==RECOVER_0,
 * this method should return the same value as pv(payments, paymentDate), inherited
 * from IDiscountCurve. */
double QuantoCDSParSpreads::annuityPV(const CashFlowArray& payments,
                                      const DateTime&      paymentDate,
                                      RecoveryType         accruedRecTyp,
                                      double               recoveryDelay,
                                      DateTime             accrueStartDate) const
{
    DayCountConventionConstSP   dcc = DayCountConventionConstSP(new Actual365F());
    CashFlowArray*              paymentsCopy = new CashFlowArray(payments);
    CashFlowArraySP             paymentsSP = CashFlowArraySP(paymentsCopy);
    double                      protectionValue = .0;
    double                      accrualValue = .0;
    double                      couponValue = .0;

    if (payments.empty())
    {
        return .0;
    }

    if (accrueStartDate.empty())
      {
	accrueStartDate = payments[0].date;
      }
    
    if (accruedRecTyp != IDiscountCurveRisky::RECOVER_0)
    {
        CDSPricer   cdsPricer 
            = CDSPricer(paymentsSP,                // paymentDates,
                        this->defaultRates(),   // defaultRates,
                        1.0,                    // notional,
                        0.0,                    // recovery,
                        true,                   // payAccruedFeed,
                        dcc,                    // swpAccrualDCC
                        paymentDate,            // valueDate,
                        accrueStartDate,            // protectionStartDate,
                        payments[payments.size()-1].date,            // protectionEndDate,
                        accrueStartDate,            // accruedStartDate,
                        getYieldCurve());
        
        if (recoveryDelay == .0)
        {
            cdsPricer.calculateDefaultPayment(  true, 
                                                protectionValue, 
                                                accrualValue,
                                                CDSPricer::NO_DELAY_DISCOUNTING);
        }
        else
        {
            cdsPricer.calculateDefaultPayment(  true, 
                                                protectionValue, 
                                                accrualValue,
                                                CDSPricer::DELAYED_DISCOUNTING,
                                                recoveryDelay);
        }

        if (accruedRecTyp == IDiscountCurveRisky::RECOVER_R)
        {
            accrualValue *= getRecovery();
        }
        else if (accruedRecTyp == IDiscountCurveRisky::RECOVER_1_MINUS_R)
        {
            accrualValue *= (1-getRecovery());
        }
    }
    
    couponValue = pv(payments, paymentDate);

    return couponValue + accrualValue;
}

/** Compute discount factor between two dates. This is the price
    * at time date1 of a zero-coupon (discount) bond that pays 1
    * at time date2 if it still exists then. 
    * Note that, if this bond is risky, this method
    * in IDiscountCurveRisky means the PV of such a bond which knocks
    * out on default with no recovery at all, and the price is 
    * contingent on no default before date1.
    * @param date1 payment settlement date (conditional on existence at date1)
    * @param date2 payment/zero-coupon bond maturity date
    * @return Discount factor between date1 & date2
    */
double QuantoCDSParSpreads::risklessPV(const DateTime& date1, 
                                       const DateTime& date2) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(date1, date2);
}
    
/** Compute price for settlement today of a zero-coupon bond
    * maturing at date. Note settlement is to TODAY, not to
    * the curve spot date. This is because some curves
    * may have ambiguous spot-dates - for example should a combined
    * credit and rates curve have spot date T+1 or T+2?
    * @param date To get discount factor/PV for
    * @return Discount factor between today & given date
    */
double QuantoCDSParSpreads::risklessPV(const DateTime& date) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(date);
}
    
/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double QuantoCDSParSpreads::risklessPV(const CashFlowArray& cashFlows,
                                       const DateTime&      baseDate) const
{
    YieldCurveConstSP dsc = getYieldCurve();
    return dsc->pv(cashFlows, baseDate);
}


double QuantoCDSParSpreads::pv(const DateTime& date1, 
                               const DateTime& date2) const 
{
    return this->defaultRates()->calcTotalDefaultPV(date1, date2) * getYieldCurve()->pv(date1, date2);
}
    
/** Compute price for settlement today of a zero-coupon bond
 * maturing at date. Note settlement is to TODAY, not to
 * the curve spot date. This is because some curves
 * may have ambiguous spot-dates - for example should a combined
 * credit and rates curve have spot date T+1 or T+2?
 * @param date To get discount factor/PV for
 * @return Discount factor between today & given date
 */
double QuantoCDSParSpreads::pv(const DateTime& date) const
{
    return pv(getValueDate(), date);
}
    
/** Calculates present value to baseDate of supplied cash flows conditional
    on continued existence (i.e. no default for a risky curve)
    Cash-flows on or before baseDate are ignored. No
    ordering of the cashflows is assumed. */
double QuantoCDSParSpreads::pv(const CashFlowArray& cashFlows,
                               const DateTime&      baseDate) const
{   
    double      value = .0;
        
    for (int i=0; i < cashFlows.size(); ++i)
        {
            const DateTime cfDate = cashFlows[i].date;

            if (cfDate.isGreater(baseDate))
                {
                    value += cashFlows[i].amount * pv(baseDate, cfDate);
                }
        }

    return value;
}

//utility for IDiscountCurveRisky methods
YieldCurveConstSP QuantoCDSParSpreads::getYieldCurve() const
{
    //we want the instrument discount curve
    return fx->getBaseCcy();
}

DateTime QuantoCDSParSpreads::getValueDate() const
{
    return fx->getToday();
}

CClassConstSP const QuantoCDSParSpreads::TYPE = CClass::registerClassLoadMethod(
    "QuantoCDSParSpreads", typeid(QuantoCDSParSpreads), load);

DRLIB_END_NAMESPACE
