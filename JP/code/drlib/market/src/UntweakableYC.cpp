//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UntweakableYC.hpp
//
//   Description : A yield curve that is built from a zero curve
//                 Motivation: for porting SRM3 tests
//                 DO NOT USE THIS FOR ANYTHING ELSE
//
//   Author      : Mark A Robson
//
//   Date        : 26 May 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UntweakableYC.hpp"
#include "edginc/CompoundBasis.hpp"
#include "edginc/RateConversion.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/Validator.hpp"

DRLIB_BEGIN_NAMESPACE

static void throwError(){
    throw ModelException("UntweakableYC", "You really can't shift this "
                         "yield curve!");
}

/** overrides CObject version to allow for easy default */
bool UntweakableYC::accept(ICollector* collector) const{
    if (!CClass::invokeAcceptMethod(this, collector)){
        // if no method registered try  vol
        return irVol->accept(collector);
    }
    return false;
}

/** Records name of iso code against yc name in market data object.
    It is invoked ONCE only
    - immediately after this object is placed in the cache. */
void UntweakableYC::initialise(MarketData* market){
    market->setYieldCurveISOCode(name, ccy);
}

/** populate from market cache */
void UntweakableYC::getMarket(const IModel* model, const MarketData* market) {
    // if specified, ask for the vol
    if (!irVol.isEmpty()){
        // NB This might result in a null vol if the model reckons we don't
        // need it
        irVol.getData(model, market);
    }
}

/** @return Yield curve's currency */
string UntweakableYC::getCcy() const{
    return ccy;
}

/** @return Yield curve's name - used to identify sensitivities */
string UntweakableYC::getName() const{
    return name;
}

/** @return Yield curve's spot date */
DateTime UntweakableYC::getSpotDate() const{
    return valueDate;
}

/** @return Yield curve's today date */
DateTime UntweakableYC::getToday() const {
    return today;
}

/** Useful accessor methods */
ExpiryArrayConstSP UntweakableYC::getExpiries() const
{
    return ExpiryArraySP(new ExpiryArray(0));
}

StringArrayConstSP UntweakableYC::getInstruments() const
{
    return StringArraySP(new StringArray(0));
}

/** Returns tradeDate + n business days where n = spot offset */
DateTime UntweakableYC::settles(const DateTime& tradeDate) const{
    // we haven't got any holidays ... (not yet anyway)
    return tradeDate.rollDate(spotOffset);
}

/** Returns rateMaturity->toDate(rateStart)  */
DateTime UntweakableYC::rateMaturity(
    const DateTime&         rateStart,
    const Expiry*           rateMaturity,
    const BadDayConvention* bdc) const { // optional
    // we have no holidays or bad day convention
    return rateMaturity->toDate(rateStart);
}


 //// switch to using growth curve
void UntweakableYC::setProjectionCurve(bool useEstimatingCurve) const {
    if (useEstimatingCurve && !growthZeroCurve){
        throw ModelException("UntweakableYC::setProjectionCurve",
                             "Growth zero curve has not been specified");
    }
    useDiscCurve = !useEstimatingCurve;
}

/** Returns name identifying yield curve for rho parallel */
string UntweakableYC::sensName(RateShift* shift) const{
    return ccy;
}

/** Shifts the object using given shift */
bool UntweakableYC::sensShift(RateShift* shift){
    throwError();
    return false;
}

/** Returns name identifying yield curve for additive or
 * multiplicative weighted shift */
string UntweakableYC::sensName(YCWeightedShift* shift) const{
    return ccy;
}

/** Shifts the object using given shift */
bool UntweakableYC::sensShift(YCWeightedShift* shift){
    throwError();
    return false;
}

/** Returns name identifying yield curve for rho parallel */
string UntweakableYC::sensName(const RateParallel* shift) const{
    return ccy;
}

/** Shifts the object using given shift */
TweakOutcome UntweakableYC::sensShift(const PropertyTweak<RateParallel>& shift){
    throwError();
    return TweakOutcome(shift.coefficient,false);
}

/** Returns the name of the yield curve - used to determine whether
    to tweak the object */
string UntweakableYC::sensName(const RatePointwise* shift) const{
    return ccy;
}

/** Return the array of expiries (ie maturities/benchmark dates) that
    need to be tweaked for this  yield curve */
ExpiryWindowArrayConstSP UntweakableYC::sensQualifiers(const RatePointwise* shift) const{
    throwError();
    return ExpiryWindow::series(ExpiryArrayConstSP());
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
TweakOutcome UntweakableYC::sensShift(const PropertyTweak<RatePointwise>& shift) {
    throwError();
    return TweakOutcome(shift.coefficient,false);
}

/** NOTE: IRRatePointwise is supported for exposure reporting only. */
/** Returns the name of the yield curve - used to determine whether
to tweak the object */
string UntweakableYC::sensName(const IRRatePointwise* shift) const{
    return name;
}

/** Return the array of expiries (ie maturities/benchmark dates) that
need to be tweaked for this yield curve for exposure reporting. */
ExpiryWindowArrayConstSP UntweakableYC::sensQualifiers(const IRRatePointwise* shift) const{
    const DateTimeArray &zcDates = zeroCurve->getDates();
    ExpiryArraySP benchmarkDatesSP(new ExpiryArray);
    for (int i = 0; i < zcDates.size(); ++i)
        benchmarkDatesSP->push_back(BenchmarkDateSP(new BenchmarkDate(zcDates[i])));
    return ExpiryWindow::series(benchmarkDatesSP);
}

/** Actual tweaking is not supported. */
TweakOutcome UntweakableYC::sensShift(const PropertyTweak<IRRatePointwise>& shift) {
    throwError(); 
    return TweakOutcome(shift.coefficient,false);
}

bool UntweakableYC::sensShift(Theta* shift){
    DateTime newDate(shift->rollDate(today));
    if (!newDate.equals(today)){
        throwError();
    }
    return true;
}

/** grab the dates used in the zero curve */
DateTimeArray UntweakableYC::zeroDates() const{
    return (useDiscCurve? zeroCurve: growthZeroCurve)->getDates();
}

/** Returns an processed vol - which combines the vol market data with the
instrument data in the volRequest */
CVolProcessed* UntweakableYC::getProcessedVol(
    const CVolRequest* volRequest) const{
    if (irVol.isEmpty()){
        throw ModelException(
            "UntweakableYC::getProcessedVol",
            "Pricing model requires IR Vol be supplied for "
            "UntweakableYC curve with name "+getName());
    }
    return irVol->getProcessedVol(volRequest, this);
}

// these to move to RiskyCurve
/** this function calculates the discount factor based on the assumption
    that the on which the recovery is based is provided externally.
    This allows to use different methodologies
    (PV, face value + accrued etc.) to be included easily */
double UntweakableYC::riskyPV(const DateTime& lodate,
                       const DateTime& hidate,
                       double          cashFlow,
                       double          recoveryNotional) const{
    return pv(lodate, hidate) * cashFlow;
}

double UntweakableYC::fwd(const DateTime&        payDate,
                       const DateTime&           refixDate, 
                       const Expiry*             rateMat,
                       const BadDayConvention*   bdc, // optional
                       const DayCountConvention* dcc,
                       const bool                isCMS) const{
    throw ModelException("UntweakableYC::fwd","Not implemented");
}


/** make a risky curve from a credit spread curve */
IYieldCurveSP UntweakableYC::makeRiskyCurve(
    const CreditSpreadCurve& spreadCurve,
    const  DateTime*    maturityDate) const
{
    throw ModelException("UntweakableYC::makeRiskyCurve","Not implemented");
}


CreditSpreadCurveSP UntweakableYC::makeCreditSpreadCurve(
	const string&        name,
	const CashFlowArray& defaultRates,
	double               recovery) const
{
    throw ModelException("UntweakableYC::makeCreditSpreadCurve","Not implemented");
}


CashFlowArraySP UntweakableYC::getRatesAndDates() const
{
    throw ModelException("UntweakableYC::getRatesAndDates","Not implemented");
}

IYieldCurveSP UntweakableYC::createForwardCurve(const DateTime& forwardDate) const
{
    throw ModelException("UntweakableYC::createForwardCurve","Not implemented");
}


IObject* UntweakableYC::defaultConstructor(){
    return new UntweakableYC();
}

/** Invoked when Class is 'loaded' */
void UntweakableYC::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(UntweakableYC, clazz);
    SUPERCLASS(YieldCurve);
    IMPLEMENTS(IDeterministicYieldCurve);
    IMPLEMENTS(ITweakableWithRespectTo<RateParallel>);
    IMPLEMENTS(ITweakableWithRespectTo<RatePointwise>);
    IMPLEMENTS(ITweakableWithRespectTo<IRRatePointwise>);
    IMPLEMENTS(YCWeightedShift::IShift);
    IMPLEMENTS(Theta::Shift);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ccy, "currency name");
    FIELD(name, "rate index");
    FIELD(today, "today");
    FIELD(spotOffset, "spot offset");
    FIELD(moneyMarketDayCount, "money market day count");
    FIELD(swapDayCount, "swap market day count");
    FIELD(swapFrequency, "swap frequency");
    FIELD(zeroCurve, "discount zero curve");
    FIELD(growthZeroCurve, "growth zero curve");
    FIELD_MAKE_OPTIONAL(growthZeroCurve);
    FIELD(valueDate, "value date");
    FIELD(irVol, "Interest rate volatility");
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_NO_DESC(useDiscCurve);
    FIELD_MAKE_TRANSIENT(useDiscCurve);

    // these 2 fields are not set if created via ZERO_CURVE_MAKE
    FIELD_MAKE_OPTIONAL(moneyMarketDayCount);
    FIELD_MAKE_OPTIONAL(swapDayCount);

    ClassSetAcceptMethod(acceptValueDateCollector);
    ClassSetAcceptMethod(acceptWrapperNameCollector);
    ClassSetAcceptMethod(acceptYieldNameCollector);
}

/** Compute discount factor between two dates
     * @param lodate Lower date
     * @param hidate Upper date
     * @return Discount factor between lodate & hidate
     */
double UntweakableYC::pv(
    const DateTime& lodate,
    const DateTime& hidate) const {
    static const string method = "UntweakableYC::pv";
    try {
        double lo = (useDiscCurve? zeroCurve: growthZeroCurve)
            ->discountFactor(lodate);
        double hi = (useDiscCurve? zeroCurve: growthZeroCurve)
            ->discountFactor(hidate);
        return (hi/lo);
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

/** Compute discount factor between value date and a date
     * @param date To get discount factor for
     * @return Discount factor between value date & given date
     */
double UntweakableYC::pv(const DateTime& date) const {
    return pv(today, date);
}

double UntweakableYC::parSwapRate(const DateTime& maturity) const
{
	MaturityPeriod mp(swapFrequency);
    return couponRate(valueDate, maturity, mp, true, swapDayCount.get());
}

DateTimeArraySP UntweakableYC::getDates() const
{
    DateTimeArraySP dates = DateTimeArraySP(new DateTimeArray(zeroDates()));
    return dates;
}

/** Interpolate zero coupon rate at a date
     * @param date Interpolate zero coupon rate here
     * @return Zero coupon rate at given date
     */
double UntweakableYC::zero(const DateTime& date) const {
    return (useDiscCurve? zeroCurve: growthZeroCurve)->zeroCouponRate(date);
}

 /** Interpolate a forward rate between two dates
  * see CompoundBasis for basis values
  */
double UntweakableYC::fwd(
    const DateTime&           lodate,
    const DateTime&           hidate,
    const DayCountConvention* dcc,
    int                       basis) const {
   static const string method = "UntweakableYC::fwd";
   try {
       double disc = pv(lodate, hidate);
       double rate = RateConversion::discountToRate(disc,
                                                    lodate,
                                                    hidate,
                                                    dcc,
                                                    basis);
       return rate;
   }
   catch (exception &e) {
       throw ModelException(e, method);
   }
}

/** Interpolate a forward rate between two dates
 * if CMS, uses swapFrequency for the compounding basis
 */
double UntweakableYC::fwd(
    const DateTime&           refixDate,
    const Expiry*             rateMaturity,
    const BadDayConvention*   bdc, // optional
    const DayCountConvention* dcc,
    const bool                isCMS) const {

    static const string method = "UntweakableYC::fwd";
    try {

        DateTime lodate = refixDate; //no holidays !
        DateTime hidate = rateMaturity->toDate(lodate); //no holidays !

        double rate;

        if (isCMS)
        {
            //this is not as sophisticated as GtoSwapRate2
            //but will match GtoSwapRate
            MaturityPeriod interval(12 / swapFrequency, "M");
            rate = couponRate(lodate,
                              hidate,
                              interval,
                              false, //stub at end
                              dcc);
        }
        else
        {
            double disc = pv(lodate, hidate);
            rate = RateConversion::discountToRate(disc,
                                                  lodate,
                                                  hidate,
                                                  dcc,
                                                  CompoundBasis::SIMPLE);
        }
        return rate;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }
}

void UntweakableYC::acceptValueDateCollector(
                        const UntweakableYC* yieldCurve,
                        CValueDateCollector* collector){
    collector->valueDateValidate(yieldCurve->today, yieldCurve->getName());
}

void UntweakableYC::acceptWrapperNameCollector(const UntweakableYC*  yc,
                                               WrapperNameCollector* collector){
    collector->addName(yc->getName());
}

void UntweakableYC::acceptYieldNameCollector(const UntweakableYC* yc,
                                             YieldNameCollector*  collector){
    collector->addName(yc->getName());
}


/** Returns a key used to optimise repeated calculations of
    discount factors */
YieldCurve::IKey* UntweakableYC::logOfDiscFactorKey() const{
    return (useDiscCurve? zeroCurve: growthZeroCurve)->logOfDiscFactorKey();
}


UntweakableYC::UntweakableYC(
    const DateTime&  pValueDate, 
    const DateTime&  pBaseDate, 
    const ZeroCurve& discounting, 
    const ZeroCurve* estimating)
: YieldCurve(TYPE), 
  today(pValueDate), 
  spotOffset(pBaseDate.daysDiff(pValueDate)), 
  zeroCurve(&discounting), 
  growthZeroCurve(estimating),
  valueDate(pBaseDate),
  useDiscCurve(true)
{
}

/* for reflection */
UntweakableYC::UntweakableYC(): YieldCurve(TYPE), spotOffset(0),
                                useDiscCurve(true) {}

CClassConstSP const UntweakableYC::TYPE = CClass::registerClassLoadMethod(
    "UntweakableYC", typeid(UntweakableYC), load);


/* external symbol to allow class to be forced to be linked in */
bool UntweakableYCLinkIn(){
    return true;
}

bool UntweakableYC::validate()
{
    if (! useDiscCurve)
        throw ModelException("UntweakableYC::validate", "Not implemented for growth zero curve");

    DateTimeArray dates = zeroCurve->getDates();
    for (int i=0; i<dates.getLength()-1; ++i)
    {
        double fwd = zeroCurve->fwd(dates[i], dates[i+1], moneyMarketDayCount.get(), swapFrequency);
        if (fwd < 0)
            return false;
    }

    return true;
}

DEFINE_VALIDATOR_TYPE(UntweakableYC);

DRLIB_END_NAMESPACE

