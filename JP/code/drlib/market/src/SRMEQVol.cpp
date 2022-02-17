#include "edginc/config.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/MaturityTimePeriod.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FlatVol.hpp"


DRLIB_BEGIN_NAMESPACE

// This is crude - perhaps someone has done something better elsewhere? XXX
ExpirySP SRMEQVol::toExpiry(string& expiryStr) 
{
    static const string routine = "SRMEQVol::toExpiry";

    ExpirySP expiry(0);

    // Standard date time format has a space separating date from time
    int  breaktime = expiryStr.find(" ");
    if (breaktime == (int)std::string::npos) {
        // try as a MaturityPeriod
        try {
            expiry = ExpirySP(new MaturityPeriod(expiryStr));
        } catch (exception) {
            expiry = ExpirySP(   );
        }
    } else {
        string dt = expiryStr.substr(0, breaktime); 
        string tm = expiryStr.substr(breaktime+1);

        try {
            int breaktimeTest = expiryStr.find("-");
            if (breaktimeTest > 0) {
                expiry = ExpirySP(new BenchmarkDate(DateTime(dt, tm)));
            } else {
                expiry = ExpirySP(new MaturityTimePeriod(dt, DateTime::timeConvert(tm)));
            }
        } catch (exception) {
            expiry = ExpirySP(   );
        }
    }
    
    if (!expiry) {
        throw ModelException(routine, 
                             "Unrecognised expiry : " + expiryStr);
    }
    return expiry;
}

/** Turn BM strings into date times */
void SRMEQVol::cacheDates()
{
    compVolDate.resize(compVolBMs.size());
    compVolExpiry = ExpiryArraySP(new ExpiryArray(compVolBMs.size()));
    compVolResetDate.resize(compVolBMs.size());
    for (int i = 0; i < compVolDate.size(); i++){
        ExpirySP expiry(toExpiry(compVolBMs[i]));
        compVolDate[i] = expiry->toDate(today);
        (*compVolExpiry)[i] = expiry;

        if (compVolResetBMs.empty()){
            compVolResetDate[i] = compVolDate[i];
        } else {
            ExpirySP expiry(toExpiry(compVolResetBMs[i]));
            compVolResetDate[i] = expiry->toDate(today);
        }
    }
    spotVolDate.resize(spotVolBMs.empty()? 0: spotVolBMs.size());
    spotVolExpiry = ExpiryArraySP(new ExpiryArray(spotVolBMs.size()));
    for (int j = 0; j < spotVolDate.size(); j++){
        ExpirySP expiry(toExpiry(spotVolBMs[j]));
        spotVolDate[j] = expiry->toDate(today);
        (*spotVolExpiry)[j] = expiry;
    }
    smileDate.resize(smileBMs.size());
    smileExpiry = ExpiryArraySP(new ExpiryArray(smileBMs.size()));
    for (int k = 0; k < smileDate.size(); k++){
        ExpirySP expiry(toExpiry(smileBMs[k]));
        smileDate[k] = expiry->toDate(today);
        (*smileExpiry)[k] = expiry;
    }
    try{
        DateTime::ensureIncreasing(compVolDate, "compVolDate", false); 
        DateTime::ensureIncreasing(compVolResetDate, "compVolResetDate", false);
        DateTime::ensureIncreasing(spotVolDate, "spotVolDate", false); 
        DateTime::ensureIncreasing(smileDate, "smileDate", false);
    } catch (exception& e){
        throw ModelException(e, "SRMEQVol::cacheDates", 
                             "For SRMEQVol with name "+name);
    }
}

/** calculates the trading time between two dates */
double SRMEQVol::calcTradingTime(
    const DateTime &date1, 
    const DateTime &date2) const
{
    // no time metric supported at the moment
    return date1.yearFrac(date2);
}

/** retrieve time measure for the vol */
TimeMetricConstSP SRMEQVol::GetTimeMetric()const
{
    HolidaySP noHols(Holiday::noHolidays());
    return TimeMetricConstSP(new TimeMetric(1.0, noHols.get()));
}

void SRMEQVol::validatePop2Object()
{
    static const string method("SRMEQVol::validatePop2Object");

    if (!compVolResetBMs.empty() && 
        compVolResetBMs.size() != compVolBMs.size()){
        throw ModelException(method, "compVolBMs (#= " + 
                             Format::toString(compVolBMs.size()) +
                             ") and compVolResetBMs (#= " + 
                             Format::toString(compVolResetBMs.size()) +
                             ") must be of the same length for "+name);
    }
    if (compVol.size() != compVolBMs.size()) {
        throw ModelException(method, "compVol (#= " + 
                             Format::toString(compVol.size()) +
                             ") and compVolBMs (#= " + 
                             Format::toString(compVolBMs.size()) +
                             ") must be of the same length for "+name);
    }
    if (!spotVolBMs.empty() &&
        (spotVol.size() != spotVolBMs.size())) {
        throw ModelException(method, "spotVol (#= " + 
                             Format::toString(spotVol.size()) +
                             ") and spotVolBMs (#= " + 
                             Format::toString(spotVolBMs.size()) +
                             ") must be of the same length for "+name);
    }
    if (compVolBMs.empty() &&
        spotVolBMs.empty()) {
        throw ModelException(method, "Require at least one compVol or spotVol for "+name);
    }
    if (!compVolBMs.empty() && smileBMs.empty()) {
        throw ModelException(method, "Need at least one smile BM for "+name);
    }
    if (!spotVolBMs.empty() && smileBMs.empty()) {
        throw ModelException(method, "Need at least one smile BM for "+name);
    }
    if (smileBMs.size() != smileA1.size() ||
        smileBMs.size() != smileA2.size() ||
        smileBMs.size() != smileA3.size()) {
        throw ModelException(method, "smileBMs (#= " + 
                             Format::toString(smileBMs.size()) +
                             "), smileA1 (#= " + 
                             Format::toString(smileA1.size()) +
                             "), smileA2 (#= " + 
                             Format::toString(smileA2.size()) +
                             "), smileA3 (#= " + 
                             Format::toString(smileA3.size()) +
                             ") must all be of the same length for "+name);
    }
    for (int i = 0; i < smileA3.size(); i++){
        if (Maths::isNegative(smileA3[i])){
            throw ModelException(method, 
                                 "Each smileA3 parameter must not be negative but #" +
                                 Format::toString(i+1) +
                                 " = " + Format::toString(smileA3[i]) + " for " +name);
        }
        if (Maths::isZero(smileA3[i])){ /* flat local vol */ 
            if (!Maths::isZero(smileA2[i]) || 
                !Maths::isZero(smileA1[i])) {
                /* can't fit if no variation allowed */
                throw ModelException(method, "If smileA3 (#" + Format::toString(i+1) +
                                     ") is zero, smileA1 (" + Format::toString(smileA1[i]) +
                                     ") and smileA2 (" + Format::toString(smileA2[i]) +
                                     ") must be zero too, for " +name);
            }
        }
    }
}

/** populate from market cache */
void SRMEQVol::getMarket(
    const IModel* /*model*/, 
    const MarketData* market) 
{
    market->GetReferenceDate(today);
    cacheDates();
}

/** Shifts the object using given shift. */
bool SRMEQVol::sensShift(Theta* shift)
{
    const DateTime& newDate = shift->rollDate(today);
    if (newDate != today){
        cacheDates();
    }
    return false; // nothing else to shift
}

/** I guess we can do something with the spot (=>fwd ) vols? */
CVolProcessed* SRMEQVol::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      /*asset*/) const
{
    if (VolRequestRaw::TYPE->isInstance(volRequest) ||
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<SRMEQVol*>(this);
    } else if (ATMVolRequest::TYPE->isInstance(volRequest) ||
               LinearStrikeVolRequest::TYPE->isInstance(volRequest)) {
        // FWD_AT_MAT request for protected assets will ask for this
        // Looking just for something that doesn't break
        HolidaySP noHols(Holiday::noHolidays());
        TimeMetricSP tm(new TimeMetric(1.0, noHols.get()));
        double flatFXVol = !compVol.empty()? compVol[0] : spotVol[0];
        FlatVolSP flatVol(new FlatVol(name, 
                                      today, // baseDate
                                      tm.get(),
                                      flatFXVol));
        return flatVol->getProcessedVol(volRequest, 0);
    }
    throw ModelException("SRMEQVol::getProcessedVol", 
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported for " +name);
}

/** Returns name of vol */
string SRMEQVol::getName() const
{
    return name;
}

/** returns merged benchmark dates (spot vol and comp vol) */
DateTimeArray SRMEQVol::getVolBmDates() const 
{
    return DateTime::merge(spotVolDate,compVolDate); 
}

ExpiryArraySP SRMEQVol::getSmileParamExpiries(const IObject* obj)
{
    const SRMEQVol* srmeq_vol = dynamic_cast<const SRMEQVol*>(obj);
    if (obj){
        if (!srmeq_vol->smileExpiry){
            throw ModelException("SRMEQVol::getSmileParamExpiries", "Internal error: null smileExpiry");
        }
        return ExpiryArraySP(copy(srmeq_vol->smileExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMEQVol::getSmileParamExpiries", "Internal error: obj is not of the desired type");
}

ExpiryArraySP SRMEQVol::getCompVolExpiries(const IObject* obj)
{
    const SRMEQVol* srmeq_vol = dynamic_cast<const SRMEQVol*>(obj);
    if (obj){
        if (!srmeq_vol->compVolExpiry){
            throw ModelException("SRMEQVol::getSmileParamExpiries", "Internal error: null compVolExpiry");
        }
        return ExpiryArraySP(copy(srmeq_vol->compVolExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMEQVol::getCompVolExpiries", "Internal error: obj is not of the desired type");
}

ExpiryArraySP SRMEQVol::getSpotVolExpiries(const IObject* obj)
{
    const SRMEQVol* srmeq_vol = dynamic_cast<const SRMEQVol*>(obj);
    if (obj){
        if (!srmeq_vol->spotVolExpiry){
            throw ModelException("SRMEQVol::getSmileParamExpiries", "Internal error: null spotVolExpiry");
        }
        return ExpiryArraySP(copy(srmeq_vol->spotVolExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMEQVol::getSpotVolExpiries", "Internal error: obj is not of the desired type");
}

namespace SRMBounds {

    // for use in calibrator and initial guess
    const double volLowerBound = 0.01;
    const double volUpperBound = 1.;
    const double A1LowerBound = -15.;
    const double A1UpperBound = 15.;
    const double A2LowerBound = 0.;
    const double A2UpperBound = 10000.;
    const double A3LowerBound = 0.;
    const double A3UpperBound = 10.;
}

void SRMEQVol::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SRMEQVol, clazz);
    SUPERCLASS(CVolBase);
    IMPLEMENTS(IVolProcessed);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(Theta::IShift);
    IMPLEMENTS(TweakableWith<CompositeVegaParallelTweak>);
    IMPLEMENTS(PointwiseTweakableWith<CompositeVegaPointwiseTweak>);
    IMPLEMENTS(TweakableWith<SpotVegaParallelTweak>);
    IMPLEMENTS(PointwiseTweakableWith<SpotVegaPointwiseTweak>);
    IMPLEMENTS(IDynamicsParameter);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "SRMEQVol identifier");
    FIELD(compVolBMs, "Composite vol benchmarks");
    FIELD(compVolResetBMs, "Composite vol maturity benchmarks");
    FIELD_MAKE_OPTIONAL(compVolResetBMs);
    FIELD(compVol, "Composite vols");
    FIELD(spotVolBMs, "Spot vol benchmarks (optional)");
    FIELD_MAKE_OPTIONAL(spotVolBMs);
    FIELD(spotVol, "Spot Vols (optional)");
    FIELD_MAKE_OPTIONAL(spotVol);
    FIELD(smileBMs, "EQ Smile benchmarks");
    FIELD(smileA1, "ATM skew");
    FIELD(smileA2, "ATM crv");
    FIELD(smileA3, "maximum variation");
    FIELD_NO_DESC(today);
    FIELD_MAKE_TRANSIENT(today);
    FIELD_NO_DESC(compVolDate);
    FIELD_MAKE_TRANSIENT(compVolDate);
    FIELD_NO_DESC(compVolExpiry);
    FIELD_MAKE_TRANSIENT(compVolExpiry);
    FIELD_NO_DESC(compVolResetDate);
    FIELD_MAKE_TRANSIENT(compVolResetDate);
    FIELD_NO_DESC(spotVolDate);
    FIELD_MAKE_TRANSIENT(spotVolDate);
    FIELD_NO_DESC(spotVolExpiry);
    FIELD_MAKE_TRANSIENT(spotVolExpiry);
    FIELD_NO_DESC(smileDate);
    FIELD_MAKE_TRANSIENT(smileDate);
    FIELD_NO_DESC(smileExpiry);
    FIELD_MAKE_TRANSIENT(smileExpiry);

    ClassSetAcceptMethod(acceptCriticalDateCollector);
    
    
    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA1",
        new Range(ClosedBoundary(SRMBounds::A1LowerBound), ClosedBoundary(SRMBounds::A1UpperBound)),
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA2",
        new Range(ClosedBoundary(SRMBounds::A2LowerBound), ClosedBoundary(SRMBounds::A2UpperBound)),
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA3",
        new Range(ClosedBoundary(SRMBounds::A3LowerBound), ClosedBoundary(SRMBounds::A3UpperBound)),  // engine can be unstable if A3 > 5
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "compVol",
        new Range(OpenBoundary(0), Infinity(Infinity::Plus)),
        getCompVolExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "spotVol",
        new Range(ClosedBoundary(SRMBounds::volLowerBound), ClosedBoundary(SRMBounds::volUpperBound)), // min 1%, max 100% 
        getSpotVolExpiries);
}

void SRMEQVol::acceptCriticalDateCollector(
    const SRMEQVol* vol,
    CriticalDateCollector* collector) 
{
    // SRM wants the smile dates
    collector->addDates(vol->smileDate, vol->getClass());
}


/****************************************
 * Composite volatility - Parallel tweak 
 ****************************************/
string SRMEQVol::sensName(CompositeVegaParallelTweak* /*shift*/) const 
{
    return name;
}

bool SRMEQVol::sensShift(CompositeVegaParallelTweak* shift) 
{
    static const string method = "SRMEQVol::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            for (int i = 0; i < compVol.size(); ++i) {
                compVol[i] += shiftSize;
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, method,
                             "CompositeVegaParallelTweak failed for " + getName());
    }
    return false; // No sub-components need tweaking		   
}



/*****************************************
 * Composite volatility - Pointwise tweak 
 *****************************************/
string SRMEQVol::sensName(CompositeVegaPointwiseTweak* /*shift*/) const 
{
    return name;
}

bool SRMEQVol::sensShift(CompositeVegaPointwiseTweak* shift) 
{
    static const string method = "SRMEQVol::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            int i = shift->getExpiry()->search(compVolExpiry.get());
            compVol[i] += shiftSize;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method,
                             "CompositeVegaPointwiseTweak failed for " + getName());
    }
    return false; // No sub-components need tweaking
}

ExpiryArrayConstSP SRMEQVol::sensExpiries(
    CompositeVegaPointwiseTweak* /*shift*/) const 
{
    return compVolExpiry;
}

/***********************************
 * Spot volatility - Parallel tweak 
 ***********************************/
string SRMEQVol::sensName(SpotVegaParallelTweak* /*shift*/) const 
{
    // Return the object's name if the optional field spotVol is present,
    // or an empty string otherwise 
    return (spotVol.size() > 0) ? name : "";
}

bool SRMEQVol::sensShift(SpotVegaParallelTweak* shift) 
{
    static const string method = "SRMEQVol::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            for (int i = 0; i < spotVol.size(); ++i) {
                spotVol[i] += shiftSize;
            }
        }
    } 
    catch (exception& e) {
        throw ModelException(e, method,
                             "SpotVegaParallelTweak failed for " + getName());
    }
    return false; // No sub-components need tweaking		   
}


/************************************
 * Spot volatility - Pointwise tweak 
 ************************************/
string SRMEQVol::sensName(SpotVegaPointwiseTweak* /*shift*/) const 
{
    // Return the object's name if the optional field spotVol is present,
    // or an empty string otherwise 
    return (spotVol.size() > 0) ? name : "";
}

bool SRMEQVol::sensShift(SpotVegaPointwiseTweak* shift) 
{
    static const string method = "SRMEQVol::sensShift";

    try {
        double shiftSize = shift->getShiftSize();

        if (!Maths::isZero(shiftSize)) {
            int i = shift->getExpiry()->search(spotVolExpiry.get());
            spotVol[i] += shiftSize;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method,
                             "SpotVegaPointwiseTweak failed for " + getName());
    }
    return false; // No sub-components need tweaking
}


ExpiryArrayConstSP SRMEQVol::sensExpiries(SpotVegaPointwiseTweak* /*shift*/) const 
{
    return spotVolExpiry;
}


CClassConstSP const SRMEQVol::TYPE = CClass::registerClassLoadMethod(
    "SRMEQ::Vol", typeid(SRMEQVol), load);

/**************************/
bool SRMEQVolLoad(void) {
    return SRMEQVol::TYPE != 0;
}

DRLIB_END_NAMESPACE
