#include "edginc/config.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FlatFXVol.hpp"

DRLIB_BEGIN_NAMESPACE

void SRMFXVol::cacheDates(void){
    compVolDate.resize(compVolExpiry->size());
    compVolMatDate.resize(compVolExpiry->size());
    for (int i = 0; i < compVolDate.size(); i++){
        compVolDate[i] = (*compVolExpiry)[i]->toDate(today);
        if (!compVolMatExpiry){
            compVolMatDate[i] = compVolDate[i];
        } else {
            compVolMatDate[i] = (*compVolMatExpiry)[i]->toDate(today);
        }
    }
    spotVolDate.resize(!spotVolExpiry? 0: spotVolExpiry->size());
    for (int j = 0; j < spotVolDate.size(); j++){
        spotVolDate[j] = (*spotVolExpiry)[j]->toDate(today);
    }
    smileDate.resize(!smileExpiry? 0: smileExpiry->size());
    for (int k = 0; k < smileDate.size(); k++){
        smileDate[k] = (*smileExpiry)[k]->toDate(today);
    }
    if (forcedTimeOfDay) {
        forceTimeOfDay();
    }
    try{
        DateTime::ensureIncreasing(compVolDate, "compVolDate", false);
        DateTime::ensureIncreasing(compVolMatDate, "compVolMatDate", false);
        DateTime::ensureIncreasing(spotVolDate, "spotVolDate", false);
        DateTime::ensureIncreasing(smileDate, "smileDate", false);
    } catch (exception& e){
        throw ModelException(e, "SRMFXVol::cacheDates",
                             "For SRMFXVol with name "+name);
    }
}

/** force all times of day to be the same as today - unfortunate hack */
void SRMFXVol::forceTimeOfDay(void){
    forcedTimeOfDay = true;
    DateTime::setTimeOfDay(compVolDate, today.getTime());
    DateTime::setTimeOfDay(compVolMatDate, today.getTime());
    DateTime::setTimeOfDay(spotVolDate, today.getTime());
    DateTime::setTimeOfDay(smileDate, today.getTime());
}

SRMFXVol::SRMFXVol(string name,
    ExpiryArraySP compVolExpiry,
    ExpiryArraySP compVolMatExpiry,
    DoubleArray   compVol,
    ExpiryArraySP spotVolExpiry,
    DoubleArray   spotVol,
    ExpiryArraySP smileExpiry,
    DoubleArray   smileA1,
    DoubleArray   smileA2,
    DoubleArray   smileA3,
    DateTime      today,
    DateTimeArray compVolDate,
    DateTimeArray compVolMatDate,
    DateTimeArray spotVolDate,
    DateTimeArray smileDate
):
    name(name),
    SRMFXVolBase(VOL_TYPE),
    compVolExpiry(compVolExpiry),
    compVolMatExpiry(compVolMatExpiry),
    compVol(compVol),
    spotVolExpiry(spotVolExpiry),
    spotVol(spotVol),
    smileExpiry(smileExpiry),
    smileA1(smileA1),
    smileA2(smileA2),
    smileA3(smileA3),
    today(today),
    compVolDate(compVolDate),
    compVolMatDate(compVolMatDate),
    spotVolDate(spotVolDate),
    smileDate(smileDate),
    forcedTimeOfDay(false)
{}

/** calculates the trading time between two dates */
double SRMFXVol::calcTradingTime(const DateTime &date1,
                               const DateTime &date2) const{
    // no time metric supported at the moment
    return date1.yearFrac(date2);
}

/** retrieve time measure for the vol */
TimeMetricConstSP SRMFXVol::GetTimeMetric(void) const {
    HolidaySP noHols(Holiday::noHolidays());
    return TimeMetricConstSP(new TimeMetric(1.0, noHols.get()));
}

void SRMFXVol::validatePop2Object(void) {
    static const string method("SRMFXVol::validatePop2Object");
    if (compVolMatExpiry.get() &&
        compVolMatExpiry->size() != compVolExpiry->size()){
        throw ModelException(method, "compVolExpiry and compVolMatExpiry"
                             " must be of the same length for "+name);
    }
    for (int j = 0; j < compVolExpiry->size(); j++){
        if (!((*compVolExpiry)[j]) ||
            compVolMatExpiry.get() && !((*compVolMatExpiry)[j])){
            throw ModelException(method, "Null compVol expiry on index "+
                                 Format::toString(j)+" on "+name);
        }
    }
    for (int k = 0; spotVolExpiry.get() && k < spotVolExpiry->size(); k++){
        if (!((*spotVolExpiry)[k])){
            throw ModelException(method, "Null spotVolExpiry on index "+
                                 Format::toString(k)+" on "+name);
        }
    }
    if (smileExpiry.get()){
        for (int m = 0; m < smileExpiry->size(); m++){
            if (!((*smileExpiry)[m])){
                throw ModelException(method, "Null smileExpiry on index "+
                                     Format::toString(m)+" on "+name);
            }
        }
    }
    for (int i = 0; i < smileA3.size(); i++){
        if (Maths::isNegative(smileA3[i])){
            throw ModelException(method,
                                 "Each smileA3 parameter must be >=0");
        }
        if (Maths::isZero(smileA3[i])){ /* flat local vol */
            if (!Maths::isZero(smileA2[i]) ||
                !Maths::isZero(smileA1[i])) {
                /* can't fit if no variation allowed */
                throw ModelException(method, "If smileA3 is zero, "
                                     "smileA1 and smileA2 must be "
                                     "zero too.");
            }
        }
    }
    if (compVolExpiry->size() != compVol.size()){
        throw ModelException(method, "compVolExpiry and compVol arrays "
                             "must be the same length");
    }
    if (spotVolExpiry.get() && spotVolExpiry->size() != spotVol.size()){
        throw ModelException(method, "spotVolExpiry and spotVol arrays "
                             "must be the same length");
    }
    int smileExpirySize = !smileExpiry? 0: smileExpiry->size();
    if (smileExpirySize != smileA1.size() ||
        smileA1.size() != smileA2.size() ||
        smileA2.size() != smileA3.size()){
        throw ModelException(method, "smileExpiry and smileA1, smileA2, "
                             "smileA3 arrays must be the same length");
    }
}

/** populate from market cache */
void SRMFXVol::getMarket(const IModel* model, const MarketData* market) {
    market->GetReferenceDate(today);
    cacheDates();
}

/** Shifts the object using given shift. */
bool SRMFXVol::sensShift(Theta* shift){
    const DateTime& newDate = shift->rollDate(today);
    if (newDate != today){
        cacheDates();
    }
    return false; // nothing else to shift
}

/** returns merged benchmark dates (spot vol and comp vol) */
DateTimeArray SRMFXVol::getVolBmDates() const {
    return DateTime::merge(spotVolDate,compVolDate);
}

/** I guess we can do something with the spot (=>fwd ) vols? */
CVolProcessed* SRMFXVol::getProcessedVol(const CVolRequest* volRequest,
                                       const CAsset*      asset) const{
    if (VolRequestRaw::TYPE->isInstance(volRequest) ||
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<SRMFXVol*>(this);
    } else if (ATMVolRequest::TYPE->isInstance(volRequest) ||
               LinearStrikeVolRequest::TYPE->isInstance(volRequest)) {
        // FWD_AT_MAT request for protected assets will ask for this
        // Looking just for something that doesn't break
        HolidaySP noHols(Holiday::noHolidays());
        TimeMetricSP tm(new TimeMetric(1.0, noHols.get()));
        FlatFXVolSP flatVol( new FlatFXVol(name,
                                           today, // baseDate
                                           tm.get(),
                                           compVol[0]/*flatFXVol*/));
        return flatVol->getProcessedVol(volRequest, 0);
    }
    throw ModelException("SRMFXVol::getProcessedVol",
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported");
}


/** Returns name of vol */
string SRMFXVol::getName(void) const{
    return name;
}


ExpiryArraySP SRMFXVol::getSmileParamExpiries(const IObject* obj){
    const SRMFXVol* srmfxVol = dynamic_cast<const SRMFXVol*>(obj);
    if (obj){
        if (!srmfxVol->smileExpiry){
            // smileExpiry is optional
            return ExpiryArraySP(new ExpiryArray());
        }
        return ExpiryArraySP(copy(srmfxVol->smileExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMFXVol::getSmileParamExpiries", "Internal error: obj is not of the desired type");
}

ExpiryArraySP SRMFXVol::getCompVolExpiries(const IObject* obj){
    const SRMFXVol* srmfx_vol = dynamic_cast<const SRMFXVol*>(obj);
    if (obj){
        return ExpiryArraySP(copy(srmfx_vol->compVolExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMFXVol::getCompVolExpiries", "Internal error: obj is not of the desired type");
}

void SRMFXVol::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SRMFXVol, clazz); // FIX - used to be Vol.  Does this break public API?
    SUPERCLASS(FXVolBase);
    IMPLEMENTS(IVolProcessed);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(Theta::IShift);
    IMPLEMENTS(TweakableWith<CompositeVegaParallelTweak>);
    IMPLEMENTS(PointwiseTweakableWith<CompositeVegaPointwiseTweak>);
    IMPLEMENTS(TweakableWith<SpotVegaParallelTweak>);
    IMPLEMENTS(PointwiseTweakableWith<SpotVegaPointwiseTweak>);
    IMPLEMENTS(IDynamicsParameter);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Vol identifier");
    FIELD(compVolExpiry, "Composite vol dates");
    FIELD(compVolMatExpiry, "Composite vol maturity dates");
    FIELD_MAKE_OPTIONAL(compVolMatExpiry);
    FIELD(compVol, "Composite vols");
    FIELD(spotVolExpiry, "Spot vol dates");
    FIELD_MAKE_OPTIONAL(spotVolExpiry);
    FIELD(spotVol, "Spot Vols");
    FIELD_MAKE_OPTIONAL(spotVol);
    FIELD(smileExpiry, "FX Smile dates");
    FIELD_MAKE_OPTIONAL(smileExpiry);
    FIELD(smileA1, "ATM skew");
    FIELD_MAKE_OPTIONAL(smileA1);
    FIELD(smileA2, "ATM crv");
    FIELD_MAKE_OPTIONAL(smileA2);
    FIELD(smileA3, "maximum variation");
    FIELD_MAKE_OPTIONAL(smileA3);
    FIELD_NO_DESC(today);
    FIELD_MAKE_TRANSIENT(today);
    FIELD_NO_DESC(compVolDate);
    FIELD_MAKE_TRANSIENT(compVolDate);
    FIELD_NO_DESC(compVolMatDate);
    FIELD_MAKE_TRANSIENT(compVolMatDate);
    FIELD_NO_DESC(spotVolDate);
    FIELD_MAKE_TRANSIENT(spotVolDate);
    FIELD_NO_DESC(smileDate);
    FIELD_MAKE_TRANSIENT(smileDate);
    FIELD_NO_DESC(forcedTimeOfDay);
    FIELD_MAKE_TRANSIENT(forcedTimeOfDay);

    ClassSetAcceptMethod(acceptCriticalDateCollector);

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA1",
        new Range(ClosedBoundary(-100.), ClosedBoundary(100.)),
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA2",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)),
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "smileA3",
        new Range(ClosedBoundary(0), ClosedBoundary(10.)),  // engine can be unstable if A3 > 5
        getSmileParamExpiries);
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "compVol",
        new Range(OpenBoundary(0), Infinity(Infinity::Plus)),
        getCompVolExpiries);

	// we should add spot vol calibration to this class
	/*
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "spotVol",
        new Range(ClosedBoundary(0.01), ClosedBoundary(1.)), // min 1%, max 100% 
        getSpotVolExpiries);
    */
}

void SRMFXVol::acceptCriticalDateCollector(const SRMFXVol* vol,
                                             CriticalDateCollector* collector) {
    // SRM wants the smile dates
    collector->addDates(vol->smileDate, vol->getClass());
}

/****************************************
 * Composite volatility - Parallel tweak
 ****************************************/
string SRMFXVol::sensName (CompositeVegaParallelTweak* shift) const {
    return name;
}

bool SRMFXVol::sensShift (CompositeVegaParallelTweak* shift) {
    static const string method = "SRMFXVol::sensShift";

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
string SRMFXVol::sensName (CompositeVegaPointwiseTweak* shift) const {
    return name;
}

bool SRMFXVol::sensShift (CompositeVegaPointwiseTweak* shift) {
    static const string method = "SRMFXVol::sensShift";

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


ExpiryArrayConstSP SRMFXVol::sensExpiries (CompositeVegaPointwiseTweak* shift) const {
    return compVolExpiry;
}


/***********************************
 * Spot volatility - Parallel tweak
 ***********************************/
string SRMFXVol::sensName (SpotVegaParallelTweak* shift) const {
    // Return the object's name if the optional field spotVol is present,
    // or an empty string otherwise
    return (spotVol.size() > 0) ? name : "";
}

bool SRMFXVol::sensShift (SpotVegaParallelTweak* shift) {
    static const string method = "SRMFXVol::sensShift";

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
string SRMFXVol::sensName (SpotVegaPointwiseTweak* shift) const {
    // Return the object's name if the optional field spotVol is present,
    // or an empty string otherwise
    return (spotVol.size() > 0) ? name : "";
}

bool SRMFXVol::sensShift (SpotVegaPointwiseTweak* shift) {
    static const string method = "SRMFXVol::sensShift";

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


ExpiryArrayConstSP SRMFXVol::sensExpiries (SpotVegaPointwiseTweak* shift) const {
    return spotVolExpiry;
}

CClassConstSP const SRMFXVol::VOL_TYPE = CClass::registerClassLoadMethod(
    "SRMFX::Vol", typeid(SRMFXVol), SRMFXVol::load);

/**************************/
bool SRMFXVolLoad(void) {
    return SRMFXVol::VOL_TYPE != 0;
}

DRLIB_END_NAMESPACE