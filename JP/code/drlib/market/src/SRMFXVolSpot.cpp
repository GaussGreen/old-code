#include "edginc/config.hpp"
#include "edginc/SRMFXVolSpot.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FlatFXVol.hpp"

DRLIB_BEGIN_NAMESPACE

void SRMFXVolSpot::forceTimeOfDay(void){
    DateTime::setTimeOfDay(spotVolDate, today.getTime());
    DateTime::setTimeOfDay(smileDate, today.getTime());
}

SRMFXVolSpot::SRMFXVolSpot(
                           string name,
                           ExpiryArraySP spotVolExpiry,
                           DoubleArray   spotVol,
                           ExpiryArraySP smileExpiry,
                           DoubleArray   smileA1,
                           DoubleArray   smileA2,
                           DoubleArray   smileA3,
                           DateTime      today,
                           DateTimeArray spotVolDate,
                           DateTimeArray smileDate
):
    name(name),
    SRMFXVolBase(VOL_TYPE),
    spotVolExpiry(spotVolExpiry),
    spotVol(spotVol),
    smileExpiry(smileExpiry),
    smileA1(smileA1),
    smileA2(smileA2),
    smileA3(smileA3),
    today(today),
    spotVolDate(spotVolDate),
    smileDate(smileDate)
{}

/** calculates the trading time between two dates */
double SRMFXVolSpot::calcTradingTime(const DateTime &date1,
                               const DateTime &date2) const{
    // no time metric supported at the moment
    return date1.yearFrac(date2);
}

/** retrieve time measure for the vol */
TimeMetricConstSP SRMFXVolSpot::GetTimeMetric(void) const {
    HolidaySP noHols(Holiday::noHolidays());
    return TimeMetricConstSP(new TimeMetric(1.0, noHols.get()));
}

void SRMFXVolSpot::validatePop2Object(void) {
    static const string method("SRMFXVolSpot::validatePop2Object");
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
void SRMFXVolSpot::getMarket(const IModel* model, const MarketData* market) {
    market->GetReferenceDate(today);
}

/** I guess we can do something with the spot (=>fwd ) vols? */
CVolProcessed* SRMFXVolSpot::getProcessedVol(const CVolRequest* volRequest,
                                       const CAsset*      asset) const{
    if (VolRequestRaw::TYPE->isInstance(volRequest) ||
        VolRequestTime::TYPE->isInstance(volRequest)){
        // it's ours or can just use this
        return const_cast<SRMFXVolSpot*>(this);
    } 
    throw ModelException("SRMFXVolSpot::getProcessedVol",
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported");
}


/** Returns name of vol */
string SRMFXVolSpot::getName(void) const{
    return name;
}


ExpiryArraySP SRMFXVolSpot::getSmileParamExpiries(const IObject* obj){
    const SRMFXVolSpot* srmfxVol = dynamic_cast<const SRMFXVolSpot*>(obj);
    if (obj){
        if (!srmfxVol->smileExpiry){
            // smileExpiry is optional
            return ExpiryArraySP(new ExpiryArray());
        }
        return ExpiryArraySP(copy(srmfxVol->smileExpiry.get()));
    }
    // should never happen
    throw ModelException("SRMFXVolSpot::getSmileParamExpiries", "Internal error: obj is not of the desired type");
}

void SRMFXVolSpot::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SRMFXVolSpot, clazz); // FIX - used to be Vol.  Does this break public API?
    SUPERCLASS(FXVolBase);
    IMPLEMENTS(IVolProcessed);
//  IMPLEMENTS(IDynamicsParameter);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Vol identifier");
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
    FIELD_NO_DESC(spotVolDate);
    FIELD_MAKE_TRANSIENT(spotVolDate);
    FIELD_NO_DESC(smileDate);
    FIELD_MAKE_TRANSIENT(smileDate);

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

	// we should add spot vol calibration to this class
	/*
    Calibrator::IAdjustable::registerBootstrappableField(
        clazz, "spotVol",
        new Range(ClosedBoundary(0.01), ClosedBoundary(1.)), // min 1%, max 100% 
        getSpotVolExpiries);
    */
}

void SRMFXVolSpot::acceptCriticalDateCollector(const SRMFXVolSpot* vol,
                                             CriticalDateCollector* collector) {
    // SRM wants the smile dates
    collector->addDates(vol->smileDate, vol->getClass());
}

CClassConstSP const SRMFXVolSpot::VOL_TYPE = CClass::registerClassLoadMethod(
    "SRMFX::VolSpot", typeid(SRMFXVolSpot), SRMFXVolSpot::load);

/**************************/
bool SRMFXVolSpotLoad(void) {
    return SRMFXVolSpot::VOL_TYPE != 0;
}

DRLIB_END_NAMESPACE