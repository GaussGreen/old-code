#include "edginc/config.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/SRMFXVolSimple.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/FlatFXVol.hpp"
#include "edginc/FXAsset.hpp" // Needed for SRMFX_VolSimple2VolConverterAddin
#include "edginc/MarketDataFetcherSRM.hpp" // Needed for SRMFX_VolSimple2VolConverterAddin
#include "edginc/NonPricingModel.hpp" // Needed for SRMFX_VolSimple2VolConverterAddin

DRLIB_BEGIN_NAMESPACE

/** Turn expiries into date times */
void SRMFXVolSimple::cacheDates(){
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
    try{
        DateTime::ensureIncreasing(compVolDate, "compVolDate", false);
        DateTime::ensureIncreasing(compVolMatDate, "compVolMatDate", false);
        DateTime::ensureIncreasing(spotVolDate, "spotVolDate", false);
        DateTime::ensureIncreasing(smileDate, "smileDate", false);
    } catch (exception& e){
        throw ModelException(e, "SRMFXVolSimple::cacheDates",
                             "For SRMFXVolSimple with name "+name);
    }
}

/** calculates the trading time between two dates */
double SRMFXVolSimple::calcTradingTime(
    const DateTime &date1,
    const DateTime &date2) const
{
    // no time metric supported at the moment
    return date1.yearFrac(date2);
}

/** retrieve time measure for the vol */
TimeMetricConstSP SRMFXVolSimple::GetTimeMetric() const
{
    HolidaySP noHols(Holiday::noHolidays());
    return TimeMetricConstSP(new TimeMetric(1.0, noHols.get()));
}

 void SRMFXVolSimple::validatePop2Object()
 {
    static const string method("SRMFXVolSimple::validatePop2Object");
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
    for (int m = 0; m < smileExpiry->size(); m++){
        if (!((*smileExpiry)[m])){
            throw ModelException(method, "Null smileExpiry on index "+
                                 Format::toString(m)+" on "+name);
        }
    }
    Calibrator::IAdjustable::checkRange(this);
    if (Maths::isNegative(smileA31Y)){
        throw ModelException(method,
                             "The 1-year smileA3 parameter must be >=0");
    }
    if (Maths::isZero(smileA31Y)){ /* flat local vol */
        if (!Maths::isZero(smileA21Y) ||
            !Maths::isZero(smileA11Y)) {
            /* can't fit if no variation allowed */
            throw ModelException(method, "If the 1-year smileA3 parameter is zero, "
                                 "the 1-year smileA1 and smileA2 parameters must be "
                                 "zero too.");
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
}

/** populate from market cache */
void SRMFXVolSimple::getMarket(
    const IModel* model,
    const MarketData* market)
{
    market->GetReferenceDate(today);
    cacheDates();
}

/** Shifts the object using given shift. */
bool SRMFXVolSimple::sensShift(Theta* shift)
{
    const DateTime& newDate = shift->rollDate(today);
    if (newDate != today){
        cacheDates();
    }
    return false; // nothing else to shift
}

/** I guess we can do something with the spot (=>fwd ) vols? */
CVolProcessed* SRMFXVolSimple::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      asset) const
{
    if (VolRequestRaw::TYPE->isInstance(volRequest) ||
        VolRequestTime::TYPE->isInstance(volRequest)){
        int nbSmileDates = smileDate.size();
        DoubleArray smileA1(nbSmileDates);
        DoubleArray smileA2(nbSmileDates);
        DoubleArray smileA3(nbSmileDates);
        for (int iSmileDate = 0; iSmileDate < nbSmileDates; ++iSmileDate){
            double yearfrac = calcTradingTime(today,
                                              smileDate[iSmileDate]);
            smileA1[iSmileDate] = smileA11Y
                * exp(- smileA1ShortTermPower * log(yearfrac)
                      - smileA1LongTermExpo * (yearfrac - 1.0));
            smileA2[iSmileDate] = smileA21Y
                * exp(- smileA2ShortTermPower * log(yearfrac)
                      - smileA2LongTermExpo * (yearfrac - 1.0));
            smileA3[iSmileDate] = smileA31Y
                * exp(- smileA3ShortTermPower * log(yearfrac)
                      - smileA3LongTermExpo * (yearfrac - 1.0));
        }
        return new SRMFXVol(name,
                       compVolExpiry,
                       compVolMatExpiry,
                       compVol,
                       spotVolExpiry,
                       spotVol,
                       smileExpiry,
                       smileA1,
                       smileA2,
                       smileA3,
                       today,
                       compVolDate,
                       compVolMatDate,
                       spotVolDate,
                       smileDate);
    }
    throw ModelException("SRMFXVolSimple::getProcessedVol",
                         "Request of type "+
                         volRequest->getClass()->getName()+
                         " not supported");
}


/** Returns name of vol */
string SRMFXVolSimple::getName() const
{
    return name;
}

/** I guess we can do something with the spot (=>fwd ) vols? */
CVolProcessed* SRMFXVolSimple::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      eqAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const
{
    throw ModelException("SRMFXVolSimple::getProcessedVol",
                         "Currency struck processed vol not supported");
}

SRMFXVolSimple::SRMFXVolSimple() :
    FXVolBase(VOL_SIMPLE_TYPE)
{
}

IObject* SRMFXVolSimple::defaultConstructor()
{
    return new SRMFXVolSimple();
}

void SRMFXVolSimple::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(SRMFXVolSimple, clazz); // FIX - used to be VolSimple.  Does this break public API?
    SUPERCLASS(FXVolBase);
    IMPLEMENTS(IVolProcessed);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(Theta::IShift);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "VolSimple identifier");
    FIELD(compVolExpiry, "Composite vol dates");
    FIELD(compVolMatExpiry, "Composite vol maturity dates");
    FIELD_MAKE_OPTIONAL(compVolMatExpiry);
    FIELD(compVol, "Composite vols");
    FIELD(spotVolExpiry, "Spot vol dates");
    FIELD_MAKE_OPTIONAL(spotVolExpiry);
    FIELD(spotVol, "Spot Vols");
    FIELD_MAKE_OPTIONAL(spotVol);
    FIELD(smileExpiry, "FX Smile dates");
    FIELD(smileA11Y, "1-year ATM skew");
    FIELD(smileA1ShortTermPower, "short term skew propagation power");
    FIELD(smileA1LongTermExpo, "long term skew propagation exponent");
    FIELD(smileA21Y, "1-year ATM crv");
    FIELD(smileA2ShortTermPower, "short term crv propagation power");
    FIELD(smileA2LongTermExpo, "long term crv propagation exponent");
    FIELD(smileA31Y, "1-year maximum variation");
    FIELD(smileA3ShortTermPower, "short term maximum variation propagation power");
    FIELD(smileA3LongTermExpo, "long term maximum variation propagation exponent");
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

    ClassSetAcceptMethod(acceptCriticalDateCollector);

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "smileA11Y",
        new Range(Infinity(Infinity::Minus), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA1ShortTermPower",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA1LongTermExpo",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA21Y",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA2ShortTermPower",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA2LongTermExpo",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA31Y",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA3ShortTermPower",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "smileA3LongTermExpo",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "compVol",
        new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
}

void SRMFXVolSimple::acceptCriticalDateCollector(
    const SRMFXVolSimple*      vol,
    CriticalDateCollector* collector)
{
    // SRM wants the smile dates
    collector->addDates(vol->smileDate, vol->getClass());
}

CClassConstSP const SRMFXVolSimple::VOL_SIMPLE_TYPE = CClass::registerClassLoadMethod(
    "SRMFX::VolSimple", typeid(SRMFXVolSimple), SRMFXVolSimple::load);

/* Addin SRMFX_VolSimple2VolConverterAddin */
class SRMFX_VolSimple2VolConverterAddin: public CObject{
public:
    static CClassConstSP const TYPE;

private:
    CMarketDataSP       market;
    FXAssetWrapper      fxAsset;

    // transients

    virtual void validatePop2Object(){
        static const string method = "SRMFX_VolSimple2VolConverterAddin::validatePop2Object";
        try{
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    IObjectSP convert(){
        static const string routine = "SRMFX_VolSimple2VolConverterAddin::convert";
        try {
            // populate asset
            StringArray fxVolType(1, "SRMFX::VolSimple");
            MarketDataFetcherSP mdf(
                new MarketDataFetcherSRM("IRCalib::Smile2Q",    // not used
                                         "IRCalib::Model1FL",   // not used
                                         false,
                                         false,                 // useIRVolPair?
                                         "SRMEQ::Vol",          // not used
                                         fxVolType,
                                         "CRCalib::Smile",      // not used
                                         "FlatCDSSpotVol"));    // not used
            NonPricingModel model(mdf);
            fxAsset.getData(&model, market);
            // process VolSimple into Vol
            VolRequestRaw request;
            CVolProcessedSP processedVol(fxAsset->getProcessedVol(&request));
            return SRMFXVolSP(dynamic_cast<SRMFXVol*>(processedVol.get()));
        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    static IObjectSP convert2(SRMFX_VolSimple2VolConverterAddin* params) {
        return params->convert();
    }

    /** for reflection */
    SRMFX_VolSimple2VolConverterAddin():
    CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SRMFX_VolSimple2VolConverterAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(market, "market");
        FIELD(fxAsset, "Asset wrapper");

        Addin::registerClassObjectMethod("SRMFX_VOLSIMPLE2VOL_CONVERT",
                                         Addin::MARKET,
                                         "Converts SRMFX::VolSimple into SRMFX::Vol",
                                         TYPE,
                                         false,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)convert2);
    }

    static IObject* defaultCtor(){
        return new SRMFX_VolSimple2VolConverterAddin();
    }
};

CClassConstSP const SRMFX_VolSimple2VolConverterAddin::TYPE = CClass::registerClassLoadMethod(
    "SRMFX_VolSimple2VolConverterAddin", typeid(SRMFX_VolSimple2VolConverterAddin), load);


DRLIB_END_NAMESPACE
