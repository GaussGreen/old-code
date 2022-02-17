//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GenericNFBase.cpp
//
//   Description : Useful base class for generic n factor equity products
//                 which have a reference level and past values. (It handles
//                 theta for you)
//
//   Author      : Mark A Robson
//
//   Date        : 26 October 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/AsMultiFactors.hpp"

DRLIB_BEGIN_NAMESPACE

GenericNFBase::~GenericNFBase(){}

/* for reflection */
GenericNFBase::GenericNFBase(CClassConstSP clazz): GenericNFactor(clazz) {
    // empty
}

GenericNFBase::GenericNFBase(CClassConstSP           clazz,
                             const DateTime&         valueDate,
                             double                  notional,         
                             InstrumentSettlementSP  instSettle,
                             IMultiMarketFactorsSP   assets, 
                             const YieldCurveWrapper& discount,
                             IRefLevelSP             refLevel,
                             IPastValuesSP           pastValues) 
    : GenericNFactor(clazz, valueDate, notional, instSettle, assets, discount),
    refLevel(refLevel), pastValues(pastValues)  
{
    validatePop2Object();
}

void GenericNFBase::GetMarket(const IModel* model, const CMarketDataSP market) {
    GenericNFactor::GetMarket(model, market);
    refLevel->GetMarket(model,market);
}

void GenericNFBase::validatePop2Object() {
    /* Since we need to be able to make refLevel and pastValues optional in FR::Instrument, 
     * but nowhere else, they are made optional here in GenericNF base, but are 
     * still checked for.
     * All derived classes which override validatePop2Object() must call this method
     * before doing anything else (except where these fields should be optional, like FR::Instrument).
     */
    if (!refLevel || !pastValues) {
        throw ModelException("refLevel and pastValues must be provided");
    }
}

/** Validate instrument having aquired market data */

void GenericNFBase::Validate(){
    GenericNFactor::Validate();

    if (!isdaDateAdjust.get()) {
        isdaDateAdjust = SamplingConventionSP(new UnadjustedConvention());
    }
    // Used to pass "refLevel->getSimStartDate(valueDate)" as start date
    // but this meant fwd starting XCBs were forbidden unnecessarily.
    // So long as the asset is well-defined from the first averaging-in
    // date we should be fine.
    assets->crossValidate(refLevel->getAllDates().front(),
                          valueDate,
                          discount.get(),
                          this);
    obsMap = ObservationMapSP(new ObservationMap());
    validatePastValues();
}

// this method is called after the market data is retrieved.
// basically it fills in the missing pieces of the past values object
// so that we can guarantee that all dates we will ever need are present
// and all dates in the past have values (if not an Exception is thrown)
void GenericNFBase::validatePastValues() {
    // get PastValues to sort itself out. Basically we give it the set of
    // required dates and the assets (which should contain historic samples)
    // we know this is going to be a MultiAsset surely?
    const IAsMultiFactors* asMulti = dynamic_cast<const IAsMultiFactors*>(assets.get());
    if (!asMulti){
        throw ModelException("GenericNFBase::validatePastValues",
                                "IMultiFactors view not supported by"
                                " object of type "+
                                assets.get()->getClass()->getName());
    }
    IMultiFactorsSP mAsset(asMulti->asMultiFactors());

    // note some instruments have sampling dates per asset
    // and some have one set of sampling dates which all assets share
    // we don't currently handle the isda asjustment for the first case
    if (assetsShareSamplingDates()) {
        // get all required dates from instrument
        DateTimeArray requiredDates = samplingDates();
        // currently sources might be in the PastValues object
        ObservationSourceArraySP obs(pastValues->getSources(mAsset.get()));
        ObservationTypeArraySP obsTypes(pastValues->getObsTypes(mAsset.get()));
        bool hasSamplingConvention(pastValues->hasSamplingConvention());
        // initialise the data relating sample dates and obs dates
        obsMap->initialise(requiredDates, refLevel->getAllDates(), mAsset.get(),
                           *obs, *obsTypes, isdaDateAdjust.get(),
                           hasSamplingConvention, valueDate);
        DateTimeArrayArray reqDates = obsMap->getAllModellingDates();
        // validate past values object
        pastValues->validate(valueDate, reqDates, mAsset.get(), isdaDateAdjust.get());
    } else {
        smartConstPtr<DateTimeArrayArray> reqDates(getSamplingDatesPerAsset());
        // validate past values object
        pastValues->validate(valueDate, *reqDates, mAsset.get(), isdaDateAdjust.get());
    }
}

// implementation of PastSamplesEvent::IEventHandler interface
void GenericNFBase::getEvents(const PastSamplesEvent* samples, IModel* model,
                              const DateTime& eventDate, EventResults* events) const {
    static const string method("GenericNFBase::getEvents");

    try {
        // currently we are not building the obsMap in the case of the SuperRainbow
        // THIS NEEDS FIXING AT SOME STAGE
        if (assetsShareSamplingDates()) {
            // we know this is going to be a MultiAsset surely?
            const IAsMultiFactors* asMulti = dynamic_cast<const IAsMultiFactors*>(assets.get());
            if (!asMulti){
                throw ModelException("GenericNFBase::getEvents",
                                        "IMultiFactors view not supported by"
                                        " object of type "+
                                        assets.get()->getClass()->getName());
            }
            IMultiFactorsSP mAsset(asMulti->asMultiFactors());
            DateTimeArrayArray reqDates = obsMap->getAllModellingDates();
            pastValues->pastSamplesEvents(eventDate, reqDates, mAsset.get(),
                                        isdaDateAdjust.get(), events);
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

//// roll through time (setting historic values)
bool GenericNFBase::sensShift(Theta* theta){
    // use valueDate before it changes
    pastValues->roll(theta->getUtil(valueDate),
                     assets.get()); // then roll our past values
    GenericNFactor::sensShift(theta); // and then call parent's method
    return true; // continue to tweak components which implement Theta
}

/** indicates whether VEGA_MATRIX is sensible for this instrument */
bool GenericNFBase::avoidVegaMatrix(const IModel* model){
    // this does raise the question of whether we should have method(s)
    // on the model rather than the instrument
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return !mc.vegaMatrixSupported(this);
    }
    return true;
}

/** returns all strikes on the vol surface to which
    this instrument is sensitive */
DoubleArraySP GenericNFBase::getSensitiveStrikes(OutputNameConstSP outputName,
                                                 const IModel*      model){
    // this does raise the question of whether we should have method(s)
    // on the model rather than the instrument
    if (MonteCarlo::TYPE->isInstance(model)){
        const MonteCarlo& mc = dynamic_cast<const MonteCarlo&>(*model);
        return mc.getSensitiveStrikes(this, outputName);
    }
    throw ModelException("GenericNFBase::getSensitiveStrikes",
                         "Sensitive Strikes not supported");
}

/** when to stop tweaking - defer to product */
DateTime GenericNFBase::endDate(const Sensitivity* sensControl) const{
    return productEndDate(this, sensControl);
}

// most products have a single set of sampling dates so this defaults to true
// some products e.g. SuperRainbow have a set per asset
bool GenericNFBase::assetsShareSamplingDates() const {
    return true;
}

const DateTimeArrayArray* GenericNFBase::getSamplingDatesPerAsset() const {
    throw ModelException("GenericNFBase::getSamplingDatesPerAsset",
                "Internal error, shouldn't call this function except "
                "on specific implementation");
}

/** Invoked when Class is 'loaded' */
void GenericNFBase::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(GenericNFBase, clazz);
    SUPERCLASS(GenericNFactor);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(LastSensDate);
    IMPLEMENTS(ISensitiveStrikes);
    IMPLEMENTS(PastSamplesEvent::IEventHandler);
    FIELD(refLevel, "How to compute reference levels");
    FIELD_MAKE_OPTIONAL(refLevel);
    FIELD(pastValues, "Past Values per Asset");
    FIELD_MAKE_OPTIONAL(pastValues);
    FIELD(isdaDateAdjust, "ISDA holiday adjustment details");
    FIELD_MAKE_OPTIONAL(isdaDateAdjust);
    FIELD(obsMap, "Handles map between ");
    FIELD_MAKE_TRANSIENT(obsMap);
}

CClassConstSP const GenericNFBase::TYPE = CClass::registerClassLoadMethod(
    "GenericNFBase", typeid(GenericNFBase), load);

DRLIB_END_NAMESPACE

