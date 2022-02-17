//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquityBase.hpp
//
//   Description : SimpleEquity class -> EquityBase class. contains basically everything that was in SimpleEquity
//                 but the registration of vol and equity is left to SimpleEquity which now inherits from this class.
//
//   Date        : 30 Nov 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/VolRequestLN.hpp"
#include "edginc/AllStrikes.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/VolTypeSensitiveStrikes.hpp"
#include "edginc/VolDeltaShiftSize.hpp"

DRLIB_BEGIN_NAMESPACE

/** Populates the object with the market data that this object
    needs.  This method is invoked as the getMarket chains down
    from the instrument to the specific instance of market
    data. The default implementation provided by MarketObject is
    to do nothing. Market data objects that require other pieces
    of market data (eg an XCB requires assets) need to override
    this method. Here we pull out the vol from the market data */
void EquityBase::getMarket(const IModel* model, const MarketData* market){
    try{
        vol.getData(model, market);
        // then ask the data (that we already have) to fill in its missing bits
        equity->getMarket(model, market);
    } catch (exception& e){
        throw ModelException(e, "EquityBase::getMarket");
    }
}
    
/** returns the spot price */
double EquityBase::getSpot() const{
    return equity->spot();
}

/** Returns fair value of stock price */
double EquityBase::fairValue() const
{
    return equity->fairValue();
}

/** Returns the name (not the ISO code) of the yield curve used to
grow the stock */
string EquityBase::getYCName() const
{
    return equity->getYCName();
}

/** Returns the ISO code (not the name) of the yield curve used to
grow the stock */
string EquityBase::getYCIsoCode() const
{
    return equity->getYCIsoCode();
}

DateTime EquityBase::getValueDate() const {
    return equity->getValueDate();
}

/** returns the asset name */
string EquityBase::getName() const{
    return equity->getName();
}

// the IMarketObservable interface for retrieving a single sample
double EquityBase::pastValue(const DateTime&            sampleDate,
                            const ObservationType*      obsType,
                            const ObservationSource*    source,
                            const FixingType*           fixType,
                            const IObservationOverride* overrides,
                            const SamplingConvention*   sampleRule) const{
    return equity->pastValue(sampleDate, obsType, source,
                             fixType, overrides, sampleRule);
}

// IMarketObservable - retrieve a single observation date
// Returns false if obs is to be omitted
bool EquityBase::observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const {
    return equity->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double EquityBase::addPastSampleEvent(const DateTime&       sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const {
    return equity->addPastSampleEvent(sampleDate, obsType, source, fixType,
                                      overrides, sampleRule, collector);
}

// the IMarketObservable interface for 
// is the given date a holiday for the relevant source
bool EquityBase::isHoliday(const DateTime&            sampleDate,
                           const ObservationSource*   source) const {
    return equity->isHoliday(sampleDate, source);
}

/** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
CVolProcessed * EquityBase::getProcessedVol(
    const CVolRequest* volRequest) const{
    return vol->getProcessedVol(volRequest, this);
}

/** Calculates the expected spot price of the asset at the given date */
double EquityBase::fwdValue(const DateTime& date) const{
    return equity->fwdValue(date);
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void EquityBase::fwdValue(const DateTimeArray& dates,
                            CDoubleArray&        result) const{
    equity->fwdValue(dates, result);
}

/** Calculates the expected spot prices of the asset at the given dates
    respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
void EquityBase::fwdValue(const DateTimeArray&     dates,
                            const FwdValueAlgorithm& algo,
                            CDoubleArray&            result) const{
    equity->fwdValue(dates, algo, result);
}

/** Calculates the expected spot price of the asset at the given date if
    the spot price had the given value spot on spotDate */
double EquityBase::fwdFwd(const DateTime& spotDate,
                            double          spot, 
                            const DateTime& fwdDate) const {
     DateTimeArraySP dateArray(new DateTimeArray(1));
     CDoubleArraySP  result(new CDoubleArray(1));
 
     (*dateArray)[0] = fwdDate;
 
     equity->fwdFwd(spotDate,
                    spotDate,
                    spot,
                    *(dateArray.get()), 
                    *(result.get()));
 
     return (*result)[0];
}

/** Array version of fwdFwd */
void EquityBase::fwdFwd(const DateTime&      spotDate,
                        double               spot, 
                        const DateTimeArray& fwdDates,
                        DoubleArray&         results) const
{
     equity->fwdFwd(spotDate,
                    spotDate,
                    spot,
                    fwdDates, 
                    results);
}

/** Calculate the settlement date associated with a given trade date */
DateTime EquityBase::settleDate(const DateTime& tradeDate) const{
    return equity->settles(tradeDate);
}

/** Create a struck equity from this simple equity */
CAssetSP EquityBase::createStruckEquity(
    const MarketData* market,
    const string&     payOutYCName){
    try{
        // need to get hold of fx asset's name
        string growYCName = equity->getYCName();
        string fxAssetName = market->getFXName(growYCName, payOutYCName);
        return CAssetSP(new StruckEquity("", // use default
                                equity.get(),
                                vol->getName(),
                                fxAssetName));
    } catch (exception& e){
        throw ModelException(e, "EquityBase::createStruckEquity", "Failed to"
                             " create struck ("+payOutYCName+")"
                             " asset for "+equity->getName());
    }
}

/** Create a protected equity from this simple equity */
CAssetSP EquityBase::createProtEquity(
    const MarketData* market,
    const string&     payOutYCName){
    try{
        // need to get hold of fx vol's name (so need fx asset's name)
        string growYCName = equity->getYCName();
        string fxAssetName = market->getFXName(growYCName, payOutYCName);
        // then get hold of fx asset and ask for fx vol name
        IObjectSP fxAsset(market->GetData(fxAssetName, FXAsset::TYPE));
        FXAsset& fx = dynamic_cast<FXAsset&>(*fxAsset);
        
        return CAssetSP(new ProtEquity("", // use default
                              equity.get(),
                              vol->getName(),
                              fx.getVolName(),
                              payOutYCName,
                              fx.useCoherentGreeks()));
    } catch (exception& e){
        throw ModelException(e, "EquityBase::createProtEquity", "Failed to"
                             " create protected ("+payOutYCName+")"
                             " asset for "+equity->getName());
    }
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string EquityBase::sensName(VegaSkewParallel* shift) const{
    return vol->getName();
}
    
bool EquityBase::sensShift(VegaSkewParallel* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot());
    return true; // continue on and do the vol
}

/** Returns the name of the vol - used to determine whether to tweak
    the object */
string EquityBase::sensName(VegaSkewPointwise* shift) const{
    return vol->getName();
}
    
ExpiryArrayConstSP EquityBase::sensExpiries(VegaSkewPointwise* shift) const{
    return ExpiryArrayConstSP(); // return null SP - the vol will do this
}

bool EquityBase::sensShift(VegaSkewPointwise* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot());
    return true; // continue on and do the vol
}

/** Returns the name of the EQUITY - used to determine whether to tweak
    the object */
string EquityBase::sensName(DeltaSurface* shift) const{
    return equity->getName();
}
    
bool EquityBase::sensShift(DeltaSurface* shift){
    // store the underlying spot price
    shift->setSpot(equity->spot(), equity->getName(), vol->getName());
    return true; // continue on and do the vol
}

// implementation of VolRelativeShift::IShift interface
string EquityBase::sensName(VolRelativeShift* shift) const{
    return VolRelativeShift::IShift::TYPE->isInstance(vol.get())
        ? vol->getName(): "";
}
    
bool EquityBase::sensShift(VolRelativeShift* shift){
    // store the underlying spot price
    shift->setSpot(equity->getValueDate(), this);
    return true; // continue on and do the vol
}

/** Returns the name of the stock/asset - used to determine
    whether to shift the object */
string EquityBase::sensName(SpotLevelProbability* shift) const {
    return equity->getName();
}

/** Shifts the object using given shift. Return true to make the
    infrastructure keep tweaking the components within the object
    which implements this interface */
bool EquityBase::sensShift(SpotLevelProbability* shift) {
    double spot = shift->spotLevel(equity->getValueDate(), this);
    SpotLevel level(spot);
    equity->sensShift(&level);  // set equity spot to our new level
    return false;  // all done;
}

/** returns sensitive strikes for a given vol request */
void EquityBase::getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const
{
    // to do: alter getSensitiveStrikes signature
    const CVolRequestLN& request = dynamic_cast<const CVolRequestLN&>(
        *(IObject*)(volRequest));
    if (outputName->equals(vol->getName()) )
    {
        if (sensStrikeDesc.forwardOnly == false) {
            request.getSensitiveStrike(equity->spot(),
                                   sensitiveStrikes);
        }
    } else if (IVolTypeSensitiveStrikes::TYPE->isInstance(vol.get())) {
        // may have to defer to the vol
        const IVolTypeSensitiveStrikes* myVol = 
                dynamic_cast<const IVolTypeSensitiveStrikes*>(vol.get());
        myVol->getSensitiveStrikes(this, volRequest, outputName, 
                                    sensStrikeDesc, sensitiveStrikes);
    }
}

void EquityBase::acceptValueDateCollector(const EquityBase*    asset, 
                                          CValueDateCollector* collector)
{
    asset->equity->accept(collector);
}

void EquityBase::acceptDeltaShift(const EquityBase*   asset, 
                                  ShiftSizeCollector* collector)
{
    const IObject* obj = asset->vol.get();
    const IVolDeltaShiftSize* volDelShift = 
                        dynamic_cast<const IVolDeltaShiftSize*>(obj);
    if (volDelShift) {
        volDelShift->adjustDeltaShiftSize(collector,
                                    asset->getName(),
                                    asset->getSpot());
    }
}

void EquityBase::acceptImntCcy(const EquityBase*  asset,
                               AssetCcyCollector* collector)
{
  collector->currencyValidate(asset->getYCIsoCode(),
                              asset->getName());

}

void EquityBase::acceptHoliday(const EquityBase* asset,
                               HolidayCollector* collector)
{
    collector->setHoliday(asset->equity->getMarketHolidays());
}

/** dollar div and yield div converted to dollar, yield converted to amount 
    but the type name is kept as PERCENET for dollar div treatment.
    better to add a new type AMOUNT_FROM_YIELD in dividend */
void EquityBase::acceptCriticalDateCollector(
    const EquityBase*      equity, 
    CriticalDateCollector* collector)
{
    equity->equity->accept((ICollector*)collector);
    equity->vol->accept((ICollector*)collector);
}

DividendListConstSP  EquityBase::getDivList() const {
    return equity->getDivList();
}

/** returns the name of the vol base object */
string EquityBase::getVolName() const {
    return vol.getName();
}

/** given a current spot level, get the next strike on the vol surface where 
    the slope is non-differentiable */
double EquityBase::getNextStrike(const double& strike,
                                   bool          isUp,
                                   bool&         offSurface) const
{
    double volStrike;
    const IObject* obj = vol.get();
    const INextStrike* nextStrike = dynamic_cast<const INextStrike*>(obj);
    if ( nextStrike )
    {
        volStrike = nextStrike->getNextStrike(strike,
                                              isUp,
                                              offSurface);
    }
    else
    {
        volStrike  = 0.0;
        offSurface = true;
    }
    return volStrike;
}

PDFCalculator* EquityBase::pdfCalculator(const PDFRequest* request) const {
    if (IPDFCalculator::TYPE->isInstance(vol.get())) {
        const IPDFCalculator* pdf = dynamic_cast< const IPDFCalculator*>(vol.get());
        return pdf->getPDFCalculator(request, this);
    }
    throw ModelException("EquityBase::pdfCalculator",
                         "vol of type (" + vol->getClass()->getName() +
                         ") has no pdf calculator");
}

/* for IEqVolNamePair */
bool EquityBase::getNamePairs(string& eqName, string& volName) const {
    eqName = getName();
    volName = getVolName();
    return false;  // no more assets inside here
}

/** constructor  */
EquityBase::EquityBase(const CClassConstSP& clazz,
                       const Equity*        equity,
                       const CVolBase*      vol): 
    CAsset(clazz), equity(copy(equity)), vol(copy(vol)) {}

/* for reflection */
EquityBase::EquityBase(): CAsset(TYPE){}

class EquityBaseHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPrivate();
        REGISTER(EquityBase, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(SpotLevelProbability::Shift);
        IMPLEMENTS(VegaSkewParallel::IShift);
        IMPLEMENTS(VegaSkewPointwise::IShift);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(INextStrike);
        IMPLEMENTS(DeltaSurface::IShift);
        IMPLEMENTS(VolRelativeShift::IShift);
        IMPLEMENTS(IEqVolNamePair);
        EMPTY_SHELL_METHOD(defaultEquityBase);
        ClassSetAcceptMethod(EquityBase::acceptValueDateCollector);
        ClassSetAcceptMethod(EquityBase::acceptDeltaShift);
        ClassSetAcceptMethod(EquityBase::acceptImntCcy);
        ClassSetAcceptMethod(EquityBase::acceptHoliday);
        ClassSetAcceptMethod(EquityBase::acceptCriticalDateCollector);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    static IObject* defaultEquityBase(){
        return new EquityBase();
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const EquityBase*  asset,
                                        DividendCollector* collector){
        collector->addDivs(asset->equity->getDivList(), 0 /* non struck*/);
    }

};

CClassConstSP const EquityBase::TYPE = CClass::registerClassLoadMethod(
    "EquityBase", typeid(EquityBase), EquityBaseHelper::load);

DRLIB_END_NAMESPACE

