//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IlliquidStock.cpp
//
//   Description : Equity asset for those Pyramid stocks with no vol
//
//   Author      : Andrew J Swain
//
//   Date        : 12 September 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/IlliquidStock.hpp"
#include "edginc/DividendCollector.hpp"


DRLIB_BEGIN_NAMESPACE
/** Populates the object with the market data that this object
    needs.  This method is invoked as the getMarket chains down
    from the instrument to the specific instance of market
    data. The default implementation provided by MarketObject is
    to do nothing. Market data objects that require other pieces
    of market data (eg an XCB requires assets) need to override
    this method. Here we pull out the vol from the market data */
void IlliquidStock::getMarket(const IModel* model, const MarketData* market){
    try {
        // then ask the data (that we already have) to fill in its missing bits
        equity->getMarket(model, market);
    } 
    catch (exception& e) {
        throw ModelException(e, "IlliquidStock::getMarket");
    }
}
    
/** returns the spot price */
double IlliquidStock::getSpot() const{
    return equity->spot();
}

/** Returns fair value of stock price */
double IlliquidStock::fairValue() const
{
    return equity->fairValue();
}

/** returns the asset name */
string IlliquidStock::getName() const{
    return equity->getName();
}

// the IMarketObservable interface for retrieving a single sample
double IlliquidStock::pastValue(const DateTime&             sampleDate,
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
bool IlliquidStock::observationDate(const DateTime&           sampleDate,
                                    const ObservationSource*  source,
                                    const SamplingConvention* sampleRule,
                                    DateTime*                 obsDate) const {
    return equity->observationDate(sampleDate, source, sampleRule, obsDate);
}

// the IMarketObservable interface for retrieving past samples events
double IlliquidStock::addPastSampleEvent(const DateTime&    sampleDate,
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
bool IlliquidStock::isHoliday(const DateTime&            sampleDate,
                              const ObservationSource*   source) const {
    return equity->isHoliday(sampleDate, source);
}

/** this obviously won't work by definition */
CVolProcessed * IlliquidStock::getProcessedVol(
    const CVolRequest* volRequest) const {
    throw ModelException("IlliquidStock::getProcessedVol",
                         "illiquid stock (" + 
                         equity->getName() + 
                         ") doesn't have a vol");
}

/** Calculates the expected spot price of the asset at the given date */
double IlliquidStock::fwdValue(const DateTime& date) const{
    return equity->fwdValue(date);
}

/** Calculates the expected spot price of the asset at each of the
    given dates */
void IlliquidStock::fwdValue(const DateTimeArray& dates,
                            CDoubleArray&        result) const {
    equity->fwdValue(dates, result);
}

/** Calculates the expected spot prices of the asset at the given dates
    respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
void IlliquidStock::fwdValue(const DateTimeArray&     dateList,
                             const FwdValueAlgorithm& algo,
                             CDoubleArray&            result) const{
    equity->fwdValue(dateList, algo, result);
}    

/** Returns the name (not the ISO code) of the asset ccy */
string IlliquidStock::getYCName() const {
    return equity->getYCName();
}


/** Calculate the settlement date associated with a given trade date */
DateTime IlliquidStock::settleDate(const DateTime& tradeDate) const{
    return equity->settles(tradeDate);
}

void IlliquidStock::acceptValueDateCollector(
    const IlliquidStock* asset, 
    CValueDateCollector* collector)
{
    asset->equity->accept(collector);
}


void IlliquidStock::acceptImntCcy(
    const IlliquidStock* asset,
    AssetCcyCollector*   collector)
{
  collector->currencyValidate(asset->equity->getYCIsoCode(),
                              asset->getName());

}

void IlliquidStock::acceptHoliday(
    const IlliquidStock*    asset,
    HolidayCollector*       collector)
{
    collector->setHoliday(asset->equity->getMarketHolidays());
}

void IlliquidStock::acceptCriticalDateCollector(
    const IlliquidStock*   asset, 
    CriticalDateCollector* collector)
{
    asset->equity->accept((ICollector*)collector);
}

/** return a pdf calculator */
PDFCalculator* IlliquidStock::pdfCalculator(const PDFRequest* request) const {
    throw ModelException("IlliquidStock::pdfCalculator",
                         "illiquid stock (" + 
                         equity->getName() + 
                         ") doesn't have a vol");
}

/** Allows asset to be currency struck. This method fails as we've
    got an illiquid stock */
CVolProcessed* IlliquidStock::getProcessedVol(
    const CVolRequest* volRequest,
    const CAsset*      struckAsset,
    const FXAsset*     fxAsset,
    const Correlation* eqFXCorr) const{
    return getProcessedVol(volRequest);
}


/* for reflection */
IlliquidStock::IlliquidStock(): CAsset(TYPE){}
    
class IlliquidStockHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IlliquidStock, clazz);
        SUPERCLASS(Asset);
        IMPLEMENTS(IAssetFairValue);
        IMPLEMENTS(Asset::IStruckable);
        EMPTY_SHELL_METHOD(defaultIlliquidStock);
        FIELD(equity, "Stock");
        ClassSetAcceptMethod(IlliquidStock::acceptValueDateCollector);
        ClassSetAcceptMethod(IlliquidStock::acceptImntCcy);
        ClassSetAcceptMethod(IlliquidStock::acceptHoliday);
        ClassSetAcceptMethod(acceptDividendCollector);
    }

    /** pass on our dividends */
    static void acceptDividendCollector(const IlliquidStock* asset,
                                        DividendCollector*   collector){
        collector->addDivs(asset->equity->getDivList(), 0 /* non struck*/);
    }
    static IObject* defaultIlliquidStock(){
        return new IlliquidStock();
    }
};

CClassConstSP const IlliquidStock::TYPE = CClass::registerClassLoadMethod(
    "IlliquidStock", typeid(IlliquidStock), IlliquidStockHelper::load);

DRLIB_END_NAMESPACE
