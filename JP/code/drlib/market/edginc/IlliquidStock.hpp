//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IlliquidStock.hpp
//
//   Description : Equity asset for those Pyramid stocks with no vol
//
//   Author      : Andrew J Swain
//
//   Date        : 12 September 2001
//
//
//----------------------------------------------------------------------------

#ifndef ILLIQUIDSTOCK_HPP
#define ILLIQUIDSTOCK_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Equity asset for those Pyramid stocks with no vol */
class MARKET_DLL IlliquidStock: public CAsset,
                     public virtual IAssetFairValue,
                     public virtual Asset::IStruckable{
public:
    static CClassConstSP const TYPE;

    /** Populates the object with the market data that this object
        needs.  This method is invoked as the getMarket chains down
        from the instrument to the specific instance of market
        data. The default implementation provided by MarketObject is
        to do nothing. Market data objects that require other pieces
        of market data (eg an XCB requires assets) need to override
        this method. Here we pull out the vol from the market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;

    /** Returns fair value of stock price */
    virtual double fairValue() const;

    /** returns the asset name */
    virtual string getName() const;

    // the IMarketObservable interface for retrieving a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const;

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const;

    // the IMarketObservable interface for 
    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const ObservationSource*   source) const;

    /** this obviously won't work by definition */
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot price of the asset at each of the
        given dates */
    virtual void fwdValue(const DateTimeArray& dates,
                          CDoubleArray&        result) const;

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;

    /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;
     
    /** Allows asset to be currency struck. This method fails as we've
        got an illiquid stock */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      struckAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;
private:
    friend class IlliquidStockHelper;

    IlliquidStock();
    IlliquidStock(const IlliquidStock& rhs);
    IlliquidStock& operator=(const IlliquidStock& rhs);

    static void acceptValueDateCollector(const IlliquidStock* asset, 
                                         CValueDateCollector* collector);

    static void acceptImntCcy(const IlliquidStock* asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const IlliquidStock* asset,
                              HolidayCollector* collector);

    static void acceptCriticalDateCollector(const IlliquidStock* asset,
                                            CriticalDateCollector* collector);

    // fields
    EquitySP equity;
};

typedef smartPtr<IlliquidStock> IlliquidStockSP;
typedef smartConstPtr<IlliquidStock> IlliquidStockConstSP;

// support for wrapper class
typedef MarketWrapper<IlliquidStock> IlliquidStockWrapper;

DRLIB_END_NAMESPACE
#endif
