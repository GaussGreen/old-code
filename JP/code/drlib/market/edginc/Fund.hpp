//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Fund.hpp
//
//   Description : Fund of assets
//
//   Author      : Andrew J Swain
//
//   Date        : 24 August 2001
//
//
//----------------------------------------------------------------------------

#ifndef FUND_HPP
#define FUND_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/StartDateCollector.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/FutureExpiryCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/Theta.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/Spot.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/Correlation.hpp"
#include "edginc/RateParallel.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/IEqVolNamePair.hpp"
#include "edginc/ProxyVol.hpp"
#include "SimpleEquity.hpp"

DRLIB_BEGIN_NAMESPACE

class Fund;
typedef smartPtr<Fund> FundSP;
typedef smartConstPtr<Fund> FundConstSP;

/** Fund of assets */
class MARKET_DLL Fund: public CAsset,
            virtual public CAsset::IStruckable,
            virtual public ObjectIteration::IOverride,
            virtual public IAssetFairValue {
public:
    static CClassConstSP const TYPE;

    /** Validation */
    virtual void validatePop2Object();

    /** Returns fair value of fund price */
    virtual double fairValue() const;

    /** Pull out the component assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;
    
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

    /** Returns the name (not the ISO code) of the fund ccy */
    virtual string getYCName() const;

   /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;

    /** returns sensitive strikes for a given vol request */
    virtual void getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                CResults*       results,
                                const DateTime& maturityDate) const;

    /** used to turn off shifting through anything but the SimpleEquity */
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        asset together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. Note that the struckAsset is indeed the struckAsset cf
        'this' which is the non struck asset */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      struckAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** Returns the Funds dividends */
    DividendListConstSP getDivList() const;

private:
    friend class FundHelper;
    friend class FundVol;
    friend class XLGetCorrelationMatrixAddin;

    Fund();
    Fund(const Fund& rhs);
    Fund& operator=(const Fund& rhs);

    static void acceptCollector(const Fund* asset, 
                                StartDateCollector* collector);

    static void acceptNameCollector(const Fund* asset, 
                                    AssetNameCollector* collector);

    static void acceptFutureCollector(const Fund* asset, 
                                      FutureExpiryCollector* collector);

    static void acceptDeltaShift(const Fund* asset,
                                 ShiftSizeCollector* collector);

    static void acceptImntCcy(const Fund* asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const Fund* asset,
                              HolidayCollector* collector);

    static void acceptCriticalDateCollector(const Fund* asset,
                                            CriticalDateCollector* collector);

    //  fields
    EquitySP           fund;        // spot, divs, borrow
    CAssetWrapperArray components;  // array of component assets
    CStringArray       ccyTreatments; // None, Struck or Prot for each asset
    DoubleArray        weights;     // percentage of each asset 
    HolidayWrapper     marketHols;  // market holidays
    DateTime           baseDate;
    DoubleMatrix       correlations; /* between assets - this should be 
                                        initialised in the getMarketData 
                                        method using Correlations */
    TimeMetricSP       timeMetric;  // for composite vol
    ExpiryArraySP      spreadDates; // ditto
    DoubleArray        volSpread;   // ditto
    DoubleArray        volSkew;     // ditto
    double             minVol;      // ditto

    // use term structure of correlation when computing composite vol
    bool               useCorrelationTerm;    

    // transient fields
    ProxyVolSP         proxyVol;
    SimpleEquitySP     realAsset;
};

DRLIB_END_NAMESPACE
#endif
