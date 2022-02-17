//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Future.hpp
//
//   Description : A Future as an asset.
//
//   Author      : Mark A Robson
//
//   Date        : 19 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FUTURE_HPP
#define EDG_FUTURE_HPP
#include "edginc/Asset.hpp"
#include "edginc/AssetNameCollector.hpp"
#include "edginc/FutureExpiryCollector.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of CAsset for a single stock with no currency treatment */
class MARKET_DLL Future: public CAsset,
              public Theta::IShift{
public:
    static CClassConstSP const TYPE;
    friend class FutureHelper;

    /** Validation */
    virtual void validatePop2Object();

    /** returns the spot price */
    virtual double getSpot() const;
    
    /** returns the asset name */
    virtual string getName() const;

    /** Pull out the underlying asset and value date */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;
    
    /** Calculates the expected spot price of the asset at the given date.*/
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot price of the asset at each of the
        given dates */
    virtual void fwdValue(const DateTimeArray&     dates,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

    /** Calculates the expected spot price of the asset at each of the
        given dates */
    virtual void fwdValue(const DateTimeArray& dates,
                          CDoubleArray&        result) const;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;


//    /** Calculates the expected spot price of the asset at the given date if
//        the spot price had the given value spot on spotDate */
//    virtual double fwdFwd(const DateTime& spotDate,
//                          double          spot, 
//                          const DateTime& fwdDate) const;

   /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;

    virtual void getSensitiveStrikes(
                    const CVolRequest* volRequest,
                    OutputNameConstSP outputName,
                    const SensitiveStrikeDescriptor& sensStrikeDesc,
                    DoubleArraySP sensitiveStrikes) const;

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                CResults*       results,
                                const DateTime& maturityDate) const;
     
    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** constructor */
    Future(const  string&     name, // optional
           const Asset*       asset,
           const DateTime&    valueDate,
           const DateTime&    maturity);

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

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

private:
    Future();
    Future(const Future& rhs);
    Future& operator=(const Future& rhs);

    static void acceptFutureCollector(const Future* asset, 
                                      FutureExpiryCollector* collector);

    static void acceptValueDateCollector(const Future* asset, 
                                         CValueDateCollector* collector);

    static void acceptImntCcy(const Future* asset,
                              AssetCcyCollector* collector);

    /* utility-- throws exception if date.isGreater(maturity)
       Error message is:
       caller: cannot do calculation for date dd-mmm-yyyy
       after futures expiration dd-mmm-yyyy */
    void throwDateAfterExpiry(const DateTime &date,
                              const string &caller) const;

    ////// fields /////
    string            name;   // optional
    CAssetWrapper     asset;
    DateTime          valueDate;
    DateTime          maturity;
    string            ccyTreatment;
    YieldCurveWrapper yc;
};

typedef smartPtr<Future> FutureSP;
typedef smartConstPtr<Future> FutureConstSP;

DRLIB_END_NAMESPACE
#endif
