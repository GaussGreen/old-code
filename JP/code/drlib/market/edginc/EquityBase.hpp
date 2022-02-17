//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EquityBase.hpp
//
//   Description : base equity class. contains basically everything that was in SimpleEquity 
//                 but the registration of vol and equity is left to SimpleEquity which now inherits from this class.
//
//   Date        : 30 Nov 2001
//
//----------------------------------------------------------------------------

#ifndef EQUITY_BASE_HPP
#define EQUITY_BASE_HPP
#include "edginc/Asset.hpp"
#include "edginc/Equity.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/VegaSkewParallel.hpp"
#include "edginc/VegaSkewPointwise.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/ShiftSizeCollector.hpp"
#include "edginc/AssetCcyCollector.hpp"
#include "edginc/HolidayCollector.hpp"
#include "edginc/NextStrike.hpp"
#include "edginc/AssetFairValue.hpp"
#include "edginc/SpotLevelProbability.hpp"
#include "edginc/HaveEquity.hpp"
#include "edginc/DeltaSurface.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/DDeltaDVol.hpp"
#include "edginc/IEqVolNamePair.hpp"

DRLIB_BEGIN_NAMESPACE
class StruckEquity;
class ProtEquity;

/** Implementation of CAsset for a single stock with no currency treatment */
class MARKET_DLL EquityBase: public CAsset,
                  virtual public SpotLevelProbability::Shift,
                  virtual public VegaSkewParallel::IShift,
                  virtual public VegaSkewPointwise::IShift,
                  virtual public IAssetFairValue,
                  virtual public INextStrike,
                  virtual public IHaveEquity,
                  virtual public DeltaSurface::IShift,
                  virtual public VolRelativeShift::IShift,
                  virtual public IEqVolNamePair
{
public:
    static CClassConstSP const TYPE;
    friend class EquityBaseHelper;
    friend class PseudoEquityBase;

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

    /** Returns valueDate */
    virtual DateTime getValueDate() const; 

    /** Returns the name (not the ISO code) of the yield curve used to
        grow the stock */
    string getYCName() const;

    /** Returns the ISO code (not the name) of the yield curve used to
        grow the stock */
    string getYCIsoCode() const;

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

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
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

    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const;

    virtual void fwdFwd(const DateTime&      spotDate,
                        double               spot, 
                        const DateTimeArray& fwdDates,
                        DoubleArray&         results) const;


    /** Calculate the settlement date associated with a given trade date */
    DateTime settleDate(const DateTime& tradeDate) const;
     
    /** Create a struck equity from this simple equity */
    virtual CAssetSP createStruckEquity(const MarketData* market,
                                     const string&     payOutYCName);

    /** Create a protected equity from this simple equity */
    virtual CAssetSP createProtEquity(const MarketData* market,
                                 const string&     payOutYCName);

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewParallel* shift) const;
    /** Shifts the object for VegaSkewParallel */
    virtual bool sensShift(VegaSkewParallel* shift);

    /** Returns the name of the vol - used to determine whether to tweak
        the object */
    virtual string sensName(VegaSkewPointwise* shift) const;
    /** Return the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    /** Shifts the object for VegaSkewPointwise */
    virtual bool sensShift(VegaSkewPointwise* shift);

    /// implementation of DeltaSurface::IShift interface
    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

    // implementation of VolRelativeShift::IShift interface
    virtual string sensName(VolRelativeShift* shift) const;
    virtual bool sensShift(VolRelativeShift* shift);

    /** Returns the name of the stock/asset - used to determine
        whether to shift the object */
    virtual string sensName(SpotLevelProbability* shift) const;

    /** Shifts the object using given shift. Return true to make the
        infrastructure keep tweaking the components within the object
        which implements this interface */
    virtual bool sensShift(SpotLevelProbability* shift);


    /** returns sensitive strikes for a given vol request */
    virtual void getSensitiveStrikes(
        const CVolRequest* volRequest,
        OutputNameConstSP outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP sensitiveStrikes) const;

    DividendListConstSP getDivList() const;

    /** could replace this by a collector */
    string getVolName() const;

    /** return the vol wrapper */
    const CVolBaseWrapper& getVol()const {return vol;}

    virtual EquitySP getEquity() const{ return equity;}

    /** given a current spot level, get the next strike on the vol surface where 
        the slope is non-differentiable */
    virtual double getNextStrike(const double& strike,
                                 bool          isUp,
                                 bool&         offSurface) const;

    /** return a pdf calculator */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /* for IEqVolNamePair */
    virtual bool getNamePairs(string& eqName, string& volName) const;

protected:
    EquityBase(const CClassConstSP& clazz): CAsset(clazz){}
    EquityBase(const CClassConstSP& clazz,
               const Equity*   equity,
               const CVolBase* vol);
    EquityBase();

    static void acceptValueDateCollector(const EquityBase* asset, 
                                         CValueDateCollector* collector);

    static void acceptDeltaShift(const EquityBase* asset, 
                                 ShiftSizeCollector* collector);

    static void acceptImntCcy(const EquityBase* asset,
                              AssetCcyCollector* collector);

    static void acceptHoliday(const EquityBase* asset,
                              HolidayCollector* collector);

    static void acceptCriticalDateCollector(const EquityBase* simpleEquity, 
                                            CriticalDateCollector* collector);

    // data fields
    EquitySP        equity; // $unregistered
    CVolBaseWrapper vol; // $unregistered
private:
    EquityBase(const EquityBase& rhs);
    EquityBase& operator=(const EquityBase& rhs);
};

typedef smartPtr<EquityBase> EquityBaseSP;
typedef smartConstPtr<EquityBase> EquityBaseConstSP;

// support for wrapper class
typedef MarketWrapper<EquityBase> EquityBaseWrapper;

DRLIB_END_NAMESPACE
#endif
