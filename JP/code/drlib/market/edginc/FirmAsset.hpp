//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FirmAsset.hpp
//
//   Description : Data for firm asset diffusion models
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : September 3, 2001
//
//
//----------------------------------------------------------------------------

#ifndef FIRMASSET_HPP
#define FIRMASSET_HPP

#include <string>
#include "edginc/config.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Class.hpp"
#include "edginc/Asset.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/LiquiditySpreadCurve.hpp"
#include "edginc/AssetVegaParallel.hpp"
#include "edginc/AssetVegaPointwise.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VegaPointwise.hpp"
#include "edginc/LiquiditySpreadRhoParallel.hpp"
#include "edginc/LiquiditySpreadPointwise.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Equity.hpp"
#include "edginc/ObjectIteration.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"


DRLIB_BEGIN_NAMESPACE

/* Forward declare the classes refered to from CDSParSpreads if we are not
 * interested in their storage properties (ie, they are used through (smart) 
 * pointers and therefore their include files are not required here) */
FORWARD_DECLARE_WRAPPER(CDSParSpreads);


class MARKET_DLL FirmAsset : public CAsset,
                  virtual public ObjectIteration::IOverride, 
                  public AssetVegaParallel::IShift,
                  public AssetVegaPointwise::IShift,
                  public LiquiditySpreadRhoParallel::IShift,
                  public ITweakableWithRespectTo<LiquiditySpreadPointwise>,
                  public Theta::IShift {
public:
    static CClassConstSP const TYPE;
    friend class FirmAssetHelper;

    virtual ~FirmAsset();

    void validatePop2Object();

    virtual string getName() const {
        return name;
    }
 
    double CalcNoDefaultProb(const DateTime& fromDate, 
                             const DateTime& toDate);

    double CalcNoDefaultProb(const DateTime& fromDate, 
                             const DateTime& toDate,
                             const double& spotPrice);


    virtual void getMarket(const IModel* model, const MarketData* market);

    void calculateProcessedVol(const DateTime& maturityDate) const;

    double getAssetSpot() const;

    void validate();

    CDSParSpreadsWrapper getLiquiditySpreadCurve() const;

    /** Returns name identifying vol for asset vega parallel */
    virtual string sensName(AssetVegaParallel* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(AssetVegaParallel* shift);
    
    /** Returns name identifying vol for vega pointwise */
    virtual string sensName(AssetVegaPointwise* shift) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryArrayConstSP sensExpiries(AssetVegaPointwise* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(AssetVegaPointwise* shift);

    /** Returns name identifying vol for asset vega parallel */
    virtual string sensName(LiquiditySpreadRhoParallel* shift) const;
    /** Shifts the object using given shift */
    virtual bool sensShift(LiquiditySpreadRhoParallel* shift);

    /** Returns name identifying vol for vega pointwise */
    virtual string sensName(const LiquiditySpreadPointwise*) const;
    /** Returns the array of expiries (ie maturities/benchmark dates) that
        need to be tweaked for this vol */
    virtual ExpiryWindowArrayConstSP sensQualifiers(const LiquiditySpreadPointwise*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<LiquiditySpreadPointwise>&);

    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    // switch off the tweaks which we don't want to do for assets
    bool recurse(const CFieldConstSP& field,
                 const CClassConstSP& targetClass) const;

    /** returns the spot price */
    virtual double getSpot() const;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const;
    

    /** Calculates the expected spot price of the asset at the given date.
        Do not use this repeatedly to calculate values over a set of
        dates (poor performance) - instead use other fwdValue method */
    virtual double fwdValue(const DateTime& date) const;

    /** Calculates the expected spot prices of the asset at the given dates
        respecting any 'algorithmic' choices set in the FwdValueAlgorithm */
    virtual void fwdValue(const DateTimeArray&     dateList,
                          const FwdValueAlgorithm& algo,
                          CDoubleArray&            result) const;

    /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const;

    /** return a pdf calculator provided our vol supports it */
    virtual PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const;

    double  getLiquiditySpread(const DateTime& valueDate, const DateTime& maturity,
                               const BadDayConvention* bdc, const Holiday* hols);

    double getLiquiditySpreadsRecovery()const;

    CAssetWrapper getEquityAsset();

    double getDefaultBarrier();

    double getLambda();

    /** record forwards at maturity*/
    void recordFwdAtMat(OutputRequest*  request,
                        CResults*       results,
                        const DateTime& maturityDate) const;

    FirmAsset(const string&     name,
              const double&     assetVol,
              const double&     lambda,
              const double&     globalRecovery,
              const double&     debtPerShare,
              EquitySP          equity,
              const double&     liquiditySpread,
              const bool&       payAccrued,
              const double&     recovery,
              const DateTime&   valueDate);

    void createStruckFirmAsset(const IModel*     model,
                               const MarketData* market, 
                               const string&     payOutYCName);

    double getFXRate();

    //// Returns the ccy treatment for the asset. Default returns 
    virtual string getCcyTreatment() const;

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

protected:
    FirmAsset(const FirmAsset& rhs);
    FirmAsset& operator=(const FirmAsset& rhs);

    FirmAsset();

    FirmAsset(const CClassConstSP& clazz);

    static void acceptWrapperNameCollector(const FirmAsset* firmAsset,
                                           WrapperNameCollector* collector);

    string                  name;
    CVolBaseWrapper         assetVol;
    double                  lambda;
    double                  lBar;
    double                  debt;
    CAssetWrapper           asset;
    CDSParSpreadsWrapper    cdsLiquiditySpreads;
    DateTime                valueDate;

    // private members
    mutable CVolProcessedBSSP processedVol;
    mutable double            processedAssetVol;
    bool                      isStruck;
    CAssetConstSP             plainAsset;
};

typedef smartConstPtr<FirmAsset> FirmAssetConstSP;
typedef smartPtr<FirmAsset>      FirmAssetSP;

// support for wrapper class
typedef MarketWrapper<FirmAsset> FirmAssetWrapper;

DRLIB_END_NAMESPACE

#endif
