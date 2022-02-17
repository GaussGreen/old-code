//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetDDE.hpp
//
//   Description : spread curve for DDE. keep functionality of par/credit spread
//                                      and provide link to DDEModule
//
//
//   Author      : Qing Hou
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ASSET_DDE_H
#define EDG_ASSET_DDE_H

#include "edginc/CDSParSpreads.hpp"
#include "edginc/CleanSpreadCurve.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/CanBeRisky.hpp"
#include "edginc/SpreadEquityFunc.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/VolBase.hpp"
#include "edginc/Spot.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/ParSpreadParallel.hpp"
#include "edginc/ParSpreadPointwise.hpp"
#include "edginc/DDEParams.hpp"

DRLIB_BEGIN_NAMESPACE

class CleanSpreadVolCurve;
class EquitySpreadCorrCurve;
class SpreadEquityFunc;

/** base class to support credit sensitivity calc for DDE. 
 ** responsible for getting appropriate clean spread curve 
 **/
class MARKET_DLL CreditCurveDDE: virtual public CreditSpreadRhoParallel::RestorableShift,
                      virtual public CreditSpreadRhoPointwise::IRestorableShift,
                      virtual public AdjCreditSpreadRhoParallel::RestorableShift,
                      virtual public AdjCreditSpreadRhoPointwise::IRestorableShift,
                      virtual public CreditSpreadLevel::Shift,
                      virtual public ITweakableWithRespectTo<ParSpreadParallel>,
                      virtual public ITweakableWithRespectTo<ParSpreadPointwise>
{
public:
    static const string DEF_CREDIT_TYPE;

    virtual ~CreditCurveDDE();

    //////////////////////////////////////////////////
    //                      CreditSpreadCurve functions                     //
    //////////////////////////////////////////////////

    /** SpreadRhoParallel support */
    string sensName(CreditSpreadRhoParallel* shift) const
    { return creditName; }
    bool sensShift(CreditSpreadRhoParallel* shift);
    void sensRestore(CreditSpreadRhoParallel* shift);

    /** SpreadRhoPointwise support */
    string sensName(CreditSpreadRhoPointwise* shift) const
    { return creditName; }
    ExpiryArrayConstSP sensExpiries(CreditSpreadRhoPointwise* shift) const
    { return creditExpiries; }
    bool sensShift(CreditSpreadRhoPointwise* shift);
    void sensRestore(CreditSpreadRhoPointwise* shift);
    
    /** SpreadRhoParallel support */
    string sensName(AdjCreditSpreadRhoParallel* shift) const
    { return creditName; }
    bool sensShift(AdjCreditSpreadRhoParallel* shift);
    void sensRestore(AdjCreditSpreadRhoParallel* shift);

    /** AdjSpreadRhoPointwise support */
    string sensName(AdjCreditSpreadRhoPointwise* shift) const
    { return creditName; }
    ExpiryArrayConstSP sensExpiries(AdjCreditSpreadRhoPointwise* shift) const
    { return creditExpiries; }
    bool sensShift(AdjCreditSpreadRhoPointwise* shift);
    void sensRestore(AdjCreditSpreadRhoPointwise* shift);
    
    /** Implements CreditSpreadLevel scenario */
    /** Shifts the object using given shift (see CreditSpreadLevel::Shift)*/
    string sensName(CreditSpreadLevel* shift) const
    { return creditName; }
    bool sensShift(CreditSpreadLevel* shift);

    //////////////////////////////////////////////////
    //                      ParSpreadCurve functions                        //
    //////////////////////////////////////////////////

    /** ParSpreadParallel support */
    string sensName(const ParSpreadParallel*) const
    { return creditName; }
    TweakOutcome sensShift(const PropertyTweak<ParSpreadParallel>&);
    
    /** ParSpreadRhoPointwise support */
    string sensName(const ParSpreadPointwise*) const;
    ExpiryWindowArrayConstSP sensQualifiers(const ParSpreadPointwise*) const;
    TweakOutcome sensShift(const PropertyTweak<ParSpreadPointwise>& shift);
    
    //////////////////////////////////////////////////
    //                      Own functions                                           //
    //////////////////////////////////////////////////

    /** mark that the asset needs to be calibrated */
    void calibCredit();
    void calibCreditRestore();
        
    void getMarket(const IModel* model, const MarketData* market, const YieldCurve *yc);

    CleanSpreadCurveSP getCleanSpreadCurve(DateTime valueDate) const;

    ExpiryArrayConstSP getExpiries() const
    { return creditExpiries; }

    bool isCreditTweaked() const { return creditTweaked; }

protected:
    CreditCurveWrapper              creditCurve;    // input credit curve. never actually used after initial setup
    string                                  creditType;             // type of credit curve

    // internal field
    bool                                    creditTweaked;
    bool                                    useParSpread;   // false if is CreditSpreadCurve
    string                                  creditName;
    YieldCurveSP                    creditYC;
    ExpiryArrayConstSP              creditExpiries;

    CreditSpreadCurveSP             creditSprds;
    CDSParSpreadsSP                 parSprds;

    mutable CleanSpreadCurveSP              cleanSprds; // to accommodate getCleanSpreadCurve() const
    CleanSpreadCurveSP              cleanSprdsBak;

protected:
    CreditCurveDDE();

    /** build and retrieve various spread curves */
    bool buildCreditSpreads();

    bool buildCDSParSpreads();

    bool buildCleanSpreadCurve(DateTime valueDate) const;

    CreditSpreadCurveSP getCreditSpreadsForShift();
    ParSpreadCurveSP getParSpreadsForShift();
private:
    CreditCurveDDE(const CreditCurveDDE& rhs);
    CreditCurveDDE& operator=(const CreditCurveDDE& rhs);
};

typedef smartConstPtr<CreditCurveDDE> CreditCurveDDEConstSP;
typedef smartPtr<CreditCurveDDE> CreditCurveDDESP;


/** base class to support equity sensitivity calc for DDE.
 **/
class MARKET_DLL EquityDDE :       public IHaveEquity,
                        virtual public ITweakableWithRespectTo<Spot>,
                        virtual public DeltaDDE::RestorableShift,
                        virtual public MuParallel::IShift,
                        virtual public MuSpecial::IShift,
                        virtual public MuPointwise::IShift,
                        virtual public SpotLevel::Shift,
                        virtual public SpotShift::Shift,
                        virtual public RhoBorrowParallel::RestorableShift,
                        virtual public RhoBorrowPointwise::IRestorableShift,
                        virtual public BorrowParallelShift::Shift,
                        virtual public BorrowLevel::Shift,
                        virtual public IRestorableWithRespectTo<VolParallel>,
                        virtual public ITweakableWithRespectTo<VolPointwise>,
                        virtual public VegaSkewParallel::IShift,
                        virtual public VegaSkewPointwise::IShift,
                        virtual public VolLevel::Shift,
                        virtual public VolRelativeShift::IShift
{
public:

    /** Delta support */
    string sensName(const Spot*) const;
    TweakOutcome sensShift(const PropertyTweak<Spot>& shift);
    void sensRestore(const PropertyTweak<Spot>& shift);

    /** DeltaDDE support */
    string sensName(DeltaDDE* shift) const;
    bool sensShift(DeltaDDE* shift);
    void sensRestore(DeltaDDE* shift);

    /** MU_PARALLEL support */
    string sensName(MuParallel* shift) const;
    bool sensShift(MuParallel* shift);

    /** MU_S support */
    string sensName(MuSpecial* shift) const;
    bool sensShift(MuSpecial* shift);

    /** MU_POINTWISE support */
    string sensName(MuPointwise* shift) const;
    bool sensShift(MuPointwise* shift);
  
    /** SpotLevel support */
    string sensName(SpotLevel* shift) const;
    bool sensShift(SpotLevel* shift);

    /** SpotShift support */
    string sensName(SpotShift* shift) const;
    bool sensShift(SpotShift* shift);

    /** RhoBorrowParallel support */
    string sensName(RhoBorrowParallel* shift)const;
    bool sensShift(RhoBorrowParallel* shift);
    void sensRestore(RhoBorrowParallel* shift);

    /** RhoBorrowPointwise support */
    string sensName(RhoBorrowPointwise* shift)const;
    ExpiryArrayConstSP sensExpiries(RhoBorrowPointwise* shift)const;
    bool sensShift(RhoBorrowPointwise* shift);
    void sensRestore(RhoBorrowPointwise* shift);  

    /** BorrowParallelShift support */
    string sensName(BorrowParallelShift* shift)const;
    bool sensShift(BorrowParallelShift* shift);

    /** BorrowLevel support */
    string sensName(BorrowLevel* shift)const;
    bool sensShift(BorrowLevel* shift);

    /** VolParallel */
    string sensName(const VolParallel*) const;
    TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);
    void sensRestore(const PropertyTweak<VolParallel>& tweak);

    /** VegaPointwise Interface */
    string sensName(const VolPointwise*) const;
    ExpiryWindowArrayConstSP sensQualifiers(const VolPointwise*) const;
    TweakOutcome sensShift(const PropertyTweak<VolPointwise>&);

    /** VegaSkewParallel support */
    string sensName(VegaSkewParallel* shift) const;
    bool sensShift(VegaSkewParallel* shift);

    /** VegaSkewPointwise support */
    string sensName(VegaSkewPointwise* shift) const;
    ExpiryArrayConstSP sensExpiries(VegaSkewPointwise* shift) const;
    bool sensShift(VegaSkewPointwise* shift);

    /** VolLevel support */
    string sensName(VolLevel* shift) const;
    bool sensShift(VolLevel* shift);

    /** VolRelativeShift support */
    string sensName(VolRelativeShift* shift) const;
    bool sensShift(VolRelativeShift* shift);

    //////////////////////////////////////////////////
    //                      Own functions                                           //
    //////////////////////////////////////////////////

    /** mark that the asset needs to be calibrated */
    void calibEquity();
    void calibEquityRestore();
    void getMarket(const IModel* model, const MarketData* market);

    EquitySP getEquity() const;

    CVolBaseConstSP getVol() const;

    string getVolName() const;

    TimeMetricConstSP getTimeMetric() const;

    bool isEquityTweaked() const { return equityTweaked; }

protected:
    CAssetWrapper       asset;
        
    // internal field
    bool                    equityTweaked;

    CVolBaseConstSP         volBase;
    TimeMetricConstSP       volTimeMetric;

protected:
    EquityDDE();
};

typedef smartPtr<EquityDDE> EquityDDESP;
typedef smartConstPtr<EquityDDE> EquityDDEConstSP;


/** Implementation of CAsset for a DDE risky equity.
 ** It's composed of EquityDDE and CreditCurveDDE */
class MARKET_DLL AssetDDE: public CAsset,
                public CreditCurveDDE,
                public EquityDDE,
                public virtual CAsset::ICanHaveDDE,
                public virtual CAsset::IStruck,
                virtual public Theta::IShift  // theta affects both equity and credit
{
public:
    static CClassConstSP const TYPE;
    friend class AssetDDEHelper;
    friend class DDEModule;

    //////////////////////////////////////////////////
    //                      Asset functions                                         //
    //////////////////////////////////////////////////


    double getSpot() const 
    { return asset->getSpot(); }

    string getName() const 
    { return name; }

    string getTrueName() const; 

    CVolProcessed* getProcessedVol(const CVolRequest* volRequest) const;

    double fwdValue(const DateTime& date) const 
    { return asset->fwdValue(date); }

    void fwdValue(const DateTimeArray& dateList,
                  CDoubleArray&        result) const
    { asset->fwdValue(dateList, result); }

    void fwdValue(const DateTimeArray&     dateList,
                  const FwdValueAlgorithm& algo,
                  CDoubleArray&            result) const
    { asset->fwdValue(dateList, algo, result); }

    string getYCName() const
    { return asset->getYCName(); }

    void recordFwdAtMat(OutputRequest*  request,
                        CResults*       results,
                        const DateTime& maturityDate) const
    { asset->recordFwdAtMat(request, results, maturityDate); }

    double fwdFwd(const DateTime& spotDate,
                  double          spot, 
                  const DateTime& fwdDate) const
    { return asset->fwdFwd(spotDate, spot, fwdDate); }

    void fwdFwd(const DateTime&      spotDate,
                double               spot, 
                const DateTimeArray& fwdDates,
                DoubleArray&         results) const
    { asset->fwdFwd(spotDate, spot, fwdDates, results); }

    double expNumberShares( const DateTime& valueDate, 
                            const DateTime& fwdDate,
                            const bool& convAdjust) const
    { return asset->expNumberShares( valueDate, fwdDate, convAdjust); }

    DateTime settleDate(const DateTime& tradeDate) const
    { return asset->settleDate(tradeDate); }

    void getSensitiveStrikes(
        const CVolRequest* volRequest,
        OutputNameConstSP outputName,
        const SensitiveStrikeDescriptor& sensStrikeDesc,
        DoubleArraySP sensitiveStrikes) const
    {
        asset->getSensitiveStrikes(volRequest, outputName, sensStrikeDesc, sensitiveStrikes);
    }

    IMultiFactors* asMultiFactors() const;

    PDFCalculator* pdfCalculator(const PDFRequest* request) const;

    void getMarket(const IModel* model, const MarketData* market);

    /** THETA support */
    bool sensShift(Theta* shift);

    //////////////////////////////////////////////////
    //                      Own functions                                           //
    //////////////////////////////////////////////////

    void validatePop2Object();

    const DDEParams *getDDEParams() const;

    /** function related to ccy struck */
    void makeStruck(const IModel*  model,                           // name for struck asset
                    const MarketData* market,
                    const string&     payOutYCName);// name for pay out ccy

    void makeProt(const IModel*  model,                                     // name for prot asset
                  const MarketData* market,
                  const string&     payOutYCName);// name for pay out ccy

    /** (IStruck interface) Returns the fx spot */
    virtual double getFXSpot() const;

    /** (IStruck interface) Returns the fx forward value */
    virtual double fxFwdValue(const DateTime& date) const;

    bool isDDE() const { return true; }

    bool hasRiskyVol() const; // return wether have risky vol overwrite

    CVolBaseSP getRiskyVol(DateTime valueDate) const; // return wether have risky vol overwrite

    //////////////////////////////////////////////////
    //                      DDEModule functions                                     //
    //////////////////////////////////////////////////

    void getPv(DateTime valueDate, const DateTimeArray &dates, DoubleArray &pvs) const;
    void getFwd(const DateTimeArray &dates, DoubleArray &fwds) const;
    CVolProcessedBS* getUnderlyerProcessedVol(
        const CVolRequestLN* volRequest) const; // get underlyer's vol
 
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
    AssetDDE();
    //AssetDDE(const AssetDDE& rhs);
    //AssetDDE& operator=(const AssetDDE& rhs);

    string          name;
    double          riskyness;                  // range (0,1) for eq risk participation

    DDEParamsSP ddeParams;

    // risky vol overwrite parameters
    ExpiryArraySP   riskyVolExps;
    DoubleArraySP   riskyVolVals;
};

// smart pointer support
typedef smartPtr<AssetDDE> AssetDDESP;
typedef smartConstPtr<AssetDDE> AssetDDEConstSP;


DRLIB_END_NAMESPACE
#endif
