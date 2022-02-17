//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiAsset.hpp
//
//   Description :
//
//   Author      : Linus Thand, November 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#ifndef QLIB_MULTIASSET_HPP
#define QLIB_MULTIASSET_HPP

#include "edginc/MultiFactors.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/LocalCorrSqueeze.hpp"

DRLIB_BEGIN_NAMESPACE


/** This is not an asset, but captures multiple factors data and
    implements a choice of views of them. */
class MARKET_DLL MultiAsset: public CObject,
                  public virtual IMultiMarketFactors,
                  public virtual IAsMultiFactors, 
                  public virtual Theta::IShift {
    class AsMultiFactors;
    friend class AsMultiFactors;
    YieldCurveWrapper  yc;              // multi's currency
    IMarketFactorWrapperArray assets;        // array of factors ('assets')
    CStringArray       ccyTreatments; // None, Struck or Prot for each asset
    DoubleMatrix       correlations;  /* between assets - this should be 
                                           initialised in the getMarketData 
                                           method using  Correlations */
    // transient & tweakable fields
    CorrelationCommonArray  corrObjects;        // array of (n^2-n)/2 correlation objects
    CorrelationTermArray    corrTermArray;      // array of (n^2-n)/2 correlation term objects
    LocalCorrSqueezeArray   localCorrSqueezeArray;
    IntArray                localCorrSqueezeIndexArray;
    DateTime                valueDate;
    BoolArray               skipFwdCorrelation; // used for bonds, commodities, nasty funds, ...

    // not used variables (but too painful to remove from interface)
    bool                    useCorrelationTermNotUsed;  // $unregistered
    DoubleMatrix            shortTermSpreadNotUsed;  // $unregistered
    DoubleMatrix            longTermSpreadNotUsed; // $unregistered
    CorrelationTermArray    corrTermArrayNotUsed; // $unregistered

    // For Sampras, many assets are diffused but there are very few
    // correlations.  Storing 0 correlations results in few Gigs of
    // memory being needed.  The following parameter was added to
    // address this.
    // When strictCorr = false, missing correlations are ignored
    // so that the associated corrObjects[pos] will be null.
    // When strictCorr = true, missing correlations result in exceptions.
    // True is the default value as this is the previous behavior in Qlib.
    bool               strictCorr;    


public:
    static CClassConstSP const TYPE;

    void checkFactorRange(const char* routine, int iFactor) const;

    /** Returns the number of factors in this collection (this does not include
        MarketFactors which are contained within top level MarketFactors) */
    virtual int numFactors() const;

    //// for backwards compatibility
    virtual int NbAssets() const;

    //// Returns the correlations between the immediate factors.
    //// This is numFactors() x numFactors() 
    CDoubleMatrixConstSP factorsCorrelationMatrix() const;

    /** Returns the name of the specified 'market factor' */
    virtual string getName(int index) const;
    
    /** Returns the MarketFactor with the specifed index */
    virtual IMarketFactorConstSP getFactor(int index) const;

    /** Pull out the component assets & correlations from the market data */
    void getMarket(const IModel* model, const MarketData* market);

    /** Validation */
    void validatePop2Object();

    // Maps onto assetCrossValidate for each constituent
    void crossValidate(const DateTime&       startDate,
                       const DateTime&       valueDate,
                       const YieldCurve*     discCcy,
                       const CInstrument*    instrument) const;
   
    
    // implementation below
    virtual IMultiFactors* asMultiFactors() const;
    
    /** Returns a list of factor indexes indicating which asset is
        sensitive to combination of supplied SensControl and
        OutputName. If crossAssetSensitivities is true then this
        only covers sensitivities such as phi which can live at
        this level otherwise these are excluded. name must not be
        NULL */
    IntArray getSensitiveFactors(
        const IPerNameSensitivity* sens,
        bool                       crossAssetSensitivities,
        const OutputName*          name) const;

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    IntArray getSensitiveFactors(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const;

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    IntArray getSensitiveAssets(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const;

    /** Shifts the object using given shift (see Theta::Shift)*/
    bool sensShift(Theta* shift);

        
    /** record forwards at maturity for each factor */
    void recordFwdAtMat(OutputRequest*  request,
                        Results*        results,
                        const DateTime& maturityDate) const;

    // returns ccy treatment
    string getCcyTreatment(int iFactor) const;

    MultiAsset(YieldCurveWrapper yc,
               IMarketFactorWrapperArray assets,
               CStringArray ccyTreatments);

private:
    /* for reflection */
    MultiAsset();
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultMultiAsset();
};


DECLARE(MultiAsset)

/** View of MultiAsset as IMultiFactors */
class MultiAsset::AsMultiFactors: public CObject, 
                  public virtual IMultiFactors {
    CAssetArray       assets;        // array of assets
    MultiAssetConstSP multiAsset;    // the original
public:
    static CClassConstSP const TYPE;

    AsMultiFactors(const CAssetArray& assets,
                   MultiAssetConstSP  multiAsset);

    // MultiFactors interface. Here this is a degenerate case with
    // NbFactors==NbAssets
    int NbAssets() const;

    int numFactors() const;

    /** Returns the name of the specified 'market factor' */
    virtual string getName(int index) const;

    /** Returns the MarketFactor with the specifed index */
    virtual IMarketFactorConstSP getFactor(int index) const;
    
    /** Pull out the component assets & correlations from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    
    /** Validate that the market data and instrument are consistent with
        each other etc */
    virtual void crossValidate(const DateTime&       startDate,
                               const DateTime&       valueDate,
                               const YieldCurve*     discCcy,
                               const CInstrument*    instrument) const;

    void checkAssetRange(const char* routine, int iAsset) const;
    
    string assetGetName(int iAsset) const;
    string assetGetTrueName(int iAsset) const;
    string assetGetCcyTreatment(int iAsset) const;
    double assetGetSpot(int iAsset) const;
    double assetFwdValue(int             iAsset,
                         const DateTime& date) const;
    void assetFwdValue(int                  iAsset,
                       const DateTimeArray& dates,
                       CDoubleArray&        result) const;
        
    int NbFactors() const;
    void checkFactorRange(const char* routine, int iFactor) const;

    string factorGetName(int iFactor) const;

    CVolProcessed * factorGetProcessedVol(
        int                iFactor,
        const CVolRequest* volRequest) const;
        
    double factorGetSpot(int iFactor) const;

    double factorFwdValue(int             iFactor,
                          const DateTime& date) const;
        
    void factorFwdValues(int                  iFactor,
                         const DateTimeArray& dates,
                         CDoubleArray&        result) const;

    CDoubleMatrixConstSP factorsCorrelationMatrix() const;
        
    virtual BoolArray getSkipFwdCorrArray() const;

    //virtual CorrelationArray getCorrObjArray() const;
    virtual CorrelationCommonArray getCorrObjArray() const;
    virtual CorrelationTermArray getCorrTermObjArray() const;
    virtual LocalCorrSqueezeArray getLocalCorrSqueezeArray() const {
        return multiAsset->localCorrSqueezeArray;
    }
    virtual IntArray getLocalCorrSqueezeIndexArray() const {
        return multiAsset->localCorrSqueezeIndexArray;
    }
    virtual TimeMetricArray getTimeMetricArray() const;

    /** record forwards at maturity for each factor */
    void recordFwdAtMat(OutputRequest*  request,
                        CResults*       results,
                        const DateTime& maturityDate) const;

    void CollapseFactors(
        int             len,                 // length of each array
        const double**  factorLevels,        // [NbFactors] arrays
        double**        assetLevels) const;
        
    /** Using supplied outputName to select the relevent asset, adds
        sensitive strikes to the sensitiveStrikes parameter */
    void assetGetSensitiveStrikes(
        int                               iAsset,
        const CVolRequest*                volRequest,
        const OutputNameConstSP&          outputName,
        const SensitiveStrikeDescriptor&  sensStrikeDesc,
        const DoubleArraySP&              sensitiveStrikes) const;

    /** Returns a reference to the i-th Asset */
    const CAsset& getAsset(int iAsset) const;
        
    virtual IntArray getSensitiveFactors(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const;

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    virtual IntArray getSensitiveAssets(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const;

    /** Returns PDF calculator for iAsset */
    PDFCalculator* assetPdfCalculator(const PDFRequest* request,
                                      int               iAsset) const;
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};




DRLIB_END_NAMESPACE


#endif
