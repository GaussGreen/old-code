//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRMGenSV.hpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------
#ifndef MCPathConfigSRMGenSV_HPP
#define MCPathConfigSRMGenSV_HPP

#include "edginc/QMCGenDiffusibleAsset.hpp"
#include "edginc/IElemStateVariableGenVisitor.hpp"

DRLIB_BEGIN_NAMESPACE 


/** Holds information about the different state variable generators */
class MCARLO_DLL SV : public IElemStateVariableGenVisitor
{
public:

    /** fields of struct SV */
    string                  domISOCode;    // of payoff currency
    YieldCurveConstSP       domYC;         // domestic yield curve
    const IMultiMarketFactors* multiFactors;  // for resolving spot requests; also gives access to dates for SV
    int                     numFX;
    int                     numEq;
    int                     numCr;
	int						numEnrg;
    const IPastValues*      pastValues; // historic values

//    int                     storedNumIRFactors;
    int                     storedNumAllFactors;

    vector<const SVGenSpot*>   spots;         // MC SV Gen for asset spot

    map<string, QMCGenDiffusibleIRSP>      mapIRGen;        // data per ccy IR, keyed by ccy name
    map<string, QMCGenDiffusibleFXSP>      mapFXGen;         // data per ccy FX, keyed by ccy name
    map<string, QMCGenDiffusibleEquitySP>  mapEquityGen;     // data per equity, keyed by equity name
    map<string, QMCGenDiffusibleCreditSP>  mapCreditGen;     // data per credit, keyed by credit name
	map<string, QMCGenDiffusibleEnergySP>  mapEnergyGen;    // data per energy, keyed by energy name
	map<string, QMCGenDiffusibleBasisSpreadSP> mapBasisGen;  // data per basis, keyed by basis name

    vector<QMCGenDiffusibleIRSP>          orderedCcys;   // mapIRGen ordered into an array
    vector<QMCGenDiffusibleEquitySP>      orderedEquities;//mapEquityGen ordered into an array
    vector<QMCGenDiffusibleFXSP>          orderedFXs;    //mapFXGen ordered into an array
    vector<QMCGenDiffusibleCreditSP>      orderedCredits;// mapCreditGen ordered into an array
	vector<QMCGenDiffusibleEnergySP>      orderedEnergys;// mapEnergyGen ordered into an array
	vector<QMCGenDiffusibleBasisSpreadSP> orderedBasis;  // mapBasisGen ordered into an array
    vector<QMCGenPureJumpSP>              orderedPureJumps;     // pure Poisson generators

    vector<IQMCGenDiffusibleAssetSP>      byAssetPosition;// all the polymorphic asset generators 

    vector<const SVGenPathWeight*>        pathWeightSVGen;

    // these two data members is a temporary solution for being somehow 
    // unable to lookup eq and fx by name (this is legacy issue, i just moved the solution)
    map<int, QMCGenDiffusibleFXSP> mFXbyId; // .first is not a size_t as it can be -1
    map<int, QMCGenDiffusibleEquitySP> mEQbyId;// .first is not a size_t as it can be -1

    // constructor
    SV(MCPathConfigSRM*       pathConfig,
       const MCProductClient* prodClient); 

    /** Counts the total number of factors in simulation */
    int numAllFactors() const { return storedNumAllFactors; }

    // this is to implement 'visitor' model of SVGen sorting
    virtual void processSVGen(const SVGenDiscFactor* df);
    virtual void processSVGen(const SVGenExpectedDiscFactor* edf);
    virtual void processSVGen(const SVGenSurvivalDiscFactor* sdf);
    virtual void processSVGen(const SVGenExpectedSurvivalDiscFactor* esdf);
    virtual void processSVGen(const SVGenAggregatedSurvivalDiscFactor* asdf);
    virtual void processSVGen(const SVGenExpectedSpot* es);
    virtual void processSVGen(const SVGenSpot* s);
    virtual void processSVGen(const SVGenExpectedEnergyFuture* eEnrgFutureGen);
    virtual void processSVGen(const SVGenExpectedBasisFwdSpread* s);
    virtual void processSVGen(const SVGenDateOfDefault* dod);
    virtual void processSVGen(const SVGenPathWeight* pw);
	// catch-all default that throws an error
    virtual void processSVGen(const IElemStateVariableGen* base);
    
    int processAllSVGens(StateVariableCollectorSP    svCollector);

private:

    /** Check that all FX asset have the domestic currency as their 'baseCcy'. */
    void validateBaseFXCurrencies(const CAssetArray& diffAssets, YieldCurveConstSP domYC);

    /** Sets ccy treatment for supporting assets. */
    void setCurrencyTreatment(IQMCGenDiffusibleAssetSP perAssetData, IMarketFactorConstSP factor);

    /** Creates a QMCGenDiffusibleIRSP according to the Rates model type */
    QMCGenDiffusibleIRSP diffusibleIRGenSPFromModelType(string ratesModelType);

    /** Saves this discount yield curve in a map against the isocode of the ccy */
    void saveDiscountYC(MCPathConfigSRM* pathConfig,
                        YieldCurveConstSP discYC);

    /** Creates a QMCGenDiffusibleCreditSP according to the Credit model type */
    QMCGenDiffusibleCreditSP diffusibleCreditGenSPFromModelType(string creditModelType, MCPathConfigSRM* pathConfig);

    /** Saves CDS curve in a map against the name of the credit */
    void buildPerCreditMap(MCPathConfigSRM*     pathConfig,
                            ICDSParSpreadsArray& diffCDSs);

    /** Saves CDS curve in a map against the name of the credit for CID model */
    void buildPerCreditMapCID(MCPathConfigSRM*     pathConfig,
                              ICDSParSpreadsArray& diffCDSs);

	/** Create map to store energy future related info */
	void buildPerEnergyMap(MCPathConfigSRM * pathConfig,
                            EnergyFuturesCurveArray & diffEnrgs);

	/** call this after we have processed energy SV gens */
	void buildPerBasisMap(
		MCPathConfigSRM * pathConfig,
		IBasisIndexCurveArray & diffBasis
		);

    // Build 'mapEquityGen' and 'mapFXGen' maps and extend 'mapIRGen' map.
    // Use compositeAssets and diffAssets to connect to multiFactors/SVs
    void extractFromMultiFactors(MCPathConfigSRM* pathConfig,
                                 string domISOCode,
                                 CAssetArray& diffAssets,
                                 vector<const SVGenSpot*> &spots);

    void setupAssetCurrencyTreatment();
    
    /**    assigns all the randomIdx inside every PerXXX structure */
    void assignRandomIdx(const MCProductClient* prod); 
    void verifyRandomIdx(void); //< sanity checks for randomIdx assignment
    
    /** arranges all assets in their diffusion order, stores it in byAssetPosition */
    void setByAssetPosition();
};
typedef refCountPtr<SV> SVSP;


DRLIB_END_NAMESPACE
#endif

