//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCGenDiffusibleAsset.hpp
//
//   Description : A generator for creating QMCDiffusibleAsset classes
//                 this file contains only generic, i.e. non diffusion model
//                 specific
//                 collections of data, that MCPathConfigSRMGenSV will use
//
//
//
//----------------------------------------------------------------------------
#ifndef QMCGenDiffusibleAsset_HPP
#define QMCGenDiffusibleAsset_HPP

#include "edginc/YieldCurve.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/SRMEquityDiffuse.hpp"  // temporarily...

#include "edginc/SRMFXUtil.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMEnergyUtil.hpp"
#include "edginc/SRMBasisSpreadHJMUtil.hpp"
#include "edginc/SRMICE.hpp"

#include "edginc/QMCPureJumps.hpp"
#include "edginc/QMCCreditCIDJumps.hpp"
#include "edginc/QMCCreditComposite.hpp"
#include "edginc/QMCDetermRates.hpp"

#include "edginc/MCPathConfig.hpp"

#include "edginc/DECLARE.hpp"
#include "edginc/ISupportPathConfigSRM.hpp"

DRLIB_BEGIN_NAMESPACE

class QMCGenDiffusibleIR;
class QMCGenDiffusibleFX;
class QMCGenDiffusibleEquity;
class QMCGenDiffusibleCredit;
class QMCGenDiffusibleBasisSpread;
class QMCGenDiffusibleEnergy;
class MCPathConfig;
class SparseDoubleMatrix;


class IQMCGenDiffusibleAsset : public VirtualDestructorBase {
public:
    virtual ~IQMCGenDiffusibleAsset() {}

    IQMCGenDiffusibleAsset() : randomIdx(-1), ccyTreatment(ccyVanilla) {}


    // these methods are numbered in the order they will be called (important!):
    // (*) denotes an optional method


    // (1)
    virtual void setDomesticIRGen(QMCGenDiffusibleIR* /*pDomesticCcy*/) {} // override only if necessary

    // (2)
    virtual void initializeDiffusibleAssets() =0;

    // (3)
    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm) =0; // legacy of MultiFactor...

    // (4)
    virtual DateTimeArray getAssetDates() { return getDiffusibleAsset()->getAssetDates(); }
    virtual DateTimeArray getCriticalDates(MCPathConfig* mcPathConfig,
                                           const DateTime& start,       // likely to be "today"
                                           const DateTime& finish); // the latest requested date

    // (5)
    void createUtil(MCPathConfig*       mcPathConfig,
                    const DateTime&        today,
                    DateTimeArrayConstSP   simDates);

    // (5.5)
    /** For multi-factor assets, the expansion of the correlations to the multi-factors
    is handles by this method.  The matrix passed in has the correct dimensions
    so just store the results. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const = 0;

    // (5.6)
    /** For 1 factor return beta.  For multi-factor, return expanded beta vector 
        of dimension equal to the number of factors */
    virtual vector<double> expandTier2Betas(double beta) const = 0;

    // (5.7)
    /** Expand the asset-asset correlation matrix.  Implementation should be
    provided at the pathConfig specific level e.g. SRMGenDiffusibleXYZ or 
    CIDGenDiffusibleXYZ.  Also note that the corr matrix output is "flattened" 
    row-by-row into a 1D double vector. */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr/*inter-asset*/) = 0;
    // expand corrs from IR to current asset
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) = 0;
    // expand corrs from Credit to current asset
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) = 0;
    // expand corrs from Equity to current asset
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) = 0;
    // expand corrs from FX to current asset
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) = 0;
    // expand corrs from BasisSpread to current asset
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) = 0;
    // expand corrs from Energy to current asset
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) = 0;


    // (6)(*) - only if necessary
    /** Create series of swaptions for this IR factor. Also create and store
        associated pricing class. Append requested state variables to SVGens arrays. */
    virtual void    pushbackICERecord(StateVarDBase&, vector<ICE>& ) {} // override in IR class. necessary for ICE

    // (7)
    virtual void setDiffusibleAsset(
            const DateTime&     today,              // base date of the run
            MCPathConfig* mcPathConfig,       // to get max/min boundaries, etc
            const IPastValues*  pastValues,         // historic values
            DependenceSP        dependence);    // factorized correlations


    // might be unnecessary, but yet pretty handy.
    virtual bool    isRates() { return false; }
    virtual bool    isFX() { return false; }
    virtual bool    isEquity() { return false; }
    virtual bool    isCredit() { return false; }
    virtual bool    isEnrg() { return false; }
    virtual bool    isBS() { return false; }
    virtual bool    isJump() { return false; }
    virtual bool    isComposite() { return false; }
    virtual bool    isDefined() { return randomIdx>=0; }

    virtual QMCRatesUtilSP getRatesUtil() { return QMCRatesUtilSP(); } // override in IR class. necessary for correlation calc

    virtual int     numFactors() const { return 1; } // only those who are !=1 would override this

    // sets the random Idx to be used for this object, advances the index by the used number of factors, returns used number of factors
    // can be overloaded if something more fancy is needed, like initializing composite objects
    virtual int     setRandomIdx(int &rndIdx) 
    { 
        QLIB_VERIFY(randomIdx < 0, "Attempting to reassign the random idx for the asset " + 
            name + ". Internal error. ");
        if (!numFactors()) return 0;
        randomIdx = rndIdx; 
        rndIdx += numFactors(); 
        return numFactors(); 
    }

    virtual void setCurrencyTreatment(CurrencyTreatment _ccyTreatment);

    virtual int     getRandomIdx() const { return randomIdx; }


    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() =0;
    string          name;      // of the asset - convenient

protected:
    int                 randomIdx; // starting position of the asset in the "noise" array
    CurrencyTreatment   ccyTreatment;


};

DECLARE(IQMCGenDiffusibleAsset);

class SVGenExpectedSpot;

class QMCGenDiffusibleFX : public IQMCGenDiffusibleAsset
{
public:
    QMCGenDiffusibleFX(): assetIdx(-1), pForQMCGenDiffusibleIR(NULL), 
        multiSpots(NULL), pDomQMCGenDiffusibleIR(NULL) {}

    virtual void setDomesticIRGen(QMCGenDiffusibleIR* pDomesticCcy)
        { pDomQMCGenDiffusibleIR = pDomesticCcy; }

    // (5.5)
    /** Currently, all implementations are one factor so do nothing at this level. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Currently, FX is one dim'l so implement at this level */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(1,beta); }

    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm); // whether moment-matching is ON for a given SV

    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleFX(); }
    virtual bool isFX() { return true; }

    virtual IQMCDiffusibleFXSP getDiffusibleFX() = 0;
    virtual SRMFXUtilSP       getFXUtil() = 0; // TODO: see if necessary...

    string              ccy;      // allows accessing of mapIRGen (map keyed by ISO codes)
    int                 assetIdx; // index of this asset into multi; can be <0 meaning no corresponding entry
    FXAssetConstSP      fxAsset;  // fx asset

    vector<const SVGenSpot*>   *multiSpots; // just a local pointer to a global multispot for now...
    vector<const SVGenExpectedSpot*>  fxExpSpots;// fx spot state variables

    // naked pointers here to avoid having to resolve circular dependencies
    QMCGenDiffusibleIR*    pForQMCGenDiffusibleIR;  // allows easy access of corresponding IRGen
    QMCGenDiffusibleIR*    pDomQMCGenDiffusibleIR;  // allows easy access of domestic IRGen
    int getAssetIdx(void)   { return assetIdx;}
};
DECLARE(QMCGenDiffusibleFX);


class SVGenDiscFactor;
class SVGenExpectedDiscFactor;

class QMCGenDiffusibleIR : public IQMCGenDiffusibleAsset { // includes domestic currency
public:
    typedef pair<vector<const SVGenDiscFactor*>,
             vector<const SVGenExpectedDiscFactor*> > DFPair;
    typedef pair<IYieldCurveConstSP, IYieldCurveConstSP> YCPair;


    QMCGenDiffusibleIR(): nFactors(-1), pQMCGenDiffusibleFX(NULL)  {}
    virtual ~QMCGenDiffusibleIR() {}

    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm);

    virtual bool isRates() { return true; }
    virtual int numFactors() const { return nFactors; }
    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleIR(); }

    virtual IQMCDiffusibleInterestRateSP getDiffusibleIR() =0;
    virtual QMCRatesUtilSP                 getRatesUtil() =0;// TODO: see if necessary...

    int                     nFactors; // number of IR factors for diffusion
    YCPair                  ycs;  // discount/diffusion
    DFPair                  dfs;  // discount factors state variables
    QMCGenDiffusibleFX*        pQMCGenDiffusibleFX; // allows accessing of perFXs
};
DECLARE(QMCGenDiffusibleIR);


class QMCGenDiffusibleEquity : public IQMCGenDiffusibleAsset {
public:
    QMCGenDiffusibleEquity(): assetIdx(-1), multiSpots(NULL) {}

    virtual void setDomesticIRGen(QMCGenDiffusibleIR* pDomesticCcy)
        { pDomQMCGenDiffusibleIR = pDomesticCcy; }

    /** Currently, all implementations are one factor so do nothing at this level. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Currently, Equity is one dim'l so implement at this level */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(1,beta); }

    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm);

    virtual bool isEquity() { return true; }
    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleEQ(); }
    virtual IQMCDiffusibleEQSP getDiffusibleEQ() = 0;
    virtual SRMEquityUtilSP getEquityUtil() = 0;// TODO: see if necessary...

    string            ccy;      // allows accessing of mapIRGen (map keyed by ISO codes)
    int               assetIdx; // index of this asset into multi; can be <0 meaning no corresponding entry

    CAssetConstSP     eqAsset;

    vector<const SVGenSpot*> *multiSpots; // just a local pointer to a global multispot for now...
    vector<const SVGenExpectedSpot*> expSpots; // eq spot state variables

    QMCGenDiffusibleIR*  pQMCGenDiffusibleIR;  // for easy link to underlying currency data
                                         // will be a foreign curve for quantoed equity
    QMCGenDiffusibleIR*  pDomQMCGenDiffusibleIR;  // allows easy access of domestic IRGen

    // TODO: see if this is necessary here
    // Important dates for timeline
    DateTimeArray criticalDates(const DateTime& start,
                                const DateTime& end);
    int getAssetIdx(void) { return assetIdx;}


};
DECLARE(QMCGenDiffusibleEquity);

class SVGenSurvivalDiscFactor;
class SVGenExpectedSurvivalDiscFactor;
class SVGenAggregatedSurvivalDiscFactor;
class SVGenDateOfDefault;

class QMCGenDiffusibleCredit : public IQMCGenDiffusibleAsset {
public:
    /*typedef pair<vector<const SVGenSurvivalDiscFactor*>,
             vector<const SVGenExpectedSurvivalDiscFactor*> > SDFPair;
    */
    struct SDFGenerators
    {
        vector<const SVGenSurvivalDiscFactor*> first;
        vector<const SVGenExpectedSurvivalDiscFactor*> second;
        vector<const SVGenAggregatedSurvivalDiscFactor*> third;
        vector<const SVGenDateOfDefault*> dod;
    };

    QMCGenDiffusibleCredit(): pQMCGenDiffusibleIR(NULL) {}

    /** Currently, all implementations are one factor so do nothing at this level. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Currently, Credit is one dim'l so implement at this level */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(1,beta); }


    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm);

    virtual bool isCredit() { return true; }
    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleCredit(); }
    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit()= 0;
    virtual SRMCreditUtilSP getCreditUtil() = 0;// TODO: see if necessary...

    string                    ccy;    //underlying currency name
    ICDSParSpreadsConstSP     cds;  // used only in createUtil
    SDFGenerators             sdf; // populated in SV::processSVGen

    QMCGenDiffusibleIR* pQMCGenDiffusibleIR;// for easy link to underlying ccy
};
DECLARE(QMCGenDiffusibleCredit);

class SVGenExpectedEnergyFuture;

class QMCGenDiffusibleEnergy : public IQMCGenDiffusibleAsset {
public:
	QMCGenDiffusibleEnergy():
        pQMCGenDiffusibleIR(NULL),
        tierLevel(1),
        pParentQMCGenDiffusibleEnergy(NULL)
        {}

	virtual void createStateVars(StateVarDBase&, 
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm);

    virtual bool isEnrg() { return true; }
	virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleEnergy(); }
	virtual IQMCDiffusibleEnergySP getDiffusibleEnergy() = 0;
	virtual int numFactors() const {return nFactors; }

	string ccy; // underlying currency name
	EnergyFuturesCurveConstSP futureCurve;
	string modelType; // which energy model to use
	int nFactors; // number of factor in the energy model

    vector<const SVGenExpectedEnergyFuture*> expFpGens;

	QMCGenDiffusibleIR * pQMCGenDiffusibleIR;// for easy link to underlying ccy
    QMCGenDiffusibleEnergy * pParentQMCGenDiffusibleEnergy; // for linking to parent energy
    bool isTier1; // true if no parent
    int tierLevel;
    void setTierLevel(); // only call this after the parent relationships are set up for all names
};
DECLARE(QMCGenDiffusibleEnergy);

// for comparing two QMCGenDiffusibleEnergySP based on their tier levels:
class DiffEnrgGenSPComparator : public binary_function<QMCGenDiffusibleEnergySP, QMCGenDiffusibleEnergySP, bool>
{
public:
    DiffEnrgGenSPComparator() {}
    bool operator() (const QMCGenDiffusibleEnergySP& a, const QMCGenDiffusibleEnergySP& b);
};

class SVGenExpectedBasisFwdSpread;

class QMCGenDiffusibleBasisSpread : public IQMCGenDiffusibleAsset {
public:
	QMCGenDiffusibleBasisSpread(): pQMCGenDiffusibleIR(NULL) {}

	virtual void createStateVars(StateVarDBase&, 
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm);
	virtual bool isBS() { return true; }
    /** Currently, all implementations are one factor so do nothing at this level. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Currently, Basis is one dim'l so implement at this level */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(1,beta); }

	virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getDiffusibleBasis(); }
	virtual IQMCDiffusibleBasisIndexSP getDiffusibleBasis() = 0;
	virtual SRMBasisSpreadUtilSP getSPUtil() = 0;

	string ccy; // underlying currency name
	IBasisIndexCurveConstSP basisCurve;
	vector<const SVGenExpectedBasisFwdSpread*> expSPGens;

	QMCGenDiffusibleIR * pQMCGenDiffusibleIR;// for easy link to underlying ccy
};
DECLARE(QMCGenDiffusibleBasisSpread);

// This class is non-abstract as pure jumps are non-asset or model specific
class QMCGenPureJump : public IQMCGenDiffusibleAsset,
                        public virtual ISupportPathConfigSRM
{
public:
    QMCGenPureJump() {} 

    // (2)
    virtual void initializeDiffusibleAssets();

    // (3) no pure jump SV so far
    virtual void createStateVars(StateVarDBase&,
                                 vector<vector<SVPathSP> > &assetPaths,
                                 bool mm){}

    // (4) no dates necessary
    virtual DateTimeArray getAssetDates() { return DateTimeArray(); }
    virtual DateTimeArray getSRMCriticalDates(const MCPathConfigSRM* mcPathConfig,
                                           const DateTime& start,  // likely to be "today"
                                           const DateTime& finish) // the latest requested date
                            { return DateTimeArray(); }

    // (5) no corresponding util
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates) {}

    /** (5.5) Currently, all implementations are one factor so do nothing at this level. */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Currently, pure jump is one dim'l so implement at this level */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(1,beta); }

    // (5.7)
    /** Not implemented yet, dummy implementations at this level.  Future 
    implementations should go to the pathConfig level e.g. SRM or CID */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    // (7)
    virtual void setSRMDiffusibleAsset(
            const DateTime&     today,              // base date of the run
            const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
            const IPastValues*  pastValues,         // historic values
            DependenceSP        dependence);        // factorized correlations


    // might be unnecessary, but yet pretty handy.
    virtual bool    isJump() { return true; }
    virtual bool    isDefined() { return true; }

    virtual int     numFactors() const { return 0; } // only those who are !=1 would override this

    virtual IQMCDiffusibleAssetBaseSP getDiffusibleAsset() { return getPureJumpAsset(); }
    QMCPureJumpsSP   getPureJumpAsset() { return pPureJumpsAsset; }
                                               // need this date as we are not to produce jumps after it

    QMCPureJumpsSP   pPureJumpsAsset;
};
DECLARE(QMCGenPureJump);

class MCARLO_DLL QMCGenDiffusibleCreditCIDJumps : public QMCGenDiffusibleCredit,
                                                  public virtual ISupportPathConfigSRM
{
public:
    // Implement ISupportPathConfigSRM
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Not implemented yet, dummy implementations at this level.  Future 
    implementations should go to the pathConfig level e.g. SRM or CID */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    virtual void setSRMDiffusibleAsset(
        const DateTime&        today,              // base date of the run
        const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
        const IPastValues*     pastValues,         // historic values
        DependenceSP           dependence);        // factorized correlations

    virtual DateTimeArray getSRMCriticalDates(
        const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, 
        const DateTime& finish) {return DateTimeArray();}

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit() { return pAssetCreditCIDJumps; }
    virtual SRMCreditUtilSP getCreditUtil() { return SRMCreditUtilSP(); }
    virtual int     numFactors() const { return 0; } // not using any correlated brownians

    QMCGenDiffusibleCreditCIDJumps() : isFullMC(true) {}

    QMCCreditCIDJumpsSP pAssetCreditCIDJumps; // "shell" of a diffusible credit asset

    QMCGenPureJumpSP    pJumpSource; // generator of the jump source
    bool                isFullMC;    // whether the monte-carlo is full MC
};
DECLARE(QMCGenDiffusibleCreditCIDJumps);

class MCARLO_DLL QMCGenDiffusibleCreditComposite : public QMCGenDiffusibleCredit,
    public virtual ISupportPathConfigSRM
{
public:
    virtual void createSRMUtil(MCPathConfigSRM*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates);

    // (5.7)
    /** Not implemented yet, dummy implementations at this level.  Future 
    implementations should go to the pathConfig level e.g. SRM or CID */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    virtual void setSRMDiffusibleAsset(
            const DateTime&        today,              // base date of the run
            const MCPathConfigSRM* mcPathConfig,       // to get max/min boundaries, etc
            const IPastValues*     pastValues,         // historic values
            DependenceSP           dependence);        // factorized correlations

    virtual DateTimeArray getSRMCriticalDates(
        const MCPathConfigSRM* mcPathConfig,
        const DateTime& start, 
        const DateTime& finish) { return DateTimeArray(); }

    virtual void initializeDiffusibleAssets();

    virtual IQMCDiffusibleCreditSpreadSP getDiffusibleCredit() 
    { return pAssetCreditComposite; }
    virtual SRMCreditUtilSP getCreditUtil() { return SRMCreditUtilSP(); } // no util

    virtual int     numFactors() const { return 0; } // not using any correlated brownians
    virtual bool    isComposite() { return true; }

    QMCCreditCompositeSP pAssetCreditComposite; // "shell" of a diffusible credit asset

    QMCGenDiffusibleCreditComposite() : catastrophicRecoveryRate(0.0), regularRecoveryRate(0.0) {}

    // components of this composite credit
    vector<double>                       weights; 
    vector<QMCGenDiffusibleCreditSP>     components;

    // here for the time being -- but should be encapsulated better
    double catastrophicRecoveryRate;
    double regularRecoveryRate;

    virtual void addComponent(QMCGenDiffusibleCreditSP component, double weight);
};
DECLARE(QMCGenDiffusibleCreditComposite);

class QMCGenDetermRates : public QMCGenDiffusibleIR,
                          public ISupportPathConfigSRM 
{ // includes domestic currency
public:

    QMCGenDetermRates() { nFactors = 0; }
    virtual ~QMCGenDetermRates() {}

    virtual IQMCDiffusibleInterestRateSP getDiffusibleIR() { return pAssetRates; }
    virtual QMCRatesUtilSP               getRatesUtil() { return QMCRatesUtilSP(); }

    virtual void createSRMUtil(MCPathConfig*       mcPathConfig,
                            const DateTime&        today,
                            DateTimeArrayConstSP   simDates) {}

    /** One factor so do nothing */
    virtual void expandFactors(SparseDoubleMatrix& outputSparseCorrMatrix ) const {}

    // (5.6)
    /** Skip deterministic assets */
    virtual vector<double> expandTier2Betas(double beta) const { return vector<double>(0); }

    // (5.7)
    /** Skip deterministic assets */
    virtual vector<double> expandAssetAssetCorrs(IQMCGenDiffusibleAsset* asset2, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToIRCorrs(QMCGenDiffusibleIR* irAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToCRCorrs(QMCGenDiffusibleCredit* crAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToEQCorrs(QMCGenDiffusibleEquity* eqAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToFXCorrs(QMCGenDiffusibleFX* fxAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToBSCorrs(QMCGenDiffusibleBasisSpread* bsAsset, double corr) { return vector<double>(0); }
    virtual vector<double> expandAssetToENCorrs(QMCGenDiffusibleEnergy* enAsset, double corr) { return vector<double>(0); }

    virtual void pushbackICERecord(StateVarDBase& svdbICE, vector<ICE>& vICE) {}
    virtual DateTimeArray getSRMCriticalDates(
        const MCPathConfig* mcPathConfig,
        const DateTime& start, 
        const DateTime& finish)
    {    return DateTimeArray(); }


    virtual void initializeDiffusibleAssets();

    virtual void setSRMDiffusibleAsset(
            const DateTime&     today,              // base date of the run
            const MCPathConfigSRM* mcPathConfig,    // to get max/min boundaries, etc
            const IPastValues*  pastValues,         // historic values
            DependenceSP        dependence);        // factorized correlations

    QMCDetermRatesSP pAssetRates; // "shell" of a diffusible IR asset

};
DECLARE(QMCGenDetermRates);

DRLIB_END_NAMESPACE
#endif // QMCGenDiffusibleAsset

