//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MCPathConfigSRMGen.hpp
//
//   Description : A generator of paths using stochastic rates
//                 SRM = stochastic rate model
//
//   Date        : 27 May 2004
//
//
//----------------------------------------------------------------------------
#ifndef MCPathConfigSRMGen_HPP
#define MCPathConfigSRMGen_HPP

#include "edginc/SRMSwaption.hpp"
#include "edginc/SRMCorrelation.hpp"
#include "edginc/FactorCorrelation.hpp"

#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SRMFXDiffuse.hpp"

#include "edginc/MCPathGenerator.hpp"
#include "edginc/MCRandom.hpp"
#include "edginc/MCPathConfigSRMGenSV.hpp"
#include "edginc/Dependence.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/QMCCreditDiffuse.hpp"
#include "edginc/SRMEnergyDiffuse.hpp"
#include "edginc/SRMICE.hpp"
#include "edginc/BetaCorrelation.hpp"

#include "edginc/DependenceGauss.hpp"
#include "edginc/DependenceGaussTerm.hpp"

DRLIB_BEGIN_NAMESPACE

class MCPathConfigSRM;
class SVGenDiscFactor;
class SVGenSurvivalDiscFactor;
class SVGenExpectedDiscFactor;
class SVGenExpectedSurvivalDiscFactor;


// placeholder for srm corr data (for corr term as well as corr mapping)
struct MCARLO_DLL SrmCorrData {
    /** Constructor */
    SrmCorrData(MCPathConfigSRM*   pathConfig);

    CorrelationCommonArray      corrObjArray;   // (n^2-n)/2 correlations
    CorrelationCommonArray      eqeqCorrObjArray;
    BetaCorrelationArray        betaCorrArray;
    CorrelationTermArray        eqeqCorrTermObjArray;
    bool                        strictCorr;
};
typedef refCountPtr<SrmCorrData> SrmCorrDataSP;


class MCARLO_DLL MCPathConfigSRMGen: public virtual MCPathGenerator, // For backward comp
                            public virtual MCStatelessPathGen,
                            public virtual IStateVariableGen::IStateGen,
                            public DependenceMakerGaussSrm::Support {

public:

    //// constructor - most of the work is done at construction
    MCPathConfigSRMGen(
        int                      numSimPaths,
        MCPathConfigSRM*         mcPathConfig,
        const MCPathGeneratorSP& pastPathGenerator,
        const MCProductClient*   prodClient,
        DateTimeArray&           simDates,
        Control*                 control,
        Results*                 results);

    /** MCpathGenerator methods */
    // Deprecated methods
    int NbSimAssets() const {
        throw ModelException("MCPathConfigSRMGen::NbSimAssets",
            "Method is retired for StateVars");
    }

    const double* Path(int /*iAsset*/, int /*iPath*/) const {
        throw ModelException("MCPathConfigSRMGen::Path",
            "Method is retired for StateVars");
    };

    double refLevel(int /*iAsset*/, int /*iPath*/) const {
        throw ModelException("MCPathConfigSRMGen::refLevel",
            "Method is retired for StateVars");
    }

    double maxDriftProduct(int /*iAsset*/) const {
        throw ModelException("MCPathConfigSRMGen::maxDriftProduct",
            "Method is retired for StateVars");
    }

    int begin(int /*iAsset*/) const {
        throw ModelException("MCPathConfigSRMGen::begin",
            "Method is retired for StateVars");
    }

    int end(int /*iAsset*/) const{
        throw ModelException("MCPathConfigSRMGen::end",
            "Method is retired for StateVars");
    }

    //// live methods
    bool hasPast() const {
        return havePast;
    }

    int getPathIndex() const {
        return pathIdx;
    }

    double getPathWeight() const {
        return pathWeight;
    }

    virtual void advance();

    virtual void reset();

    virtual void generatePath(int pathIdx);

    /** Returns true if this path generator is being used to 'simulate'
        the past */
    virtual bool doingPast() const{
        return false;
    }

    /** Returns the state variable corresponding to generator.
        Part of the IStateVariableGen::IStateGen IFace */
    virtual IStateVariableSP create(const IStateVariableGen* svGen){
        return svDBase.find(svGen);
    }

private:    
    
    virtual SparseDoubleMatrixSP getBetaCorrelations(
        const DependenceMakerGaussSrm* dependenceMaker) const;

    virtual double getMaxBetaCorr(
        const DependenceMakerGaussSrm* dependenceMaker) const {
            return betaCorrMax;
        }

    /** general helper method */
    int nbEqEqAssets() const { return sv->numEq;}

    /** helper method to be used by dependence maker */
    virtual DateTimeArray getSimDates() const { return simDates; }

    /** DEPENDENCE MAKER GAUSS SRM && NO CORR MAPPING */
    virtual vector<SparseDoubleMatrixSP> createSparseGaussMatrixArray(
        const DependenceMakerGaussSrm* dependenceMaker) const;

    /** DEPENDENCE MAKER GAUSS SRM && YES CORR MAPPING */
    virtual vector<SparseDoubleMatrixSP> createSparseGaussMatrixArray(
        const DependenceMakerGaussSrm*  dependenceMaker,
        const IntArray&                 fwdCorrIndexArray,
        const DateTimeArray&            fwdCorrDatesArray) const;

    /** DEPENDENCE MAKER GAUSS TERM SRM && NO/YES CORR MAPPING */
    virtual DoubleMatrix getFwdVarAtDates() const;

    virtual void getCorrelationData(CorrelationCommonArray& corrObjArray,
                                    CorrelationTermArray&   corrTermObjArray) const;

    virtual vector<SparseDoubleMatrixSP> createSparseGaussTermMatrixArray(
        const DependenceMakerGaussSrm*  dependenceMaker,
        DoubleMatrixArraySP                 fwdCorrelations,
        const IntArray&                     fwdCorrIndexArray,
        const DateTimeArray&                fwdCorrDatesArray) const;

    /** SrmEqVolData needed for CorrTS and/or CorrMapping */
    void populateSrmEqVolData( SRMRatesHJMUtil* baseRatesUtil );

    /** creation of master sparse collection - NO extention for multiple IR factors !!!
        done only once in order to pull out correlations from market
        output request for initial corr massaging is done here
        => CALLED ONLY ONCE, since highly inefficient */
    SparseCollectionSP createMasterSparseCollection() const;

    /** check pos definiteness, potentially modify and extend for multiplie IR factors
        return (potentially) modified & extended sparse matrix */
    SparseDoubleMatrixSP checkAndExtendSparseCorrMatrix(
                SparseDoubleMatrixSP    inputSparseCorrMatrix,
                double                  eigenValueFloor,
                double                  maxSqError,
                DateTime                startDate,
                DateTime                endDate) const ;

    /** retrieve raw EqEq corrmatrix from corr obj array via string search */
    CDoubleMatrixSP getOrderedEqEqCorrelations() const;

    /** use raw const EqEq or time-dep fwd EqEq corrs and do corr mapping */
    void computeAdjEqEqCorrs(DoubleMatrixArraySP&   fwdCorrelations,
                             const IntArray&        fwdCorrIndexArray,
                             const DateTimeArray&   fwdCorrDatesArray) const;    

    /** for CorrTS as well as CorrMap
        create sparse corr matrices using original IR corrs and modified EQ corrs */
    vector<SparseDoubleMatrixSP> replaceOrderedEqEqCorrelations(
                CDoubleMatrixSP     origEqEqCorrelations,
                DoubleMatrixArraySP adjEqEqCorrelations,
                double              eigenValueFloor,
                double              maxSqError) const;

    //SparseDoubleMatrixSP retrieveBetaCorrelations() const;

    void initializeRandGen(int numSimPaths, MCPathConfigSRM* mcPathConfig);
    void enforceFirstMomentMatching(int numSimPaths);

    // placeholder for srm vol data (for corr term as well as corr mapping)
    struct MCARLO_DLL SrmEqVolData {
        /** constructor */
        SrmEqVolData(const DoubleArrayArray&  assetSpotVols,
                     const DoubleArrayArray&  fwdAssetCompVars,
                     const DoubleArrayArray&  totalAssetCompVars,
                     const DoubleArrayArray&  irAssetIntegral,
                     const DoubleArray&       irIrIntegral);

        DoubleArrayArray    assetSpotVols;      // sigma(t)
        DoubleArrayArray    fwdAssetCompVars; // int_0^t compSigma^2(u) du
        DoubleArrayArray    totalAssetCompVars;
        DoubleArrayArray    irAssetIntegral;    // 2*rho*int_0^t sigma(u) sigma_P(u,t) du
        DoubleArray         irIrIntegral;       // int_0^t sigma_P^2(u,t) du
    };

    typedef refCountPtr<SrmEqVolData> SrmEqVolDataSP;

    /** fields of MCPathConfigSRM::Gen */
    Control*        control;
    Results*        results;
    int             pathIdx;
    vector<IQMCDiffusibleAssetBase*> allAssets;

    MCRandomSP      randomGen;  // Manages random number generation
    IQMCRNGManagerSP rngMgr;    // provides correlated numbers + auxilary ones

    int             iceCalibFreq; // how often to invoke ICE calib
    int             iceCalibIdx;  // next idx to invoke ICE calib
    vector<ICE>     icePerCcy;   // per currency
    StateVarDBase   svDBase;        //!< Collection of Generators + statevars
    vector<IAdvanceableStateVariable*> advanceSvSet; // for stateless advance
    bool            havePast;
    vector<vector<SVPathSP> > assetPaths; // temporary holding area
    DependenceSP    dependence;
    SVSP            sv;
    SrmCorrDataSP   corrData;
    SrmEqVolDataSP  thisSrmEqVolData;
    DateTimeArray&  simDates;
    double          betaCorrMax;

    vector<QMCStrataConstSP> strataByPath;
    vector<double>           fractionByPath;
    double          pathWeight;
};

DRLIB_END_NAMESPACE
#endif

