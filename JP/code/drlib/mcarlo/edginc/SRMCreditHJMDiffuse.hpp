//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditHJMDiffuse.hpp
//
//   Description : Markovian HJM / Ritchken-Sankarasubramanian credit path generation
//
//
//----------------------------------------------------------------------------

#ifndef SRMCreditHJMDIFFUSE_HPP
#define SRMCreditHJMDIFFUSE_HPP

#include "edginc/IQMCDiffusibleAsset.hpp"
#include "edginc/SVGenSurvivalDiscFactor.hpp"
#include "edginc/SVGenExpectedSurvivalDiscFactor.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesHJMDiffuse.hpp"
#include "edginc/SRMCreditDiffuse.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/DECLARE.hpp"
#include <set>

#include "edginc/SRMCreditHJMUtil.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"

DRLIB_BEGIN_NAMESPACE

class IQMCHelperTimeLogic;

class SRMCreditHJM3F;

/** base class for Markovian HJM credit path generation */
class SRMCreditHJMDiffuse : public SRMCreditDiffuse, public DiffusionDriver {
public:

    /** constructor */
    SRMCreditHJMDiffuse(QMCRatesDiffuseSP srmRatesDiffuse);  // NULL must not happen

    /** destructor */
    virtual ~SRMCreditHJMDiffuse();

    /** initialization */
    void setSRMCreditHJMDiffuse(
        int                    _randomIndex,
        const DateTime&        _today,
        SRMCreditHJMUtilSP     _srmCreditHJMUtil,
        double                 _NbSigmasMax,
        double                 _NbSigmasMin,
        const vector<double>&  _corrCRIR,
        double                 _crFxCorr,
        const vector<double>&  _prob); // historic dates

    /** set the credit-rates correlation */
    virtual void setCrIrCorr(const vector<double>& corrCRIR) = 0;

    /** finalizes the timeline, allocates memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);


        /** add all dates that will be needed later, extends base class definition */
    virtual void addAggregatedDates(                
        const DateTimeArray& _sdfDates,
        const DateTimeArray& _esdfRequestedDates,
        const DateTimeArray& _esdfForwardDates);


protected:

    /** compute kT and gT factors for CR as well as for IR */
    virtual void computeKGT(int             dateIdx,
                            double          crKT,
                            double          crGT,
                            vector<double>& irKT,
                            vector<double>& irGT) = 0;

    /* important difference to SRMRatesHJMDiffuse (at least for the time being)
        origSVol, which is populated via SRMRatesUtil::getSpotVols(), is not
        needed here and spotVol lives on SRMCreditHJMUtil */

    /** populates fields required for calculating sigmaL during simulation
        a subset of crdiffuse::CalcSigmaL -- do all possible precalcs here ... */
    void calcSigmaLParams(const DateTimeArray&);

    /** returns credit gfactor */
    double getGFactorCR(FwdIdx i, FwdIdx j);

    /** returns deterministic ratio of log surv probs */
    double getLnProbRatio(FwdIdx i, FwdIdx j);

    const vector<double>&   sigmaR;         // IR vols accross dates, passed in when created in MCPathConfigSRM

    SRMCreditHJMUtilSP srmCreditHJMUtil;
    SRMRatesHJMUtilSP ratesHJMUtil;      // Stored for convenience

    vector<double>  fwdLnProbRatios; // [esdfForwardDates.size()]
    vector<double>  partialIntegral; // [esdfForwardDates.size()]
    vector<double>  zeta;   // [esdfRequestedDates.size()]



    double          qLeft;
    double          qRight;
    double          pivotRatio;         // invariably zero it seems


    // cached variables
    vector<double>  logFwdProbSimple;   // deterministic from discount zero curve. [simDates-1]
    vector<double>  svol;               // contains spot vols [simDates-1]
    vector<double>  MaxEffRateCR;       // the cutoff limit at each ZDate time point
    vector<double>  MinEffRateCR;       // the cutoff limit at each ZDate time point
    vector<double>  lbar;               // lbar at each ZDate time point
    vector<double>  lbarT;              // lbarT at each ZDate time point

    // computed discount factors at probIndexes
    bool            zeroQ;              // true, if qLeft = qRight = 0

    /** returns qPivot */
    double          getQPivot(size_t i) const;

    double          NbSigmasMax;
    double          NbSigmasMin;

    /** trim internal arrays to the actual diffusion */
    void            trimToDiffusion(void);


};

DECLARE(SRMCreditHJMDiffuse);


/** credit Markovian HJM diffusion with one factor CR and one factor IR */
class SRMCreditHJM1F: public SRMCreditHJMDiffuse {
    enum {nfact=1};
public:

    /** constructor */
    SRMCreditHJM1F(QMCRatesDiffuseSP srmRatesDiffuse);

    /** set the credit-rates correlation */
    virtual void setCrIrCorr(const vector<double>& corrCRIR);

    /** finalizes the timeline, allocates memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

    /** generates path across all dates */
    virtual void generatePathAndSurvivalRates(
        IQMCRNGManagerSP rngMgr);

    /** for capturing state of diffusion */
    struct DiffusedState {
        double GAMMA3, THETA0, PHI03, PHI33;
    };

    /** for capturing diffused state on each date on which calculations are performed */
    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGFactors {
        double kFactorCR, gFactorCR;
        double kFactor0IR, gFactor0IR;
    };
    struct KGTFactors {
        double kTfactorCR, gTfactorCR;
        double kTfactor0IR, gTfactor0IR;
    };

   /** accesses the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedSurvivalDiscFactor(FwdIdx measurementDateIdx, FwdIdx futureDateIdx);

    /** accesses the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

protected:

    /** compute kT and gT factors for CR as well as for IR */
    virtual void computeKGT(int             dateIdx,
                            double          crKT,
                            double          crGT,
                            vector<double>& irKT,
                            vector<double>& irGT);
private:

    double                  alpha0;         // for rates bit, backed out from ratesHJMUtil
    double                  rho03;          // correlation between CR and IR

    vector<DiffusedState>   expProb;        /* for computing expected survival probability
                                               stores the values of GAMMA3, THETA0, PHI03, PHI33
                                               length is equal to required nb where expProb has to be stored */
    vector<KGFactors>       kgFactors;      /* for computing k and g factors for CR as well as for IR
                                               length is equal to nb of simulation dates */
    vector<KGTFactors>      kgtFactors;     /* for computing kT and gT factors for CR as well as for IR
                                               length is equal to nb of simulation dates */
};

DECLARE(SRMCreditHJM1F);


/** credit Markovian HJM diffusion with one factor CR and two factor IR */
class SRMCreditHJM2F : public SRMCreditHJMDiffuse {
    enum {nfact=2};
public:

    /** constructor */
    SRMCreditHJM2F(QMCRatesDiffuseSP srmRatesDiffuse);

    /** set the credit-rates correlation */
    virtual void setCrIrCorr(const vector<double>& corrCRIR);

    /** finalizes the timeline, allocates memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

    /** generates path across all dates */
    virtual void generatePathAndSurvivalRates(
        IQMCRNGManagerSP rngMgr);

    /** for capturing state of diffusion */
    struct DiffusedState {
        double GAMMA3, THETA0, THETA1, PHI03, PHI13, PHI33;
    };

    /** for capturing diffused state on each date on which calculations are performed */
    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGFactors {
        double kFactorCR, gFactorCR;
        double kFactor0IR, gFactor0IR, kFactor1IR, gFactor1IR;
    };
    struct KGTFactors {
        double kTfactorCR, gTfactorCR;
        double kTfactor0IR, gTfactor0IR, kTfactor1IR, gTfactor1IR;
    };

    /** accesses the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

    /** accesses the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

protected:

    /** compute kT and gT factors for CR as well as for IR */
    virtual void computeKGT(int             dateIdx,
                            double          crKT,
                            double          crGT,
                            vector<double>& irKT,
                            vector<double>& irGT);
private:

    double                  alpha0;         // for rates bit, backed out from ratesHJMUtil
    double                  alpha1;
    double                  rho03;          // correlation between CR and factor zero IR
    double                  rho13;          // correlation between CR and factor one IR

    vector<DiffusedState>   expProb;        /* for computing expected survival probability
                                               stores the values of GAMMA3, THETA0, THETA1, PHI03, PHI01, PHI33
                                               length is equal to required nb where expProb has to be stored */
    vector<KGFactors>       kgFactors;      /* for computing k and g factors for CR as well as for IR
                                               length is equal to nb of simulation dates */
    vector<KGTFactors>      kgtFactors;     /* for computing kT and gT factors for CR as well as for IR
                                               length is equal to nb of simulation dates */
};

DECLARE(SRMCreditHJM2F);


/** credit Markovian HJM diffusion with one factor CR and three factor IR */
class SRMCreditHJM3F: public SRMCreditHJMDiffuse {
    enum {nfact=3};
public:

    /** constructor */
    SRMCreditHJM3F(QMCRatesDiffuseSP srmRatesDiffuse);

    /** set the credit-rates correlation */
    virtual void setCrIrCorr(const vector<double>& corrCRIR);

    /** finalizes the timeline, allocates memory */
    virtual void finalizePathGenerator(DateTimeArrayConstSP allDates);

    /** generates path across all dates */
    virtual void generatePathAndSurvivalRates(
        IQMCRNGManagerSP rngMgr);

    /** for capturing state of diffusion */
    struct DiffusedState {
        double GAMMA3, THETA0, THETA1, THETA2, PHI03, PHI13, PHI23, PHI33;
    };

    /** for capturing diffused state on each date on which calculations are performed */
    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGFactors {
        double kFactorCR, gFactorCR;
/*        double kFactor0IR, gFactor0IR, kFactor1IR, gFactor1IR, kFactor2IR, gFactor2IR;*/
    };
    struct KGTFactors {
        double kTfactorCR, gTfactorCR;
        double kTfactor0IR, gTfactor0IR, kTfactor1IR, gTfactor1IR, kTfactor2IR, gTfactor2IR;
    };

    /** accesses the expected value ExpSDF(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

    /** accesses the natural log of the expected value ExpSDF(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
    virtual double getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j);

protected:

    /** compute kT and gT factors for CR as well as for IR */
    virtual void computeKGT(int             dateIdx,
                            double          crKT,
                            double          crGT,
                            vector<double>& irKT,
                            vector<double>& irGT);

private:

    double                  alpha0;         // for rates bit, backed out from ratesHJMUtil
    double                  alpha1;
    double                  alpha2;
    double                  rho03;          // correlation between CR and factor zero IR
    double                  rho13;          // correlation between CR and factor one IR
    double                  rho23;          // correlation between CR and factor two IR

    vector<DiffusedState>   expProb;        /* for computing expected survival probability
                                               stores the values of GAMMA3, THETA0, THETA1, THETA2, PHI03, PHI01, PHI23, PHI33
                                               length is equal to required nb where expProb has to be stored */
    vector<KGFactors>       kgFactors;      /* for computing k and g factors for CR as well as for IR
                                               length is equal to nb of simulation dates */
    vector<KGTFactors>      kgtFactors;     /* for computing kT and gT factors for CR as well as for IR
                                               length is equal to nb of simulation dates */

    struct                  State;
    struct                  AuxArgs;

    class Observer;
    friend class Observer; // Argh!FIXME: fix priveleges

    void advanceDiffusion(SRMCreditHJM3F::State&, const SRMCreditHJM3F::AuxArgs&, int i);

    /** trim internal arrays to the actual diffusion; non-virtual! */
    void trimToDiffusion(void);
};

//declare smartPtr versions
DECLARE(SRMCreditHJM3F);

DRLIB_END_NAMESPACE

#endif // SRMCreditHJMDIFFUSE_HPP
