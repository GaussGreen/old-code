//----------------------------------------------------------------------------
//
//   Group       : Cross Asset QR
//
//   Filename    : SRMBasisSpreadHJMDiffuse.hpp
//
//   Description : Rates SP model diffusion
//
//   Date        :
//
//----------------------------------------------------------------------------

#ifndef SRMBasisSpreadHJMDiffuse_HPP
#define SRMBasisSpreadHJMDiffuse_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/QMCBasisSpreadDiffuse.hpp"
#include "edginc/SRMBasisSpreadHJMUtil.hpp"
#include "edginc/SRMRatesHJMDiffuse.hpp"
#include "edginc/BasisIndexCurve.hpp"


#include <set>

DRLIB_BEGIN_NAMESPACE

//class SRMBasisSpreadUtil;
//class SRMIRUtil;

//// implementation of Basis Spread path generator
class SRMBasisSpreadHJMDiffuse : public QMCBasisSpreadDiffuse
{
public:
    virtual ~SRMBasisSpreadHJMDiffuse(){}
    //// constructor
    SRMBasisSpreadHJMDiffuse(
        SRMRatesHJMDiffuseSP srmIRBase, // NULL must not happen
        IBasisIndexCurveConstSP _basisIndexCurve);  // unfortunately we need curve before "setSP" command

    void setSRMBasisSpreadHJMDiffuse(
        int                    randomIndex,
        const DateTime&        today,
        SRMBasisSpreadUtilSP   srmSPUtil,
        double                 NbSigmasMax,
        double                 NbSigmasMin,
        const vector<double>&  corrSPIR,
        double                 spFxCorr );

        virtual void setSpIrCorr(
            const vector<double>& corrSPIR,
            const vector<double>& alphaIR) = 0;

   /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    /** getting the simulation start date */
    virtual DateTime getBaseDate() { return today; }


    virtual void addAggregatedDates(
                            const DateTimeArray& spot,
                            const DateTimeArray& measurementDates,
                            const DateTimeArray& resetDates);

    virtual IQMCDiffusibleInterestRateSP getUnderlyingIRAsset() { return irAsset; }

    /** Returns pointer to internal storage for MaxDiff(usion)Date / MaxCurveMaturity */
    virtual QMCHelperBoundedDiffusion*    getDiffusionBound() { return &diffBound;}

protected:

    /* important difference to SRMIR::Base (at least for the time being)
        origSVol, which is populated via SRMIRUtil::getSpotVols(), is not
        needed here and spotVol lives on SRMCRUtil */

    // does a collar
    static double COLLAR(double amt, double cap, double flr) {
        return Maths::max(Maths::min(amt,cap),flr);
    }

    /** populates fields required for calculating sigmaL during simulation
        is called in constructor of SRMSP::SP1F, SRMSP::SP2F or SRMSP::SP3F
        a subset of crdiffuse::CalcSigmaL -- do all possible precalcs here ... */
    void calcSigmaSPParams();

    size_t getNumSimDates(void) {return diffusionDates.size();}
    size_t getNbRequestedDates(void) {return getForwardDates().size(); /*requestedDates.size();*/
    }


    IBasisIndexCurveConstSP basisIndexCurve;
    SRMRatesHJMDiffuseSP    irAsset;
    const vector<double>&   sigmaR;         // IR vols accross dates, passed in when created in MCPathConfigSRM
    const vector<double>*   sigmaFX;        // FX vols accross dates, passed in when created in MCPathConfigSRM

    int             randomIndex;        // index into random numbers

    double          qSpread;            // normal/lognormal switch

    double          spFxCorr;           // correlation(CR,FX), always one single nb

    //// cached variables
    // expProbIndexes -> expSpreadIndexes
    vector<int>     expSpreadIndexes;     // indexes for when we compute expected survival probabilities [Nesdf]
    vector<double>  spVol;               // contains spot vols [simDates-1]
    vector<double>  MaxEffRateSP;       // the cutoff limit at each ZDate time point
    vector<double>  MinEffRateSP;       // the cutoff limit at each ZDate time point
    vector<double>  sqrtYearFrac;       // sqrt of yearFrac between sim dates
    vector<double>  spBar;               // deterministic from discount zero curve. [simDates-1]

    // computed discount factors at probIndexes
    bool            zeroQ;          // true, if qLeft = qRight = 0
    SRMRatesHJMUtilSP srmIRUtil;      // Stored for convenience
    SRMBasisSpreadUtilSP srmSPUtil;      // will be modified at the finalize() time

    int             todayIdx;

    vector<int>     fwdIdx2RequestedIdx; // projection of ForwardDates to RequestedDates
    vector<int>     fwdIdx2IRedfIdx; // projection onto EDFs
    vector<int>     fwdPlusIdx2IRedfIdx;  // projection of resetDates to EDFs
    // DateTimeArray   forwardDates;   // requested dates and reset dates
    DateTimeArray   diffusionDates;
    double          getFwdSpread(FwdIdx i) { return fwdSpread[i]; }
    vector<double>  partialHFactor; // [forwardDates.size()]


private:

    double         NbSigmasMax;
    double         NbSigmasMin;

    DateTime        today;

    vector<double>  fwdSpread; // [forwardDates.size()]

    QMCHelperBoundedDiffusion    diffBound;
};

DECLARE(SRMBasisSpreadHJMDiffuse);

/***********************************************************************************************************/
/***********************************************************************************************************/
/***********************************************************************************************************/

// one factor SP and one factor IR
class SRMBasisSpreadHJM1F: public SRMBasisSpreadHJMDiffuse {
    enum {nfact=1};  // to do: How is this used?
public:
    //// constructor
    SRMBasisSpreadHJM1F(
        SRMRatesHJMDiffuseSP srmIRBase,
        IBasisIndexCurveConstSP _basisIndexCurve)
            : SRMBasisSpreadHJMDiffuse(srmIRBase, _basisIndexCurve),alpha0(-999.),rho0(-999.) {}

    virtual void setSpIrCorr(const vector<double>& corrSPIR, const vector<double>& alphaIR)
    {
        alpha0=alphaIR[0];
        rho0=corrSPIR[0];
    }

    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates ...
        essentially crdiffuse::DiffuseCR_1F plus a part of crdiffuse::CalcSigmaLCR */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    // for capturing state of diffusion
    struct DiffusedState{
        double OMEGA, CHI0;
    };

    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGHFactors{
        double hFactor;
        double kFactor0IR, gFactor0IR;
        double gFactor0TIR; // for lognormal models -- g(T_i, T_i+1Libor)
    };

   /** Accessing the expected value ExpSpread(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedBasisFwdSpread(
        FwdIdx measurementDateIdx,  // Idx of measDate
        FwdIdx resetDatesIdx );


private:
    double                  alpha0;         // for rates bit, backed out from srmIRUtil
    double                  rho0;          // correlation between SP and IR

    vector<DiffusedState>   expSP;        /* for computing expected survival probability
                                               stores the values of OMEGA and PHI0
                                               length is equal to required nb where expProb has to be stored */
    vector<KGHFactors>       kghFactors;      /* for computing k and g factors for SP as well as for IR
                                               length is equal to nb of simulation dates */

};


// one factor SP and two-factor underlying IR
class SRMBasisSpreadHJM2F: public SRMBasisSpreadHJMDiffuse {
    enum {nfact=2};  // to do: How is this used?
public:
    //// constructor
    SRMBasisSpreadHJM2F(
        SRMRatesHJMDiffuseSP srmIRBase,
                IBasisIndexCurveConstSP _basisIndexCurve)
            : SRMBasisSpreadHJMDiffuse(srmIRBase, _basisIndexCurve),
            alpha0(-999.),alpha1(-999.),rho0(-999.),rho1(-999.) {}

    virtual void setSpIrCorr(const vector<double>& corrSPIR, const vector<double>& alphaIR)
    {
        alpha0=alphaIR[0];
        alpha1=alphaIR[1];
        rho0=corrSPIR[0];
        rho1=corrSPIR[1];
    }

    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates ...
        essentially spdiffuse::DiffuseSP_2F plus a part of spdiffuse::CalcSigmaLCR */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    // for capturing state of diffusion
    struct DiffusedState{
        double OMEGA, CHI0, CHI1;
    };

    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGHFactors{
        double hFactor;
        double kFactor0IR, gFactor0IR;
        double kFactor1IR, gFactor1IR;
        double gFactor0TIR, gFactor1TIR; // for lognormal models -- g(T_i, T_i+1Libor)
    };

   /** Accessing the expected value ExpSpread(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedBasisFwdSpread(
        FwdIdx measurementDateIdx,  // Idx of measDate
        FwdIdx resetDatesIdx );


private:
    double                  alpha0;         // for rates bit, backed out from srmIRUtil
    double                  alpha1;         // for rates bit, backed out from srmIRUtil
    double                  rho0;           // correlation between SP and factor-zero of IR
    double                  rho1;           // correlation between SP and factor-one of IR

    vector<DiffusedState>   expSP;        /* for computing expected survival probability
                                               stores the values of OMEGA and PHI0, PHI1
                                               length is equal to required nb where expProb has to be stored */
    vector<KGHFactors>       kghFactors;      /* for computing k and g factors for SP as well as for IR
                                               length is equal to nb of simulation dates */

};





// one factor SP and three-factor underlying IR
class SRMBasisSpreadHJM3F: public SRMBasisSpreadHJMDiffuse {
    enum {nfact=3};  // to do: How is this used?
public:
    //// constructor
    SRMBasisSpreadHJM3F(
        SRMRatesHJMDiffuseSP srmIRBase,
        IBasisIndexCurveConstSP _basisIndexCurve)
            : SRMBasisSpreadHJMDiffuse(srmIRBase, _basisIndexCurve),
            alpha0(-999.),alpha1(-999.),alpha2(-999.),
            rho0(-999.),rho1(-999.),rho2(-999.) {}

    virtual void setSpIrCorr(const vector<double>& corrSPIR, const vector<double>& alphaIR)
    {
        alpha0=alphaIR[0];
        alpha1=alphaIR[1];
        alpha2=alphaIR[2];
        rho0=corrSPIR[0];
        rho1=corrSPIR[1];
        rho2=corrSPIR[2];
    }

    virtual void finalize(DateTimeArrayConstSP allDates);

    /** generate path across all dates ...
        essentially spdiffuse::DiffuseSP_3F plus a part of spdiffuse::CalcSigmaLCR */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    // for capturing state of diffusion
    struct DiffusedState{
        double OMEGA, CHI0, CHI1, CHI2;
    };

    typedef vector<DiffusedState>::const_iterator DiffusedStateIter;

    struct KGHFactors{
        double hFactor;
        double kFactor0IR, gFactor0IR;
        double kFactor1IR, gFactor1IR;
        double kFactor2IR, gFactor2IR;
        double gFactor0TIR, gFactor1TIR, gFactor2TIR; // for lognormal models -- g(T_i, T_i+1Libor)
    };

   /** Accessing the expected value ExpSpread(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
    virtual double getExpectedBasisFwdSpread(
        FwdIdx measurementDateIdx,  // Idx of measDate
        FwdIdx resetDatesIdx );


private:
    double                  alpha0;         // for rates bit, backed out from srmIRUtil
    double                  alpha1;         // for rates bit, backed out from srmIRUtil
    double                  alpha2;         // for rates bit, backed out from srmIRUtil
    double                  rho0;           // correlation between SP and factor-zero of IR
    double                  rho1;           // correlation between SP and factor-one of IR
    double                  rho2;           // correlation between SP and factor-two of IR

    vector<DiffusedState>   expSP;        /* for computing expected survival probability
                                               stores the values of OMEGA and CHI0, CHI1, CHI2
                                               length is equal to required nb where expProb has to be stored */
    vector<KGHFactors>       kghFactors;      /* for computing k and g factors for SP as well as for IR
                                               length is equal to nb of simulation dates */

};


DRLIB_END_NAMESPACE

#endif // QLIB_SRM_SP_MODEL_HPP
