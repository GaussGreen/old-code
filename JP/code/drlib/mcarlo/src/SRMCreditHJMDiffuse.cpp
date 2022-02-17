//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditHJMDiffuse.cpp
//
//   Description : Markovian HJM / Ritchken-Sankarasubramanian path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditHJMDiffuse.hpp"
#include "edginc/SRMCreditUtil.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMUtil.hpp"

#include "edginc/SVQmcImplemented.hpp"
#include "edginc/MemoryProfiler.hpp"
#include "edginc/Format.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/************************************************************************/
/* Implementation of SRMCreditHJMDiffuse class                          */
/************************************************************************/

/** constructor */
SRMCreditHJMDiffuse::SRMCreditHJMDiffuse(QMCRatesDiffuseSP srmRatesDiffuse) :
    SRMCreditDiffuse(srmRatesDiffuse),
    sigmaR(srmRatesDiffuse->getSigmaR()),
    srmCreditHJMUtil(NULL),
    ratesHJMUtil(0),
    qLeft(-999.),
    qRight(-999.),
    pivotRatio(-999.),
    zeroQ(false),
    NbSigmasMax(-999.),
    NbSigmasMin(-999.)
{}

/** destructor */
SRMCreditHJMDiffuse::~SRMCreditHJMDiffuse() {
}

/** add all dates that will be needed later, extends base class definition */
void SRMCreditHJMDiffuse::addAggregatedDates(                
    const DateTimeArray& _sdfDates,
    const DateTimeArray& _esdfRequestedDates,
    const DateTimeArray& _esdfForwardDates)
{
    QMCCreditDiffuse::addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
    // essential for HJM models as they access diffused data from HJM directly
    irAsset->addAggregatedDates(_sdfDates, _esdfRequestedDates, _esdfForwardDates);
}


/** initialization */
void SRMCreditHJMDiffuse::setSRMCreditHJMDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMCreditHJMUtilSP     _srmCreditHJMUtil,
    double                 _NbSigmasMax,
    double                 _NbSigmasMin,
    const vector<double>&  _corrCRIR,
    const double           _crFxCorr,
    const vector<double>&  _prob)    // for historic dates
{
    static const string method("SRMCreditHJMDiffuse::setSRMCreditHJMDiffuse");

    srmCreditHJMUtil = _srmCreditHJMUtil;
    qLeft=srmCreditHJMUtil->getQLeft();
    qRight=srmCreditHJMUtil->getQRight();
    randomIndex=_randomIndex;
    crFxCorr=_crFxCorr;
    NbSigmasMax=_NbSigmasMax;
    NbSigmasMin=_NbSigmasMin;
    today=_today;
    probStart=-1; // FIXME probStart : see QMCRatesDiffuse how to calc it
    prob=_prob;
    zeroQ = Maths::isZero(qLeft) && Maths::isZero(qRight);

    setRecoveryRate(srmCreditHJMUtil->getCdsCurve()->getRecovery());

    // must be available now
    SRMRatesHJMDiffuse* hjmAsset = dynamic_cast<SRMRatesHJMDiffuse*>(irAsset.get());
    if (hjmAsset == 0) {
        throw ModelException(method, "Failed to cast IR model to HJM");
    }
    ratesHJMUtil = hjmAsset->getRatesHJMUtil();
    setCrIrCorr(_corrCRIR);
}

/** returns qPivot */
double SRMCreditHJMDiffuse::getQPivot(size_t i) const {
    return lbar[i] * pivotRatio;
}

/** trim internal arrays to the actual diffusion */
void SRMCreditHJMDiffuse::trimToDiffusion(void) {
    svol.resize(lastDiffusionIdx+1);
    //    QPivot.resize(lastDiffusionIdx+1);
    if (zeroQ)
    {
        MaxEffRateCR.resize(lastDiffusionIdx+1);
        MinEffRateCR.resize(lastDiffusionIdx+1);
    }
    //    sqrtYearFrac.resize(lastDiffusionIdx);
    logFwdProbSimple.resize(lastDiffusionIdx+1);
}


/** finalize the timelines, allocate necessary memory */
/* By now SRMCreditUtil is already initialized */
void SRMCreditHJMDiffuse::finalizePathGenerator(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMCreditHJMDiffuse::finalize");

    const DateTimeArray& simDates = srmCreditHJMUtil->getSimDates();
    const int numSimDates = simDates.size();

    ASSERT(numSimDates == SRMUtil::getNumSimDates(today, *allDatesSP));

    calcFirstAndLastDiffusionIdx(simDates, *allDatesSP);
    processAllDates(allDatesSP);
    calcRemapToIRAssetIdx(getForwardForwardDates());


    if (isWholeTimelineSurvivalRequested())
        wholeTimelineLogSurvProb.resize(lastDiffusionIdx-todayIdx+1,0.0);
    svol = srmCreditHJMUtil->getSpotVols();

    srmCreditHJMUtil->computeLogFwdProbSimple(logFwdProbSimple);

    if (srmCreditHJMUtil->getMomentMatchingFlag())
    {
        srmCreditHJMUtil->computeLogFwdProbSimple(getSpotDates() /*sdfRequestedDates*/, originalProbs);
        srmCreditHJMUtil->computeLogFwdProbSimple(getForwardForwardDates() /*esdfForwardDates*/, originalExpProbs); // is reductant, see fwdLnProbRatios
    }

    srmCreditHJMUtil->calcEffRateLimit(NbSigmasMax, NbSigmasMin,
        MaxEffRateCR, MinEffRateCR);

    //    sqrtYearFrac = SRMUtil::computeSqrtYearFrac(simDates);

    srmCreditHJMUtil->computeLogFwdProbSimple(getForwardForwardDates(), fwdLnProbRatios);
    srmCreditHJMUtil->populatePartialIntCR(today, getForwardForwardDates(), partialIntegral);
    srmCreditHJMUtil->populatePartialZeta(today, getForwardDates(), zeta);

    trimToDiffusion();
}

/** populates fields required for calculating sigmaL during simulation
    a subset of crdiffuse::CalcSigmaL -- do all possible precalcs here ... */
void SRMCreditHJMDiffuse::calcSigmaLParams(const DateTimeArray& dates) // dates are srmCreditHJMUtil->getSimDates()
{
    static const string method("SRMCreditHJMDiffuse::calcSigmaLParams");
    pivotRatio = 1.0 / (1.0 + srmCreditHJMUtil->getFwdShift());
    // loop across sim dates and not across extended timeline
//     const DateTimeArray& dates = srmCreditHJMUtil->getSimDates();

    //QPivot.resize(dates.size()-1);
    lbar.resize(dates.size()-1);
    lbarT.resize(dates.size()-1);
    svol.resize(dates.size()-1);

    for (int i = 0; i < dates.size()-1; i++){
        const DateTime& CurrDate = dates[i];        // current sim date
        const DateTime& LastDate = dates.back();    // last sim date

        // Vol-driver rate's maturity  T = DateTime(currDate + offset, time)
        DateTime T0 = CurrDate.rollDate(driver > 0 ? driver : LastDate.daysDiff(CurrDate)/2);
        if (T0 < CurrDate){
            throw ModelException(method, "Internal error");
        }
        this->lbar[i] = srmCreditHJMUtil->lBar(i);                  // save lbar across dates
        this->lbarT[i] = srmCreditHJMUtil->lBar(T0);     // save lbarT across dates
  //      this->QPivot[i] = lbar[i] * pivotRatio;    // save QPivot across dates

        if (Maths::isZero(this->lbarT[i])){
            throw ModelException(method, "Zero forward intensity at " + T0.toString());
        }

        /* adjust for the effect of the forward shift. This is required because the value of SVol(.,.)
            assumes it is applied to an effR which Vladimirises to rbar. This adj is not stored in
            svol because the unadjusted vol is used elsewhere (e.g. in FX spot vol bootstrapping)  */
        if (pivotRatio <= 1.0) {
            svol[i] /= (qRight + (1-qRight)*pivotRatio);
        } else {
            svol[i] /= (qLeft + (1-qLeft)*pivotRatio);
        }
        /* calculate the vol factors and constants */
        if(!zeroQ) {
            vector<double> irKT(ratesHJMUtil->numFactors());
            vector<double> irGT(ratesHJMUtil->numFactors());
            double crKT = srmCreditHJMUtil->kFactor(CurrDate, T0);
            double crGT = srmCreditHJMUtil->gFactor(CurrDate, T0);

            ratesHJMUtil->kFactor(CurrDate, T0, irKT);
            ratesHJMUtil->gFactor(CurrDate, T0, irGT);
            computeKGT(i, crKT, crGT, irKT, irGT);
        }
    }
}

/** returns credit gfactor */
double SRMCreditHJMDiffuse::getGFactorCR(FwdIdx i, FwdIdx j) {
    SpotIdx iEDF = getTimeLogic()->getReqEDFIdx(i);
#if !defined (_MSC_VER) || (_MSC_VER >= 1300)
    ASSERT(iEDF != SpotIdx(SpotIdx::npos));
#else
    ASSERT(iEDF != -1);
#endif
    return zeta[iEDF]*(partialIntegral[j] - partialIntegral[i]);
}

/** returns deterministic ratio of log surv probs */
double SRMCreditHJMDiffuse::getLnProbRatio(FwdIdx i, FwdIdx j) {
    return fwdLnProbRatios[j] - fwdLnProbRatios[i];
}

/************************************************************************/
/* Implementation of SRMCreditHJM1F class                               */
/************************************************************************/

/** constructor */
SRMCreditHJM1F::SRMCreditHJM1F(QMCRatesDiffuseSP srmRatesDiffuse) :
    SRMCreditHJMDiffuse(srmRatesDiffuse),
    alpha0(-999.),
    rho03(-999.)
{}

/** set the credit-rates correlation */
void SRMCreditHJM1F::setCrIrCorr(const vector<double>& corrCRIR) {
    rho03=corrCRIR[0];
}

/** finalize the timelines, allocate necessary memory */
void SRMCreditHJM1F::finalizePathGenerator(DateTimeArrayConstSP allDatesSP) {
    static const string method("SRMCreditHJM1F::finalizePathGenerator");
    try {

        alpha0=ratesHJMUtil->getAlpha(0);

        SRMCreditHJMDiffuse::finalizePathGenerator(allDatesSP); // will notify srmCreditHJMUtil about new timeline
        expProb.resize(getNumESDFDates());

        // initialization of quite a few parameters
        calcSigmaLParams(srmCreditHJMUtil->getSimDates());

        const DateTimeArray& dates = srmCreditHJMUtil->getSimDates(); // same as diffusionDates ?
        const int numSimDates = dates.size();

        kgFactors.resize(dates.size()-1);

        for (int i = 0; i < numSimDates - 1; i++) {
            kgFactors[i].kFactorCR = srmCreditHJMUtil->kFactor(dates[i], dates[i+1]);
            kgFactors[i].gFactorCR = srmCreditHJMUtil->gFactor(dates[i], dates[i+1]);
            kgFactors[i].kFactor0IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kgFactors[i].gFactor0IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
        }
        srmCreditHJMUtil = SRMCreditHJMUtilSP(0);
        ratesHJMUtil = SRMRatesHJMUtilSP(0);
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** generate path across all dates */
void SRMCreditHJM1F::generatePathAndSurvivalRates(IQMCRNGManagerSP rngMgr) 
{
    SRMRatesHJM1F* srmRates = static_cast<SRMRatesHJM1F *> (getUnderlyingIRAsset().get());

    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = logFwdProbSimple.size();

    int probDatePos = probStart; // position in probIndexes/prob.
    /* set up position in expProbIndexes/expSvob. If expProbIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of probStart and subtracting it in expDatePos formula
    int expProbDatePos = expProbIndexes.front() == 0? 1: 0;
    int probDateIdx = probIndexes[probDatePos]; // date when we save survival probability
    int expProbDateIdx = expProbIndexes[expProbDatePos]; // date when we save expected survival probability
    int stopIdx = Maths::min(probDateIdx, expProbDateIdx); // when to do something

    double LnQ = 0.0;       // exp(int lambda(u,u) du), ie 1 / survival probability
    double GAMMA3 = 0.0;    // mixed term int [sigmaL(sigmaL* dt + sigmaF* dt + dWL)]
    double THETA0 = 0.0;    // mixed term int [sigmaF * sigmaL* dt]
    double PHI03 = 0.0;     // mixed term int [sigmaF * sigmaL dt]
    double PHI33 = 0.0;     // pure credit term int [sigmaL * sigmaL dt]

    // loop over all dates
    for (int i = 0; i < lastDiffusionIdx-todayIdx /*numDates*/; /* increment in body */){
        assert(i < numDates);
        double W3 = randoms[i];

        /* Retrieve precomputed k-factors & g-factors as well as kT-factors & gT-factors
           (precomputation takes place in constructor of SRMCreditHJM1F) */
        double k0 = kgFactors[i].kFactor0IR;
        double g0 = kgFactors[i].gFactor0IR;
        double k3 = kgFactors[i].kFactorCR;
        double g3 = kgFactors[i].gFactorCR;

        /* Step1a: update sigmaL (see FwdCrSprd_s and CalcSigmaLCR in crdiffuse) */
        double QPivot = getQPivot(i); // this->QPivot[i]; // (this->QPivot[i] is already lbar[i] * pivotRatio);
        double effL = 0.0;
        if(zeroQ) {
            effL = Maths::collar(QPivot, MaxEffRateCR[i], MinEffRateCR[i]);
        } else { // see FwdCrSprd_s in crdiffuse -- one factor case
            double k0T = kgtFactors[i].kTfactor0IR;
            double g0T = kgtFactors[i].gTfactor0IR;
            double k3T = kgtFactors[i].kTfactorCR;
            double g3T = kgtFactors[i].gTfactorCR;
            double aRate = lbarT[i] + k3T*GAMMA3 + k3T*g3T*PHI33
                            + k0T*THETA0 + (k0T*g3T+k3T*g0T)*PHI03;
            double lbarRatio = lbar[i] / lbarT[i];
            aRate *= lbarRatio;

            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);
            effL = (aRate >= QPivot) ?
                    (QPivot_x_OneMinusQR + (aRate * qRight)):
                    (QPivot_x_OneMinusQL + (aRate * qLeft));
                    effL = Maths::collar(effL, MaxEffRateCR[i], MinEffRateCR[i]);
        }
        double sigmaL = effL * svol[i];

        /* Step 1b: update sigmaR -- saved when IR simulation bit was done ... */
        double sigmaR0 = sigmaR[i] * alpha0;

        /* Step 1c: compute modified k-factors and modified g-factors */
        double RootDelT = srmRates->getSqrtYearFrac(i);

        double sg3 = sigmaL*g3*RootDelT;
        double sk3 = sigmaL*k3*RootDelT;

        double sg0 = sigmaR0*g0*RootDelT;
        double sk0 = sigmaR0*k0*RootDelT;

        /* Step 2: update LnQ */
        double lnFwdProb = logFwdProbSimple[i]; // ratio S(0,T) / S(0,t) = deterministic term of S(t,T)
        double lnProbTerm = -g3*(g3*PHI33/2.0 + GAMMA3) // stoch term of S(t,T)
                    - g0*(g3*PHI03 + THETA0);
        LnQ += -(lnFwdProb+lnProbTerm) + sg3*(sg3/2.0 + sg0*rho03 + W3);

        /* Step 3: update GAMMA */
        GAMMA3 = k3*(g3*PHI33 + g0*PHI03 + GAMMA3)
                    + sk3*(sg3 + sg0*rho03 + W3);

        /* Step 4: update THETA0 */
        THETA0 = k0*(g3*PHI03 + THETA0) + sk0*(sg3*rho03);

        /* Step 5: update PHI33 and PHI03 */
        PHI33 = k3*k3*PHI33 + sk3*sk3;
        PHI03 = k0*k3*PHI03 + sk0*sk3*rho03;

        /* Step 6: add cups component, if necessary ... this is related to FX */
        if (sigmaFX) {
            double sX = (*sigmaFX)[i]*RootDelT;
            GAMMA3 -= (sk3*sX*crFxCorr);
            LnQ -= (sg3*sX*crFxCorr);
        }

        i++; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i] = LnQ;

        if (i + todayIdx == stopIdx) {
            // hit an "event", ie record sthg
            if (stopIdx == probDateIdx) { // save survival probability
                prob[probDatePos] = LnQ;
                probDatePos++;
                probDateIdx = probIndexes[probDatePos];
            }
            if (stopIdx == expProbDateIdx) {
                expProb[expProbDatePos].GAMMA3 = GAMMA3;
                expProb[expProbDatePos].THETA0 = THETA0;
                expProb[expProbDatePos].PHI33 = PHI33;
                expProb[expProbDatePos].PHI03 = PHI03;
                expProbDatePos++;
                expProbDateIdx = expProbIndexes[expProbDatePos];
            }
            stopIdx = Maths::min(probDateIdx, expProbDateIdx); // refresh
        }
    }
    assert(prob.size() == probDatePos);
    assert(expProb.size() == expProbDatePos);
    // then do exp on prob vector
    for (size_t j = probStart; j < prob.size(); j++){
        prob[j] = exp(-prob[j]);
    }
}

/** accesses the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMCreditHJM1F::getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    // to avoid one virtual call specify the needed class explicitly
    return exp(SRMCreditHJM1F::getLnExpectedSurvivalDiscFactor(i, j));
}

/** accesses the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double SRMCreditHJM1F::getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    double g0 = static_cast<SRMRatesHJM1F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j));
    double g3 = getGFactorCR(i, j);

    const int iEDF = getTimeLogic()->getReqEDFIdx(i);
    const DiffusedState& state = expProb[iEDF];

    double result = getLnProbRatio(i, j) - g3*state.GAMMA3
                                            - 0.5*g3*g3*state.PHI33
                                            - g0*state.THETA0
                                            - g0*g3*state.PHI03;
    return result;
}

/** compute kT and gT factors for CR as well as for IR */
void SRMCreditHJM1F::computeKGT(
                        int             dateIdx,
                        double          crKT,
                        double          crGT,
                        vector<double>& irKT,
                        vector<double>& irGT)
{
    if ((int)kgtFactors.size() <= dateIdx)
        kgtFactors.resize(dateIdx+1);
    kgtFactors[dateIdx].kTfactorCR = crKT;
    kgtFactors[dateIdx].gTfactorCR = crGT;
    kgtFactors[dateIdx].kTfactor0IR = irKT[0];
    kgtFactors[dateIdx].gTfactor0IR = irGT[0];
}

/************************************************************************/
/* Implementation of SRMCreditHJM2F class                               */
/************************************************************************/

/** constructor */
SRMCreditHJM2F::SRMCreditHJM2F(QMCRatesDiffuseSP srmRatesDiffuse) :
    SRMCreditHJMDiffuse(srmRatesDiffuse),
    alpha0(-999.),
    alpha1(-999.),
    rho03(-999.),
    rho13(-999.)
{}

/** set the credit-rates correlation */
void SRMCreditHJM2F::setCrIrCorr(const vector<double>& corrCRIR) {
    rho03=corrCRIR[0];
    rho13=corrCRIR[1];
}

/** finalize the timelines, allocate necessary memory */
void SRMCreditHJM2F::finalizePathGenerator(DateTimeArrayConstSP allDatesSP) {
    static const string method("SRMCreditHJM2F::finalize");
    try {
        alpha0=ratesHJMUtil->getAlpha(0);
        alpha1=ratesHJMUtil->getAlpha(1);

        const int numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);
        kgtFactors.resize(zeroQ ? 0 : numSimDates-1);

        SRMCreditHJMDiffuse::finalizePathGenerator(allDatesSP); // will notify srmCreditHJMUtil about new timeline
        expProb.resize(getNumESDFDates());

        // initialization of quite a few parameters
        calcSigmaLParams(srmCreditHJMUtil->getSimDates());

        const DateTimeArray& dates = srmCreditHJMUtil->getSimDates(); // same as diffusionDates ?
        kgFactors.resize(dates.size()-1);
        assert(dates.size() == numSimDates);

        for (int i = 0; i < numSimDates - 1; i++) {
            kgFactors[i].kFactorCR = srmCreditHJMUtil->kFactor(dates[i], dates[i+1]);
            kgFactors[i].gFactorCR = srmCreditHJMUtil->gFactor(dates[i], dates[i+1]);
            kgFactors[i].kFactor0IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kgFactors[i].gFactor0IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kgFactors[i].kFactor1IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
            kgFactors[i].gFactor1IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
        }
        srmCreditHJMUtil = SRMCreditHJMUtilSP(0); // FIXME: see if it's needed
        ratesHJMUtil = SRMRatesHJMUtilSP(0);
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** generate path across all dates */
void SRMCreditHJM2F::generatePathAndSurvivalRates(IQMCRNGManagerSP rngMgr)  
{
    SRMRatesHJM2F* srmRates = static_cast<SRMRatesHJM2F *> (getUnderlyingIRAsset().get());

    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = logFwdProbSimple.size();

    int probDatePos = probStart; // position in probIndexes/prob.
    /* set up position in expProbIndexes/expSvob. If expProbIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of probStart and subtracting it in expDatePos formula
    int expProbDatePos = expProbIndexes.front() == 0? 1: 0;
    int probDateIdx = probIndexes[probDatePos]; // date when we save survival probability
    int expProbDateIdx = expProbIndexes[expProbDatePos]; // date when we save expected survival probability
    int stopIdx = Maths::min(probDateIdx, expProbDateIdx); // when to do something

    double LnQ = 0.0;      // exp(int lambda(u,u) du), ie 1 / survival probability
    double GAMMA3 = 0.0;
    double THETA0 = 0.0;
    double THETA1 = 0.0;
    double PHI03 = 0.0;
    double PHI13 = 0.0;
    double PHI33 = 0.0;

    // loop over all dates
    for (int i = 0; i < lastDiffusionIdx-todayIdx /*numDates*/; /* increment in body */){
        assert(i < numDates);
        double W3 = randoms[i];

        /* Retrieve precomputed k-factors & g-factors as well as kT-factors & gT-factors
           (precomputation takes place in constructor of SRMCreditHJM2F) */
        double k0 = kgFactors[i].kFactor0IR;
        double g0 = kgFactors[i].gFactor0IR;
        double k1 = kgFactors[i].kFactor1IR;
        double g1 = kgFactors[i].gFactor1IR;
        double k3 = kgFactors[i].kFactorCR;
        double g3 = kgFactors[i].gFactorCR;

        /* Step1a: update sigmaL (see FwdCrSprd_s and CalcSigmaLCR in crdiffuse) */
        double QPivot = getQPivot(i); //this->QPivot[i]; // (this->QPivot[i] is already lbar[i] * pivotRatio);
        double effL = 0.0;
        if(zeroQ) {
            effL = Maths::collar(QPivot, MaxEffRateCR[i], MinEffRateCR[i]);
        } else { // see FwdCrSprd_s in crdiffuse -- two factor case
            double k0T = kgtFactors[i].kTfactor0IR;
            double k1T = kgtFactors[i].kTfactor1IR;
            double g0T = kgtFactors[i].gTfactor0IR;
            double g1T = kgtFactors[i].gTfactor1IR;
            double k3T = kgtFactors[i].kTfactorCR;
            double g3T = kgtFactors[i].gTfactorCR;
            double aRate = lbarT[i] + k3T*GAMMA3 + k3T*g3T*PHI33
                            + k0T*THETA0 + (k0T*g3T+k3T*g0T)*PHI03
                            + k1T*THETA1 + (k1T*g3T+k3T*g1T)*PHI13;
            double lbarRatio = lbar[i] / lbarT[i];
            aRate *= lbarRatio;

            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);
            effL = (aRate >= QPivot) ?
                    (QPivot_x_OneMinusQR + (aRate * qRight)):
                    (QPivot_x_OneMinusQL + (aRate * qLeft));
                    effL = Maths::collar(effL, MaxEffRateCR[i], MinEffRateCR[i]);
        }
        double sigmaL = effL * svol[i];

        /* Step 1b: update sigmaR -- saved when IR simulation bit was done ... */
        double sigmaR0 = sigmaR[i] * alpha0;
        double sigmaR1 = sigmaR[i] * alpha1;

        /* Step 1c: compute modified k-factors and modified g-factors */
        double RootDelT = srmRates->getSqrtYearFrac(i);

        double sg3 = sigmaL*g3*RootDelT;
        double sk3 = sigmaL*k3*RootDelT;

        double sg0 = sigmaR0*g0*RootDelT;
        double sk0 = sigmaR0*k0*RootDelT;

        double sg1 = sigmaR1*g1*RootDelT;
        double sk1 = sigmaR1*k1*RootDelT;

        /* Step 2: update LnQ */
        double lnFwdProb = logFwdProbSimple[i]; // ratio S(0,T) / S(0,t) = deterministic term of S(t,T)
        double lnProbTerm = -g3*(g3*PHI33/2.0 + GAMMA3) // stoch term of S(t,T)
                    - g0*(g3*PHI03 + THETA0)
                    - g1*(g3*PHI13 + THETA1);

        LnQ += -(lnFwdProb+lnProbTerm) + sg3*(sg3/2.0 + sg0*rho03 + sg1*rho13 + W3);

        /* Step 3: update GAMMA */
        GAMMA3 = k3*(g3*PHI33 + g0*PHI03 + g1*PHI13 + GAMMA3)
                    + sk3*(sg3 + sg0*rho03 + sg1*rho13 + W3);

        /* Step 4: update THETA0 */
        THETA0 = k0*(g3*PHI03 + THETA0) + sk0*(sg3*rho03);
        THETA1 = k1*(g3*PHI13 + THETA1) + sk1*(sg3*rho13);

        /* Step 5: update PHI33 and PHI03 */
        PHI33 = k3*k3*PHI33 + sk3*sk3;
        PHI03 = k0*k3*PHI03 + sk0*sk3*rho03;
        PHI13 = k1*k3*PHI13 + sk1*sk3*rho13;

        /* Step 6: add cups component, if necessary ... this is related to FX */
        if (sigmaFX) {
            double sX = (*sigmaFX)[i]*RootDelT;
            GAMMA3 -= (sk3*sX*crFxCorr);
            LnQ -= (sg3*sX*crFxCorr);
        }

        i++; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i] = LnQ;


        if (i + todayIdx == stopIdx) {
            // hit an "event", ie record sthg
            if (stopIdx == probDateIdx) { // save survival probability
                prob[probDatePos] = LnQ;
                probDatePos++;
                probDateIdx = probIndexes[probDatePos];
            }
            if (stopIdx == expProbDateIdx) {
                expProb[expProbDatePos].GAMMA3 = GAMMA3;
                expProb[expProbDatePos].THETA0 = THETA0;
                expProb[expProbDatePos].THETA1 = THETA1;
                expProb[expProbDatePos].PHI03 = PHI03;
                expProb[expProbDatePos].PHI13 = PHI13;
                expProb[expProbDatePos].PHI33 = PHI33;
                expProbDatePos++;
                expProbDateIdx = expProbIndexes[expProbDatePos];
            }
            stopIdx = Maths::min(probDateIdx, expProbDateIdx); // refresh
        }
    }

    assert(prob.size() == probDatePos);
    assert(expProb.size() == expProbDatePos);

    // then do exp on prob vector
    for (size_t j = probStart; j < prob.size(); j++){
        prob[j] = exp(-prob[j]);
    }
}

/** accesses the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMCreditHJM2F::getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    return exp(SRMCreditHJM2F::getLnExpectedSurvivalDiscFactor(i, j));
}

/** accesses the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double SRMCreditHJM2F::getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    double g0 = static_cast<SRMRatesHJM2F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j), 0);
    double g1 = static_cast<SRMRatesHJM2F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j), 1);
    double g3 = getGFactorCR(i, j);

    const int iEDF = getTimeLogic()->getReqEDFIdx(i);
    const DiffusedState& state = expProb[iEDF];

    double result = getLnProbRatio(i,j) - g3*state.GAMMA3
                                        - 0.5*g3*g3*state.PHI33
                    - g0*state.THETA0 - g0*g3*state.PHI03
                    - g1*state.THETA1 - g1*g3*state.PHI13;
    return result;
}

/** compute kT and gT factors for CR as well as for IR */
void SRMCreditHJM2F::computeKGT(
                                int             dateIdx,
                                double          crKT,
                                double          crGT,
                                vector<double>& irKT,
                                vector<double>& irGT)
{
    if ((int)kgtFactors.size() <= dateIdx)
        kgtFactors.resize(dateIdx+1);
    kgtFactors[dateIdx].kTfactorCR = crKT;
    kgtFactors[dateIdx].gTfactorCR = crGT;
    kgtFactors[dateIdx].kTfactor0IR = irKT[0];
    kgtFactors[dateIdx].kTfactor1IR = irKT[1];
    kgtFactors[dateIdx].gTfactor0IR = irGT[0];
    kgtFactors[dateIdx].gTfactor1IR = irGT[1];
}

/************************************************************************/
/* Implementation of SRMCreditHJM3F class                               */
/************************************************************************/

/** constructor */
SRMCreditHJM3F::SRMCreditHJM3F(QMCRatesDiffuseSP srmRatesDiffuse) :
    SRMCreditHJMDiffuse(srmRatesDiffuse),
    alpha0(-999.),
    alpha1(-999.),
    alpha2(-999.),
    rho03(-999.),
    rho13(-999.),
    rho23(-999.)
{}

/** set the credit-rates correlation */
void SRMCreditHJM3F::setCrIrCorr(const vector<double>& corrCRIR) {
    rho03=corrCRIR[0];
    rho13=corrCRIR[1];
    rho23=corrCRIR[2];
}

/** finalize the timelines, allocate necessary memory */
void SRMCreditHJM3F::finalizePathGenerator(DateTimeArrayConstSP allDatesSP) {
    static const string method("SRMCreditHJM3F::finalize");
    try {
        alpha0=ratesHJMUtil->getAlpha(0);
        alpha1=ratesHJMUtil->getAlpha(1);
        alpha2=ratesHJMUtil->getAlpha(2);

        SRMCreditHJMDiffuse::finalizePathGenerator(allDatesSP); // srmCreditHJMUtil is already initialized

        const DateTimeArray& dates = srmCreditHJMUtil->getSimDates(); // same as diffusionDates ?
        const int numSimDates = dates.size();
        const int lastDiffIdx = getDiffusionBound()->getMaxDiffDate().findLower(dates);
//        cerr << "Asset " << this << " numSimDates= " << numSimDates << " lastDiffIdx= " << lastDiffIdx << endl;

        expProb.resize(getNumESDFDates() ); // TODO: trim (e)DFs to maxDiff/maxCurve dates

        // initialization of quite a few parameters
        calcSigmaLParams(dates);

        kgFactors.resize(numSimDates-1); // OK

        for (int i = 0; i < numSimDates - 1; i++) {
            kgFactors[i].kFactorCR = srmCreditHJMUtil->kFactor(dates[i], dates[i+1]);
            kgFactors[i].gFactorCR = srmCreditHJMUtil->gFactor(dates[i], dates[i+1]);
            // The kgFactors[i] are aliased with those from IR
//             kgFactors[i].kFactor0IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
//             kgFactors[i].gFactor0IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
//             kgFactors[i].kFactor1IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
//             kgFactors[i].gFactor1IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
//             kgFactors[i].kFactor2IR = ratesHJMUtil->kFactor(dates[i], dates[i+1], 2 /* nb IR factor */);
//             kgFactors[i].gFactor2IR = ratesHJMUtil->gFactor(dates[i], dates[i+1], 2 /* nb IR factor */);
        }

        srmCreditHJMUtil = SRMCreditHJMUtilSP(0);
        ratesHJMUtil = SRMRatesHJMUtilSP(0);
        trimToDiffusion();
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

template <class T>
    void trim(T& obj, size_t maxSize) {
    if (obj.size() > maxSize)
        obj.resize(maxSize);
}

/** compute kT and gT factors for CR as well as for IR */
void SRMCreditHJM3F::computeKGT(
    int             dateIdx,
    double          crKT,
    double          crGT,
    vector<double>& irKT,
    vector<double>& irGT)
{
    if ((int)kgtFactors.size() <= dateIdx)
        kgtFactors.resize(dateIdx+1);
    kgtFactors[dateIdx].kTfactorCR = crKT;
    kgtFactors[dateIdx].gTfactorCR = crGT;
    kgtFactors[dateIdx].kTfactor0IR = irKT[0];
    kgtFactors[dateIdx].kTfactor1IR = irKT[1];
    kgtFactors[dateIdx].kTfactor2IR = irKT[2];
    kgtFactors[dateIdx].gTfactor0IR = irGT[0];
    kgtFactors[dateIdx].gTfactor1IR = irGT[1];
    kgtFactors[dateIdx].gTfactor2IR = irGT[2];
}

/** trim internal arrays to the actual diffusion; non-virtual! */
void SRMCreditHJM3F::trimToDiffusion(void)
{
    size_t newSize = lastDiffusionIdx +1;
    trim(kgFactors, newSize);
    trim(kgtFactors, newSize);
    trim(lbar, newSize);
    trim(lbarT, newSize);
    trim(MaxEffRateCR, newSize);
    trim(MinEffRateCR, newSize);
/*    fwdIdx2IRedfIdx.clear();
    probIndexes.clear();
    expProbIndexes.clear();*/
}

struct SRMCreditHJM3F::State {
    double LnQ;      // exp(int lambda(u,u) du), ie 1 / survival probability
    double GAMMA3;
    double THETA0;
    double THETA1;
    double THETA2;
    double PHI03;
    double PHI13;
    double PHI23;
    double PHI33;
};

struct SRMCreditHJM3F::AuxArgs {
    const double* randoms;
    double        sigmaFX; // at step i:  0.0 or from FX: sigmaFX[i]
};

class SRMCreditHJM3F::Observer {
    SRMCreditHJM3F&      cr;
    public:
        Observer( SRMCreditHJM3F& _cr) : cr(_cr)
        {


        }
        void observe(const SRMCreditHJM3F::State& state, int i);
};

void SRMCreditHJM3F::advanceDiffusion(SRMCreditHJM3F::State& state,
                                const SRMCreditHJM3F::AuxArgs& aux,
                                int i)
{
    SRMRatesHJM3F* srmRates = static_cast<SRMRatesHJM3F *> (getUnderlyingIRAsset().get());

    ASSERT(srmRates != NULL);

    const KGFactors& kgFactor = kgFactors[i];

    double& LnQ = state.LnQ;      // exp(int lambda(u,u) du), ie 1 / survival probability
    double& GAMMA3 = state.GAMMA3;
    double& THETA0 = state.THETA0;
    double& THETA1 = state.THETA1;
    double& THETA2 = state.THETA2;
    double& PHI03 = state.PHI03;
    double& PHI13 = state.PHI13;
    double& PHI23 = state.PHI23;
    double& PHI33 = state.PHI33;

    const double& W3 = aux.randoms[i];

        /* Retrieve precomputed k-factors & g-factors as well as kT-factors & gT-factors
    (precomputation takes place in constructor of SRMCreditHJM3F) */
/*    const double k0 = kgFactor.kFactor0IR; assert(k0 == srmRates->varsPerTimePoint[i].k0);
    const double g0 = kgFactor.gFactor0IR; assert(g0 == srmRates->varsPerTimePoint[i].g0);
    const double k1 = kgFactor.kFactor1IR; assert(k1 == srmRates->varsPerTimePoint[i].k1);
    const double g1 = kgFactor.gFactor1IR; assert(g1 == srmRates->varsPerTimePoint[i].g1);
    const double k2 = kgFactor.kFactor2IR; assert(k2 == srmRates->varsPerTimePoint[i].k2);
    const double g2 = kgFactor.gFactor2IR; assert(g2 == srmRates->varsPerTimePoint[i].g2);
    */
    const double k0 = srmRates->varsPerTimePoint[i].k0;
    const double g0 = srmRates->varsPerTimePoint[i].g0;
    const double k1 = srmRates->varsPerTimePoint[i].k1;
    const double g1 = srmRates->varsPerTimePoint[i].g1;
    const double k2 = srmRates->varsPerTimePoint[i].k2;
    const double g2 = srmRates->varsPerTimePoint[i].g2;

    const double k3 = kgFactor.kFactorCR;
    const double g3 = kgFactor.gFactorCR;

    /* Step1a: update sigmaL (see FwdCrSprd_s and CalcSigmaLCR in crdiffuse) */
    double QPivot = getQPivot(i); //this->QPivot[i]; // (this->QPivot[i] is already lbar[i] * pivotRatio);
    double effL;
    if(zeroQ) {
        effL = QPivot;
    } else { // see FwdCrSprd_s in crdiffuse -- three factor case
        const KGTFactors& kgtFactor = kgtFactors[i];
        double k0T = kgtFactor.kTfactor0IR;
        double k1T = kgtFactor.kTfactor1IR;
        double k2T = kgtFactor.kTfactor2IR;
        double g0T = kgtFactor.gTfactor0IR;
        double g1T = kgtFactor.gTfactor1IR;
        double g2T = kgtFactor.gTfactor2IR;

        double k3T = kgtFactor.kTfactorCR;
        double g3T = kgtFactor.gTfactorCR;
        double aRate = lbarT[i] + k3T*GAMMA3 + k3T*g3T*PHI33
                    + k0T*THETA0 + (k0T*g3T+k3T*g0T)*PHI03
                    + k1T*THETA1 + (k1T*g3T+k3T*g1T)*PHI13
                    + k2T*THETA2 + (k2T*g3T+k3T*g2T)*PHI23;
        double lbarRatio = lbar[i] / lbarT[i];
        aRate *= lbarRatio;

        double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
        double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);
        effL = (aRate >= QPivot) ?
                    (QPivot_x_OneMinusQR + (aRate * qRight)):
                (QPivot_x_OneMinusQL + (aRate * qLeft));
    }
    effL = Maths::collar(effL, MaxEffRateCR[i], MinEffRateCR[i]);
    double sigmaL = effL * svol[i];

    /* Step 1b: update sigmaR -- saved when IR simulation bit was done ... */
    double sigmaR0 = sigmaR[i] * alpha0;
    double sigmaR1 = sigmaR[i] * alpha1;
    double sigmaR2 = sigmaR[i] * alpha2;

    /* Step 1c: compute modified k-factors and modified g-factors */
    double RootDelT = srmRates->getSqrtYearFrac(i);

    double sg3 = sigmaL*g3*RootDelT;
    double sk3 = sigmaL*k3*RootDelT;

    double sg0 = sigmaR0*g0*RootDelT;
    double sk0 = sigmaR0*k0*RootDelT;

    double sg1 = sigmaR1*g1*RootDelT;
    double sk1 = sigmaR1*k1*RootDelT;

    double sg2 = sigmaR2*g2*RootDelT;
    double sk2 = sigmaR2*k2*RootDelT;

    /* Step 2: update LnQ */
    double lnFwdProb = logFwdProbSimple[i]; // ratio S(0,T) / S(0,t) = deterministic term of S(t,T)
    double lnProbTerm = -g3*(g3*PHI33/2.0 + GAMMA3) // stoch term of S(t,T)
                - g0*(g3*PHI03 + THETA0)
                - g1*(g3*PHI13 + THETA1)
                - g2*(g3*PHI23 + THETA2);

    LnQ += -(lnFwdProb+lnProbTerm) + sg3*(sg3/2.0 + sg0*rho03 + sg1*rho13 + sg2*rho23 + W3);

    /* Step 3: update GAMMA */
    GAMMA3 = k3*(g3*PHI33 + g0*PHI03 + g1*PHI13 + g2*PHI23 + GAMMA3)
            + sk3*(sg3 + sg0*rho03 + sg1*rho13 + sg2*rho23 + W3);

    /* Step 4: update THETA0 */
    THETA0 = k0*(g3*PHI03 + THETA0) + sk0*(sg3*rho03);
    THETA1 = k1*(g3*PHI13 + THETA1) + sk1*(sg3*rho13);
    THETA2 = k2*(g3*PHI23 + THETA2) + sk2*(sg3*rho23);

    /* Step 5: update PHI33 and PHI03 */
    PHI33 = k3*k3*PHI33 + sk3*sk3;
    PHI03 = k0*k3*PHI03 + sk0*sk3*rho03;
    PHI13 = k1*k3*PHI13 + sk1*sk3*rho13;
    PHI23 = k2*k3*PHI23 + sk2*sk3*rho23;

    /* Step 6: add cups component, if necessary ... this is related to FX */
    double sX = aux.sigmaFX*RootDelT; // if no FX, pass sigmaFX = 0.0
    GAMMA3 -= (sk3*sX*crFxCorr);
    LnQ -= (sg3*sX*crFxCorr);

   // return state; // modified in place
}

void SRMCreditHJM3F::Observer::observe(const SRMCreditHJM3F::State& state, int i) {
    // We check if stop is requested at a given date by quering cr.timeLogic
    // If sdf was requested a !npos index will be returned., same for ESDF
    DateTime stopDate = cr.getTimeLineDate(i + cr.todayIdx);
    SpotIdx probDatePos = cr.getSpotIndex(stopDate);
    FwdIdx  fwdIdx = cr.getForwardForwardIndex(stopDate);
#if !defined (_MSC_VER) || (_MSC_VER >= 1300)
    SpotIdx expProbDatePos = (fwdIdx == FwdIdx::npos) ? SpotIdx(SpotIdx::npos) : cr.getTimeLogic()->getReqEDFIdx(fwdIdx);
    SpotIdx SpotIdxNpos = SpotIdx::npos;
#else
    SpotIdx expProbDatePos = (fwdIdx == -1) ? -1 : cr.getTimeLogic()->getReqEDFIdx(fwdIdx);
    SpotIdx SpotIdxNpos = -1;
#endif

    // Check if have to store SDF state
    if (probDatePos != SpotIdxNpos) { // save survival probability
            cr.prob[probDatePos] = exp(- state.LnQ);
    }
    // Check if we are to store ESDF state
    if (expProbDatePos  != SpotIdxNpos) {
        cr.expProb[expProbDatePos].GAMMA3 = state.GAMMA3;
        cr.expProb[expProbDatePos].THETA0 = state.THETA0;
        cr.expProb[expProbDatePos].THETA1 = state.THETA1;
        cr.expProb[expProbDatePos].THETA2 = state.THETA2;
        cr.expProb[expProbDatePos].PHI03 = state.PHI03;
        cr.expProb[expProbDatePos].PHI13 = state.PHI13;
        cr.expProb[expProbDatePos].PHI23 = state.PHI23;
        cr.expProb[expProbDatePos].PHI33 = state.PHI33;
    }
}

/** generate path across all dates */
void SRMCreditHJM3F::generatePathAndSurvivalRates(IQMCRNGManagerSP rngMgr)  
{
    SRMCreditHJM3F::AuxArgs auxArgs = {NULL, 0.0};
    SRMCreditHJM3F::State    state = {0,0,0,0,0,0,0,0,0};

    auxArgs.randoms = rngMgr->getCorrelatedRandoms(randomIndex);
    auxArgs.sigmaFX = 0;

    Observer    observer(*this);

    const int numDates = logFwdProbSimple.size();

    assert(lastDiffusionIdx - todayIdx <= numDates);

    for (int i = 0; i < lastDiffusionIdx - todayIdx; ++i)
    {
        if (sigmaFX)
            auxArgs.sigmaFX = (*sigmaFX)[i]; // update external state: FIXME: eventually we will not be able to read FX state by index!
        advanceDiffusion(state, auxArgs, i); // move from state[i] to [i+1]; no side-effects

        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i+1] = state.LnQ;

        observer.observe(state, i+1);     // fills in diffused state arrays; starts with 1 for compatib. with the old SRMCredit
    }
}

#if 0
/** generate path across all dates */
void SRMCreditHJM3F::generatePathAndSurvivalRates(IQMCRNGManagerSP rngMgr)  
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = logFwdProbSimple.size();

    int probDatePos = probStart; // position in probIndexes/prob.
    /* set up position in expProbIndexes/expSvob. If expProbIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of probStart and subtracting it in expDatePos formula
    int expProbDatePos = expProbIndexes.front() == 0? 1: 0;
    int probDateIdx = probIndexes[probDatePos]; // date when we save survival probability
    int expProbDateIdx = expProbIndexes[expProbDatePos]; // date when we save expected survival probability
    int stopIdx = Maths::min(probDateIdx, expProbDateIdx); // when to do something

    double LnQ = 0.0;      // exp(int lambda(u,u) du), ie 1 / survival probability
    double GAMMA3 = 0.0;
    double THETA0 = 0.0;
    double THETA1 = 0.0;
    double THETA2 = 0.0;
    double PHI03 = 0.0;
    double PHI13 = 0.0;
    double PHI23 = 0.0;
    double PHI33 = 0.0;

    // loop over all dates
    for (int i = 0; i < lastDiffusionIdx - todayIdx /*numDates*/; /* increment in body */){
        assert(i < numDates);
        double W3 = randoms[i];

        /* Retrieve precomputed k-factors & g-factors as well as kT-factors & gT-factors
           (precomputation takes place in constructor of SRMCreditHJM3F) */
        double k0 = kgFactors[i].kFactor0IR;
        double g0 = kgFactors[i].gFactor0IR;
        double k1 = kgFactors[i].kFactor1IR;
        double g1 = kgFactors[i].gFactor1IR;
        double k2 = kgFactors[i].kFactor2IR;
        double g2 = kgFactors[i].gFactor2IR;
        double k3 = kgFactors[i].kFactorCR;
        double g3 = kgFactors[i].gFactorCR;

        /* Step1a: update sigmaL (see FwdCrSprd_s and CalcSigmaLCR in crdiffuse) */
        double QPivot = this->QPivot[i]; // (this->QPivot[i] is already lbar[i] * pivotRatio);
        double effL = 0.0;
        if(zeroQ) {
            effL = QPivot;
        } else { // see FwdCrSprd_s in crdiffuse -- three factor case
            double k0T = kgtFactors[i].kTfactor0IR;
            double k1T = kgtFactors[i].kTfactor1IR;
            double k2T = kgtFactors[i].kTfactor2IR;
            double g0T = kgtFactors[i].gTfactor0IR;
            double g1T = kgtFactors[i].gTfactor1IR;
            double g2T = kgtFactors[i].gTfactor2IR;
            double k3T = kgtFactors[i].kTfactorCR;
            double g3T = kgtFactors[i].gTfactorCR;
            double aRate = lbarT[i] + k3T*GAMMA3 + k3T*g3T*PHI33
                            + k0T*THETA0 + (k0T*g3T+k3T*g0T)*PHI03
                            + k1T*THETA1 + (k1T*g3T+k3T*g1T)*PHI13
                            + k2T*THETA2 + (k2T*g3T+k3T*g2T)*PHI23;
            double lbarRatio = lbar[i] / lbarT[i];
            aRate *= lbarRatio;

            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);
            effL = (aRate >= QPivot) ?
                    (QPivot_x_OneMinusQR + (aRate * qRight)):
                    (QPivot_x_OneMinusQL + (aRate * qLeft));
        }
        effL = Maths::collar(effL, MaxEffRateCR[i], MinEffRateCR[i]);
        double sigmaL = effL * svol[i];

        /* Step 1b: update sigmaR -- saved when IR simulation bit was done ... */
        double sigmaR0 = sigmaR[i] * alpha0;
        double sigmaR1 = sigmaR[i] * alpha1;
        double sigmaR2 = sigmaR[i] * alpha2;

        /* Step 1c: compute modified k-factors and modified g-factors */
        double RootDelT = sqrtYearFrac[i];

        double sg3 = sigmaL*g3*RootDelT;
        double sk3 = sigmaL*k3*RootDelT;

        double sg0 = sigmaR0*g0*RootDelT;
        double sk0 = sigmaR0*k0*RootDelT;

        double sg1 = sigmaR1*g1*RootDelT;
        double sk1 = sigmaR1*k1*RootDelT;

        double sg2 = sigmaR2*g2*RootDelT;
        double sk2 = sigmaR2*k2*RootDelT;

        /* Step 2: update LnQ */
        double lnFwdProb = logFwdProbSimple[i]; // ratio S(0,T) / S(0,t) = deterministic term of S(t,T)
        double lnProbTerm = -g3*(g3*PHI33/2.0 + GAMMA3) // stoch term of S(t,T)
                    - g0*(g3*PHI03 + THETA0)
                    - g1*(g3*PHI13 + THETA1)
                    - g2*(g3*PHI23 + THETA2);

        LnQ += -(lnFwdProb+lnProbTerm) + sg3*(sg3/2.0 + sg0*rho03 + sg1*rho13 + sg2*rho23 + W3);

        /* Step 3: update GAMMA */
        GAMMA3 = k3*(g3*PHI33 + g0*PHI03 + g1*PHI13 + g2*PHI23 + GAMMA3)
                    + sk3*(sg3 + sg0*rho03 + sg1*rho13 + sg2*rho23 + W3);

        /* Step 4: update THETA0 */
        THETA0 = k0*(g3*PHI03 + THETA0) + sk0*(sg3*rho03);
        THETA1 = k1*(g3*PHI13 + THETA1) + sk1*(sg3*rho13);
        THETA2 = k2*(g3*PHI23 + THETA2) + sk2*(sg3*rho23);

        /* Step 5: update PHI33 and PHI03 */
        PHI33 = k3*k3*PHI33 + sk3*sk3;
        PHI03 = k0*k3*PHI03 + sk0*sk3*rho03;
        PHI13 = k1*k3*PHI13 + sk1*sk3*rho13;
        PHI23 = k2*k3*PHI23 + sk2*sk3*rho23;

        /* Step 6: add cups component, if necessary ... this is related to FX */
        if (sigmaFX) {
            double sX = (*sigmaFX)[i]*RootDelT;
            GAMMA3 -= (sk3*sX*crFxCorr);
            LnQ -= (sg3*sX*crFxCorr);
        }

        i++; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i] = LnQ;


        if (i + todayIdx == stopIdx) {
            // hit an "event", ie record sthg
            if (stopIdx == probDateIdx) { // save survival probability
                prob[probDatePos] = LnQ;
                probDatePos++;
                probDateIdx = probIndexes[probDatePos];
            }
            if (stopIdx == expProbDateIdx) {
                expProb[expProbDatePos].GAMMA3 = GAMMA3;
                expProb[expProbDatePos].THETA0 = THETA0;
                expProb[expProbDatePos].THETA1 = THETA1;
                expProb[expProbDatePos].THETA2 = THETA2;
                expProb[expProbDatePos].PHI03 = PHI03;
                expProb[expProbDatePos].PHI13 = PHI13;
                expProb[expProbDatePos].PHI23 = PHI23;
                expProb[expProbDatePos].PHI33 = PHI33;
                expProbDatePos++;
                expProbDateIdx = expProbIndexes[expProbDatePos];
            }
            stopIdx = Maths::min(probDateIdx, expProbDateIdx); // refresh
        }
    }
    assert(prob.size() == probDatePos);
    assert(expProb.size() == expProbDatePos);
    // then do exp on prob vector
    for (size_t j = probStart; j < prob.size(); j++){
        prob[j] = exp(-prob[j]);
    }
}
#endif

/** accesses the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMCreditHJM3F::getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    return exp(SRMCreditHJM3F::getLnExpectedSurvivalDiscFactor(i, j));
}

/** accesses the natural log of the expected value ExpSDF(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made. */
double SRMCreditHJM3F::getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j) {
    double g0 = static_cast<SRMRatesHJM3F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j), 0);
    double g1 = static_cast<SRMRatesHJM3F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j), 1);
    double g2 = static_cast<SRMRatesHJM3F*>(irAsset.get())->getGFactor(getIRfwdIdx(i),
                                                getIRfwdIdx(j), 2);
    double g3 = getGFactorCR(i, j);

    const int iEDF = getTimeLogic()->getReqEDFIdx(i);
    const DiffusedState& state = expProb[iEDF];

    double result = getLnProbRatio(i,j)  - g3*state.GAMMA3 - 0.5*g3*g3*state.PHI33
                    - g0*state.THETA0 - g0*g3*state.PHI03
                    - g1*state.THETA1 - g1*g3*state.PHI13
                    - g2*state.THETA2 - g2*g3*state.PHI23;
    return result;
}

DRLIB_END_NAMESPACE
