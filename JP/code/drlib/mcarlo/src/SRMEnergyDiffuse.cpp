//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEnergyDiffuse.cpp
//
//   Description : An implementation for a template (new) asset diffusion
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCEnergyDiffuse.hpp"
#include "edginc/SRMEnergyDiffuse.hpp"
#include "edginc/SRMEnergyUtil.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/IQMCRNGManager.hpp"
#include <cassert>

#if 0
#include <fstream> // for logging to file
#endif

DRLIB_BEGIN_NAMESPACE

/************************************************************************
 * These classes form the interface for all the diffusible assets       *
 * they are convenient for a few generic methods acting on all the      *
 * different assets uniformly and for model-independent access to the   *
 * data the models produce                                              *
 ************************************************************************/



/**********************************************************************************/
/******************************* Oil Model Stuffs *********************************/
/**********************************************************************************/
SRMEnergyOilDiffuse::SRMEnergyOilDiffuse(QMCRatesDiffuseSP srmRatesDiffuse) :
    QMCEnergyDiffuse(srmRatesDiffuse),
    rhoFxEnrg0(-999.),
    rhoFxEnrg1(-999.)
    {}

void SRMEnergyOilDiffuse::setEnrgFxCorr(const vector<double> & _enrgFxCorr)
{
    ASSERT(!_enrgFxCorr.empty());
    rhoFxEnrg0 = _enrgFxCorr[0];
    rhoFxEnrg1 = _enrgFxCorr[1];
}

void SRMEnergyOilDiffuse::finalize(DateTimeArrayConstSP allDates)
{
    // call the underlying finalize method
    QMCEnergyDiffuse::finalize(allDates);

    static const string method("SRMEnergyOilDiffuse::finalize");
    DateTimeArray efpRequestedDates = getForwardDates();
    DateTimeArray efpForwardDates = getForwardForwardDates();

    // get initial future price(s) for the requested maturity date(s)
    initPrices.resize(efpForwardDates.size());
    for (size_t i = 0; i < initPrices.size(); ++i) {
        initPrices[i] = log( futureCurve->interpolate(efpForwardDates[i]) );
    }

    // allocate mem for state vector
    states.resize(efpRequestedDates.size());

    // compute quanto adjustment if necessary:
    if (sigmasFx) {
        MUs.resize(nDiffusionSteps);
        calcQuantoAdjustments();
    }

    // pre-diffuse the Phi's: k-factors and states vector should be ready:
    preDiffusePhi();

    return;
}

// compute quanto adjustment for the 2-factor oil model
void SRMEnergyOilDiffuse::calcQuantoAdjustments()
{
    ASSERT(MUs.size() == dtSqrts.size());
    ASSERT(ks.size() == dtSqrts.size());
    // TODO: shall we check the size of the sigmasFx array?

    // get calibrated quantities from util class:
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();
    double sigma1, sigma1Bar, sigma2Bar;

    for (size_t i = 0; i < nDiffusionSteps; ++i)
    {
        sigma1 = sigmas[i][0];
        sigma1Bar = sigmas[i][1];
        sigma2Bar = sigmas[i][2];

        // compute and store quanto adjustments in the mu vector:
        MUs[i].mu0 = -rhoFxEnrg0 * sigma1 / alpha * (1.0 - ks[i]);
        MUs[i].mu1 = -(rhoFxEnrg0*sigma1Bar + rhoFxEnrg1*sigma2Bar) *
                     dtSqrts[i] * dtSqrts[i];
    }
}

// only call this at the end of model-dependent finalize() method!
void SRMEnergyOilDiffuse::preDiffusePhi()
{
    double k;
    double dtSqrt;
    double sigma1, sigma1Bar, sigma2Bar;

    // get calibrated model vols
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();

    // set up position in efpIndexes. If efpIndexes[0] is 0 then
    // we just store 0 for the GAMMA/PHI and start at efpIndexes[1]
    size_t expFpIdx = 0;
    if (efpIndexes[0] == todayIdx) // first diffusion date == today?
    {
        ++expFpIdx; // bypass the first expected future date
        states[0].PHI00 = 0.; // and set the phi's to 0
        states[0].PHI01 = 0.;
        states[0].PHI11 = 0.;
    }

    double PHI00 = 0.;
    double PHI01 = 0.;
    double PHI11 = 0.;

    for (size_t i = 0, dateIdx = todayIdx + 1; // first date to diffuse is the sim date immediately follows today
        i < nDiffusionSteps;
        ++i, ++dateIdx
        )
    {
        k = ks[i];
        dtSqrt = dtSqrts[i];
        sigma1 = sigmas[i][0];
        sigma1Bar = sigmas[i][1];
        sigma2Bar = sigmas[i][2];

        // diffuse PHI
        PHI00 = k*k*PHI00 + sigma1*sigma1*(1 - k*k)*0.5/alpha;
        PHI11 = PHI11 + (sigma1Bar*sigma1Bar + sigma2Bar*sigma2Bar)*dtSqrt*dtSqrt;
        PHI01 = k*PHI01 + sigma1*sigma1Bar*(1 - k)/alpha;

        if (dateIdx == efpIndexes[expFpIdx])
        {
            // store internal diffusion states
            states[expFpIdx].PHI00 = PHI00;
            states[expFpIdx].PHI01 = PHI01;
            states[expFpIdx].PHI11 = PHI11;
            ++expFpIdx;
        }
    }
    return;
}

void SRMEnergyOilDiffuse::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    // get random numbers
    const double* randoms = rngMgr->getCorrelatedRandoms(rndIdx);

    // get calibrated model vols
    //const DoubleMatrix & sigmas = enrgUtil->getSigmas();
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();

    // set up position in efpIndexes. If efpIndexes[0] is 0 then
    // we just store 0 for the GAMMA/PHI and start at efpIndexes[1]
    size_t expFpIdx = 0;
    if (efpIndexes[0] == todayIdx) // first diffusion date == today?
    {
        ++expFpIdx; // bypass the first expected future date
        states[0].GAMMA0 = 0.;  // and set the states to zero
        states[0].GAMMA1 = 0.;
    }

    double W0, W1;
    double k;
    double dtSqrt;
    double sigma1, sigma1Bar, sigma2Bar, sigmaFx;

    // setup pointers to random numbers
    const double* Z0 = rngMgr->getCorrelatedRandoms(rndIdx);
    const double* Z1 = rngMgr->getCorrelatedRandoms(rndIdx, 1);

    // initialze diffusion states
    double GAMMA0 = 0.;
    double GAMMA1 = 0.;

    for (size_t i = 0, dateIdx = todayIdx + 1; // first date to diffuse is the sim date immediately follows today
        i < nDiffusionSteps;
        ++i, ++dateIdx)
    {
        k = ks[i];
        W0 = Z0[i];
        W1 = Z1[i];
        dtSqrt = dtSqrts[i];
        sigma1 = sigmas[i][0];
        sigma1Bar = sigmas[i][1];
        sigma2Bar = sigmas[i][2];

        // diffuse Gamma
        GAMMA0 = k*GAMMA0 + sigma1*sqrt((1.0 - k*k)*0.5/alpha)*W0;
        GAMMA1 = GAMMA1 + sigma1Bar*dtSqrt*W0 + sigma2Bar*dtSqrt*W1;

        // add quanto adjustments if necessary
        if (sigmasFx) {
            sigmaFx = (*sigmasFx)[i];
            GAMMA0 += MUs[i].mu0 * sigmaFx;
            GAMMA1 += MUs[i].mu1 * sigmaFx;
        }

        //++i;
        //if (i + todayIdx == efpIndexes[expFpIdx])
        if (dateIdx == efpIndexes[expFpIdx])
        {
            // store internal diffusion states
            states[expFpIdx].GAMMA0 = GAMMA0;
            states[expFpIdx].GAMMA1 = GAMMA1;
            ++expFpIdx;

            #if 0
            // dump random numbers at the measurement date on the path
            ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
            file << W0 << "\t" << W1 << endl;
            #endif
        }
    }

#if 0
    // dump all random numbers along the path:
    ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
    if (file.is_open()) {
        for (size_t i = 0; i < nDiffusionSteps; ++i) {
            file << i+1 << "\t" << Z0[i] << "\t" << Z1[i] << endl;
        }
    }
#endif

#if 0
    // dump random numbers at a particular point on the path:
    ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
    if (file.is_open()) {
        file << Z0[1] << "\t" << Z1[1] << endl;
    }
#endif

}

/** Accessing the expected value ExpFP(md, fp) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMEnergyOilDiffuse::getExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                       FwdIdx j/*futureDateIdx*/)
{
    return exp(getLnExpectedPrice(i, j));
}

/** Accessing the natural log of the ExpectedFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made.
    Note: both i and j should be obtained by calling getExpFuturePriceIndex(date),
    in other words, i and j should be quoted as idx of forward maturity dates **/
double SRMEnergyOilDiffuse::getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                         FwdIdx j/*futureDateIdx*/)
{
    ASSERT(0<= i && i < (int)fwdIdx2RequestedIdx.size()); // FIXME: remove
    int fwd2ReqIdx = fwdIdx2RequestedIdx[i];
    ASSERT(fwd2ReqIdx >= 0);
    const DiffusedState & state = states[fwd2ReqIdx]; // get the diffused state
    double lnF0 = initPrices[j]; // get the log initial future price
    double kij0 = kTs[j] / kts[fwd2ReqIdx]; // get the k-factors for interval [t, T]
    double kij1 = 1.0;

    // compute fwd prices based on remembered internal diffusion states
    double shock = kij0*state.GAMMA0 + kij1*state.GAMMA1 -
                    0.5 * (kij0*(kij0*state.PHI00 + kij1*state.PHI01) +
                           kij1*(kij0*state.PHI01 + kij1*state.PHI11));

    return lnF0 + shock;
}


/**********************************************************************************/
/******************************* Gas Model Stuffs *********************************/
/**********************************************************************************/
SRMEnergyGasDiffuse::SRMEnergyGasDiffuse(QMCRatesDiffuseSP srmRatesDiffuse) :
    QMCEnergyDiffuse(srmRatesDiffuse),
    rhoFxEnrg0(-999.),
    rhoFxEnrg1(-999.)
    {}

void SRMEnergyGasDiffuse::setEnrgFxCorr(const vector<double> & _enrgFxCorr)
{
    ASSERT(!_enrgFxCorr.empty());
    rhoFxEnrg0 = _enrgFxCorr[0];
    rhoFxEnrg1 = _enrgFxCorr[1];
}

void SRMEnergyGasDiffuse::finalize(DateTimeArrayConstSP allDates)
{
    // call the underlying finalize method
    QMCEnergyDiffuse::finalize(allDates);

    static const string method("SRMEnergyGasDiffuse::finalize");

    DateTimeArray efpRequestedDates = getForwardDates();
    DateTimeArray efpForwardDates = getForwardForwardDates();

    // get initial future price(s) for the requested maturity date(s)
    initPrices.resize(efpForwardDates.size());
    for (size_t i = 0; i < initPrices.size(); ++i) {
        initPrices[i] = log( futureCurve->interpolate(efpForwardDates[i]) );
    }

    // Gas model also need k-factors for beta:
    double beta = enrgUtil->getBeta();
    k2s.assign(nDiffusionSteps, 0.0);
    kt2s.assign(efpRequestedDates.size(), 0.0);
    kT2s.assign(efpForwardDates.size(), 0.0);
    SRMEnergyUtilBase::computeKFactors(ks, beta, dtSqrts);
    SRMEnergyUtilBase::computeCumKFactors(kts, beta, today, efpRequestedDates);
    SRMEnergyUtilBase::computeCumKFactors(kTs, beta, today, efpForwardDates);

    // allocate mem for state vector
    states.resize(efpRequestedDates.size());

    // compute quanto adjustment if necessary:
    if (sigmasFx) {
        MUs.resize(nDiffusionSteps);
        calcQuantoAdjustments();
    }

    // pre-diffuse the Phi's: k-factors and states vector should be ready:
    preDiffusePhi();

    return;
}

// only call this at the end of model-dependent finalize() method!
void SRMEnergyGasDiffuse::preDiffusePhi()
{
    double k, k2;
    double dtSqrt;
    double sigma1, sigma1Bar, sigma2Bar;

    // get calibrated model vols
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();
    double beta = enrgUtil->getBeta();

    // set up position in efpIndexes. If efpIndexes[0] is 0 then
    // we just store 0 for the GAMMA/PHI and start at efpIndexes[1]
    size_t expFpIdx = 0;
    if (efpIndexes[0] == todayIdx) // first diffusion date == today?
    {
        ++expFpIdx; // bypass the first expected future date
        states[0].PHI11 = 0.; // and set the phi's to 0
        states[0].PHI22 = 0.;
        states[0].PHI33 = 0.;
        states[0].PHI21 = 0.;
        states[0].PHI31 = 0.;
        states[0].PHI32 = 0.;
    }

    double PHI11 = 0.;
    double PHI22 = 0.;
    double PHI33 = 0.;
    double PHI21 = 0.;
    double PHI31 = 0.;
    double PHI32 = 0.;

    for (size_t i = 0, dateIdx = todayIdx + 1; // first date to diffuse is the sim date immediately follows today
        i < nDiffusionSteps;
        ++i, ++dateIdx
        )
    {
        k = ks[i];
        k2 = k2s[i];
        dtSqrt = dtSqrts[i];
        sigma1 = sigmas[i][0];
        sigma1Bar = sigmas[i][1];
        sigma2Bar = sigmas[i][2];

        // diffuse PHI
        PHI11 = k*k*PHI11 + sigma1*sigma1*(1 - k*k)*0.5/alpha;
        PHI22 = k2*k2*PHI22 + sigma1*sigma1*(1 - k2*k2)*0.5/beta;
        PHI33 = PHI33 + (sigma1Bar*sigma1Bar + sigma2Bar*sigma2Bar)*dtSqrt*dtSqrt;
        PHI21 = k*k2*PHI21 + sigma1*sigma1Bar*(1 - k*k2)/(alpha + beta);
        PHI31 = k*PHI31 + sigma1*sigma1*(1 - k)/alpha;
        PHI32 = k2*PHI32 + sigma1*sigma1Bar*(1 - k2)/beta;

        if (dateIdx == efpIndexes[expFpIdx])
        {
            // store internal diffusion states
            states[expFpIdx].PHI11 = PHI11;
            states[expFpIdx].PHI22 = PHI22;
            states[expFpIdx].PHI33 = PHI33;
            states[expFpIdx].PHI21 = PHI21;
            states[expFpIdx].PHI31 = PHI31;
            states[expFpIdx].PHI32 = PHI32;
            ++expFpIdx;
        }
    }
    return;
}

/** Accessing the expected value ExpFP(md, fp) where md is a
simulated measurement date and fd is some future date after the
measurement is made. */
double SRMEnergyGasDiffuse::getExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                      FwdIdx j/*futureDateIdx*/)
{
    return exp(getLnExpectedPrice(i, j));
}

/** Accessing the natural log of the ExpectedFP(md, fd)
where md is a simulated measurement date and fd is some future
date after the measurement is made.
Note: both i and j should be obtained by calling getExpFuturePriceIndex(date),
in other words, i and j should be quoted as idx of forward maturity dates **/
double SRMEnergyGasDiffuse::getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                        FwdIdx j/*futureDateIdx*/)
{
    ASSERT(0<= i && i < (int)fwdIdx2RequestedIdx.size()); // FIXME: remove
    int fwd2ReqIdx = fwdIdx2RequestedIdx[i];
    ASSERT(fwd2ReqIdx >= 0); // FIXME: remove
    const DiffusedState & state = states[fwd2ReqIdx]; // get the diffused state
    double lnF0 = initPrices[j]; // get the log initial future price
    double kij1 = kTs[j] / kts[fwd2ReqIdx]; // alpha-k-factors for interval [t, T]
    double kij2 = kT2s[j] / kt2s[fwd2ReqIdx]; // beta-k-factors for interval [t, T]
    double kij3 = 1.0;

    // compute fwd prices based on remembered internal diffusion states
    double shock = kij1*state.GAMMA1 + kij2*state.GAMMA2 + kij3*state.GAMMA3
        - 0.5 * (kij1*(kij1*state.PHI11 + kij2*state.PHI21 + kij3*state.PHI31) +
                 kij2*(kij1*state.PHI21 + kij2*state.PHI22 + kij3*state.PHI32) +
                 kij3*(kij1*state.PHI31 + kij2*state.PHI32 + kij3*state.PHI33));

    return lnF0 + shock;
}


////////////////////////////////////////////////////////////////////////
/////////////////////// Tier 2 Energy Stuffs ///////////////////////////
////////////////////////////////////////////////////////////////////////
SRMEnergyTier2Diffuse::SRMEnergyTier2Diffuse(
    QMCRatesDiffuseSP srmRatesDiffuse,
    QMCEnergyDiffuseSP _parent)
    :
    QMCEnergyDiffuse(srmRatesDiffuse),
    rhoFxEnrg(-999.),
    parent(_parent)
{
    getDiffusionBound()->add(_parent->getDiffusionBound());
}

/*
SRMEnergyTier2Diffuse::SRMEnergyTier2Diffuse(QMCRatesDiffuseSP srmRatesDiffuse) :
    QMCEnergyDiffuse(srmRatesDiffuse),
    rhoFxEnrg(-999.)
    {}

void SRMEnergyTier2Diffuse::setParent(QMCEnergyDiffuseSP _parent)
{
    ASSERT(_parent.get());
    parent = _parent;
}*/

void SRMEnergyTier2Diffuse::setEnrgFxCorr(const vector<double> & _enrgFxCorr)
{
    ASSERT(!_enrgFxCorr.empty());
    rhoFxEnrg = _enrgFxCorr[0];
}

/** finalize the timelines, allocate necessary memory */
void SRMEnergyTier2Diffuse::finalize(DateTimeArrayConstSP allDates)
{
    // make sure parent is initialized
    ASSERT(parent.get());

    // call the underlying finalize method
    QMCEnergyDiffuse::finalize(allDates);

    static const string method("SRMEnergyTier2Diffuse::finalize");

    DateTimeArray efpRequestedDates = getForwardDates();
    DateTimeArray efpForwardDates = getForwardForwardDates();

    // get the initial tier2 spread(s) for the requested maturity date(s)
    initPrices.resize(efpForwardDates.size());
    for (size_t i = 0; i < initPrices.size(); ++i) {
        initPrices[i] = futureCurve->interpolate(efpForwardDates[i]);
    }

    // allocate memory for internal state variable vector
    GAMMAs.assign(efpRequestedDates.size(), 0.0);

    // compute quanto adjustment if necessary:
    if (sigmasFx) {
        MUs.resize(nDiffusionSteps);
        calcQuantoAdjustments();
    }

    // create tier2-tier1 measurement (requested) dates mapping
    DateTimeArray parentFwdDates = parent->getForwardForwardDates(); /*getTimeLogic()->getFwdEDFDates();*/
    assert(DateTime::isSubset(parentFwdDates, efpRequestedDates));
    child2ParentReqIdxs = DateTime::getIndexes(parentFwdDates, efpRequestedDates);

    // create maturity mapping indices
    const DateTimeArray & parentForwardDates = parent->getForwardForwardDates();
    child2ParentLowerIdxs = DateTime::getFloorProjection(efpForwardDates, parentForwardDates);
    child2ParentUpperIdxs = DateTime::getCeilingProjection(efpForwardDates, parentForwardDates);

    // compute interpolation table
    EnergyUnderlyerConstSP energyUnderlyer = futureCurve->getEnergyUnderlyer();
    string modelType = energyUnderlyer->getUnderlyerType();
    HolidayConstSP holiday = energyUnderlyer->getHoliday();

#if 0
    // get ready for dumping nBizDays:
    ofstream file("c:/debugEnergySRM3.txt", ios_base::app);
#endif

    if (CString::equalsIgnoreCase(modelType, EnrgOilStr))
    {
        // Oil model interpolation:
        parentLowerInterpRatios.clear();
        parentUpperInterpRatios.clear();
        for (int i = 0; i < efpForwardDates.size(); ++i)
        {
            int nBizDayLower = 0;
            int nBizDayUpper = 0;
            int nBizDayTotal = 0;
            double lowerInterpRatio;
            double upperInterpRatio;
            DateTime parentFwdDateLower = parentForwardDates[child2ParentLowerIdxs[i]];
            DateTime parentFwdDateUpper = parentForwardDates[child2ParentUpperIdxs[i]];

            /*
            // use the number of business days to compute interpolation ratios:
            for (int offset = 1; offset <= parentFwdDateUpper.daysDiff(parentFwdDateLower); ++offset)
            {
                DateTime interDate = parentFwdDateLower.rollDate(offset);
                if (energyUnderlyer->isBusinessDay(interDate)) {
                    ++nBizDayTotal;
                    if (interDate <= efpForwardDates[i])
                        ++nBizDayLower;
                    else
                        ++nBizDayUpper;
                }
            }
#if 0
            // dump nBizDays:
            if (file.is_open()) {
                file << nBizDayLower << "\t" << nBizDayUpper << "\t";
            }
#endif
            */

            // use the number of business days to compute interpolation ratios:
            if (parentFwdDateLower != parentFwdDateUpper)
            {
                nBizDayLower = holiday->businessDaysDiff(parentFwdDateLower, efpForwardDates[i]);
                nBizDayUpper = holiday->businessDaysDiff(efpForwardDates[i], parentFwdDateUpper);
                nBizDayTotal = nBizDayLower + nBizDayUpper;
            }
#if 0
            // dump nBizDays:
            if (file.is_open()) {
                file << nBizDayLower << "\t" << nBizDayUpper << endl;
            }
#endif

            if (nBizDayTotal == 0) {// either the parent fwd maturities are the same, or they are both over the weekend (or holidays)!!
                lowerInterpRatio = 0.5;
                upperInterpRatio = 0.5;
            }
            else {
                lowerInterpRatio = (double)nBizDayLower / (double)nBizDayTotal;
                upperInterpRatio = 1.0 - lowerInterpRatio;
            }

            parentLowerInterpRatios.push_back(lowerInterpRatio);
            parentUpperInterpRatios.push_back(upperInterpRatio);
        }
    }
    else if (CString::equalsIgnoreCase(modelType, EnrgGasStr))
    {
        // Gas model interpolation:
        parentLowerInterpRatios.assign(efpForwardDates.size(), 1.0);
        parentUpperInterpRatios.assign(efpForwardDates.size(), 0.0);
    }
    else if (CString::equalsIgnoreCase(modelType, EnrgPowStr))
    {
        // Power model interpolation:
        parentLowerInterpRatios.assign(efpForwardDates.size(), 1.0);
        parentUpperInterpRatios.assign(efpForwardDates.size(), 0.0);
    }
}

// compute quanto adjustment for the 2-factor oil model
void SRMEnergyTier2Diffuse::calcQuantoAdjustments()
{
    ASSERT(MUs.size() == dtSqrts.size());
    ASSERT(ks.size() == dtSqrts.size());
    // TODO: shall we check the size of the sigmasFx array?

    // get calibrated quantities from util class:
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();

    for (size_t i = 0; i < nDiffusionSteps; ++i)
    {
        // compute and store quanto adjustments in the mu vector:
        MUs[i] = -rhoFxEnrg * sigmas[i][0] / alpha * (1.0 - ks[i]);
    }
}

void SRMEnergyTier2Diffuse::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    // get random numbers
    const double* randoms = rngMgr->getCorrelatedRandoms(rndIdx);

    // get calibrated model vols
    const DoubleMatrix & sigmas = enrgUtil->getExtendedSigmas();
    double alpha = enrgUtil->getAlpha();

    // set up position in eSpIndexes. If eSpIndexes[0] is 0 then
    // we just store 0 for the GAMMA/PHI and start at efpIndexes[1]
    size_t expFpIdx = 0;
    if (efpIndexes[0] == todayIdx) // first diffusion date == today?
    {
        ++expFpIdx; // bypass the first expected future date
        GAMMAs[0] = 0.0;    // and set the states to zero
    }

    double W;
    double k;
    double dtSqrt;
    double sigma;

    // setup pointers to random numbers
    const double* Z = rngMgr->getCorrelatedRandoms(rndIdx);

    // initialze diffusion states
    double GAMMA = 0.;

    for (size_t i = 0, dateIdx = todayIdx + 1; // first date to diffuse is the sim date immediately follows today
        i < nDiffusionSteps;
        ++i, ++dateIdx)
    {
        k = ks[i];
        W = Z[i];
        dtSqrt = dtSqrts[i];
        sigma = sigmas[i][0];

        // diffuse Gamma
        GAMMA = k*GAMMA + sigma*sqrt((1.0 - k*k)*0.5/alpha)*W;

        // add quanto adjustments if necessary
        if (sigmasFx) {
            GAMMA += MUs[i] * (*sigmasFx)[i];
        }

        //if (i + todayIdx == efpIndexes[expFpIdx])
        if (dateIdx == efpIndexes[expFpIdx])
        {
            // store internal diffusion states
            GAMMAs[expFpIdx] = GAMMA;
            ++expFpIdx;
        }
    }
}

/** Accessing the expected value ExpFP(md, fp) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMEnergyTier2Diffuse::getLnExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                          FwdIdx j/*futureDateIdx*/)
{
    return log(getExpectedPrice(i, j));
}

/** Accessing the natural log of the ExpectedFP(md, fd)
    where md is a simulated measurement date and fd is some future
    date after the measurement is made.
    Note: both i and j should be obtained by calling getExpFuturePriceIndex(date),
    in other words, i and j should be quoted as idx of forward maturity dates **/
double SRMEnergyTier2Diffuse::getExpectedPrice( FwdIdx i/*measurementDateIdx*/,
                                          FwdIdx j/*futureDateIdx*/)
{
    ASSERT(i >= 0 && i < (int)fwdIdx2RequestedIdx.size()); // FIXME: remove
    int fwd2ReqIdx = fwdIdx2RequestedIdx[i];
    
    ASSERT(fwd2ReqIdx >= 0); // FIXME: remove
    double GAMMA = GAMMAs[fwd2ReqIdx]; // get the diffused state
    double S0 = initPrices[j]; // get time 0 spread
    double kij = kTs[j] / kts[fwd2ReqIdx]; // get the k-factors for interval [t, T]

    // compute fwd price based on remembered internal diffusion states
    double F = S0 + kij*GAMMA;

    // get parent price and do interpolation if necessary
    int parentReqDateIdx = child2ParentReqIdxs[fwd2ReqIdx];
    if (child2ParentLowerIdxs[j] == child2ParentUpperIdxs[j]) { // no interpolation
        F += parent->getExpectedPrice(parentReqDateIdx, child2ParentLowerIdxs[j]);
    }
    else { // need interpolation
        F += parent->getExpectedPrice(parentReqDateIdx, child2ParentLowerIdxs[j]) * parentLowerInterpRatios[j] +
            parent->getExpectedPrice(parentReqDateIdx, child2ParentUpperIdxs[j]) * parentUpperInterpRatios[j];
    }

    return F;
}

void SRMEnergyTier2Diffuse::addAggregatedDates(
    const DateTimeArray& spot,
    const DateTimeArray& measurementDates,
    const DateTimeArray& resetDates)
{
    IQMCDiffusibleEnergy::addAggregatedDates(spot, measurementDates, resetDates);
    parent->addAggregatedDates(DateTimeArray(), measurementDates, DateTimeArray()); // add measurement dates to parent
    
    updateMaxMaturity(spot, measurementDates, resetDates);
}



/** Create the appropriate tier 1 energy model */
QMCEnergyDiffuseSP SRMEnergyDiffuseCreate(
                                     string modelType,
                                     QMCRatesDiffuseSP srmRatesDiffuse)
{
    const string & method = "SRMEnergyDiffuseCreate";
    QMCEnergyDiffuseSP pAssetEnrg;

    if (CString::equalsIgnoreCase(modelType, EnrgOilStr)) {
        pAssetEnrg = QMCEnergyDiffuseSP(new SRMEnergyOilDiffuse(srmRatesDiffuse));
    }
    else if (CString::equalsIgnoreCase(modelType, EnrgGasStr)) {
        throw ModelException(method, "Energy gas model is not available!");
    }
    else if (CString::equalsIgnoreCase(modelType, EnrgPowStr)) {
        throw ModelException(method, "Energy power model is not available!");
    }
    else
        throw ModelException(method, "Invalid energy model type: " + modelType);

    return pAssetEnrg;
}

DRLIB_END_NAMESPACE
