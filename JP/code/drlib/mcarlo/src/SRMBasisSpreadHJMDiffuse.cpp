//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMBasisSpreadHJMDiffuse.hpp"
#include "edginc/SRMBasisSpreadHJMUtil.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/SVQmcImplemented.hpp"
#include "edginc/IQMCRNGManager.hpp"

#include <algorithm>

USING_DRLIB_NAMESPACE

void SRMBasisSpreadHJMDiffuse::setSRMBasisSpreadHJMDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMBasisSpreadUtilSP   _srmSPUtil,
    double                 _NbSigmasMax,
    double                 _NbSigmasMin,
    const vector<double>&  _corrSPIR,
    double                 _spFxCorr )
{
    srmSPUtil=_srmSPUtil;
    //basisIndexCurve = srmSPUtil->getBasisCurve();

    qSpread=srmSPUtil->getQ();
    randomIndex=_randomIndex;
    spFxCorr=_spFxCorr;
    NbSigmasMax=_NbSigmasMax;
    NbSigmasMin=_NbSigmasMin;
    today=_today;
    zeroQ = Maths::isZero(qSpread);

    srmIRUtil=irAsset->getRatesHJMUtil(); // must be available for the type we use
    vector<double> alphas(3);
    srmIRUtil->getAlpha(alphas);

    setSpIrCorr(_corrSPIR, alphas);

}


/** finalize the timelines, allocate necessary memory */
void SRMBasisSpreadHJMDiffuse::finalize(DateTimeArrayConstSP pallDates)
{
    static const string method("SRMBasisSpreadHJMDiffuse::finalize");

    const DateTimeArray& allDates = *pallDates;
    diffusionDates = DateTime::merge(
                            DateTimeArray(1, today),
                            today.getFutureDates(allDates));

    todayIdx = today.find(allDates);

    // InitDiffusionDates(allDates); Already done in derived class
    const int numSimDates = diffusionDates.size();

    DateTimeArray requestedDates = getForwardDates();
    DateTimeArray forwardDates = getForwardForwardDates();

    // translation of ESDF dates
    assert(DateTime::isSubset(allDates, requestedDates));
    expSpreadIndexes = DateTime::getIndexes(allDates, requestedDates); // Nb of requested dates


    // Projection of union(t_i, T_j) onto {t_i}
    assert(DateTime::isSubset(forwardDates, requestedDates));
    fwdIdx2RequestedIdx = DateTime::getProjection(forwardDates, requestedDates); // created map (projection)  esdfForwardDates->esdfRequestedDates

    // aux mapping to map between CR fwd dates to IR fwd dates
    fwdIdx2IRedfIdx.resize(forwardDates.size());
    fwdPlusIdx2IRedfIdx.resize(forwardDates.size());

    // TODO:  The array size of fwdPlusIdx2IRedfIdx reflects both
    // resetDates and measurement dates.  However, -1 is returned
    // by irAsset when measurementDate+liborOffset is not found.
    for(int i=0; i < forwardDates.size(); ++i) 
    {
        fwdIdx2IRedfIdx[i] = irAsset->getForwardForwardIndex(forwardDates[i]);
        fwdPlusIdx2IRedfIdx[i] = irAsset->getForwardForwardIndex(basisIndexCurve->getRefPaymentDate(forwardDates[i]));
    }

    // make life easier - add request for index off the end
    expSpreadIndexes.push_back(todayIdx+numSimDates+1);

    if (!zeroQ) // only need this one if not purely normal diffusion
        spBar = srmSPUtil->getSpotSpreads();

    spVol = srmSPUtil->getSpotVols();

    srmSPUtil->calcEffRateLimit(NbSigmasMax, NbSigmasMin,
                               MaxEffRateSP, MinEffRateSP);
    sqrtYearFrac = srmSPUtil->computeSqrtYearFrac();
    partialHFactor = srmSPUtil->computePartialH(today, forwardDates);
    fwdSpread = srmSPUtil->computeFwdParSpreads(forwardDates);
    if (!zeroQ) // lognormal works only if all the fwd spreads are strictly positive
        if (! Maths::isPositive(*min_element(fwdSpread.begin(), fwdSpread.end())))
            throw ModelException(method, (string)"The basis index spread " +
                basisIndexCurve->getName()+ " has negative values but the requested diffusion is not purely normal.");

}

/** get all the dates important for asset */
//DateTimeArray SRMBasisSpreadHJMDiffuse::getAssetDates()
// {
//     DateTime::doSortUniq(requestedDates);
//     return requestedDates;
// }

/** Getting the integer index corresponding to a given measurement
    date (md) or a future date (fd) for accessing Expected SDF. */
//FwdIdx SRMBasisSpreadHJMDiffuse::getDateIndex(const DateTime& date)
// {
//     FwdIdx idx = date.find(forwardDates);
//     assert(forwardDates[idx] == date);
//     return idx;
// }

/** This is a declaration that ExpDF(md, fd[i]) should be produced
    by diffusion and so it might be requested later. */
//void SRMBasisSpreadHJMDiffuse::addExpectedFwdSpreadDates(
//     const DateTime& measurementDate,
//     const DateTimeArray& resetDates)
// {
//     requestedDates.push_back(measurementDate);
//     forwardDatesSet.insert(measurementDate);
//     forwardDatesSet.insert(resetDates.begin(), resetDates.end());
//     DateTimeArray payDates( resetDates.size() );
//     for ( int i = 0; i < resetDates.size(); ++i ) {
//         payDates[i] = basisIndexCurve->getRefPaymentDate( resetDates[i] );
//     }
//     irAsset->addExpectedDFDates(measurementDate, payDates);
// }

void SRMBasisSpreadHJMDiffuse::addAggregatedDates(
                const DateTimeArray& spot,
                const DateTimeArray& measurementDates,
                const DateTimeArray& resetDates)
{
    IQMCDiffusibleBasisIndex::addAggregatedDates(spot, measurementDates, resetDates);
    DateTimeArray payDates( resetDates.size() );
    for ( int i = 0; i < resetDates.size(); ++i ) {
        payDates[i] = basisIndexCurve->getRefPaymentDate( resetDates[i] );
    }
    irAsset->addAggregatedDates(DateTimeArray(), measurementDates, payDates);
    
    updateMaxMaturity(spot, measurementDates, resetDates);
}

//// constructor
SRMBasisSpreadHJMDiffuse::SRMBasisSpreadHJMDiffuse(
    SRMRatesHJMDiffuseSP srmIRBase,IBasisIndexCurveConstSP _basisIndexCurve) :
        // meaningful initialization
        irAsset(srmIRBase),
        basisIndexCurve(_basisIndexCurve),
        sigmaR(srmIRBase->getSigmaR()),
        sigmaFX(srmIRBase->getSigmaFX()),
        // default initialization with 'bad' values
        randomIndex(-999),
        qSpread(-999.),
        spFxCorr(-999.),
        zeroQ(false),
        srmIRUtil(NULL),
        srmSPUtil(NULL),
        todayIdx(-999),
        NbSigmasMax(-999.),
        NbSigmasMin(-999.),
        diffBound(srmIRBase->getDiffusionBound())
{}



///////////  Spread on top of a 1-factor IR  ///////////////////////////

void SRMBasisSpreadHJM1F::finalize(DateTimeArrayConstSP allDates)
{
    static const string method("SRMBasisSpreadHJM1F::finalize");
    try {

        SRMBasisSpreadHJMDiffuse::finalize(allDates); // will notify srmSPUtil about new timeline

        const int numSimDates = getNumSimDates();
        expSP.resize(getNbRequestedDates());

        const DateTimeArray& dates = diffusionDates; // same as diffusionDates ?
        kghFactors.resize(dates.size()-1);
        assert(dates.size() == numSimDates);

        for (int i = 0; i < numSimDates - 1; i++) {
            kghFactors[i].hFactor = srmSPUtil->hFactor(dates[i], dates[i+1]);
            kghFactors[i].kFactor0IR = srmIRUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kghFactors[i].gFactor0IR = srmIRUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            if (!zeroQ)
                kghFactors[i].gFactor0TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        0 /* nb IR factor */);
            else
                kghFactors[i].gFactor0TIR = 0;
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }

}

/** generate path across all dates ...
    essentially spdiffuse::DiffuseSP_1F plus a part of spdiffuse::CalcSigmaLSP */
void SRMBasisSpreadHJM1F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = getNumSimDates()-1; // Nb of intervals

    int expSpreadPos = expSpreadIndexes.front() == 0? 1: 0;
    int stopIdx = expSpreadIndexes[expSpreadPos]; // date when we save expected survival probability

    double OMEGA = 0;
    double CHI0  = 0;

    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){
        double Z = randoms[i];
        double del_t = sqrtYearFrac[i]*sqrtYearFrac[i];
        double RootDelT = sqrtYearFrac[i];

        /* Retrieve precomputed k-factors & g-factors
           (precomputation takes place in constructor of SP1F) */
        double k0 = kghFactors[i].kFactor0IR;
        double g0 = kghFactors[i].gFactor0IR; // has g(T_i, T_{i+1}
        double h = kghFactors[i].hFactor;


        /* main calculations */

        double S0   = rho0 * alpha0 * k0 * del_t;
        double sumV = rho0 * alpha0 * g0 * del_t;

        /* Step1: update spread vol */
        /*--------------------------*/


        double effSP =  1.0;
        if (!zeroQ)
        {
            // h=1 in this formula, as we take h(T,T)
            // g is for g(T, T+1Libor)
            double g0T =kghFactors[i].gFactor0TIR; //has g(T_i, T_i + 1LiborTenor)
            double spNow = spBar[i] + g0T*CHI0 + OMEGA;

            effSP = 1.0 - qSpread + qSpread * spNow / spBar[i];
            effSP = COLLAR(effSP, MaxEffRateSP[i], MinEffRateSP[i]);
        }

        double sigmaSp = effSP * spVol[i];
        double sigmaIrSp = sigmaR[i] * sigmaSp;

        /* Step2: update OMEGA */
        /*---------------------*/

        OMEGA = h*(OMEGA+g0*CHI0 +
                   sigmaIrSp * sumV +
                   sigmaSp*RootDelT*Z);

        /* Step3: update CHI */
        /*-------------------*/

        CHI0 = h*(k0*CHI0 + sigmaIrSp*S0);

        /* Step4: add the CUPS component if necessary */
        /*--------------------------------------------*/

        if (sigmaFX != NULL )
        {
            OMEGA -= spFxCorr * (*sigmaFX)[i] * sigmaSp * h * del_t;
        }

        i++; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (i + todayIdx == stopIdx) {
                expSP[expSpreadPos].OMEGA = OMEGA;
                expSP[expSpreadPos].CHI0 = CHI0;
                expSpreadPos++;
                stopIdx= expSpreadIndexes[expSpreadPos];
        }
    } // for each date i
}


/** Accessing the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMBasisSpreadHJM1F::getExpectedBasisFwdSpread(
    FwdIdx measurementDateIdx,  // Idx of measDate
    FwdIdx resetDatesIdx )
{
    double FwdSp = getFwdSpread( resetDatesIdx );

    /* calc other deterministic quantities */

    double g0 = static_cast<SRMRatesHJM1F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            0);

    double h = partialHFactor[resetDatesIdx]/partialHFactor[measurementDateIdx];

    /* main calculations */
    int expSpreadPos = fwdIdx2RequestedIdx[measurementDateIdx];
    double OMEGA = expSP[expSpreadPos].OMEGA;
    double CHI0  = expSP[expSpreadPos].CHI0;

    double spread = FwdSp + h*(g0*CHI0 + OMEGA);
    return spread;

}




///////  Spread on top of a 2-factor IR  /////////////////

void SRMBasisSpreadHJM2F::finalize(DateTimeArrayConstSP allDates)
{
    static const string method("SRMBasisSpreadHJM1F::finalize");
    try {

        SRMBasisSpreadHJMDiffuse::finalize(allDates); // will notify srmSPUtil about new timeline

        const int numSimDates = getNumSimDates();
        expSP.resize(getNbRequestedDates());

        const DateTimeArray& dates = diffusionDates; // same as diffusionDates ?
        kghFactors.resize(dates.size()-1);
        assert(dates.size() == numSimDates);

        for (int i = 0; i < numSimDates - 1; ++i) {
            kghFactors[i].hFactor = srmSPUtil->hFactor(dates[i], dates[i+1]);
            kghFactors[i].kFactor0IR = srmIRUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kghFactors[i].gFactor0IR = srmIRUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kghFactors[i].kFactor1IR = srmIRUtil->kFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
            kghFactors[i].gFactor1IR = srmIRUtil->gFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
            if (!zeroQ)
            {
                kghFactors[i].gFactor0TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        0 /* nb IR factor */);
                kghFactors[i].gFactor1TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        1 /* nb IR factor */);
            }
            else
                kghFactors[i].gFactor0TIR = kghFactors[i].gFactor1TIR = 0;
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }

}

/** generate path across all dates ...
    essentially spdiffuse::DiffuseSP_2F plus a part of spdiffuse::CalcSigmaLSP */
void SRMBasisSpreadHJM2F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = getNumSimDates()-1; // Nb of intervals

    int expSpreadPos = expSpreadIndexes.front() == 0? 1: 0;
    int stopIdx = expSpreadIndexes[expSpreadPos]; // date when we save expected survival probability

    double OMEGA = 0;
    double CHI0  = 0;
    double CHI1  = 0;

    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){
        double Z = randoms[i];
        double del_t = sqrtYearFrac[i]*sqrtYearFrac[i];
        double RootDelT = sqrtYearFrac[i];

        /* Retrieve precomputed k-factors & g-factors
           (precomputation takes place in constructor of SP1F) */
        double k0 = kghFactors[i].kFactor0IR;
        double g0 = kghFactors[i].gFactor0IR; // has g(T_i, T_{i+1}
        double k1 = kghFactors[i].kFactor1IR;
        double g1 = kghFactors[i].gFactor1IR; // has g(T_i, T_{i+1}
        double h = kghFactors[i].hFactor;


        /* main calculations */

        double S0   = rho0 * alpha0 * k0 * del_t;
        double S1   = rho1 * alpha1 * k1 * del_t;
        double sumV = rho0 * alpha0 * g0 * del_t + rho1 * alpha1 * g1 * del_t;

        /* Step1: update spread vol */
        /*--------------------------*/


        double effSP =  1.0;
        if (!zeroQ)
        {
            // h=1 in this formula, as we take h(T,T)
            // g is for g(T, T+1Libor)
            double g0T =kghFactors[i].gFactor0TIR; //has g(T_i, T_i + 1LiborTenor)
            double g1T =kghFactors[i].gFactor1TIR; //has g(T_i, T_i + 1LiborTenor)
            double spNow = spBar[i] + g0T*CHI0 + g1T*CHI1 + OMEGA;

            effSP = 1.0 - qSpread + qSpread * spNow / spBar[i];
            effSP = COLLAR(effSP, MaxEffRateSP[i], MinEffRateSP[i]);
        }

        double sigmaSp = effSP * spVol[i];
        double sigmaIrSp = sigmaR[i] * sigmaSp;

        /* Step2: update OMEGA */
        /*---------------------*/

        OMEGA = h*(OMEGA + g0*CHI0 + g1*CHI1 +
                   sigmaIrSp * sumV +
                   sigmaSp*RootDelT*Z);

        /* Step3: update CHI */
        /*-------------------*/

        CHI0 = h*(k0*CHI0 + sigmaIrSp*S0);
        CHI1 = h*(k1*CHI1 + sigmaIrSp*S1);

        /* Step4: add the CUPS component if necessary */
        /*--------------------------------------------*/

        if (sigmaFX != NULL )
        {
            OMEGA -= spFxCorr * (*sigmaFX)[i] * sigmaSp * h * del_t;
        }

        ++i; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (i + todayIdx == stopIdx) {
                expSP[expSpreadPos].OMEGA = OMEGA;
                expSP[expSpreadPos].CHI0 = CHI0;
                expSP[expSpreadPos].CHI1 = CHI1;
                expSpreadPos++;
                stopIdx= expSpreadIndexes[expSpreadPos];
        }
    } // for each date i
}


/** Accessing the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMBasisSpreadHJM2F::getExpectedBasisFwdSpread(
    FwdIdx measurementDateIdx,  // Idx of measDate
    FwdIdx resetDatesIdx )
{
    double FwdSp = getFwdSpread( resetDatesIdx );

    /* calc other deterministic quantities */

    double g0 = static_cast<SRMRatesHJM2F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            0);

    double g1 = static_cast<SRMRatesHJM2F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            1);

    double h = partialHFactor[resetDatesIdx]/partialHFactor[measurementDateIdx];

    /* main calculations */
    int expSpreadPos = fwdIdx2RequestedIdx[measurementDateIdx];
    double OMEGA = expSP[expSpreadPos].OMEGA;
    double CHI0  = expSP[expSpreadPos].CHI0;
    double CHI1  = expSP[expSpreadPos].CHI1;

    double spread = FwdSp + h*(g0*CHI0 + g1*CHI1 + OMEGA);
    return spread;

}





///////  Spread on top of a 3-factor IR  /////////////////

void SRMBasisSpreadHJM3F::finalize(DateTimeArrayConstSP allDates)
{
    static const string method("SRMBasisSpreadHJM1F::finalize");
    try {

        SRMBasisSpreadHJMDiffuse::finalize(allDates); // will notify srmSPUtil about new timeline

        const int numSimDates = getNumSimDates();
        expSP.resize(getNbRequestedDates());

        const DateTimeArray& dates = diffusionDates; // same as diffusionDates ?
        kghFactors.resize(dates.size()-1);
        assert(dates.size() == numSimDates);

        for (int i = 0; i < numSimDates - 1; ++i) {
            kghFactors[i].hFactor = srmSPUtil->hFactor(dates[i], dates[i+1]);
            kghFactors[i].kFactor0IR = srmIRUtil->kFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kghFactors[i].gFactor0IR = srmIRUtil->gFactor(dates[i], dates[i+1], 0 /* nb IR factor */);
            kghFactors[i].kFactor1IR = srmIRUtil->kFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
            kghFactors[i].gFactor1IR = srmIRUtil->gFactor(dates[i], dates[i+1], 1 /* nb IR factor */);
            kghFactors[i].kFactor2IR = srmIRUtil->kFactor(dates[i], dates[i+1], 2 /* nb IR factor */);
            kghFactors[i].gFactor2IR = srmIRUtil->gFactor(dates[i], dates[i+1], 2 /* nb IR factor */);
            if (!zeroQ)
            {
                kghFactors[i].gFactor0TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        0 /* nb IR factor */);
                kghFactors[i].gFactor1TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        1 /* nb IR factor */);
                kghFactors[i].gFactor2TIR = srmIRUtil->gFactor(
                        dates[i],
                        basisIndexCurve->getRefPaymentDate( dates[i] ),
                        2 /* nb IR factor */);
            }
            else
                kghFactors[i].gFactor0TIR =
                kghFactors[i].gFactor1TIR =
                kghFactors[i].gFactor2TIR = 0;
        }
    } catch (exception& e){
        throw ModelException(e, method);
    }

}

/** generate path across all dates ...
    essentially spdiffuse::DiffuseSP_3F plus a part of spdiffuse::CalcSigmaLSP */
void SRMBasisSpreadHJM3F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    int numDates = getNumSimDates()-1; // Nb of intervals

    int expSpreadPos = expSpreadIndexes.front() == 0? 1: 0;
    int stopIdx = expSpreadIndexes[expSpreadPos]; // date when we save expected survival probability

    double OMEGA = 0;
    double CHI0  = 0;
    double CHI1  = 0;
    double CHI2  = 0;

    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){
        double Z = randoms[i];
        double del_t = sqrtYearFrac[i]*sqrtYearFrac[i];
        double RootDelT = sqrtYearFrac[i];

        /* Retrieve precomputed k-factors & g-factors
           (precomputation takes place in constructor of SP1F) */
        double k0 = kghFactors[i].kFactor0IR;
        double g0 = kghFactors[i].gFactor0IR; // has g(T_i, T_{i+1}
        double k1 = kghFactors[i].kFactor1IR;
        double g1 = kghFactors[i].gFactor1IR; // has g(T_i, T_{i+1}
        double k2 = kghFactors[i].kFactor2IR;
        double g2 = kghFactors[i].gFactor2IR; // has g(T_i, T_{i+1}
        double h = kghFactors[i].hFactor;


        /* main calculations */

        double S0   = rho0 * alpha0 * k0 * del_t;
        double S1   = rho1 * alpha1 * k1 * del_t;
        double S2   = rho2 * alpha2 * k2 * del_t;
        double sumV =
                        rho0 * alpha0 * g0 * del_t +
                        rho1 * alpha1 * g1 * del_t +
                        rho2 * alpha2 * g2 * del_t;

        /* Step1: update spread vol */
        /*--------------------------*/


        double effSP =  1.0;
        if (!zeroQ)
        {
            // h=1 in this formula, as we take h(T,T)
            // g is for g(T, T+1Libor)
            double g0T =kghFactors[i].gFactor0TIR; //has g(T_i, T_i + 1LiborTenor)
            double g1T =kghFactors[i].gFactor1TIR; //has g(T_i, T_i + 1LiborTenor)
            double g2T =kghFactors[i].gFactor2TIR; //has g(T_i, T_i + 1LiborTenor)
            double spNow = spBar[i] + g0T*CHI0 + g1T*CHI1 + g2T*CHI2 + OMEGA;

            effSP = 1.0 - qSpread + qSpread * spNow / spBar[i];
            effSP = COLLAR(effSP, MaxEffRateSP[i], MinEffRateSP[i]);
        }

        double sigmaSp = effSP * spVol[i];
        double sigmaIrSp = sigmaR[i] * sigmaSp;

        /* Step2: update OMEGA */
        /*---------------------*/

        OMEGA = h*(OMEGA + g0*CHI0 + g1*CHI1 + g2*CHI2 +
                   sigmaIrSp * sumV +
                   sigmaSp*RootDelT*Z);

        /* Step3: update CHI */
        /*-------------------*/

        CHI0 = h*(k0*CHI0 + sigmaIrSp*S0);
        CHI1 = h*(k1*CHI1 + sigmaIrSp*S1);
        CHI1 = h*(k2*CHI2 + sigmaIrSp*S2);

        /* Step4: add the CUPS component if necessary */
        /*--------------------------------------------*/

        if (sigmaFX != NULL )
        {
            OMEGA -= spFxCorr * (*sigmaFX)[i] * sigmaSp * h * del_t;
        }

        ++i; /* up the loop counter ...
             this comes back to the fact with eg 5 dates and 4 diffusion steps; it's easier to identify the index
             of the date you want rather than the index of the following one */
        if (i + todayIdx == stopIdx) {
                expSP[expSpreadPos].OMEGA = OMEGA;
                expSP[expSpreadPos].CHI0 = CHI0;
                expSP[expSpreadPos].CHI1 = CHI1;
                expSP[expSpreadPos].CHI2 = CHI2;
                expSpreadPos++;
                stopIdx= expSpreadIndexes[expSpreadPos];
        }
    } // for each date i
}


/** Accessing the expected value ExpSDF(md, fd) where md is a
    simulated measurement date and fd is some future date after the
    measurement is made. */
double SRMBasisSpreadHJM3F::getExpectedBasisFwdSpread(
    FwdIdx measurementDateIdx,  // Idx of measDate
    FwdIdx resetDatesIdx )
{
    double FwdSp = getFwdSpread( resetDatesIdx );

    /* calc other deterministic quantities */

    double g0 = static_cast<SRMRatesHJM3F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            0);

    double g1 = static_cast<SRMRatesHJM3F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            1);

    double g2 = static_cast<SRMRatesHJM3F*>(irAsset.get())->
        getGFactor(
            fwdIdx2IRedfIdx[measurementDateIdx],
            fwdPlusIdx2IRedfIdx[resetDatesIdx], // should be idx of j+tenor of libor
            2);


    double h = partialHFactor[resetDatesIdx]/partialHFactor[measurementDateIdx];

    /* main calculations */
    int expSpreadPos = fwdIdx2RequestedIdx[measurementDateIdx];
    double OMEGA = expSP[expSpreadPos].OMEGA;
    double CHI0  = expSP[expSpreadPos].CHI0;
    double CHI1  = expSP[expSpreadPos].CHI1;
    double CHI2  = expSP[expSpreadPos].CHI2;

    double spread = FwdSp + h*(g0*CHI0 + g1*CHI1 + g2*CHI2 + OMEGA);
    return spread;

}


