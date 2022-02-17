//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMEquityDiffuse.cpp
//
//   Description : A generator of paths using stochastic rates
//                 for equity Assets
//
//   Date        : 13 Aug 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DividendCollector.hpp"

DRLIB_BEGIN_NAMESPACE  


const string SRMEquityDiffuse::CONSTANT_SPOT_VOL("CONSTANT_SPOT_VOL");
const string SRMEquityDiffuse::NO_FAILURE_ALLOWED("NO_FAILURE_ALLOWED");
const string SRMEquityDiffuse::USE_LAST_LEVEL("USE_LAST_LEVEL");


SRMEquityDiffuse::SRMEquityDiffuse(IQMCDiffusibleInterestRateSP _domIR) :
    // non-trivial initialization
    // trivial ini (with 'bad' values)
    randomIndex(-999),
    randoms(NULL),
    LnE(-999.),
    todayIdx(-999),
    stopIdx(-999),
    spotEqPos(-999),
    expSpotEqPos(-999),
    ccyTreatment(ccyVanilla),
    corrEqFX(-999.),
    spotEqStart(-999),
    origLnE(-999.),
    divPos(-999),
    domDiscYCIdx(size_t(-1)),
    initialized(false),
    diffBound(_domIR.get() ? _domIR->getDiffusionBound() : NULL),
    smileIsFlatEverywhere(true),
    spotFXFullPath(NULL)
{
    domIR = QMCRatesDiffuseSP(dynamic_cast<QMCRatesDiffuse*>(_domIR.get()));
    irLnMONEY = &(_domIR->getDomLnMONEY());
    sigmaFX = _domIR->getSigmaFX();
        
    //spotFXFullPath set in setSRMEquityDiffuse
}


/* SRMEquityDiffuse is constructed in several steps: 
 * 1) Construct SRMEquityDiffuse                [ see SRMGenDiffusibleEquity::initializeDiffusibleAssets() ]
 * 2) Call setSRMEquityDiffuse (this method)    [ see SRMGenDiffusibleEquity::setDiffusibleAsset() ]
 * 3) Call finalize                             [ see MCPathConfigSRMGen::MCPathConfigSRMGen() ]
 */
void  SRMEquityDiffuse::setSRMEquityDiffuse(
                        int                   _randomIndex,
                        const DateTime &      _today,
                        SRMEquityUtilSP       _equityUtil,
                        const vector<double>& _past,
                        int _timePointsPerYear)
{
    timePointsPerYear = _timePointsPerYear;
    randomIndex = _randomIndex;
    today = _today;
    equityUtil = _equityUtil;
    spotEq = _past;
    spotEqStart = spotEq.size(); // will be re-checked in finalize
    ccyTreatment = equityUtil->getCurrencyTreatment();
    if (ccyTreatment == ccyStruck) {
        spotFXFullPath = domIR->getSpotFXFullPath();
    }
    domDiscYCIdx = domIR->getDiscYCIdx(); // get Idx of disc YC
} 

/** finalize the timelines, allocate necessary memory */
void SRMEquityDiffuse::finalize(DateTimeArrayConstSP allDatesSP) // FIXME
{
    static string method = (char *)__FUNCTION__;
    DateTimeArray spotRequestedDates = getSpotDates();
    DateTimeArray expSpotRequestedDates = getForwardDates();
    DateTimeArray expSpotForwardDates = getForwardForwardDates();

    QLIB_VERIFY(DateTime::isSubset(expSpotForwardDates, expSpotRequestedDates), "expSpotForwardDates must contain expSpotRequestedDates");

    size_t i;
    const size_t Nall = SRMUtil::getNumSimDates(today, *allDatesSP);

    const size_t Nefx = expSpotRequestedDates.size();
    const size_t Nefwd = expSpotForwardDates.size();
    const size_t Nfx  =    spotRequestedDates.size();

    DateTimeArray   diffusionDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);

    QLIB_VERIFY(diffusionDates.size() == Nall,
        Format::toString("diffusionDates has size %l, should have size %l", diffusionDates.size(), Nall));

    yearFrac.assign(Nall-1, 0.0);
    sqrtYearFrac.assign(Nall-1, 0.0);

    for (i = 0; i < yearFrac.size(); i++){
        yearFrac[i] = diffusionDates[i].yearFrac(diffusionDates[i+1]);
        sqrtYearFrac[i] = sqrt(yearFrac[i]);
    }
    todayIdx = today.find((*allDatesSP));

    spotEqStart = today.findUpper(spotRequestedDates);
    // skip today as diffusion results are always in future
    if (spotEqStart >= 0 && spotEqStart != spotRequestedDates.size() &&
        spotRequestedDates[spotEqStart] == today)
            ++spotEqStart;

    const DateTimeArray& simDates = equityUtil->getBaseIR()->getSimDates();

    QLIB_VERIFY(simDates.size() == Nall,
        Format::toString("simDates has size %l, should have size %l", diffusionDates.size(), Nall));
    QLIB_VERIFY(DateTime::isSubset(simDates, diffusionDates),
        "simDates should contain diffusionDates");
    QLIB_VERIFY(DateTime::isSubset((*allDatesSP), spotRequestedDates),
        "allDates should contain spotRequestedDates");

    spotEqIndexes = DateTime::getIndexes((*allDatesSP), spotRequestedDates);
    spotEqIndexes.push_back(todayIdx+simDates.size()+1);//make life easier

    QLIB_VERIFY(DateTime::isSubset((*allDatesSP), expSpotRequestedDates),
        "allDates should contain expSpotRequestedDates");
    expSpotEqIndexes = DateTime::getIndexes((*allDatesSP), expSpotRequestedDates);
    expSpotEqIndexes.push_back(todayIdx+simDates.size()+1); //make life easier

//  vector<int> expSpotForwardDateIndexes = DateTime::getIndexes((*allDatesSP), expSpotForwardDates);
//  fwd2requested = DateTime::getProjection(expSpotForwardDates, expSpotRequestedDates);

    fwd2DomIRFwd.resize(Nefwd);
    for(i=0; i < Nefwd; ++i)
    {
        fwd2DomIRFwd[i] = domIR->getForwardForwardIndex(expSpotForwardDates[i]);
    }

    // spotEQ contains history EQ rates only. So must resize
    spotEq.resize(Nfx);
    expSpotEq.resize(Nefx);

    SpotVol = equityUtil->getSpotVol();
    sigmaEQ.resize(simDates.size()-1);
    
    // calculate sigmaEQ.front() - in fact it's just SpotVol[0]
    sigmaEQ[0] = SpotVol[0];

    const DoubleArray& tmpLnFwdEQ = equityUtil->getLnFwdEQ();
    LnFwdEq = vector<double>(tmpLnFwdEQ.begin(), tmpLnFwdEQ.end());
   
    EqSmile_a1 = equityUtil->getEQSmile_a1();
    EqSmile_a2 = equityUtil->getEQSmile_a2();
    EqSmile_a3 = equityUtil->getEQSmile_a3();

    smileIsFlatEverywhere = true;
    for (i = 0;i<EqSmile_a1.size();i++)
    {
        if (EqSmile_a1[i] != 0.0 || EqSmile_a2[i] != 0.0)
            smileIsFlatEverywhere = false;
    }

    origLnE = equityUtil->getLogSpotEQ();

    if (equityUtil->getMomentMatchingFlag())
    {  // for moment matching: initialize originalLnFwdSpot and originalLnExpFwdSpot
       // I think LnFwdEq can be used
        // assume all dates are > valuedate
        equityUtil->computeOriginalSpotPrices(spotRequestedDates, originalLnFwdSpot);
        equityUtil->computeOriginalSpotPrices(expSpotForwardDates, originalLnExpFwdSpot);
    }
   
    // Dividends - need to know each div date so can stop and pay it
    // This is common across all SVs but is per asset (implicit - we have one
    // of these per asset). CHECK!!
    exDivIndexes = DateTime::getIndexes(simDates, equityUtil->getDivExDates());
    exDivIndexes.push_back(simDates.size()+1);//make life easier - allows safe reference after the last (which will do nothing)

    // This may contain continuous divs so treat carefully
    // Don't simply remove them since the treatment elsewhere for
    // ex-dates may be affected. Actually expect at most cts divs
    // after all discrete divs so shouldn't really be an issue.
    DividendListSP divList(equityUtil->getDivs());
    const DividendArray& divArray = divList->getArray();

    lnDivYield.resize(divArray.size());
    for(int j=0; j<divArray.size(); j++) {
        double amt = divArray[j].getDivAmount();
        if (!Maths::isPositive(1.0-amt)) {
            throw ModelException(method,
                                 "Bad divAmount!");
        }
        Dividend::TDividendType divType = divArray[j].getDivType();
        if (divType == Dividend::PERCENT) {
            lnDivYield[j] = log(1.0-amt);
        } else if (divType == Dividend::CONTINUOUS) {
            lnDivYield[j] = 0.0; // cts div handled elsewhere
        } else {
            throw ModelException(method,
                                 "Unsupported dividend type detected!");
        }
    }

    
    // aka PopulateCfactorArrayEQ
    // CfactorArrayEQ[i] is the contribution from simDates[j] to simDates[j+1]
    // Implements borrow and cts divs
    CfactorArrayEQ = vector<double>(yearFrac.size(),0.0);
    BorrowCurveConstSP bc = equityUtil->getBorrowCurve();
    DateTime firstContDivDate;
    bool hasCtsDiv = divList->hasContinuousDividend(firstContDivDate);
    for(size_t j=0; j<CfactorArrayEQ.size(); j++) {
        if ( bc.get() && !bc->isZero() ) {
            CfactorArrayEQ[j] = -bc->continuousFwdRate(simDates[j], simDates[j+1]) * yearFrac[j];
        }
        // simDates include all ex div dates so we can be confident of not skipping
        // a cts div. Hence each interval between simDates will have a well-defined
        // single cts div rate (if any). This arrangement is meant to be reasonably efficient
        // for the likely case (single cts div at the end of all other discrete divs
        // should mean this clause entered rarely).
        if (hasCtsDiv && firstContDivDate <= simDates[j]) {
            double ctsDiv = 0.0;
            hasCtsDiv = divList->ctsDivPayingAtDate(simDates[j], ctsDiv);
            // note this is now a continuously compounded rate
            CfactorArrayEQ[j] -= ctsDiv * yearFrac[j];
        }
    }

    // Ccy prot adjustment : eq/fx corr
    corrEqFX = equityUtil->getInpCorrEqFX();

    initialized = true;
    equityUtil = SRMEquityUtilSP(0);
}

// Convenient in several places - may have a more natural place
void SRMEquityDiffuse::getDivData(AssetConstSP       asset,
                       const DateTime&    start,
                       const DateTime&    end,
                       DateTimeArray&     divExDates, // time adjusted to fit timeline (->SOD)
                       bool               withDivList, // if withDivList must supply yc too
                       DividendListSP&    divList, // optional
                       IYieldCurveConstSP yc)      // optional
{
    static const string method("SRMEquityDiffuse::getDivData");
    // A rather ugly piece of validation - the code below does not capture
    // the case where a continuous dividend is before the start, so here
    // we capture all dividends and simply forbid that case. It is not
    // a practical situation in any case.
    DateTime beforeAllTime(0,0);
    DividendListSP checkDivList =
        DividendCollector::divsBetweenDates(asset.get(),
                                            beforeAllTime,
                                            beforeAllTime,
                                            end,
                                            DividendCollector::NONE);
    DateTime firstContDivDate;
    if (checkDivList->hasContinuousDividend(firstContDivDate) &&
        firstContDivDate <= start) {
            throw ModelException(method,
                                 "Continuous dividends paying before start are not supported");
    }

    // Extract discrete divs. We ultimately need yield divs but
    // do not exploit the DividendCollector conversion here since we subsequently
    // merge duplicate divs on the same date. If we converted to yields too early
    // the merge would introduce an error (a methodology issue - see comment
    // in DividendList::aggregate).
    DividendListSP localDivList =
        DividendCollector::divsBetweenDates(asset.get(),
                                            start,
                                            start,
                                            end,
                                            DividendCollector::NONE);
    if (withDivList) {
        // need a yc
        if (!yc.get()) {
            throw ModelException(method, "Internal error!");
        }
        // With a YieldCurve we can build an aggregated div list with payment
        // dates accounted for by a PV from payment date to ex-date
        // done deterministically using provided yc.
        // Note that this does more than aggregate - it transforms divs to
        // be ones without payment delay (even unrepeated divs).
        // It does not however aggregate yield & dollar divs on the same date, hence
        // the later call (see below).
        divList = DividendListSP(DividendList::aggregate(*localDivList.get(),
                                                         yc.get()));
        // While we ultimately wish to have yield divs, the aggregation is not
        // consistent with aggregation of dollar. So, through some hoops first...
        divList->convertYieldDivs(asset.get()); // all into dollar first
        // Once converted we try the merging again.
        divList = DividendListSP(DividendList::aggregate(*divList.get(),
                                                         yc.get()));
        // And finally we consider it safe to convert to yield divs
        divList->convertDollarDivs(asset.get());
        divExDates = *divList->getExDivDates();
    } else {
        // Otherwise we're just interested in the dates,
        // so we combine duplicates manually
        divExDates = *localDivList->getExDivDates();
        DateTime::removeDuplicates(divExDates,
                                   true/*ignoreTimeOfDay*/);
    }
    // Ex-div dates are forced (by EDR lib itself) to have EXDIV_TIME so need to
    // extract the date and use a more appropriate time. SOD seems to be
    // a "standard" within SRM3.
    int aTime = DateTime::START_OF_DAY_TIME;
    for(int i=0; i<divExDates.size(); i++) {
        divExDates[i] = DateTime(divExDates[i].getDate(),
                                 aTime);
    }
}


// Here mainly for visibility of SRMEQVol::TYPE
DateTimeArray SRMEquityDiffuse::getCriticalDates(AssetConstSP    asset,
                                      const DateTime& start,
                                      const DateTime& end)
{
    // Vol dates, and dividend dates
    DateTimeArray critDates(
        CriticalDateCollector::collectVolDates(
            asset, SRMEQVol::TYPE, start, end)); // smile
    DateTimeArray divDates;
    DividendListSP divList(0);
    IYieldCurveConstSP yc(0);
    getDivData(asset, start, end,
               divDates,
               false, // withDivList - only want dates here
               divList,
               yc);
    critDates = DateTime::merge(critDates, divDates);
    return critDates;
}

       /** get all the dates important for asset */
// DateTimeArray SRMEquityDiffuse::getAssetDates()
// {
//
//     DateTime::doSortUniq(spotRequestedDates);
//     DateTime::doSortUniq(expSpotRequestedDates);
//
//     // TODO: add critical dates retrieval
//     return DateTime::merge(spotRequestedDates, expSpotRequestedDates);
// }

    /** Getting an integer index corresponding to a given date
        for accessing SpotFX. */
// SpotIdx SRMEquityDiffuse::getSpotPriceIndex(const DateTime& date)
// {
//     SpotIdx idx =  date.find(spotRequestedDates);
//     assert(idx >= 0 && idx < spotRequestedDates.size());
//     assert(spotRequestedDates[idx] == date);
//     return idx;
// }

    /** This is a declaration that SpotFXDF(md[i]) should be produced by
        diffusion and so it can be requested later. */
// void    SRMEquityDiffuse::addSpotPriceDates(const DateTimeArray& measurementDates)
// {
//     spotRequestedDates.insert(spotRequestedDates.end(), measurementDates.begin(), measurementDates.end());
// }
//


    /** (3) Getting the simulated ExpectedFX between [measDate, futureDate] */

    /** Accessing the ExpectedFX(md, fd) where md is a
        simulated measurement date and fd is some future date after the
        measurement is made. */
double  SRMEquityDiffuse::getExpectedPrice(FwdIdx measurementDateIdx,
                                    FwdIdx futureDateIdx)
{
    return exp(getLnExpectedPrice(measurementDateIdx, futureDateIdx));
}

    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
double SRMEquityDiffuse::getLnExpectedPrice(FwdIdx i,
                        FwdIdx j)
{
    int    idx = getTimeLogic()->getReqEDFIdx(i); //fwd2requested[i];
    ASSERT(idx >= 0 && (size_t)idx < expSpotEq.size());
    
    double lnSpotEQ = expSpotEq[idx];
    double dom =  domIR->getLnExpectedDiscFactor(domDiscYCIdx, fwd2DomIRFwd[i], fwd2DomIRFwd[j]);
    double lnRes =  lnSpotEQ - dom;
#if defined(DEBUGSRM)
    fprintf(stderr, "SRMEquityDiffuse::getLnExpectedPrice(%d, %d): idx= %d lnSpotEQ= %g dom= %g res= %g\n", i, j, idx,
                    lnSpotEQ, dom, lnRes);
#endif
    // previously lnRes was passed to SRMexp(), we implement same logic here.
    lnRes = Maths::isPositive(lnRes - SRMConstants::EXP_CUTOFF) ?
        SRMConstants::EXP_CUTOFF : lnRes;
    return lnRes;
}

    /** Getting the integer index corresponding to a given measurement
        date (md) or a future date (fd) for accessing ExpectedFX. */
// FwdIdx SRMEquityDiffuse::getExpectedPriceIndex(const DateTime& date)
// {
//     FwdIdx idx =  date.find(expSpotForwardDates);
//     assert(idx >= 0 && idx < expSpotForwardDates.size());
//     assert(expSpotForwardDates[idx] == date);
//     return idx;
// }

/** SRMEquityDiffuse::addAggregatedDates has to notify domIR about dates */
void   SRMEquityDiffuse::addAggregatedDates(
                        const DateTimeArray& spot,
                        const DateTimeArray& forward,
                        const DateTimeArray& forwardForward)
{
    QMCEquityDiffuse::addAggregatedDates(spot, forward, forwardForward);
    domIR->addAggregatedDates(spot, forward, forwardForward);
    updateMaxMaturity(spot, forward, forwardForward);
}

double SRMEquityDiffuse::getOriginalLnFwdSpot(SpotIdx dateIdx) 
{
    if (originalLnFwdSpot.empty())
        throw ModelException("SRMEquityDiffuse::getOriginalLnFwdSpot", "originalLnFwdSpot is not initialized for moment matching");
    return originalLnFwdSpot[dateIdx];
}

double SRMEquityDiffuse::getOriginalLnExpectedFwdSpot(FwdIdx dateIdx) 
{
    if (originalLnExpFwdSpot.empty())
        throw ModelException("SRMEquityDiffuse::getOriginalLnExpectedFwdSpot", "originalLnExpFwdSpot is not initialized for moment matching");
    return originalLnExpFwdSpot[dateIdx];
}


void SRMEquityDiffuse::setupBigStepIndicesAndYearFracs(const vector<int> &bigStepBools, vector<int> &bsIndices, vector<double> &bsYearFracs, vector<double> &bsYearFracsSqrt)
{
    //The bigStepYearFrac must be the time till the next bigStep.    
    double yearFracTotal = 0.0;
    double cfactorTotal = 0.0;
    size_t numDates = sigmaEQ.size();
    
    for (size_t simDateIdx = 0; simDateIdx < numDates+1; ++simDateIdx)
    {
        if (bigStepBools[simDateIdx])
        {
            bsIndices.push_back(simDateIdx);
            bsYearFracs.push_back(yearFracTotal);
            bsYearFracsSqrt.push_back(sqrt(yearFracTotal));
            yearFracTotal = 0.0;
        }

        if (simDateIdx < numDates)
        {
            yearFracTotal += yearFrac[simDateIdx];
        }        
    }

    size_t numBig = bsYearFracs.size();
    for (size_t i = 0;i<numBig - 1;i++)
    {
        bsYearFracs[i] = bsYearFracs[i+1];
        bsYearFracsSqrt[i] = bsYearFracsSqrt[i+1];
    }

    bsYearFracs[numBig-1] = yearFrac[numDates-1];
    bsYearFracsSqrt[numBig-1] = sqrtYearFrac[numDates-1];
}

size_t SRMEquityDiffuse::setupCriticalBigStepBools(vector<int> &bigStepBools)
{
    divPos = 0; // position in exDivIndexes
    spotEqPos = spotEqStart; // position in spotEQ/spotEQIndexes
    expSpotEqPos = 0;
    stopIdx = Maths::min(spotEqIndexes[spotEqPos],
        expSpotEqIndexes[expSpotEqPos]);

    //'critical' dates such as dates where the smile parameters change must be simulated.
    //After that we add more dates, up to 'timePointsPerYear'.
    double atmp = 0.0;
    double btmp = 0.0;
    double ctmp = 0.0;
    size_t numDates = sigmaEQ.size();
    size_t numBigSteps = 0;
    
    bigStepBools.resize(numDates+1);

    //Day after the last simDateIdx is a special 'big step': we diffuse up to it, but not past it.
    bigStepBools[numDates] = true;
    for (size_t simDateIdx = 0; simDateIdx < numDates; ++simDateIdx)
    {
        bool isBigStep = false;
        if (atmp != EqSmile_a1[simDateIdx] || btmp != EqSmile_a2[simDateIdx] || ctmp != EqSmile_a3[simDateIdx])
        {
            atmp = EqSmile_a1[simDateIdx];
            btmp = EqSmile_a2[simDateIdx];
            ctmp = EqSmile_a3[simDateIdx];
            isBigStep = true;
        }

        if (simDateIdx == 0)
            isBigStep = true;

        if (simDateIdx + todayIdx == stopIdx){
            isBigStep = true;

            if (stopIdx == spotEqIndexes[spotEqPos]){
                spotEqPos++;
            }
            if (stopIdx == expSpotEqIndexes[expSpotEqPos]){
                expSpotEqPos++;
            }
            stopIdx = Maths::min(spotEqIndexes[spotEqPos],
                expSpotEqIndexes[expSpotEqPos]);
        }

        if (simDateIdx + 1 == exDivIndexes[divPos]) {
            isBigStep = true;
            divPos++;
        }

       // if (simDateIdx % 4 == 0)
            isBigStep = true;

        if (isBigStep)
        {
            numBigSteps++;
        }

        bigStepBools[simDateIdx] = isBigStep;
    }

    return numBigSteps;
}

// symbol (referenced by MonteCarloLib.cpp) to ensure file gets linked in
bool SRMEQLinked = true;

DRLIB_END_NAMESPACE
