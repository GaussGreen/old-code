//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMFXDiffuse.cpp
//
//   Description : A generator of paths using stochastic rates
//                 for FX Assets
//
//   Date        : 13 Aug 2004
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/SRMFXUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/IQMCRNGManager.hpp"

DRLIB_BEGIN_NAMESPACE


const string SRMFXDiffuse::CONSTANT_SPOT_VOL("CONSTANT_SPOT_VOL");
const string SRMFXDiffuse::NO_FAILURE_ALLOWED("NO_FAILURE_ALLOWED");
const string SRMFXDiffuse::USE_LAST_LEVEL("USE_LAST_LEVEL");


/** Prepare to simulate a new path. Returns the sigmaFX for the first
    step */
double SRMFXDiffuse::begin(IQMCRNGManagerSP rngMgr){
    randoms = rngMgr->getCorrelatedRandoms(randomIndex); // save pointer to randoms
#if defined (DEBUGSRM)
    fprintf(stderr, "totalAssets.numcols=%d totalAssets.numrows=%d randomIdx=%d rnd[0]=%g rnd[1]=%g\n", randomNums.numCols(), randomNums.numRows(), randomIndex, randoms[0], randoms[1]);

#endif
    LnQ = origLnQ; // reset
    if (completePathRequested) {
        spotFXComplete[0] = LnQ;
    }
    spotFXPos = spotFXStart; // position in spotFX/spotFXIndexes
    expSpotFXPos = 0;
    stopIdx = Maths::min(spotFXIndexes[spotFXPos],
                         expSpotFXIndexes[expSpotFXPos]);
    // reset explosion booleans
    isIrExploded    = false;
    explosionHelper = 0.0;

    return sigmaFX.front(); // computed up front
}

/** Diffuses the fx to the next date. Returns the sigmaFX for the next
    date. idx will be 0 to move from the first date to the second. */
double SRMFXDiffuse::moveToNextDate(double forLnMONEY, // for new date
                             int    simDateIdx){
    // update LnQ
    double sigmaFX = this->sigmaFX[simDateIdx];
    LnQ += sigmaFX * (sqrtYearFrac[simDateIdx] * randoms[simDateIdx] -
                      0.5 * sigmaFX * yearFrac[simDateIdx]);
#if defined (DEBUGSRM)
        string s = diffusionDates[simDateIdx].toString();
        fprintf(stderr, "simDateIdx= %d [%s]\nLnQ= %g domLnMONEY= %g forLnMONEY= %g sigmaFX= %g sYF= %g yF= %g rnd= %g\n",
           simDateIdx, s.c_str(), LnQ, (*domLnMONEY)[simDateIdx], forLnMONEY, sigmaFX, sqrtYearFrac[simDateIdx],
           yearFrac[simDateIdx], randoms[simDateIdx]);
#endif

    // store LnQ if needed
    simDateIdx++; // makes indexing easier for dates
    if (completePathRequested) {
        double lnSpotFX = LnQ + (*domLnMONEY)[simDateIdx-1] - forLnMONEY;
        spotFXComplete[simDateIdx] = lnSpotFX;
    }



    if (simDateIdx + todayIdx == stopIdx){
        // must save the current spot fx (or actually log of it here)
        double lnSpotFX = LnQ + (*domLnMONEY)[simDateIdx-1] - forLnMONEY;

        // check whether an explosion took place

        if ( SRMUtil::SRMExplosion((*domLnMONEY)[simDateIdx-1]) && !isIrExploded ) {
            isIrExploded = true;
            explosionHelper = SRMConstants::EXP_CUTOFF;
        }
        if ( SRMUtil::SRMExplosion(forLnMONEY) && !isIrExploded ) {
            isIrExploded = true;
            explosionHelper = -SRMConstants::EXP_CUTOFF;
        }
        // spot fx indexing includes historic dates
        if (stopIdx == spotFXIndexes[spotFXPos]){
            spotFX[spotFXPos] = lnSpotFX;
            spotFXisIrExp[spotFXPos] = isIrExploded;
            spotFXCutoff[spotFXPos] = explosionHelper;
            spotFXPos++;
        }
        if (stopIdx == expSpotFXIndexes[expSpotFXPos]){
            expSpotFX[expSpotFXPos] = lnSpotFX;
            expSpotFXisIrExp[expSpotFXPos] = isIrExploded;
            expSpotFXCutoff[expSpotFXPos] = explosionHelper;
            expSpotFXPos++;
        }
        stopIdx = Maths::min(spotFXIndexes[spotFXPos],
                             expSpotFXIndexes[expSpotFXPos]);
    }
    if (simDateIdx == (int) SpotVol.size()){
        return 0.0; // we're on the last simulation date - value is irrelevant
    }
    // compute new value for sigmaFX
    double A_Smile = FXSmile_a1[simDateIdx]; /* skew  */
    double B_Smile = FXSmile_a2[simDateIdx]; /* smile */
    double spotvol = SpotVol[simDateIdx];
    if (A_Smile == 0 && B_Smile == 0){ /* no smile, can save time */
        sigmaFX = spotvol;
    } else {
        double x = LnQ + (*domLnMONEY)[simDateIdx-1] - forLnMONEY -
            LnFwdFx[simDateIdx-1];
        double C_Smile = FXSmile_a3[simDateIdx];
        double w = C_Smile * x;
        const double SMILE_CUTOFF = 100.0;
        if (w < SMILE_CUTOFF && w > - SMILE_CUTOFF) {
            double y = exp(w);
            double z = 1.0/y;
            /* loc vol = spotvol * (1 + A * Tanh(C * x) +
             *                   B * ( 1 - 1/Cosh(C * x)) ) */
            sigmaFX = spotvol * (1.0 + A_Smile * (y-z)/(y+z)
                                 + B_Smile * (1.0 - 2.0/(y+z)));
        } else if (w <= - SMILE_CUTOFF) {
            /* y = 0 */
            sigmaFX = spotvol * (1.0 - A_Smile + B_Smile);
        } else {
            /* z = 0 */
            sigmaFX = spotvol * (1.0 + A_Smile + B_Smile);
        }
    }
    this->sigmaFX[simDateIdx] = sigmaFX; // save it
    return sigmaFX;
}

/** Called at the end of the path */
void SRMFXDiffuse::end(){
    // turn log of spotFX into the real thing
    for (unsigned int i = spotFXStart; i < spotFX.size(); i++){
        if (spotFXisIrExp[i]) {
            spotFX[i] = exp(spotFXCutoff[i]);
        } else {
            spotFX[i] = exp(spotFX[i]);
        }
    }
}

/** Allows quanto equity/credit adjustment */
const vector<double>& SRMFXDiffuse::getSigmaFX() const {
    return sigmaFX;
}

SRMFXDiffuse::~SRMFXDiffuse()
{
// /*    delete timeLogic;
//     timeLogic = NULL;*/
}

/** New constructor of an empty shell */
SRMFXDiffuse::SRMFXDiffuse(
        IQMCDiffusibleInterestRateSP _domesticIRAsset,
        IQMCDiffusibleInterestRateSP _foreignIRAsset
        ) :
    //meaningful initializations
    //default ini with 'bad' values
    randomIndex(-999),
    randoms(NULL),
    LnQ(-999.),
    todayIdx(-999),
    stopIdx(-999),
    spotFXPos(-999),
    expSpotFXPos(-999),
    spotFXStart(-999),
    origLnQ(-999.),
    isIrExploded(false),
    explosionHelper(0.0),
    domDiscYCIdx(size_t(-1)),
    forDiscYCIdx(size_t(-1)),
    initialized(false),
    datesFixed(false),
//     timeLogic(new QMCHelperCachingTimeLogic),
    diffBound(_domesticIRAsset->getDiffusionBound(),
              _foreignIRAsset->getDiffusionBound()),
              completePathRequested(false)

{
    domesticIRAsset = QMCRatesDiffuseSP(dynamic_cast<QMCRatesDiffuse*>(_domesticIRAsset.get()));
    foreignIRAsset = QMCRatesDiffuseSP(dynamic_cast<QMCRatesDiffuse*>(_foreignIRAsset.get()));

    if (foreignIRAsset->getFXBase().get())
        throw ModelException("SRMFXDiffuse::SRMFXDiffuse()",
            "The currency ... attempted to attach more than one SRMFXDiffuse process.");


    domLnMONEY = &( _domesticIRAsset->getDomLnMONEY() );
    foreignIRAsset->setFXBase(QMCFXBaseDiffuseSP(this));
}


/** setting all the values */
void SRMFXDiffuse::setSRMFX(
             int                   _randomIndex,
             const DateTime&       _today,
             SRMFXUtilSP           _fxUtil,
             const vector<double>& _spotFX)       // historic spot FX
{
    randomIndex = _randomIndex;
    spotFX = _spotFX;
    today = _today;
    fxUtil = _fxUtil;


    domDiscYCIdx = domesticIRAsset->getDiscYCIdx(); // get Idx of disc YC
    forDiscYCIdx = foreignIRAsset->getDiscYCIdx();
}

/** finalize the timelines, allocate necessary memory */
void SRMFXDiffuse::finalize(DateTimeArrayConstSP allDatesSP)
{

    DateTimeArray spotRequestedDates = getSpotDates();
    DateTimeArray expSpotRequestedDates = getForwardDates();
    DateTimeArray expSpotForwardDates = getForwardForwardDates();

    assert(DateTime::isSubset(expSpotForwardDates, expSpotRequestedDates));

    size_t i;
    const size_t Nall = SRMUtil::getNumSimDates(today, *allDatesSP);
    const size_t Nefx = expSpotRequestedDates.size();
    const size_t Nefwd = expSpotForwardDates.size();
    const size_t Nfx  =    spotRequestedDates.size();

    DateTimeArray   diffusionDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);
    assert(Nall == diffusionDates.size());
    yearFrac.assign(Nall-1, 0.0);
    sqrtYearFrac.assign(Nall-1, 0.0);

    for (i = 0; i < yearFrac.size(); i++){
        yearFrac[i] = diffusionDates[i].yearFrac(diffusionDates[i+1]);
        sqrtYearFrac[i] = sqrt(yearFrac[i]);
    }



    // TODO getSubIndexes should throw an error if input date is not found in allDates
    assert(DateTime::isSubset((*allDatesSP), spotRequestedDates)); // throws if dfDates is not a subset
    spotFXIndexes = DateTime::getIndexes((*allDatesSP), spotRequestedDates); // indexes for when we save discount factors [Nfx]

    assert(DateTime::isSubset((*allDatesSP), expSpotRequestedDates));
    expSpotFXIndexes = DateTime::getIndexes((*allDatesSP), expSpotRequestedDates);// indexes for when we compute expected  [Nefx]

    todayIdx = today.find((*allDatesSP));

    // df contains history disc factors only. So must resize // i.e. increase the size

    spotFXStart = today.findUpper(spotRequestedDates);
    if (spotFXStart >= 0 && spotFXStart != spotRequestedDates.size() &&
        spotRequestedDates[spotFXStart] == today)
            ++spotFXStart; // skip today as diffusion results are always in futureDates
    //assert(spotFXStart == 0);

    spotFX.resize(Nfx);
    spotFXisIrExp.resize(Nfx);
    spotFXCutoff.resize (Nfx);

    expSpotFX.resize(Nefx);
    expSpotFXisIrExp.resize(Nefx);
    expSpotFXCutoff.resize(Nefx);


//    fwd2requested = DateTime::getProjection(expSpotForwardDates, expSpotRequestedDates);

    // Leave it here
    fwd2DomIRFwd.resize(Nefwd);
    fwd2ForIRFwd.resize(Nefwd);
    for(i=0; i < Nefwd; ++i)
    {
        fwd2DomIRFwd[i] = domesticIRAsset->getForwardForwardIndex(expSpotForwardDates[i]);
        fwd2ForIRFwd[i] = foreignIRAsset->getForwardForwardIndex(expSpotForwardDates[i]);
    }

    spotFXIndexes.push_back(todayIdx+Nall+1);//make life easier
    expSpotFXIndexes.push_back(todayIdx+Nall+1); // ditto

    assert(diffusionDates == fxUtil->getDomIR()->getSimDates());

    SpotVol = fxUtil->getSpotVol();
    sigmaFX.resize(Nall-1);
    spotFXComplete.resize(Nall);
    // calculate sigmaFX.front() - in fact it's just SpotVol[0]
    sigmaFX.front() = SpotVol.front();
    const DoubleArray& tmpLnFwdFx = fxUtil->getLnFwdFx();
    LnFwdFx.assign(tmpLnFwdFx.begin(), tmpLnFwdFx.end());
    FXSmile_a1 = fxUtil->getFXSmile_a1();
    FXSmile_a2 = fxUtil->getFXSmile_a2();
    FXSmile_a3 = fxUtil->getFXSmile_a3();
    origLnQ = fxUtil->getLogSpotFX();

    initialized = true;
    datesFixed = true;

    if (fxUtil->getMomentMatchingFlag())
    { // for moment matching: initialize originalLnFwdSpot and originalLnExpFwdSpot
        // I think I can use LnFwdEq
        // assume all dates > valuedate
        fxUtil->computeOriginalSpotPrices(spotRequestedDates, originalLnFwdSpot);
        fxUtil->computeOriginalSpotPrices(expSpotForwardDates, originalLnExpFwdSpot);
    }



    fxUtil = SRMFXUtilSP(0);
}



///////////////////////////////////////////// IQMC ///////////////////////

/** generate path across all dates. */
/** not done this way in SRMFXDiffuse model -- we do it one time-step at a time
being called from QMCRatesDiffuse */

void SRMFXDiffuse::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{}
    /** Accessing the natural log of the ExpectedFX(md, fd)
        where md is a simulated measurement date and fd is some future
        date after the measurement is made. */
double SRMFXDiffuse::getLnExpectedFX(
                                FwdIdx i,
                                FwdIdx j)
{
        SpotIdx    idx = getTimeLogic()->getReqEDFIdx(i); //fwd2requested[i];
        double dom =  domesticIRAsset->getLnExpectedDiscFactor(domDiscYCIdx, fwd2DomIRFwd[i], fwd2DomIRFwd[j]);
        double frn =  foreignIRAsset ->getLnExpectedDiscFactor(forDiscYCIdx, fwd2ForIRFwd[i], fwd2ForIRFwd[j]);
        double fx;

        if(expSpotFXisIrExp[idx]) {
            fx  = expSpotFXCutoff[idx];
        } else if (SRMUtil::SRMExplosion(frn)) {
            fx = -SRMConstants::EXP_CUTOFF;
        } else if (SRMUtil::SRMExplosion(dom)) {
            fx = SRMConstants::EXP_CUTOFF;
        } else {
            fx = expSpotFX[idx] + frn - dom;
        }
#if defined(DEBUGSRM)
        fprintf(stderr, "#SRMFXDiffuse::getLnExpectedFX i= %d j= %d fx= %g frn= %g dom= %g res= %g\n", int(i), int(j), expSpotFX[idx], frn, dom, fx);
#endif
        return fx;
}

    /** This is a declaration that ExpectedFX(md, fd[i]) should be produced
        by diffusion and so it might be requested later. */
void   SRMFXDiffuse::addAggregatedDates(
                        const DateTimeArray&    spot,
                        const DateTimeArray&    fwdDates,
                        const DateTimeArray&    fwdFwdDates)
{
    
    QMCFXBaseDiffuse::addAggregatedDates(spot, fwdDates, fwdFwdDates); // call base class
//    getTimeLogic()->addAggregatedDates(spot, fwdDates, fwdFwdDates);
//    getTimeLogic()->addAggregatedDates(DateTimeArray(), DateTimeArray(), fwdDates);
    // FIXME: do we need to add these dates to IR?
    domesticIRAsset->addAggregatedDates (spot, fwdDates, fwdFwdDates);
    foreignIRAsset->addAggregatedDates  (spot, fwdDates, fwdFwdDates);
        
    updateMaxMaturity(spot, fwdDates, fwdFwdDates);
    
}


double SRMFXDiffuse::getOriginalLnFwdSpot(SpotIdx dateIdx) 
{
    if (originalLnFwdSpot.empty())
        throw ModelException("SRMEquityDiffuse::getOriginalLnFwdSpot", "originalLnFwdSpot is not initialized");
    return originalLnFwdSpot[dateIdx];
}

double SRMFXDiffuse::getOriginalLnExpectedFwdSpot(FwdIdx dateIdx) 
{
    if (originalLnExpFwdSpot.empty())
        throw ModelException("SRMEquityDiffuse::getOriginalLnExpectedFwdSpot", "originalLnExpFwdSpot is not initialized");
    return originalLnExpFwdSpot[dateIdx];
}

const vector<double>& SRMFXDiffuse::requestFullSpotPath() {
    completePathRequested = true;
    return spotFXComplete;
}


// symbol (referenced by MonteCarloLib.cpp) to ensure file gets linked in
bool SRMFXLinked = true;

DRLIB_END_NAMESPACE
