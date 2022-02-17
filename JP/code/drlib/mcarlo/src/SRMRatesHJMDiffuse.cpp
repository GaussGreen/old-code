//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : SRMRatesHJMDiffuse.cpp (code previously in MCPathConfigSRM.cpp)
//
//   Description :
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesHJMDiffuse.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMFXDiffuse.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/QMCHelperCachingTimeLogic.hpp"
#include "edginc/IQMCRNGManager.hpp"
#include <cassert>
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** populates fields required for calculating sigmaR0 during simulation.
    A subset of irdiffuse::CalcSigmaR */
//// Depends on dates
void SRMRatesHJMDiffuse::calcSigmaRParams()
{
    static const string method("SRMRatesHJMDiffuse::calcSigmaRParams");
    double pivotRatio = 1.0 / (1.0 + ratesHJMUtilSP->getFwdShift());
    const DateTimeArray& dates = ratesHJMUtilSP->getSimDates();
    const DateTime& LastExpiry = ratesHJMUtilSP->getLastExpiry();

    for (int i = 0; i < dates.size()-1; i++)
    {
        const DateTime& CurrDate = dates[i];
        const DateTime& LastDate = dates.back(); // last sim date
        double svol = this->origSVol[i];
        double rbar = ratesHJMUtilSP->rBar(i);
        double QPivot = rbar * pivotRatio;

        // the case with normal IR is done separately for better efficiency
        if (zeroQ)
        {
            double MaxEffRate = MaxEffRateIR[i];
            double MinEffRate = MinEffRateIR[i];
            double effR = Maths::collar(QPivot, MaxEffRate, MinEffRate);
            svol /= pivotRatio;
            this->svol[i] = svol; // save
            this->effRSVol[i] = effR * svol;
        }
        else
        {
            // save QPivot
            this->QPivot[i] = QPivot;
            this->rbar[i] = rbar; // and save rbar
            DateTime T0;

            if (driver > 0)
            {
                T0 = CurrDate.rollDate(driver);
            }
            else
            if (!calibrateAgainstSwaptionVols)
            {
                /* Vol-driver rate's maturity */
                T0 = CurrDate.rollDate(LastDate.daysDiff(CurrDate)/2);
            }
            else
            { /*  BM calibration routine */
                if (ratesHJMUtilSP->getProcessedVol().get())
                {
                    if (CurrDate <= LastExpiry)
                    {
                        /* Non-parametric version of T0 */
                        T0 = ratesHJMUtilSP->sTau(CurrDate);
                    } else {
                        DateTime firstT0(ratesHJMUtilSP->sTau(CurrDate));
                        DateTime lastT0(firstT0>LastDate? firstT0: LastDate);
                        double slope = firstT0.yearFrac(lastT0)
                            / LastExpiry.yearFrac(LastDate);
                        int stepSize = (int) (slope *
                                            CurrDate.daysDiff(LastExpiry));
                        T0 = firstT0.rollDate(stepSize);
                        if (T0.getDate() == CurrDate.getDate())
                        {
                            T0 = T0.rollDate(1);
                        }
                    }
                }
            }

            if (T0 < CurrDate)
            {
                throw ModelException(method, "Internal error");
            }

            /* adjust for the effect of the forward shift. This is required
               because the value of SVol(.,.) assumes it is applied to an
               effR which Vladimirises to rbar
               This adj is not stored in svol because the unadjusted vol
               is used elsewhere (e.g. in FX spot vol bootstrapping)  */

            if (pivotRatio <= 1.0)
            {
                svol /= qRight + (1-qRight)*pivotRatio;
            } else {
                svol /= qLeft + (1-qLeft)*pivotRatio;
            }
            this->svol[i] = svol; // save

            double rbarT  = ratesHJMUtilSP->rBar(T0);

            if (Maths::isZero(rbarT))
            {
                throw ModelException(method, "Zero forward rate at "+
                                     T0.toString());
            }
            double rbarRatio = rbar/rbarT;
            /* calculate the vol factors and constants */
            vector<double> kT(ratesHJMUtilSP->numFactors()); // reserve some space
            vector<double> gT(ratesHJMUtilSP->numFactors()); // reserve some space

            ratesHJMUtilSP->kFactor(CurrDate, T0, kT);
            ratesHJMUtilSP->gFactor(CurrDate, T0, gT);

            computeXAndY(i, kT, gT, rbarRatio);
        }
    }
}

// make sure all elementary types are initialized to some unreal values
SRMRatesHJMDiffuse::SRMRatesHJMDiffuse() :
    QMCRatesDiffuse(),
    ratesHJMUtilSP(NULL),
    qLeft(-999),
    qRight(-999),
    zeroQ(true),
    pivotRatio(-999),
    discYCIdx(-1),
    diffYCIdx(-1)
{
}

SRMRatesHJMDiffuse::~SRMRatesHJMDiffuse()
{

}

void SRMRatesHJMDiffuse::setSRMRatesHJMDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMRatesHJMUtilSP      _srmRatesHJMUtil,
    bool                   _saveDomLnMoney,
    bool                   _saveSigmaR,
    double                 _NbSigmasMax,
    double                 _NbSigmasMin,
    const vector<double>&  _df,    // for historic dates - not used ?
    const double*          _irFxCorr,
    bool                   _calibrateAgainstSwaptionVols)
{
    static const string method("SRMRatesHJMDiffuse::setSRMRatesHJMDiffuse");

    randomIndex=_randomIndex;
    saveDomLnMoney=_saveDomLnMoney;
    saveSigmaR=_saveSigmaR;
    NbSigmasMax=_NbSigmasMax;
    NbSigmasMin=_NbSigmasMin;
    today=_today;
    calibrateAgainstSwaptionVols=_calibrateAgainstSwaptionVols;
    ratesHJMUtilSP = _srmRatesHJMUtil;
    qLeft=ratesHJMUtilSP->getQLeft();
    qRight=ratesHJMUtilSP->getQRight();
    zeroQ = Maths::isZero(qLeft) && Maths::isZero(qRight);
    origSVol=ratesHJMUtilSP->getSpotVols();
    // setting these named parameters in the derived classes
    setIrFxCorr(_irFxCorr);
    setAlphaAndRho();
}

// TODO: should allDates be the same as the futureDates ? (i.e. do we throw away past dates?)
// TODO2: prev implementation included past dates in allDates. For now we assume allDates do not have past dates.

void SRMRatesHJMDiffuse::finalize(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMRatesHJMDiffuse::finalize");

    DateTimeArray dfRequestedDates = getSpotDates();
    DateTimeArray edfRequestedDates = getForwardDates();
    DateTimeArray  edfForwardDates = getForwardForwardDates();

    SRM_VERIFY(DateTime::isSubset(edfForwardDates, edfRequestedDates),
        "Internal error in expected discount factor dates",
        method);

    SRM_VERIFY(ratesHJMUtilSP.get() != NULL,
        "Internal error : ratesHJMUtil void",
        method);

    DateTimeArray   diffusionDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);

    const size_t Nall =  SRMUtil::getNumSimDates(today, *allDatesSP); // rename later
    const size_t Nedf =  edfRequestedDates.size(); // rename later
    const size_t Ndf  =  dfRequestedDates.size();

    const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);

    SRM_VERIFY(numSimDates == diffusionDates.size(),
        "Internal error : problem in timeline",
        method);

    domLnMONEY = vector<double>(saveDomLnMoney? numSimDates-1: 0);
    sigmaR =     vector<double>(saveSigmaR? numSimDates-1: 0);

    effRSVol = vector<double>(numSimDates-1); // for zero q, effR * svol[Nall]
    QPivot =   vector<double>(numSimDates-1); // [Nall]
    svol =     vector<double>(numSimDates-1); /* contains SpotVol initially. If zeroQ is false
                             then processed in calcSigmaRParams */ // [Nall]
    rbar =     vector<double>(numSimDates-1); // [Nall]

    //// cached variables
    // TODO getSubIndexes should throw an error if input date is not found in allDates
    // throws if dfDates is not a subset
    SRM_VERIFY(DateTime::isSubset((*allDatesSP), dfRequestedDates),
        "Internal error : dfRequestDates not in AllDates",
        method);

    dfIndexes = DateTime::getIndexes((*allDatesSP),  dfRequestedDates); // indexes for when we save discount factors [Ndf]

    SRM_VERIFY(DateTime::isSubset((*allDatesSP), edfRequestedDates),
        "Internal error : edfRequestDates not in AllDates",
        method);

    expDFIndexes = DateTime::getIndexes((*allDatesSP), edfRequestedDates); // indexes for when we compute expected  [Nedf]

    if (!fx && expDFIndexes.empty() && dfIndexes.empty()){
        throw ModelException(method, "Internal error - no "
                                     "results required from diffused path");
    }

    todayIdx = today.find((*allDatesSP));
    // df contains history disc factors only. So must resize // i.e. increase the size

    //Initialise to zero: some products (e.g. TARN) use zero discount factor in order to ignore coupon payments for past dates.
    df.resize(Ndf, 0.0);

    getDiffYCIdx();
    getDiscYCIdx();

    todayIndex = today.findUpper(dfRequestedDates);
    if (todayIndex >= 0 && todayIndex != dfRequestedDates.size() && dfRequestedDates[todayIndex] == today)
        ++todayIndex; // skip today as diffusion results are always in futureDates

    SRM_VERIFY(todayIndex == dfRequestedDates.size() || dfRequestedDates[todayIndex] > today,
        "Internal error : problem in dfRequestedDates",
        method);

    // now make life easier - add request for index off the end
    size_t errIndex = (*allDatesSP).size()+numSimDates+1;
    dfIndexes.push_back(errIndex);
    expDFIndexes.push_back(errIndex);

    // assert that the diffusion dates have no past dates
    // TODO relax this requirement
    SRM_VERIFY(diffusionDates[0] >= today,
        "Internal error : diffusion can't have past dates",
        method);

    ratesHJMUtilSP->computeLogDiscFactor(logDiscFactor);

    ratesHJMUtilSP->calcEffRateLimit(NbSigmasMax, NbSigmasMin,
                                     MaxEffRateIR, MinEffRateIR);

    sqrtYearFrac = SRMUtil::computeSqrtYearFrac(diffusionDates);  // discount factors

    for(size_t i=0; i < ycForwardsDB.size(); ++i)
    {
        IYieldCurveConstSP  yc = ycForwardsDB[i].first;
        vector<double> &   fwd = ycForwardsDB[i].second;
        ratesHJMUtilSP->computeLogDiscFactor(edfForwardDates, fwd, yc); // initialize fwd between t0 and T
    }

    if (ratesHJMUtilSP->getMomentMatchingFlag())
    {
        ratesHJMUtilSP->computeLogDiscFactor(dfRequestedDates, originalDFs, ratesHJMUtilSP->getDiscYC());
    }

    assetDatesFixed = true; // no more requested dates can be added
    initialized = true; // fully initialized
}

/** called after origSVol externally recalibrated via ICE */
void SRMRatesHJMDiffuse::recalibrate(QMCRatesUtilSP thisQMCRatesUtilSP)
{
    SRMRatesHJMUtil* thisSRMRatesUtilSP = dynamic_cast<SRMRatesHJMUtil*>( 
        thisQMCRatesUtilSP.get() );
    QLIB_VERIFY( thisSRMRatesUtilSP != 0, string("Improper util class of type '") 
        + typeid(*thisQMCRatesUtilSP).name() + string("' when SRMRatesHJMUtil was expected") );

    double fwdShift = thisSRMRatesUtilSP->getFwdShift();
    if (!zeroQ && Maths::isZero(fwdShift)){
        return; // nothing to do
    }
    double pivotRatio = 1.0 / (1.0 + fwdShift);
    for (size_t i = 0; i < effRSVol.size(); i++){
        double svol = this->origSVol[i];
        if (zeroQ){
            double effR = effRSVol[i]/this->svol[i];
            svol /= pivotRatio;
            this->svol[i] = svol; // save
            effRSVol[i] = effR * svol; // and save this too
        } else {
            if (pivotRatio <= 1.0) {
                svol /= qRight + (1-qRight)*pivotRatio;
            } else {
                svol /= qLeft + (1-qLeft)*pivotRatio;
            }
            this->svol[i] = svol; // save
        }
    }
}

/////////////////////////////////////////// SRMRatesHJM1F ////////////////////
    /** Precompute X and Y variables at specified dateIdx using supplied
        parameters */
void SRMRatesHJM1F::computeXAndY(
                           int                   dateIdx,
                           const vector<double>& kT,
                           const vector<double>& gT,
                           double                rbarRatio)
{
    /* determine the Xij's and Yi's */
    X00[dateIdx] = gT[0] * kT[0] * rbarRatio;
    Y0[dateIdx]  = kT[0] * rbarRatio;
}

void SRMRatesHJM1F::setIrFxCorr(const double* irFxCorr)
{
    IrFxCorr0 = !fx? 0.0: irFxCorr[0];
}
void SRMRatesHJM1F::setAlphaAndRho()
{
    alpha=ratesHJMUtilSP->getAlpha(0);
}

void SRMRatesHJM1F::finalize(DateTimeArrayConstSP allDatesSP)
{
    try {
        const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);

        X00.resize(zeroQ? 0: numSimDates-1);
        Y0.resize(zeroQ? 0: numSimDates-1);

        SRMRatesHJMDiffuse::finalize(allDatesSP); // notify base class about allDates; may modify X00,Y0
        expDF.resize(expDFIndexes.size());

        calcSigmaRParams();
        ratesHJMUtilSP->computeKFactor(kFactor, 0);
        ratesHJMUtilSP->computeGFactor(gFactor, 0);
        ratesHJMUtilSP->populatePartialIntIR(today, getForwardForwardDates(), partialIntegral);
        ratesHJMUtilSP->populatePartialZeta(today, getForwardDates(), zeta);

        ratesHJMUtilSP = SRMRatesHJMUtilSP(0);
    } catch (exception& e){
        throw ModelException(e, "SRMRatesHJM1F::SRMRatesHJM1F");
    }
}
/** generate path across all dates. Essentially irdiffuse::DiffuseIR_1F plus
    a part of irdiffuse::CalcSigmaR */
void SRMRatesHJM1F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    // An example on how to get access to additional stream of RNGs:
    // IUniformRNGGenSP uni = rngMgr->getExtraGen(randomIndex);

    /* diffused quantities. Apologies for the upper case variables but
       names have been kept the same as in SRM3 */
    int numDates = logDiscFactor.size();
    int dfDatePos = todayIndex; // position in dfIndexes/df.
    /* set up position in expDFIndexes/expDF. If expDFIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of todayIndex and subtracting it in expDatePos formula
    int expDatePos = expDFIndexes.front() == 0? 1: 0;
    int dfDateIdx = dfIndexes[dfDatePos];
    int expDateIdx = expDFIndexes[expDatePos];
    int stopIdx = Maths::min(dfDateIdx, expDateIdx); // when to do something
    double   LnMONEY = 0.0;
    double   GAMMA0 = 0.0;
    double   PHI00 = 0.0;
    double   sigmaFX = 0; // FX data
    if (fx.get())
    {
        sigmaFX = fx->begin(rngMgr);
    }

    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){

        /* Step1: update sigmaR. START: From SRM3::CalcSigmaR */
        double sigmaR;
        if (zeroQ){
            sigmaR = effRSVol[i];

        } else {
            double aRate  = rbar[i] + X00[i] * PHI00 + Y0[i] * GAMMA0;
            double QPivot = this->QPivot[i];
            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);

            double effR = (aRate >= QPivot) ?
                (QPivot_x_OneMinusQR + (aRate * qRight)):
                (QPivot_x_OneMinusQL + (aRate * qLeft));

                effR = Maths::collar(effR, MaxEffRateIR[i], MinEffRateIR[i]);
                sigmaR = effR * svol[i];
        }

        // END: from SRM3::CalcSigmaR
        double sigmaR0 = sigmaR * alpha;
        double RootDelT = sqrtYearFrac[i];
        double g0 = gFactor[i];
        double V0 = sigmaR0 * g0 * RootDelT;   /* zero bond  spot vol    */
        double k0 = kFactor[i];
        double S0 = sigmaR0 * k0 * RootDelT;  /* inst fwd rate spot vol */

        //
        /* Step2: update LnMONEY */
        /*-----------------------*/
        double lnZeroTerm = -0.5*PHI00*g0*g0 - g0*GAMMA0;
        double lnDiscFwdZero = logDiscFactor[i];
        double W0 = randoms[i]; // random number

        LnMONEY += 0.5*V0*V0 + V0*W0 - (lnDiscFwdZero + lnZeroTerm);
        /* Step3: update GAMMA */
        /*---------------------*/
        GAMMA0 = k0*GAMMA0 + k0*g0*PHI00 + S0*V0 + S0*W0;
        /* Step4: update PHI */
        PHI00 = k0 * k0 * PHI00 + S0 * S0;

        /* Step5: add the CUPS component if necessary */
        if (fx.get()){
            double  F0 = IrFxCorr0 * sigmaFX * RootDelT;
            GAMMA0  -= (F0 * S0);
            LnMONEY -= (F0 * V0);
            sigmaFX = fx->moveToNextDate(LnMONEY, i); // get new sigmaFX
        }
        // We need this for each equity ccy
        if (!domLnMONEY.empty()) {
            domLnMONEY[i] = LnMONEY; // needed by fx's linked to dom
        }
        if(!this->sigmaR.empty()) {
            this->sigmaR[i] = sigmaR;
        }
        i++; /* up the loop counter. This comes back to the fact with eg
                5 dates there are 4 diffusion steps. It's easier to
                identify the index of the date you want rather than the
                index of the following one */
        if (i + todayIdx == stopIdx){ // TODO double check later with past dates
            // hit an 'event' - record something
            if (stopIdx == dfDateIdx){
                // save the value
                df[dfDatePos] = LnMONEY;
                dfDatePos++;
                dfDateIdx = dfIndexes[dfDatePos];
            }
            if (stopIdx == expDateIdx){
                // save gamma, and phi
                expDF[expDatePos].GAMMA0 = GAMMA0;
                expDF[expDatePos].PHI00 = PHI00;
                expDatePos++;
                expDateIdx = expDFIndexes[expDatePos]; // that's why we pushed an extra element to expDFIndexes;

            }
            stopIdx = Maths::min(dfDateIdx, expDateIdx); // refresh
        }
    }
    // then do exp on df vector
    for (size_t j = todayIndex; j < df.size(); j++){
        df[j] = exp(-df[j]);
    }
    if (fx.get()){
        fx->end();
    }
}

double SRMRatesHJM1F::getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx j) {
        const int iEDF = getFwdIdx2EdfIdx(i);
        double g = getGFactor(i, j);
        double result = getYCForward(idx, i, j) - expDF[iEDF].GAMMA0 * g - 0.5 * g * g * expDF[iEDF].PHI00;
        return result;
}


//// two factor IR
/** Precompute X and Y variables at specified dateIdx using supplied
    parameters */
void SRMRatesHJM2F::computeXAndY(
    int                   dateIdx,
    const vector<double>& kT,
    const vector<double>& gT,
    double                rbarRatio)
{
    /* determine the Xij's and Yi's */
    Vars& vars = varsPerTimePoint[dateIdx]; // for ease
    vars.X00 = gT[0] * kT[0]* rbarRatio;
    vars.X11 = gT[1] * kT[1]* rbarRatio;
    vars.X01 = (gT[0] * kT[1] + gT[1] * kT[0])* rbarRatio;
    vars.Y0  = kT[0]* rbarRatio;
    vars.Y1  = kT[1]* rbarRatio;
}

void SRMRatesHJM2F::setIrFxCorr(const double* irFxCorr)
{
    IrFxCorr0 = !fx? 0.0: irFxCorr[0];
    IrFxCorr1 = !fx? 0.0: irFxCorr[1];
}
void SRMRatesHJM2F::setAlphaAndRho()
{
    alp0=ratesHJMUtilSP->getAlpha(0);
    alp1=ratesHJMUtilSP->getAlpha(1);
    rho01=ratesHJMUtilSP->getRho(0);
}

//// constructor is trivial and in the header
//SRMRatesHJM2F::SRMRatesHJM2F(){}

void SRMRatesHJM2F::finalize(DateTimeArrayConstSP allDatesSP)
{
    const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);
    varsPerTimePoint.resize(numSimDates-1);

    SRMRatesHJMDiffuse::finalize(allDatesSP); // notify base class about allDates : may modify varsPerTimePoint!
    expDF.resize(expDFIndexes.size());

    try {
        calcSigmaRParams();
        vector<vector<double> > tmpSpace(4);
        ratesHJMUtilSP->computeKFactor(tmpSpace[0], 0);
        ratesHJMUtilSP->computeKFactor(tmpSpace[1], 1);
        ratesHJMUtilSP->computeGFactor(tmpSpace[2], 0);
        ratesHJMUtilSP->computeGFactor(tmpSpace[3], 1);

        for (size_t i = 0; i < varsPerTimePoint.size(); i++){
            varsPerTimePoint[i].k0 = tmpSpace[0][i];
            varsPerTimePoint[i].k1 = tmpSpace[1][i];
            varsPerTimePoint[i].g0 = tmpSpace[2][i];
            varsPerTimePoint[i].g1 = tmpSpace[3][i];
        }
        ratesHJMUtilSP->populatePartialIntIR(today, getForwardForwardDates(), partialIntegral);
        ratesHJMUtilSP->populatePartialZeta(today, getForwardDates(), zeta);

        ratesHJMUtilSP = SRMRatesHJMUtilSP(0);

    } catch (exception& e){
            throw ModelException(e, "SRMRatesHJM2F::finalize");
    }
}

/** generate path across all dates. Essentially irdiffuse::DiffuseIR_1F plus
    a part of irdiffuse::CalcSigmaR */
void SRMRatesHJM2F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    int numDates = logDiscFactor.size();
    int dfDatePos = todayIndex; // position in dfIndexes/df.
    /* set up position in expDFIndexes/expDF. If expDFIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of todayIndex and subtracting it in expDatePos formula
    int expDatePos = expDFIndexes.front() == 0? 1: 0;
    int dfDateIdx = dfIndexes[dfDatePos];
    int expDateIdx = expDFIndexes[expDatePos];
    int stopIdx = Maths::min(dfDateIdx, expDateIdx); // when to do something
    /* diffused quantities. Apologies for the upper case variables but
       names have been kept the same as in SRM3 */
    double   LnMONEY = 0.0;
    double   GAMMA0 = 0.0;
    double   GAMMA1 = 0.0;
    double   PHI00 = 0.0;
    double   PHI01 = 0.0;
    double   PHI11 = 0.0;
    double   sigmaFX = 0; // FX data
    if (fx.get()){
        sigmaFX = fx->begin(rngMgr);
    }
    // random numbers
    const double* Z0 = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    const double* Z1 = rngMgr->getCorrelatedRandoms(randomIndex, 1); // for ease
    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){
        const Vars& vars = varsPerTimePoint[i]; // for ease
        /* Step1: update sigmaR. START: From SRM3::CalcSigmaR */
        double sigmaR;
        if (zeroQ){
            sigmaR = effRSVol[i];
        } else {
            double aRate = rbar[i] +
                (vars.X00*PHI00 + vars.X01*PHI01 + vars.X11*PHI11)
                + (vars.Y0*GAMMA0 + vars.Y1*GAMMA1);
            double QPivot = this->QPivot[i];
            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);

            double effR = (aRate >= QPivot) ?
                (QPivot_x_OneMinusQR + (aRate * qRight)):
                (QPivot_x_OneMinusQL + (aRate * qLeft));

                effR = Maths::collar(effR, MaxEffRateIR[i], MinEffRateIR[i]);
            sigmaR = effR * svol[i];
        }
        // END: from SRM3::CalcSigmaR
        double sigmaR0 = sigmaR * alp0;
        double sigmaR1 = sigmaR * alp1;
        double RootDelT = sqrtYearFrac[i];
        double g0 = vars.g0;
        double V0 = sigmaR0 * g0 * RootDelT;   /* zero bond  spot vol    */
        double g1 = vars.g1;
        double V1 = sigmaR1 * g1 * RootDelT;   /* zero bond  spot vol    */
        double k0 = vars.k0;
        double S0 = sigmaR0 * k0 * RootDelT;  /* inst fwd rate spot vol */
        double k1 = vars.k1;
        double S1 = sigmaR1 * k1 * RootDelT;  /* inst fwd rate spot vol */

        /* Step2: update LnMONEY */
        /*-----------------------*/
        double lnZeroTerm  = -(0.5*g0*g0*PHI00 + g0*g1*PHI01 +
                               0.5*g1*g1*PHI11) -
            (g0*GAMMA0 + g1*GAMMA1);
        double lnDiscFwdZero = logDiscFactor[i];
        double W0 = Z0[i]; // random number
        double W1 = Z1[i]; // random number
        LnMONEY += (0.5*V0*V0 + V0*V1*rho01 + 0.5*V1*V1) + (V0*W0 + V1*W1)
            - (lnDiscFwdZero + lnZeroTerm);
        /* Step3: update GAMMA */
        /*---------------------*/
        GAMMA0 = k0*GAMMA0 + k0*(g0*PHI00 + g1*PHI01) +
            S0*(V0 + V1*rho01) + S0*W0;
        GAMMA1 = k1*GAMMA1 +
            k1*(g0*PHI01 + g1*PHI11) + S1*(V0*rho01 + V1) + S1*W1;
        /* Step4: update PHI */
        PHI00 =  k0*k0*PHI00 + S0*S0;
        PHI01 =  k0*k1*PHI01 + S0*S1 * rho01;
        PHI11 =  k1*k1*PHI11 + S1*S1;
        /* Step5: add the CUPS component if necessary */

        if (fx.get()){
            double  F0 = IrFxCorr0 * sigmaFX * RootDelT;
            double  F1 = IrFxCorr1 * sigmaFX * RootDelT;
            GAMMA0  -= (F0 * S0);
            GAMMA1  -= (F1 * S1);
            LnMONEY -= (F0 * V0 + F1 * V1);
            sigmaFX = fx->moveToNextDate(LnMONEY, i); // get new sigmaFX
        }
        // We need this for each equity ccy
        if (!domLnMONEY.empty()) {
            domLnMONEY[i] = LnMONEY; // needed by fx's linked to dom
        }
        if (!this->sigmaR.empty()) {
            this->sigmaR[i] = sigmaR;
        }
        i++; /* up the loop counter. This comes back to the fact with eg
                5 dates there are 4 diffusion steps. It's easier to
                identify the index of the date you want rather than the
                index of the following one */
        if (i + todayIdx == stopIdx){
            // hit an 'event' - record something
            if (stopIdx == dfDateIdx){
                // save the value
                df[dfDatePos] = LnMONEY;
                dfDatePos++;
                dfDateIdx = dfIndexes[dfDatePos];
            }
            if (stopIdx == expDateIdx){
                // save gamma, and phi
                expDF[expDatePos].GAMMA0 = GAMMA0;
                expDF[expDatePos].GAMMA1 = GAMMA1;
                expDF[expDatePos].PHI00 = PHI00;
                expDF[expDatePos].PHI01 = PHI01;
                expDF[expDatePos].PHI11 = PHI11;
                expDatePos++;
                expDateIdx = expDFIndexes[expDatePos];
            }
            stopIdx = Maths::min(dfDateIdx, expDateIdx); // refresh
        }
    }
    // then do exp on df vector
    for (size_t j = todayIndex; j < df.size(); j++){
        df[j] = exp(-df[j]);
    }
    if (fx.get()){
        fx->end();
    }
}


double SRMRatesHJM2F::getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx j) {
        const int iEDF = getFwdIdx2EdfIdx(i);
        const DiffusedState & state = expDF[iEDF];
        const double g0 = getGFactor(i, j, 0);
        const double g1 = getGFactor(i, j, 1);

        double result = getYCForward(idx, i, j) - (0.5*g0*g0*state.PHI00+ g0*g1*state.PHI01+
                        0.5*g1*g1*state.PHI11) -
                        (g0*state.GAMMA0 + g1*state.GAMMA1);
        return result;
}


//// three factor IR

/** Precompute X and Y variables at specified dateIdx using supplied
    parameters */
void SRMRatesHJM3F::computeXAndY(int                   dateIdx,
                                 const vector<double>& kT,
                                 const vector<double>& gT,
                                 double                rbarRatio)
{
    /* determine the Xij's and Yi's */
    Vars& vars = varsPerTimePoint[dateIdx]; // for ease
    vars.X00 = gT[0] * kT[0]* rbarRatio;
    vars.X11 = gT[1] * kT[1]* rbarRatio;
    vars.X22 = gT[2] * kT[2]* rbarRatio;
    vars.X01 = (gT[0] * kT[1] + gT[1] * kT[0])* rbarRatio;
    vars.X02 = (gT[0] * kT[2] + gT[2] * kT[0])* rbarRatio;
    vars.X12 = (gT[1] * kT[2] + gT[2] * kT[1])* rbarRatio;
    vars.Y0  = kT[0]* rbarRatio;
    vars.Y1  = kT[1]* rbarRatio;
    vars.Y2  = kT[2]* rbarRatio;
}

void SRMRatesHJM3F::setIrFxCorr(const double* irFxCorr)
{
    IrFxCorr0 = !fx? 0.0: irFxCorr[0];
    IrFxCorr1 = !fx? 0.0: irFxCorr[1];
    IrFxCorr2 = !fx? 0.0: irFxCorr[2];
}
void SRMRatesHJM3F::setAlphaAndRho()
{
    alp0=ratesHJMUtilSP->getAlpha(0);
    alp1=ratesHJMUtilSP->getAlpha(1);
    alp2=ratesHJMUtilSP->getAlpha(2);
    rho01=ratesHJMUtilSP->getRho(0);
    rho02=ratesHJMUtilSP->getRho(1);
    rho12=ratesHJMUtilSP->getRho(2);
}

//// constructor is trivial now
//SRMRatesHJM3F::SRMRatesHJM3F() {}

void SRMRatesHJM3F::finalize(DateTimeArrayConstSP allDatesSP)
{
    const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);
    varsPerTimePoint.resize(numSimDates-1);

    SRMRatesHJMDiffuse::finalize(allDatesSP); // notify base class about allDates; may modify varsPerTimePoint[]
    expDF.resize(expDFIndexes.size());

    try {
        calcSigmaRParams();
        vector<vector<double> > tmpSpace(6);
        ratesHJMUtilSP->computeKFactor(tmpSpace[0], 0);
        ratesHJMUtilSP->computeKFactor(tmpSpace[1], 1);
        ratesHJMUtilSP->computeKFactor(tmpSpace[2], 2);
        ratesHJMUtilSP->computeGFactor(tmpSpace[3], 0);
        ratesHJMUtilSP->computeGFactor(tmpSpace[4], 1);
        ratesHJMUtilSP->computeGFactor(tmpSpace[5], 2);

        for (size_t i = 0; i < varsPerTimePoint.size(); i++){
            varsPerTimePoint[i].k0 = tmpSpace[0][i];
            varsPerTimePoint[i].k1 = tmpSpace[1][i];
            varsPerTimePoint[i].k2 = tmpSpace[2][i];
            varsPerTimePoint[i].g0 = tmpSpace[3][i];
            varsPerTimePoint[i].g1 = tmpSpace[4][i];
            varsPerTimePoint[i].g2 = tmpSpace[5][i];
        }

        ratesHJMUtilSP->populatePartialIntIR(today, getForwardForwardDates(), partialIntegral);
        ratesHJMUtilSP->populatePartialZeta(today, getForwardDates(), zeta);

        ratesHJMUtilSP = SRMRatesHJMUtilSP(0);

    } catch (exception& e){
        throw ModelException(e, "SRMRatesHJM3F::finalize");
    }
}
/** generate path across all dates. Essentially irdiffuse::DiffuseIR_1F plus
    a part of irdiffuse::CalcSigmaR */
void SRMRatesHJM3F::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    int numDates = logDiscFactor.size();
    int dfDatePos = todayIndex; // position in dfIndexes/df.
    /* set up position in expDFIndexes/expDF. If expDFIndexes[0] is 0 then
       we just store 0 for the GAMMA/PHI and start at expDFIndexes[1] */
    // to do: review use of todayIndex and subtracting it in expDatePos formula
    int expDatePos = expDFIndexes.front() == 0? 1: 0;
    //int todayIdx = today.find(diffusionDates);
    int dfDateIdx = dfIndexes[dfDatePos];
    int expDateIdx = expDFIndexes[expDatePos];
    int stopIdx = Maths::min(dfDateIdx, expDateIdx); // when to do something
    /* diffused quantities. Apologies for the upper case variables but
       names have been kept the same as in SRM3 */
    double   LnMONEY = 0.0;
    double   GAMMA0 = 0.0;
    double   GAMMA1 = 0.0;
    double   GAMMA2 = 0.0;
    double   PHI00 = 0.0;
    double   PHI01 = 0.0;
    double   PHI02 = 0.0;
    double   PHI11 = 0.0;
    double   PHI12 = 0.0;
    double   PHI22 = 0.0;
    double   sigmaFX = 0; // FX data
    if (fx.get()){
        sigmaFX = fx->begin(rngMgr);
    }
    // random numbers
    const double* Z0 = rngMgr->getCorrelatedRandoms(randomIndex); // for ease
    const double* Z1 = rngMgr->getCorrelatedRandoms(randomIndex, 1); // for ease
    const double* Z2 = rngMgr->getCorrelatedRandoms(randomIndex, 2); // for ease
    // loop over all dates
    for (int i = 0; i < numDates; /* increment in body */){
        const Vars& vars = varsPerTimePoint[i]; // for ease
        /* Step1: update sigmaR. START: From SRM3::CalcSigmaR */
        double sigmaR;
        if (zeroQ){
            sigmaR = effRSVol[i];
        } else {
            double aRate = rbar[i] +
                (vars.X00*PHI00 + vars.X01*PHI01 + vars.X02*PHI02 +
                 vars.X11*PHI11 + vars.X12*PHI12 + vars.X22*PHI22)
                +(vars.Y0*GAMMA0 + vars.Y1*GAMMA1 + vars.Y2*GAMMA2);
            double QPivot = this->QPivot[i];
            double QPivot_x_OneMinusQL = QPivot * (1.0 - qLeft);
            double QPivot_x_OneMinusQR = QPivot * (1.0 - qRight);

            double effR = (aRate >= QPivot) ?
                (QPivot_x_OneMinusQR + (aRate * qRight)):
                (QPivot_x_OneMinusQL + (aRate * qLeft));

            effR = Maths::collar(effR, MaxEffRateIR[i], MinEffRateIR[i]);
            sigmaR = effR * svol[i];
        }
        // END: from SRM3::CalcSigmaR
        double sigmaR0 = sigmaR * alp0;
        double sigmaR1 = sigmaR * alp1;
        double sigmaR2 = sigmaR * alp2;
        double RootDelT = sqrtYearFrac[i];
        double g0 = vars.g0;
        double V0 = sigmaR0 * g0 * RootDelT;   /* zero bond  spot vol    */
        double g1 = vars.g1;
        double V1 = sigmaR1 * g1 * RootDelT;   /* zero bond  spot vol    */
        double g2 = vars.g2;
        double V2 = sigmaR2 * g2 * RootDelT;   /* zero bond  spot vol    */
        double k0 = vars.k0;
        double S0 = sigmaR0 * k0 * RootDelT;  /* inst fwd rate spot vol */
        double k1 = vars.k1;
        double S1 = sigmaR1 * k1 * RootDelT;  /* inst fwd rate spot vol */
        double k2 = vars.k2;
        double S2 = sigmaR2 * k2 * RootDelT;  /* inst fwd rate spot vol */

        /* Step2: update LnMONEY */
        /*-----------------------*/
        double lnZeroTerm = -(0.5*g0*g0*PHI00 + g0*g1*PHI01 +
                              g0*g2*PHI02 + 0.5*g1*g1*PHI11 +
                              g1*g2*PHI12 + 0.5*g2*g2*PHI22) -
            (g0*GAMMA0 + g1*GAMMA1 + g2*GAMMA2);
        double lnDiscFwdZero = logDiscFactor[i];
        double W0 = Z0[i]; // random number
        double W1 = Z1[i]; // random number
        double W2 = Z2[i]; // random number
        LnMONEY += (0.5*V0*V0 + V0*V1*rho01 + V0*V2*rho02 +
                    0.5*V1*V1 + V1*V2*rho12 + 0.5*V2*V2) +
            (V0*W0 + V1*W1 + V2*W2) - (lnDiscFwdZero + lnZeroTerm);
        /* Step3: update GAMMA */
        /*---------------------*/
        GAMMA0 = k0*GAMMA0 +
            k0*(g0*PHI00 + g1*PHI01 + g2*PHI02) +
            S0*(V0 + V1*rho01 + V2*rho02) + S0*W0;

        GAMMA1 = k1*GAMMA1 +
            k1*(g0*PHI01 + g1*PHI11 + g2*PHI12) +
            S1*(V0*rho01 + V1 + V2*rho12) + S1*W1;

        GAMMA2 = k2*GAMMA2 +
            k2*(g0*PHI02 + g1*PHI12 + g2*PHI22) +
            S2*(V0*rho02 + V1*rho12 + V2) + S2*W2;

        /* Step4: update PHI */
        PHI00 =  k0*k0*PHI00 + S0*S0;
        PHI01 =  k0*k1*PHI01 + S0*S1 * rho01;
        PHI02 =  k0*k2*PHI02 + S0*S2 * rho02;
        PHI11 =  k1*k1*PHI11 + S1*S1;
        PHI12 =  k1*k2*PHI12 + S1*S2 * rho12;
        PHI22 =  k2*k2*PHI22 + S2*S2;
        /* Step5: add the CUPS component if necessary */
        if (fx.get()){
            double  F0 = IrFxCorr0 * sigmaFX * RootDelT;
            double  F1 = IrFxCorr1 * sigmaFX * RootDelT;
            double  F2 = IrFxCorr2 * sigmaFX * RootDelT;
            GAMMA0  -= (F0 * S0);
            GAMMA1  -= (F1 * S1);
            GAMMA2  -= (F2 * S2);
            LnMONEY -= (F0 * V0 + F1 * V1 + F2 * V2);
            sigmaFX = fx->moveToNextDate(LnMONEY, i); // get new sigmaFX
        }
        // We need this for each equity ccy
        if (!domLnMONEY.empty()) {
            domLnMONEY[i] = LnMONEY; // needed by fx's linked to dom
        }
        if (!this->sigmaR.empty()) {
            this->sigmaR[i] = sigmaR;
        }
        i++; /* up the loop counter. This comes back to the fact with eg
                5 dates there are 4 diffusion steps. It's easier to
                identify the index of the date you want rather than the
                index of the following one */
        if (i + todayIdx == stopIdx){
            // hit an 'event' - record something
            if (stopIdx == dfDateIdx){
                // save the value
                df[dfDatePos] = LnMONEY;
                dfDatePos++;
                dfDateIdx = dfIndexes[dfDatePos];
            }
            if (stopIdx == expDateIdx){
                // save gamma, and phi
                expDF[expDatePos].GAMMA0 = GAMMA0;
                expDF[expDatePos].GAMMA1 = GAMMA1;
                expDF[expDatePos].GAMMA2 = GAMMA2;
                expDF[expDatePos].PHI00 = PHI00;
                expDF[expDatePos].PHI01 = PHI01;
                expDF[expDatePos].PHI02 = PHI02;
                expDF[expDatePos].PHI11 = PHI11;
                expDF[expDatePos].PHI12 = PHI12;
                expDF[expDatePos].PHI22 = PHI22;
                expDatePos++;
                expDateIdx = expDFIndexes[expDatePos];
            }
            stopIdx = Maths::min(dfDateIdx, expDateIdx); // refresh
        }
    }
    // then do exp on df vector
    for (size_t j = todayIndex; j < df.size(); j++){
        df[j] = exp(-df[j]);
    }
    if (fx.get()){
        fx->end();
    }
}

double SRMRatesHJM3F::getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx j) {
        const int iEDF = getFwdIdx2EdfIdx(i);
        const DiffusedState & state = expDF[iEDF];
        const double g0 = getGFactor(i, j, 0);
        const double g1 = getGFactor(i, j, 1);
        const double g2 = getGFactor(i, j, 2);

        double result = getYCForward(idx, i, j)
            -(0.5*g0*g0*state.PHI00 + g0*g1*state.PHI01 +
            g0*g2*state.PHI02 +
            0.5*g1*g1*state.PHI11 + g1*g2*state.PHI12 +
            0.5*g2*g2*state.PHI22) -
            (g0*state.GAMMA0 + g1*state.GAMMA1 + g2*state.GAMMA2);
        return result;
}



DRLIB_END_NAMESPACE



