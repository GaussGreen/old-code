//----------------------------------------------------------------------------
//
//   Group       : QR Cross Asset
//
//   Filename    : SRMCreditLiborDiffuse.cpp
//
//   Description : CreditLibor path generation
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SRMCreditLiborDiffuse.hpp"
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

void SRMCreditLiborDiffuse::setSRMCreditLiborDiffuse(
    int                    _randomIndex,
    const DateTime&        _today,
    SRMCreditLiborUtilSP   _srmCreditLiborUtil,
    const double           _crFxCorr,
    const vector<double>&  _prob)    // for historic dates
{
    static const string method("SRMCreditLiborDiffuse::setSRMCreditLiborDiffuse");

    srmCreditLiborUtil = _srmCreditLiborUtil;
    //qLeft		 = srmCreditLiborUtil->getQLeft();
    //qRight       = srmCreditLiborUtil->getQRight();
    randomIndex  = _randomIndex;
    crFxCorr     = _crFxCorr;
    today        = _today;
    probStart    = -1;              // FIXME probStart : see SRMIRModel how to calc it
    prob         = _prob;
    //zeroQ        = Maths::isZero(qLeft) && Maths::isZero(qRight);
    setRecoveryRate(srmCreditLiborUtil->getCdsCurve()->getRecovery());

}

/*
void SRMCreditLiborDiffuse::setCrIrCorr(const vector<double>& corrCRIR)
{
    static const string method("SRMCreditLiborDiffuse::setCrIrCorr");
    rho03 = corrCRIR[0];
}
*/

void SRMCreditLiborDiffuse::storeT(const DateTime& today, const DateTimeArray& dates) {
    static const string method("SRMCreditLiborDiffuse::storeT");
    const int numDates = dates.size();
    esdfForwardYearFracs.resize(numDates);
    for (int i = 0; i < numDates; i++) {
        double yearFrac = SRMYearFrac(today, dates[i]);
        esdfForwardYearFracs[i] = yearFrac;
    }
}


     // all the parameter initialization is done inside finalize
void SRMCreditLiborDiffuse::finalizePathGenerator(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMCreditLiborDiffuse::finalize");
    try {

    const DateTimeArray& simDates = srmCreditLiborUtil->getSimDates();
    const int numSimDates         = simDates.size();
    assert(numSimDates == SRMUtil::getNumSimDates(today, *allDatesSP));

    calcFirstAndLastDiffusionIdx(simDates, *allDatesSP);
    processAllDates(allDatesSP);
    calcRemapToIRAssetIdx(getForwardForwardDates());

    if (isWholeTimelineSurvivalRequested())
        wholeTimelineLogSurvProb.resize(lastDiffusionIdx-todayIdx+1,0.0);
    // initialize model-specific dates and parameters
    mResetDates           = srmCreditLiborUtil->getResetDates();
    mAccruedFrac          = srmCreditLiborUtil->getAccruedFrac();
    mIntensityVols        = srmCreditLiborUtil->getSpotVols();
    mInitialIntensities   = srmCreditLiborUtil->getInitialIntensities();
    mSpotSurvProb         = srmCreditLiborUtil->getSpotSurvProb();
    mYearFracToFirstReset = srmCreditLiborUtil->getYearFracToFirstReset();
    mNumInt               = srmCreditLiborUtil->getNumInt();

    //srmCreditLiborUtil->computeLogFwdProbSimple(logFwdProbSimple);


    if (srmCreditLiborUtil->getMomentMatchingFlag())
    {
        srmCreditLiborUtil->computeLogFwdProbSimple(getSpotDates() /*sdfRequestedDates*/, originalProbs);
        srmCreditLiborUtil->computeLogFwdProbSimple(getForwardForwardDates() /*esdfForwardDates*/, originalExpProbs);
    }


    sqrtYearFrac = SRMUtil::computeSqrtYearFrac(simDates);

//    srmCreditLiborUtil->computeLogFwdProbSimple(today, esdfForwardDates, fwdLnProbRatios);

//    trimToDiffusion();

    expProb.resize(getNumESDFDates());
    storeT(today, getForwardForwardDates() /*esdfForwardDates*/);

    } catch (exception& e){
        throw ModelException(e, method);
    }

}


void SRMCreditLiborDiffuse::generatePathAndSurvivalRates(IQMCRNGManagerSP rngMgr) 
{
    const double* randoms = rngMgr->getCorrelatedRandoms(randomIndex);                // for ease

	int i,j,
		firstIntAliveIdx,										    // first discrete intensity still alive
		numDates = sqrtYearFrac.size(),//logFwdProbSimple.size(),
		numInt   = mNumInt,
		frozenIntensityIdx = -1.0;											//number of state variables, i.e. discrete intensities (hazard rates)

    // set up position in expProbIndexes/expSvob. If expProbIndexes[0] is 0 then
    // we just store 0 for the state variables and start at expDFIndexes[1]
    // to do: review use of probStart and subtracting it in expDatePos formula

    int probDatePos    = probStart,                                 // position in probIndexes/prob.
        expProbDatePos = expProbIndexes.front() == 0? 1: 0,
        probDateIdx    = probIndexes[probDatePos],                  // date when we save survival probability
        expProbDateIdx = expProbIndexes[expProbDatePos],            // date when we save expected survival probability
        stopIdx        = Maths::min(probDateIdx, expProbDateIdx);   // when to do something

    double LnQ = 0.0;                                               // int lambda(u) du , lambda is the the piecewise constant implicit intensitiy process (equivalent to short rate in LMM)

    double mu               =  0.0,
           ti               =  0.0,                                 // current time in the simulation
           Tj               =  0.0,                                 // next reset date
           temp             =  0.0,
           drift            =  0.0,
           expnt            =  0.0,
           discountSurvProb =  1.,                                   // Q(t,T_i(t)) = surv prob from today to the next rest date
           lambda           =  0.0,
           frozenIntensity  = -1.0,                                  // implicit intensity, corresponds to the short rate in the standard libor market model
           alpha            =  0.0,                                  // interval fraction between event date and next reset date
           rootDelT,
           delT,
           W3;

    vector<double> simulInt(numInt,1),                              // contains the simulated H_k(t), where t = i is the simulation date
                   intensityAtResetDate(numInt+1,1);                // contains the lambda(0,T_1), lambda(T_1 ,T_2), ..., lambda(T_{n-1}, T_n)

    intensityAtResetDate[0] = - (1./mYearFracToFirstReset)*log(mSpotSurvProb[0]);     //The indexing of mAccruedFrac is such that mAccruedFrac[0] = [0, T_1];
    lambda                  =    intensityAtResetDate[0];                             //initial value relative to the period [0, t_1]

    frozenIntensity         = (1./mYearFracToFirstReset)*(1.0/mSpotSurvProb[0]-1);


    //Initialize forward intensities at time zero
    for (j=0; j < numInt; j++ )
    {
        simulInt[j] = mInitialIntensities[j];
    }

    for(i=0; i<numDates;/*increment in the body */)
    {
        W3       = randoms[i];
        rootDelT = sqrtYearFrac[i];
        delT     = rootDelT*rootDelT;
        ti      += delT;


        firstIntAliveIdx = -1;                                      // if t_j > T_n, i.e  there are no discrete intensities H_k alive to simulate


        for (j=0; j < numInt; j++)
        {
            Tj = mResetDates[j];
            if (ti <= Tj)
            {
                firstIntAliveIdx = j;
                break;
            }
            else
            {

            }
        }

        LnQ += lambda*rootDelT*rootDelT;                 // update LnQ


        if (firstIntAliveIdx == -1)                     //no H_k alive
        {
        // Do nothing: use flat lambda for lnQ

        }
        else
        {

            // evolve all state variables still alive
            mu = 0.0;
            // Simulates the intensities still alive, eq (5) in CreditMarketModelForSRM3
            for (j=firstIntAliveIdx; j < numInt; j++)
            {
                 temp                         = mAccruedFrac[j]*simulInt[j]*mIntensityVols[j]/(1+mAccruedFrac[j]*simulInt[j]);
                 mu                          += temp;
                 drift                        = mIntensityVols[j]*(mu); //change to adjust for risk neutral measure change!!
                 expnt                        = exp( (drift - 0.5*pow(mIntensityVols[j],2))*rootDelT*rootDelT+ mIntensityVols[j]* rootDelT * W3);
                 simulInt[j]                  = simulInt[j]* expnt;
            }

            //Check if the current date is a reset date. If so, freeze intensity
            if(ti<= mResetDates[firstIntAliveIdx] && ti+delT> mResetDates[firstIntAliveIdx])
            {
                frozenIntensity    = simulInt[firstIntAliveIdx];
                frozenIntensityIdx = firstIntAliveIdx;
                lambda = log(1.+mAccruedFrac[firstIntAliveIdx]*simulInt[firstIntAliveIdx])/mAccruedFrac[firstIntAliveIdx];
            }
        }




        i++;

        if (isWholeTimelineSurvivalRequested())
            wholeTimelineLogSurvProb[i] = LnQ;

        if(i + todayIdx == stopIdx)                                  // hit an "event", ie record sthg
        {
            if (stopIdx == probDateIdx)                              // save survival probability
            {
                prob[probDatePos] = LnQ;
                probDatePos++;
                probDateIdx = probIndexes[probDatePos];
            }

            if (stopIdx == expProbDateIdx)
            {
                if (firstIntAliveIdx == -1)
                {
                    //DO NOTHING
                }
                else
                {
                    alpha = (mResetDates[firstIntAliveIdx]-ti)/mAccruedFrac[frozenIntensityIdx];    // is the indexing correct????

                    temp = 1./(1 + alpha*mAccruedFrac[frozenIntensityIdx]*frozenIntensity);

                    expProb[expProbDatePos].survProb.resize(numInt,0.0);
                    for (j=firstIntAliveIdx; j < numInt; j++)
                    {
                        temp *= 1./(1.+mAccruedFrac[j]*simulInt[j]);                    //??????
                        expProb[expProbDatePos].survProb[j] = temp;
                    }


                    expProb[expProbDatePos].simulIntensity.resize(numInt,0.0);
                    expProb[expProbDatePos].firstIntAliveIdx   = firstIntAliveIdx;
                    expProb[expProbDatePos].frozenIntensityIdx = frozenIntensityIdx;
                    expProb[expProbDatePos].frozenIntensity    = frozenIntensity;

                    for (j=firstIntAliveIdx; j < numInt; j++)
                    {
                        expProb[expProbDatePos].simulIntensity[j] = simulInt[j];
                    }
                    expProbDatePos++;
                    expProbDateIdx = expProbIndexes[expProbDatePos];
                }
             }

         stopIdx = Maths::min(probDateIdx, expProbDateIdx); // refresh


    }

}
    temp = prob.size();
    for (unsigned int j = probStart; j < prob.size(); j++)
    {
        prob[j] = exp(-prob[j]);                        // then do exp on prob vector
    }

}

double SRMCreditLiborDiffuse::getLnExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j)
{
    // to avoid one virtual call specify the needed class explicitly
    return log(SRMCreditLiborDiffuse::getExpectedSurvivalDiscFactor(i, j));
}



double SRMCreditLiborDiffuse::getExpectedSurvivalDiscFactor(FwdIdx i, FwdIdx j)
{
    int k,
        fIAIdx,
        frozIdx,
        supResetIdx;

    const int iEDF        = getTimeLogic()->getReqEDFIdx(i),
            numEvent      = expProb.size();

    double ti = esdfForwardYearFracs[i],
           tj = esdfForwardYearFracs[j],
           dt = tj-ti,
           nextResetDate = -1.0,
           lastResetDate = -1.0,
           temp          = -1.0,
           survProb      = -1.0,
           alphaI        =  0.0,
           alphaK        =  0.0;

    lastResetDate = mResetDates.back();
    fIAIdx        = expProb[iEDF].firstIntAliveIdx;
    frozIdx       = expProb[iEDF].frozenIntensityIdx;
    nextResetDate = mResetDates[fIAIdx];


    // if ti > lasetResetDate, use the implied "short" intensity relative to the last accrual period
    if (ti >= lastResetDate)
    {
        // TO DO
        survProb = 1./(1+dt*expProb[iEDF].simulIntensity.back());   // CHECK!!
    }
    else if (tj <= nextResetDate)
    {
        if (frozIdx == -1)                      // ti < T_0
        {
            // TO DO
            // when tj is less the than the first rest date.
        }

        survProb = 1./(1+dt*expProb[iEDF].frozenIntensity);
    }
    else
    {

        // find the highest reset date such that tj > Tk
        for (k=fIAIdx; k<mNumInt; k++)
        {
            if (mResetDates[k]<tj)
            {

            }
            else
            {
                supResetIdx = k-1;
                break;
            }
        }

        temp     = expProb[iEDF].survProb[supResetIdx];
        alphaK   = (mResetDates[supResetIdx+1]-tj)/mAccruedFrac[supResetIdx];
        survProb = temp*(1+alphaK*mAccruedFrac[supResetIdx]*expProb[iEDF].simulIntensity[supResetIdx]);
    }

        return survProb;
}

DRLIB_END_NAMESPACE
