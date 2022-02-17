#include "edginc/config.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/SRMEquityDiffuseMappingMethod.hpp"
#include "edginc/IQMCRNGManager.hpp"

DRLIB_BEGIN_NAMESPACE  

double SRMEquityDiffuseMappingMethod::computeMappingMethodDrift(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, const vector<double> &cumulativeSigmaFX, bool isFirstDateIdx) const
{
    
    double shortRateDt = 0.0;
    double lnFwdEqDiff = 0.0;
    double spotVolSquareIntegral = 0.0;
    double cfactor = 0.0;
    double discreteDividend = 0.0;
    double cupsAdjust = 0.0;

    //*irLnMONEY[n] approximates the integral from 0 to n+1 of r_t (note irLnMoney[0] != 0 in general).
    //We may thus approximate r_n\delta_n as (*irLnMONEY)[n] - (*irLnMONEY)[n-1], where \delta_n is the
    //time interval from n to n+1.

    //LnFwdEq[n] is the log fwd at time n+1 (LnFwdEq[0] != 0 in general), so we approximate
    //\int d/dt LnFwdEq between times n and n+1 is LnFwdEq[n] - LnFwdEq[n-1]
    if (isFirstDateIdx)
    {
        shortRateDt = (*irLnMONEY)[0];
        lnFwdEqDiff = LnFwdEq[0] - origLnE;
        spotVolSquareIntegral = cumulativeSpotVol[0];
        cfactor = cumulativeCfactorArrayEQ[0];
        discreteDividend = cumulativeLnDivYield[0];

        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * cumulativeSigmaFX[0];
        }
    }
    else
    {
        shortRateDt = (*irLnMONEY)[simDateIdxB] - (*irLnMONEY)[simDateIdxA];
        lnFwdEqDiff = LnFwdEq[simDateIdxB] - LnFwdEq[simDateIdxA];
        spotVolSquareIntegral = cumulativeSpotVol[simDateIdxB] - cumulativeSpotVol[simDateIdxA];
        cfactor = cumulativeCfactorArrayEQ[simDateIdxB] - cumulativeCfactorArrayEQ[simDateIdxA];
        discreteDividend = cumulativeLnDivYield[simDateIdxB] - cumulativeLnDivYield[simDateIdxA];

        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * (cumulativeSigmaFX[simDateIdxB] - cumulativeSigmaFX[simDateIdxA]);
        }
    }
    
    double drift = 0.0;    
    drift += (shortRateDt - lnFwdEqDiff + cfactor + discreteDividend) / smileEQ;
    drift -= 0.5 * spotVolSquareIntegral * smileEQ; 
    drift += cupsAdjust;
    drift -= 0.5 * spotVolSquareIntegral * smileDashEQ;
    return drift;
}


double SRMEquityDiffuseMappingMethod::computeMappingMethodDrift3(int bigDateIdx, double smileEQ, double smileDashEQ) const
{
    const SRMEquityDiffuseMappingMethod *pThis = this;

    double spotVolSquareIntegral = cumulativeSpotVol2[bigDateIdx];
    double shortRateDt = 0.0;
    double lnFwdEqDiff = 0.0;
    double cfactor = 0.0;
    double discreteDividend = 0.0;
    double cupsAdjust = 0.0;

    double drift = 0.0;    
    drift += (shortRateDt - lnFwdEqDiff + cfactor + discreteDividend) / smileEQ;
    drift -= 0.5 * spotVolSquareIntegral * smileEQ; 
    drift += cupsAdjust;
    drift -= 0.5 * spotVolSquareIntegral * smileDashEQ;
    return drift;  
}



double SRMEquityDiffuseMappingMethod::computeCorrectedMappingMethodDrift3(int bigDateIdx, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac) const
{
    double t0 = timeSoFar;
    double t1 = timeSoFar + stepYearFrac;

    double shortRateDt = 0.0;
    double lnFwdEqDiff = 0.0;
    double spotVolSquareIntegral = cumulativeSpotVol2[bigDateIdx];
    double cfactor = 0.0;
    double discreteDividend = 0.0;
    double cupsAdjust = 0.0;

    double momShortRateDt = 0.0;
    double momLnFwdEqDiff = 0.0;
    double momSpotVolSquareIntegral = firstMomentSpotVol2[bigDateIdx];
    double momCfactor = 0.0;
    double momDiscreteDividend = 0.0;
    double momCupsAdjust = 0.0;


    
    /*
    Recall the Euler drift:

    double drift = 0.0;
    drift += (shortRateDt - lnFwdEqDiff + cfactor + discreteDividend) / smileEQ;         //LINE A
    drift -= 0.5 * spotVolSquareIntegral * smileEQ;                                      //LINE B
    drift += cupsAdjust;                                                                 //LINE C
    drift -= 0.5 * spotVolSquareIntegral * smileDashEQ;                                  //LINE D
    return drift;
    */


    double eulerDrift = computeMappingMethodDrift3(bigDateIdx, smileEQ, smileDashEQ);

    //Correction to Line A
    double mu = 1.0 / smileEQest - 1.0 / smileEQ;

    double a1 = shortRateDt - lnFwdEqDiff + cfactor + discreteDividend;
    double a2 = momShortRateDt - momLnFwdEqDiff + momCfactor + momDiscreteDividend;

    double driftA = mu * (a2 - t0 * a1);

    //Correction to Line B
    mu = smileEQest - smileEQ;

    double b1 = spotVolSquareIntegral;
    double b2 = momSpotVolSquareIntegral;

    double driftB = mu * (b2 - t0 * b1);

    //Line C adjustment simplifies, so add cupsAdjust later.

    //Correction to Line D
    mu = smileDashEQest - smileDashEQ;
    double driftD = mu * (b2 - t0 * b1);    //spotVolSquareIntegral as in Line B.

    double driftCorrection = (driftA - 0.5 * driftB - 0.5 * driftD) / stepYearFrac;

    //cups is fine in the Euler case
    return driftCorrection + eulerDrift;
}




double SRMEquityDiffuseMappingMethod::computeCorrectedMappingMethodDrift(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac, const vector<double> &cumulativeSigmaFX, const vector<double> &firstMomentSigmaFX, const vector<double> &firstMomentRt, bool isFirstDateIdx) const
{
    double t0 = timeSoFar;
    double t1 = timeSoFar + stepYearFrac;

    double shortRateDt = 0.0;
    double lnFwdEqDiff = 0.0;
    double spotVolSquareIntegral = 0.0;
    double cfactor = 0.0;
    double discreteDividend = 0.0;
    double cupsAdjust = 0.0;

    double momShortRateDt = 0.0;
    double momLnFwdEqDiff = 0.0;
    double momSpotVolSquareIntegral = 0.0;
    double momCfactor = 0.0;
    double momDiscreteDividend = 0.0;
    double momCupsAdjust = 0.0;


    //*irLnMONEY[n] approximates the integral from 0 to n+1 of r_t (note irLnMoney[0] != 0 in general).
    //We may thus approximate r_n\delta_n as (*irLnMONEY)[n] - (*irLnMONEY)[n-1], where \delta_n is the
    //time interval from n to n+1.

    //LnFwdEq[n] is the log fwd at time n+1 (LnFwdEq[0] != 0 in general), so we approximate
    //\int d/dt LnFwdEq between times n and n+1 is LnFwdEq[n] - LnFwdEq[n-1]
    if (isFirstDateIdx)
    {
        shortRateDt = (*irLnMONEY)[0];
        lnFwdEqDiff = LnFwdEq[0] - origLnE;
        spotVolSquareIntegral = cumulativeSpotVol[0];
        cfactor = cumulativeCfactorArrayEQ[0];
        discreteDividend = cumulativeLnDivYield[0];
        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * cumulativeSigmaFX[0];
        }

        momShortRateDt = firstMomentRt[0];
        momLnFwdEqDiff = firstMomentLnFwdEq[0];
        momSpotVolSquareIntegral = firstMomentSpotVol[0];
        momCfactor = firstMomentCFactor[0];
        momDiscreteDividend = firstMomentLnDiv[0];

        if (ccyTreatment != ccyVanilla) {
            momCupsAdjust = firstMomentSigmaFX[0];
        }
    }
    else
    {
        shortRateDt = (*irLnMONEY)[simDateIdxB] - (*irLnMONEY)[simDateIdxA];
        lnFwdEqDiff = LnFwdEq[simDateIdxB] - LnFwdEq[simDateIdxA];
        spotVolSquareIntegral = cumulativeSpotVol[simDateIdxB] - cumulativeSpotVol[simDateIdxA];
        cfactor = cumulativeCfactorArrayEQ[simDateIdxB] - cumulativeCfactorArrayEQ[simDateIdxA];
        discreteDividend = cumulativeLnDivYield[simDateIdxB] - cumulativeLnDivYield[simDateIdxA];

        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * (cumulativeSigmaFX[simDateIdxB] - cumulativeSigmaFX[simDateIdxA]);
        }

        momShortRateDt = firstMomentRt[simDateIdxB] - firstMomentRt[simDateIdxA];
        momLnFwdEqDiff = firstMomentLnFwdEq[simDateIdxB] - firstMomentLnFwdEq[simDateIdxA];
        momSpotVolSquareIntegral = firstMomentSpotVol[simDateIdxB] - firstMomentSpotVol[simDateIdxA];
        momCfactor = firstMomentCFactor[simDateIdxB] - firstMomentCFactor[simDateIdxA];
        momDiscreteDividend = firstMomentLnDiv[simDateIdxB] - firstMomentLnDiv[simDateIdxA];

        if (ccyTreatment != ccyVanilla) {
            momCupsAdjust = firstMomentSigmaFX[simDateIdxB] - firstMomentSigmaFX[simDateIdxA];
        }

    }
 /*
    Recall the Euler drift:

    double drift = 0.0;
    drift += (shortRateDt - lnFwdEqDiff + cfactor + discreteDividend) / smileEQ;         //LINE A
    drift -= 0.5 * spotVolSquareIntegral * smileEQ;                                      //LINE B
    drift += cupsAdjust;                                                                 //LINE C
    drift -= 0.5 * spotVolSquareIntegral * smileDashEQ;                                  //LINE D
    return drift;
*/

    //Correction to Line A
    double mu1 = t1 / smileEQ - t0 / smileEQest;
    double mu2 = 1.0 /smileEQ - 1.0 / smileEQest;

    double a1 = shortRateDt - lnFwdEqDiff + cfactor + discreteDividend;
    double a2 = momShortRateDt - momLnFwdEqDiff + momCfactor + momDiscreteDividend;

    double driftA = mu1 * a1 - mu2 * a2;

    //Correction to Line B
    mu1 = t1 * smileEQ - t0 * smileEQest;
    mu2 = smileEQ - smileEQest;

    double b1 = spotVolSquareIntegral;
    double b2 = momSpotVolSquareIntegral;

    double driftB = mu1 * b1 - mu2 * b2;

    //Line C adjustment simplifies, so add cupsAdjust later.

    //Correction to Line D
    mu1 = t1 * smileDashEQ - t0 * smileDashEQest;
    mu2 = smileDashEQ - smileDashEQest;

    double driftD = mu1 * b1 - mu2 * b2;    //spotVolSquareIntegral as in Line B.

    double drift = (driftA - 0.5 * driftB - 0.5 * driftD) / stepYearFrac;
    drift += cupsAdjust;
    return drift;
}



double SRMEquityDiffuseMappingMethod::computeCorrectedMappingMethodDrift2(int simDateIdxA, int simDateIdxB, double smileEQ, double smileDashEQ, double smileEQest, double smileDashEQest, double timeSoFar, double stepYearFrac, const vector<double> &cumulativeSigmaFX, const vector<double> &firstMomentSigmaFX, const vector<double> &firstMomentRt, bool isFirstDateIdx) const
{
    double t0 = timeSoFar;
    double t1 = timeSoFar + stepYearFrac;

    double shortRateDt = 0.0;
    double lnFwdEqDiff = 0.0;
    double spotVolSquareIntegral = 0.0;
    double cfactor = 0.0;
    double discreteDividend = 0.0;
    double cupsAdjust = 0.0;

    double momShortRateDt = 0.0;
    double momLnFwdEqDiff = 0.0;
    double momSpotVolSquareIntegral = 0.0;
    double momCfactor = 0.0;
    double momDiscreteDividend = 0.0;
    double momCupsAdjust = 0.0;


    //*irLnMONEY[n] approximates the integral from 0 to n+1 of r_t (note irLnMoney[0] != 0 in general).
    //We may thus approximate r_n\delta_n as (*irLnMONEY)[n] - (*irLnMONEY)[n-1], where \delta_n is the
    //time interval from n to n+1.

    //LnFwdEq[n] is the log fwd at time n+1 (LnFwdEq[0] != 0 in general), so we approximate
    //\int d/dt LnFwdEq between times n and n+1 is LnFwdEq[n] - LnFwdEq[n-1]
    if (isFirstDateIdx)
    {
        shortRateDt = (*irLnMONEY)[0];
        lnFwdEqDiff = LnFwdEq[0] - origLnE;
        spotVolSquareIntegral = cumulativeSpotVol[0];
        cfactor = cumulativeCfactorArrayEQ[0];
        discreteDividend = cumulativeLnDivYield[0];
        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * cumulativeSigmaFX[0];
        }

        momShortRateDt = firstMomentRt[0];
        momLnFwdEqDiff = firstMomentLnFwdEq[0];
        momSpotVolSquareIntegral = firstMomentSpotVol[0];
        momCfactor = firstMomentCFactor[0];
        momDiscreteDividend = firstMomentLnDiv[0];

        if (ccyTreatment != ccyVanilla) {
            momCupsAdjust = firstMomentSigmaFX[0];
        }
    }
    else
    {
        shortRateDt = (*irLnMONEY)[simDateIdxB] - (*irLnMONEY)[simDateIdxA];
        lnFwdEqDiff = LnFwdEq[simDateIdxB] - LnFwdEq[simDateIdxA];
        spotVolSquareIntegral = cumulativeSpotVol[simDateIdxB] - cumulativeSpotVol[simDateIdxA];
        cfactor = cumulativeCfactorArrayEQ[simDateIdxB] - cumulativeCfactorArrayEQ[simDateIdxA];
        discreteDividend = cumulativeLnDivYield[simDateIdxB] - cumulativeLnDivYield[simDateIdxA];

        if (ccyTreatment != ccyVanilla) {
            cupsAdjust = - corrEqFX * (cumulativeSigmaFX[simDateIdxB] - cumulativeSigmaFX[simDateIdxA]);
        }

        momShortRateDt = firstMomentRt[simDateIdxB] - firstMomentRt[simDateIdxA];
        momLnFwdEqDiff = firstMomentLnFwdEq[simDateIdxB] - firstMomentLnFwdEq[simDateIdxA];
        momSpotVolSquareIntegral = firstMomentSpotVol[simDateIdxB] - firstMomentSpotVol[simDateIdxA];
        momCfactor = firstMomentCFactor[simDateIdxB] - firstMomentCFactor[simDateIdxA];
        momDiscreteDividend = firstMomentLnDiv[simDateIdxB] - firstMomentLnDiv[simDateIdxA];

        if (ccyTreatment != ccyVanilla) {
            momCupsAdjust = firstMomentSigmaFX[simDateIdxB] - firstMomentSigmaFX[simDateIdxA];
        }

    }
    /*
    Recall the Euler drift:

    double drift = 0.0;
    drift += (shortRateDt - lnFwdEqDiff + cfactor + discreteDividend) / smileEQ;         //LINE A
    drift -= 0.5 * spotVolSquareIntegral * smileEQ;                                      //LINE B
    drift += cupsAdjust;                                                                 //LINE C
    drift -= 0.5 * spotVolSquareIntegral * smileDashEQ;                                  //LINE D
    return drift;
    */


    double eulerDrift = computeMappingMethodDrift(simDateIdxA, simDateIdxB, smileEQ, smileDashEQ, cumulativeSigmaFX, isFirstDateIdx);

    //Correction to Line A
    double mu = 1.0 / smileEQest - 1.0 / smileEQ;
    
    double a1 = shortRateDt - lnFwdEqDiff + cfactor + discreteDividend;
    double a2 = momShortRateDt - momLnFwdEqDiff + momCfactor + momDiscreteDividend;

    double driftA = mu * (a2 - t0 * a1);

    //Correction to Line B
    mu = smileEQest - smileEQ;
   
    double b1 = spotVolSquareIntegral;
    double b2 = momSpotVolSquareIntegral;

    double driftB = mu * (b2 - t0 * b1);

    //Line C adjustment simplifies, so add cupsAdjust later.

    //Correction to Line D
    mu = smileDashEQest - smileDashEQ;
    double driftD = mu * (b2 - t0 * b1);    //spotVolSquareIntegral as in Line B.

    double driftCorrection = (driftA - 0.5 * driftB - 0.5 * driftD) / stepYearFrac;

    //cups is fine in the Euler case
    return driftCorrection + eulerDrift;
}


/** We transform the equity SDE by:
1) Z_t = log(E_t / E_0^t) = log(E_t) - log(E_0^t)
2) K_t = K(t, Z_t)
Step 2 makes the volatility equal to 1, for suitable K (depending on the smile)
*/

/** We diffuse Y_t = K(t, Z_t)
Note that K is piecewise constant as a function of t.
We diffuse with d/dt K = 0 during these constant periods, then apply a
continuity constraint on Y_t across the discontinuity of K (as a function of t). */

void SRMEquityDiffuseMappingMethod::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    int numDates = sigmaEQ.size();
    randoms = rngMgr->getCorrelatedRandoms(randomIndex); // save pointer to randoms
    double Zt = 0.0;
    double Kt = 0.0;    //We diffuse the K transformed log moneyness...

    spotEqPos = spotEqStart; // position in spotEQ/spotEQIndexes
    expSpotEqPos = 0;
    stopIdx = Maths::min(spotEqIndexes[spotEqPos],
        expSpotEqIndexes[expSpotEqPos]);
    divPos = 0; // position in exDivIndexes

    //TODO do this in the rates diffusion:
    //integral of (short rate times t)
    vector<double> firstMomentRt(numDates);
    double timeSoFar = 0.0;
    double doubleIntegralSoFar = 0.0;
    for (int simDateIdx = 0; simDateIdx < numDates;simDateIdx++)
    {
        double delta_t = yearFrac[simDateIdx];
        timeSoFar += delta_t;

        double prev = 0;
        if (simDateIdx > 0)
            prev = (*irLnMONEY)[simDateIdx-1];

        doubleIntegralSoFar += 0.5 * delta_t * (prev + (*irLnMONEY)[simDateIdx]);

        //integration by parts:
        firstMomentRt[simDateIdx] = (*irLnMONEY)[simDateIdx] * timeSoFar - doubleIntegralSoFar;
    }


    //sigmaFX is path dependent, so unfortunately we cannot pre-compute the integral in finalize.
    vector<double> cumulativeSigmaFX(numDates);
    vector<double> firstMomentSigmaFX(numDates);

    //TODO: perform this integral inside the FX diffusion
    if (ccyTreatment != ccyVanilla) {

        //Loop over ALL dates, not just the big steps
        double totalFX = 0.0;
        timeSoFar = 0.0;
        double firstMomentSoFar = 0.0;
        for (int simDateIdx = 0; simDateIdx < numDates;simDateIdx++)
        {
            double delta_t = yearFrac[simDateIdx];
            
            double inc = (*sigmaFX)[simDateIdx] * SpotVol[simDateIdx] * delta_t;
            totalFX += inc;
            cumulativeSigmaFX[simDateIdx] = totalFX;

            //evaluate at left hand side....
            firstMomentSoFar += timeSoFar * inc;
            firstMomentSigmaFX[simDateIdx] = firstMomentSoFar;

            timeSoFar += delta_t;
        }
    }

    int smileChangePos = 0;
    int nextSmileDateIdx = -1;

    //Note 0 is a smile change idx...
    if (smileChanges.size() > 1)
        nextSmileDateIdx = smileChanges[1];

    timeSoFar = 0.0;
    //Recall the last big step marks the end of the diffusion: don't diffuse past it (hence bigStepIdxs.size() - 1)
    for (size_t bigStep = 0;bigStep<bigStepIdxs.size()-1;bigStep++)
    {
        int simDateIdx = bigStepIdxs[bigStep];
        int simDateIdxNextBigStep = bigStepIdxs[bigStep+1];

        int prevBigStepIdx = 0;
        if (bigStep > 0)
            prevBigStepIdx = bigStepIdxs[bigStep-1];
        double bigYearFrac = bigStepYearFracs[bigStep];  
        double bigYearFracSqrt = bigStepYearFracsSqrt[bigStep];

        //We add up the random numbers over the interval we have 'skipped'. This is good for two reasons:
        //1) It ensures the correlation structure (with FX, rates, etc) is preserved,
        //2) It means the results are closer to the Euler method results in a 'strong' sense (i.e. pathwise). Good for testing.

        double diffusionTerm = 0.0;
        for (int j = bigStepIdxs[bigStep];j<bigStepIdxs[bigStep+1];j++)
        {
            diffusionTerm += sqrtYearFrac[j] * randoms[j] * SpotVol[j];
        }

        
        double spotVolSquareIntegral = 0.0;
   //     if (simDateIdx == 0)
     //       spotVolSquareIntegral = cumulativeSpotVol[0];
     //   else
      //      spotVolSquareIntegral = cumulativeSpotVol[simDateIdx] - cumulativeSpotVol[prevBigStepIdx];

        spotVolSquareIntegral = cumulativeSpotVol2[bigStep];

        double sqrtSpotVolSquareIntegral = sqrt(spotVolSquareIntegral);

     
        //Kt +=  d/dt K * delta_t - 0.5 * s' * spotVol^2 * delta_t + 1/s * (ir - q_t - dlog FwdEq - 0.5 * s * s * spotVol * spotVol) * delta_t 
        //      + xrandom * spotVol  * sqrt_delta_t; 

        double smileEQ, smileDashEQ, smileEQest, smileDashEQest;
        if (simDateIdx == 0)
        {
            smileEQ = 1.0;
            smileDashEQ = EqSmile_a1[0] * EqSmile_a3[0];
        }
        else
        {
            smileStates[smileChangePos].evaluateLocalVol(Zt, &smileEQ, &smileDashEQ);            
        }

        double driftEstimate = computeMappingMethodDrift3(bigStep, smileEQ, smileDashEQ);

        //Predictor-corrector: estimate next value of K, then feed into drift...
        double Kestimate = Kt + driftEstimate + diffusionTerm; 
        double Zestimate;
        smileStates[smileChangePos].evaluateInverseKExp(Kestimate, &Zestimate);
        smileStates[smileChangePos].evaluateLocalVol(Zestimate, &smileEQest, &smileDashEQest);

       
        
        //1) Old flawed mapping method + pc
        //drift = computeMappingMethodDrift(prevBigStepIdx, simDateIdx, smileEQest, smileDashEQest, cumulativeSigmaFX, (simDateIdx == 0));
        //Kt += 0.5 * (drift + driftEstimate) + x_random * sqrtSpotVolSquareIntegral; 
        
        
        //2) no pc
        //Kt += driftEstimate + diffusionTerm; //no pc
        //double tmp = Kt + LnFwdEq[simDateIdxNextBigStep-1];


        //3) New mapping method (fixed):
        //drift = computeCorrectedMappingMethodDrift(prevBigStepIdx, simDateIdx, smileEQ, smileDashEQ, smileEQest, smileDashEQest, timeSoFar, bigYearFrac, cumulativeSigmaFX, firstMomentSigmaFX, firstMomentRt, (simDateIdx == 0));
        //double driftt = computeCorrectedMappingMethodDrift2(prevBigStepIdx, simDateIdx, smileEQ, smileDashEQ, smileEQest, smileDashEQest, timeSoFar, bigYearFrac, cumulativeSigmaFX, firstMomentSigmaFX, firstMomentRt, (simDateIdx == 0));
        //Kt += driftt + x_random * sqrtSpotVolSquareIntegral; 


        //4) Totally corrected indices
        double driftt = computeCorrectedMappingMethodDrift3(bigStep, smileEQ, smileDashEQ, smileEQest, smileDashEQest, timeSoFar, bigYearFrac);
        Kt += driftt + diffusionTerm; 

        //5) Corrected indices, old flawed method....
        //double drift = computeMappingMethodDrift3(bigStep, smileEQest, smileDashEQest);
        //Kt += 0.5 * (drift + driftEstimate) + diffusionTerm; 



        //K_t now up to date... invert to get Z_t
        smileStates[smileChangePos].evaluateInverseKExp(Kt, &Zt);

        if (simDateIdxNextBigStep + todayIdx == stopIdx){
            // must save the current spot EQ (or actually log of it here)

            //We are diffusing the log moneyness, so need to add the log forward
            //Note: The LnFwdEq array is offset by 1: see SRMEquityUtil::setTimeLine
            double lnSpotEQ = Zt + LnFwdEq[simDateIdxNextBigStep-1];

            if (ccyTreatment == ccyStruck) {
                lnSpotEQ += (*spotFXFullPath)[simDateIdxNextBigStep];
            }

            if (stopIdx == spotEqIndexes[spotEqPos]){
                spotEq[spotEqPos] = SRMExp(lnSpotEQ);
                spotEqPos++;
            }
            if (stopIdx == expSpotEqIndexes[expSpotEqPos]){
                // Inconsistent treatment versus above but having lnSpotEQ
                // is more useful - then exp() is called only once in
                // SRMEquityDiffuse::ExpSpotSV::calculatePath
                expSpotEq[expSpotEqPos] = lnSpotEQ;
                expSpotEqPos++;
            }
            stopIdx = Maths::min(spotEqIndexes[spotEqPos],
                expSpotEqIndexes[expSpotEqPos]);
        }

        if (simDateIdxNextBigStep < (int)SpotVol.size())
        {
            //Does the smile change shape at the next date?
            if (simDateIdxNextBigStep < numDates && simDateIdxNextBigStep == nextSmileDateIdx)  
            {
                if (smileChangePos < (int)smileChanges.size())
                {
                    ++smileChangePos;
                    nextSmileDateIdx = smileChanges[smileChangePos+1];
                }
                else
                    nextSmileDateIdx = -1;

                //Need to diffuse Kt across the smile discontinuity.
                //We have already inverted Kt using the old smile to obtain Zt.
                //Update Kt using Zt and new smile values
                smileStates[smileChangePos].evaluateKExp(Zt, &Kt);
            }
        }

        timeSoFar += bigYearFrac;
    } //for()
}

void SRMEquityDiffuseMappingMethod::finalize(DateTimeArrayConstSP allDates)
{
    //base class finalization
    SRMEquityDiffuse::finalize(allDates);

    double timeSoFar = 0.0;
    double firstMomentSoFar = 0.0;
    double doubleIntegralSoFar = 0.0;
    double integralSoFar = 0.0;
    double totalSoFar = 0.0;
    size_t i = 0;

    //compute first moment of LnFwdEq, using integration by parts:
    firstMomentLnFwdEq.resize(SpotVol.size());
    for (i = 0;i<SpotVol.size();i++)
    {
        timeSoFar += yearFrac[i];
        doubleIntegralSoFar += yearFrac[i] * LnFwdEq[i];
        firstMomentLnFwdEq[i] = timeSoFar * LnFwdEq[i] - doubleIntegralSoFar;
    }

    //Compute cumulative spot vol (integral of spot vol squared)
    cumulativeSpotVol.resize(SpotVol.size());
    firstMomentSpotVol.resize(SpotVol.size());
    timeSoFar = 0.0;
    firstMomentSoFar = 0.0;
    for (i = 0;i<SpotVol.size();i++)
    {    
        integralSoFar += yearFrac[i] * SpotVol[i] * SpotVol[i];
        cumulativeSpotVol[i] = integralSoFar;

        double left = SpotVol[i] * SpotVol[i];
        double right = left;
        if (i < SpotVol.size()-1)
            right = SpotVol[i+1] * SpotVol[i+1];

        firstMomentSoFar +=  0.5 * (timeSoFar * left + (timeSoFar + yearFrac[i])* right) * yearFrac[i];
        firstMomentSpotVol[i] = firstMomentSoFar;
        timeSoFar += yearFrac[i];
    }

    
    //sum the dividend yields...
    cumulativeLnDivYield.resize(sigmaEQ.size());
    firstMomentLnDiv.resize(sigmaEQ.size());
    int divPos = 0;
    for (i = 0;i<sigmaEQ.size();i++)
    {
        double delta_t = yearFrac[i];

        //We use i+1 here: in the old Euler loop the ex-dividend was taken at date simDateIdx + 1
        if (i+1 == exDivIndexes[divPos]) {
            // no dollar divs as such, no pv around div pmt date vs div ex date
            // so the formula is a bit simpler and no need to go via actual fwd eq
            totalSoFar += lnDivYield[divPos];
            firstMomentSoFar += lnDivYield[divPos] * timeSoFar;
            divPos++;
        }
        cumulativeLnDivYield[i] = totalSoFar;
        firstMomentLnDiv[i] = firstMomentSoFar;
        timeSoFar += delta_t;
    }

    timeSoFar = 0.0;
    firstMomentSoFar = 0.0;

    //integral of CfactorEQ
    cumulativeCfactorArrayEQ.resize(CfactorArrayEQ.size());
    firstMomentCFactor.resize(CfactorArrayEQ.size());
    double cumulativeCFactorSoFar = 0.0;
    for (i = 0;i<CfactorArrayEQ.size();i++)
    {
        //Note: yearFrac already included in CFactorArrayEQ
        cumulativeCFactorSoFar += CfactorArrayEQ[i];
        cumulativeCfactorArrayEQ[i] = cumulativeCFactorSoFar;

        firstMomentSoFar += timeSoFar * CfactorArrayEQ[i];
        firstMomentCFactor[i] = firstMomentSoFar;
        timeSoFar += yearFrac[i];
    }

    smileChanges.clear();
    double atmp = 0.0;
    double btmp = 0.0;
    double ctmp = 0.0;

    for (int i = 0;i<(int)EqSmile_a1.size();i++)
    {
        if (i==0 || atmp != EqSmile_a1[i] || btmp != EqSmile_a2[i] || ctmp != EqSmile_a3[i])
        {
            smileChanges.push_back(i);
            atmp = EqSmile_a1[i];
            btmp = EqSmile_a2[i];
            ctmp = EqSmile_a3[i];

            HyperLocalVolState localVol(EqSmile_a1[i], EqSmile_a2[i], EqSmile_a3[i]);
            smileStates.push_back(localVol);
        }
    }



    vector<int> bigStepBools;
    size_t numBigSteps = setupCriticalBigStepBools(bigStepBools);

    //Finally we can assemble the list of 'big steps'
    bigStepYearFracs.clear();
    bigStepYearFracsSqrt.clear();
    setupBigStepIndicesAndYearFracs(bigStepBools, bigStepIdxs, bigStepYearFracs, bigStepYearFracsSqrt);

    calculateBigStepIntegral(SpotVol, true, bigStepBools, cumulativeSpotVol2);
    calculateBigStepFirstMoment(SpotVol, true, bigStepBools, firstMomentSpotVol2);

    /*

    //work out where the big steps are
    bigStepIdxs.clear();
    int smileChangePos = 0;
    int nextSmileDateIdx = -1;

    spotEqPos = spotEqStart; // position in spotEQ/spotEQIndexes
    expSpotEqPos = 0;
    stopIdx = Maths::min(spotEqIndexes[spotEqPos],
        expSpotEqIndexes[expSpotEqPos]);
    
    //Note 0 is a smile change idx...
    if (smileChanges.size() > 1)
        nextSmileDateIdx = smileChanges[1];


    
    double yearFracTotal = 0.0;
    int numDates = sigmaEQ.size();
    for (int simDateIdx = 0; simDateIdx < numDates; ++simDateIdx)
    {
        
        bool isBigStep = false;
        bool smileChange = false;   
        if (simDateIdx + 1 < numDates && simDateIdx + 1 == nextSmileDateIdx)
            smileChange = true;

        //Various conditions may force us to update K:
        if (smileChange)
            isBigStep = true;

        if (simDateIdx + 1 + todayIdx == stopIdx)
            isBigStep = true;

        if (simDateIdx == 0)
            isBigStep = true;

        if (isBigStep)
        {
            bigStepIdxs.push_back(simDateIdx);
            bigStepYearFracs.push_back(yearFracTotal);
            yearFracTotal = 0.0;
        }
        double yf = yearFrac[simDateIdx];
        yearFracTotal += yf;

        //Does the smile change shape at the next date?
        if (smileChange)  
        {
            if (smileChangePos < (int)smileChanges.size())
            {
                ++smileChangePos;
                nextSmileDateIdx = smileChanges[smileChangePos+1];
            }
            else
                nextSmileDateIdx = -1;
        }

        if (simDateIdx + 1 + todayIdx == stopIdx){
            if (stopIdx == spotEqIndexes[spotEqPos])
                spotEqPos++;
            if (stopIdx == expSpotEqIndexes[expSpotEqPos]) 
                expSpotEqPos++;
            stopIdx = Maths::min(spotEqIndexes[spotEqPos],
                expSpotEqIndexes[expSpotEqPos]);
        }
    }

    size_t numBig = bigStepYearFracs.size();
    for (size_t i = 0;i<numBig - 1;i++)
    {
        bigStepYearFracs[i] = bigStepYearFracs[i+1];
//        bigStepYearFracsSqrt[i] = bigStepYearFracsSqrt[i+1];
    }
    bigStepYearFracs[numBig-1] = yearFrac[numDates-1];
  //  bigStepYearFracsSqrt[numBig-1] = sqrtYearFrac[numDates-1];

  */

}

void SRMEquityDiffuseMappingMethod::calculateBigStepFirstMoment(const vector<double> &f, bool square, const vector<int> &bigStepBools, vector<double> &bsFirstMoment)
{

    size_t numDates = sigmaEQ.size();
    double firstMomentSoFar = 0.0;
    double timeSoFar = 0.0;
    for (size_t simDateIdx = 0; simDateIdx < numDates+1; ++simDateIdx)
    {
        //cfactor integral at time n should include CfactorArrayEQ[n]
        if (bigStepBools[simDateIdx])
        {
            bsFirstMoment.push_back(firstMomentSoFar);
            firstMomentSoFar = 0.0;
        }
        if (simDateIdx < numDates)
        {
            double left = 0.0;
            double right = 0.0;

            if (square)
            {
                left = f[simDateIdx] * f[simDateIdx];
                right = left;
                if (simDateIdx < numDates - 1)
                    right = f[simDateIdx+1] * f[simDateIdx+1];
            }
            else
            {
                left = f[simDateIdx];
                right = left;
                if (simDateIdx < numDates - 1)
                    right = f[simDateIdx+1];
            }

            firstMomentSoFar +=  0.5 * (timeSoFar * left + (timeSoFar + yearFrac[simDateIdx])* right) * yearFrac[simDateIdx];
            timeSoFar += yearFrac[simDateIdx];
        }        
    }

    size_t numBig = bigStepYearFracs.size();
    for (size_t i = 0;i<numBig - 1;i++)
    {
        bsFirstMoment[i] = bsFirstMoment[i+1];
    }

}

void SRMEquityDiffuseMappingMethod::calculateBigStepIntegral(const vector<double> &f, bool square, const vector<int> &bigStepBools, vector<double> &bsIntegral)
{
    size_t numDates = sigmaEQ.size();
    double integralTotal = 0.0;
    for (size_t simDateIdx = 0; simDateIdx < numDates+1; ++simDateIdx)
    {
        //cfactor integral at time n should include CfactorArrayEQ[n]
        if (bigStepBools[simDateIdx])
        {
            bsIntegral.push_back(integralTotal);
            integralTotal = 0.0;
        }
        if (simDateIdx < numDates)
        {
            if (square)
            {
                integralTotal += f[simDateIdx] * f[simDateIdx] * yearFrac[simDateIdx];
            }
            else
            {
                integralTotal += f[simDateIdx] * yearFrac[simDateIdx];
            }
        }        
    }

    size_t numBig = bigStepYearFracs.size();
    for (size_t i = 0;i<numBig - 1;i++)
    {
        bsIntegral[i] = bsIntegral[i+1];
    }
}

DRLIB_END_NAMESPACE  
