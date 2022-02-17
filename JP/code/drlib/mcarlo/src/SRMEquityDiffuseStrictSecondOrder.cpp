#include "edginc/config.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEquityDiffuseSecondOrder.hpp"
#include "edginc/SRMEquityDiffuseStrictSecondOrder.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/IQMCRNGManager.hpp"

DRLIB_BEGIN_NAMESPACE 

//JDee: I _think_ there is a bug in the maths of the 2nd order approx: BUGFIX_FACTOR should be 0.5
//Need to lookup in Kloeden and Platen.
#define BUGFIX_FACTOR 0.5

/** Simulate a new path = CalcStateQty_EQ + DiffuseEQ */
void SRMEquityDiffuseStrictSecondOrder::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    int numDates = this->sigmaEQ.size();

    randoms = rngMgr->getCorrelatedRandoms(randomIndex); // save pointer to randoms
    LnE = origLnE; // reset
    spotEqPos = spotEqStart; // position in spotEQ/spotEQIndexes
    expSpotEqPos = 0;
    stopIdx = Maths::min(spotEqIndexes[spotEqPos],
        expSpotEqIndexes[expSpotEqPos]);
    divPos = 0; // position in exDivIndexes

    double sigmaEQ = this->sigmaEQ[0];
    double volDashEQ = SpotVol[0] * EqSmile_a1[0] * EqSmile_a3[0];
    double volDashDashEQ = SpotVol[0] * EqSmile_a2[0] * EqSmile_a3[0] * EqSmile_a3[0];

    // loop over all dates
    // loop over interesting dates
    // carry out integration of quantities between interesting dates

    //Recall the last big step marks the end of the diffusion: don't diffuse past it (hence bigStepIdxs.size() - 1)
    for (size_t i = 0; i < bigStepIdxs.size() - 1; ++i/* DO NOT increment in body ;)*/){

        int simDateIdx = bigStepIdxs[i];
        int simDateIdxNextBigStep = bigStepIdxs[i+1];

        double delta_t = bigStepYearFracs[i];
        double sqrt_delta_t = bigStepYearFracsSqrt[i];

        //We add up the random numbers over the interval we have 'skipped'. This is good for two reasons:
        //1) It ensures the correlation structure (with FX, rates, etc) is preserved,
        //2) It means the results are closer to the Euler method results in a 'strong' sense (i.e. pathwise). Good for testing.

        double x_random_integral = 0.0;
        for (int j = bigStepIdxs[i];j<bigStepIdxs[i+1];j++)
        {
            x_random_integral += sqrtYearFrac[j] * randoms[j];
        }

        
        //Compute new values for vol, voldash, voldashdash-++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        double spotvol = SpotVol[simDateIdx];
        sigmaEQ = spotvol;
        volDashEQ = 0.0;
        volDashDashEQ = 0.0;

        bool bNontrivialSmile = false;
        if (!smileIsFlatEverywhere && simDateIdx > 0){
            double A_Smile = EqSmile_a1[simDateIdx]; /* skew  */
            double B_Smile = EqSmile_a2[simDateIdx]; /* smile */
            if (A_Smile != 0.0 || B_Smile != 0)
            {
                bNontrivialSmile = true;

                //irLnMoney and LnFwdEq are offset by 1.
                double x = LnE + (*irLnMONEY)[simDateIdx-1] - LnFwdEq[simDateIdx-1];
                double C_Smile = EqSmile_a3[simDateIdx];
                double w = C_Smile * x;
                const double SMILE_CUTOFF = 100.0;
                if (w < SMILE_CUTOFF && w > - SMILE_CUTOFF) {
                    double y = exp(w);
                    double z = 1.0/y;
                    double tanh_cx = (y-z)/(y+z);
                    double tanh_cx_sqr = tanh_cx * tanh_cx;
                    double sech_cx = 2.0/(y+z);
                    /* loc vol = spotvol * (1 + A * Tanh(C * x) +
                    *                   B * ( 1 - 1/Cosh(C * x)) ) */
                    sigmaEQ = spotvol * (1.0 + A_Smile * tanh_cx
                        + B_Smile * (1.0 - sech_cx));

                    volDashEQ = spotvol * C_Smile * ( A_Smile * (1.0 - tanh_cx_sqr) + B_Smile * tanh_cx * sech_cx );
                    volDashDashEQ = spotvol * C_Smile * C_Smile
                        * ( -2.0 * A_Smile * tanh_cx * (1.0 - tanh_cx_sqr)
                        + B_Smile * ( (1.0 - 2.0 * tanh_cx_sqr) * sech_cx) );

                } else if (w <= - SMILE_CUTOFF) {
                    /* y = 0 */
                    sigmaEQ = spotvol * (1.0 - A_Smile + B_Smile);
                    volDashEQ = 0.0;
                    volDashDashEQ = 0.0;
                } else {
                    /* z = 0 */
                    sigmaEQ = spotvol * (1.0 + A_Smile + B_Smile);
                    volDashEQ = 0.0;
                    volDashDashEQ = 0.0;
                }
            }

            sigmaEQ = Maths::max(sigmaEQ, 0.01 * spotvol);
        }

        double mu_1 = - 0.5 * sigmaEQ * sigmaEQ;
        double mu_2 = cfactorBigIntegrals[i];
        
        /* CUPS adjust if necessary
        actually could use sigmaFX as a flag here */
        if (ccyTreatment != ccyVanilla) {
            //  mu_2 -= corrEqFX * (*sigmaFX)[simDateIdx] * sigmaEQ * yearFrac[simDateIdx];
            mu_1 -= corrEqFX * (*sigmaFX)[simDateIdx] * sigmaEQ;
        }


        // Some things about dividends - we use simDateIdx+1
        // so that a div on this date is treated as ex- and not cum-
        // eqdiffuse.c lines 445 - 492
        if (simDateIdx+1 == exDivIndexes[divPos]) {
            // no dollar divs as such, no pv around div pmt date vs div ex date
            // so the formula is a bit simpler and no need to go via actual fwd eq
            mu_2 += lnDivYield[divPos];
            divPos++;
        }

        //Note: mu_dt has delta_t factor
        double mu_dt = mu_1 * delta_t + mu_2;

        /* Main formula for 2nd order discretisation */

        //1st order terms
        //double first_order = sigmaEQ * x_random * sqrt_delta_t + mu_dt;
        double first_order = sigmaEQ * x_random_integral + mu_dt;

        //TODO: only compute if smile is nonzero
        double second_order = 0;
        if (bNontrivialSmile)
        {
            double mu_dash = - sigmaEQ * volDashEQ; //first derivative of 'drift' with respect to VOL
            double mu_dashdash = - volDashEQ * volDashEQ - sigmaEQ * volDashDashEQ;

            /*
            * Natural arrangement of terms (we shall reorder the terms to make a small optimisation)
            second_order += sigmaEQ * volDashEQ * (x_random * x_random - 1) * delta_t * 0.5;
            second_order += mu_dt * volDashEQ * x_random * sqrt_delta_t * 0.5       //mu already has a delta_t factor
            + (sigmaEQ * sigmaEQ * volDashDashEQ * 0.5
            - sigmaEQ * sigmaEQ * volDashEQ) * x_random * sqrt_delta_t * delta_t * 0.5;
            second_order -= mu_dt * sigmaEQ * volDashEQ * delta_t * 0.5;
            second_order -= (volDashEQ * volDashEQ + sigmaEQ * volDashDashEQ)  * sigmaEQ * sigmaEQ * delta_t * delta_t * BUGFIX_FACTOR * 0.5;
            */

            double sigmaEQSquare = sigmaEQ * sigmaEQ;

            /*
            second_order += sigmaEQ * volDashEQ * 0.5 * (delta_t * x_random * x_random - delta_t - delta_t * mu_dt);
            second_order += mu_dt * volDashEQ * x_random * sqrt_delta_t * 0.5       //mu already has a delta_t factor
                + sigmaEQSquare * (volDashDashEQ * 0.5
                - volDashEQ) * x_random * sqrt_delta_t * delta_t * 0.5;
            second_order -= (volDashEQ * volDashEQ + sigmaEQ * volDashDashEQ)  * sigmaEQSquare * delta_t * delta_t * BUGFIX_FACTOR * 0.5;
            */

            second_order += sigmaEQ * volDashEQ * 0.5 * (x_random_integral * x_random_integral - delta_t - delta_t * mu_dt);
            second_order += mu_dt * volDashEQ * x_random_integral * 0.5       //mu already has a delta_t factor
                + sigmaEQSquare * (volDashDashEQ * 0.5
                - volDashEQ) * x_random_integral * delta_t * 0.5;
            second_order -= (volDashEQ * volDashEQ + sigmaEQ * volDashDashEQ)  * sigmaEQSquare * delta_t * delta_t * BUGFIX_FACTOR * 0.5;

        }

        LnE += first_order;
        LnE += second_order;

        /* End of main formula for 2nd order discretisation */

        // Store LnE if needed
        if (simDateIdxNextBigStep + todayIdx == stopIdx){
            // Must save the current spot EQ (or actually log of it here).
            // The irLnMoney array is offset by 1. We need to discount LnE, which is diffused up to the date at simDateIdxNextBigStep.
            double lnSpotEQ = LnE + (*irLnMONEY)[simDateIdxNextBigStep - 1];
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
    }
}

void SRMEquityDiffuseStrictSecondOrder::finalize(DateTimeArrayConstSP allDates)
{

    SRMEquityDiffuseSecondOrder::finalize(allDates);

    size_t numDates = sigmaEQ.size();
    
    vector<int> bigStepBools;
    size_t numBigSteps = setupCriticalBigStepBools(bigStepBools);

    //We may wish to concentrate the time points in regions where the volatility is greatest.
    vector<double> volSize(SpotVol.size());
    for (size_t i = 0;i<volSize.size();i++)
    {
        double a1 = EqSmile_a1[i];
        double a2 = EqSmile_a2[i];
        double a3 = EqSmile_a3[i];
        double sv = SpotVol[i];

        double alpha = (sv * sv * sv) * (sv * sv * sv);
        double beta = sv * sv * a1 * a3;
        beta = beta * beta * beta;
        double gamma = 0.25 * sv * sv * sv * a2 * a3 * a3;
        gamma = gamma * gamma;

        volSize[i] = Maths::max(Maths::max(alpha, beta), gamma);

    }

    vector<int> initialBigStepLocations;
    for (size_t simDateIdx = 0; simDateIdx < numDates+1; ++simDateIdx)
    {
        if (bigStepBools[simDateIdx])
        {
            initialBigStepLocations.push_back(simDateIdx);
        }
    }
    
    double numYears = allDates->front().yearFrac(allDates->back());
    size_t totalGapSize = numDates - numBigSteps;
    if (totalGapSize)
    {
        int toAllocate = Maths::min((int)(timePointsPerYear * numYears - numBigSteps), (int)numDates);
        //toAllocate = addTimePointsRoundDown(totalGapSize, toAllocate, initialBigStepLocations, volSize, bigStepBools);       
        //toAllocate = addTimePointsInOrder(toAllocate, bigStepBools);
        toAllocate = addTimePointsDivideLargestGap(toAllocate, volSize, bigStepBools);
        toAllocate = addTimePointsInOrder(toAllocate, bigStepBools);

        size_t addedPoints = addTimePointsFillLargeGaps(SRMConstants::MAX_DIFFUSION_STEP_SIZE, bigStepBools);
    }

    //Finally we can assemble the list of 'big steps'
    bigStepYearFracs.clear();
    bigStepYearFracsSqrt.clear();
    cfactorBigIntegrals.clear();
    setupBigStepIndicesAndYearFracs(bigStepBools, bigStepIdxs, bigStepYearFracs, bigStepYearFracsSqrt);

    double cfactorTotal = 0.0;
    for (size_t simDateIdx = 0; simDateIdx < numDates+1; ++simDateIdx)
    {
        //cfactor integral at time n should include CfactorArrayEQ[n]
        if (bigStepBools[simDateIdx])
        {
            cfactorBigIntegrals.push_back(cfactorTotal);
            cfactorTotal = 0.0;
        }
        if (simDateIdx < numDates)
        {
            cfactorTotal += CfactorArrayEQ[simDateIdx];
        }        
    }
    
    size_t numBig = bigStepYearFracs.size();
    for (size_t i = 0;i<numBig - 1;i++)
    {
        cfactorBigIntegrals[i] = cfactorBigIntegrals[i+1];
    }
}

double SRMEquityDiffuseStrictSecondOrder::getGapWeight(size_t nBegin, size_t nEnd, const vector<double> &volSize)
{
    double weight = 0.0;
    const double volCoeff = 5.0;
    for (size_t j = nBegin;j<nEnd;j++)
    {
        weight += (1 + volCoeff * volSize[j])* yearFrac[j];
        //weight += volSize[j]* yearFrac[j];
    }
    return weight;
}

size_t SRMEquityDiffuseStrictSecondOrder::addTimePointsRoundDown(size_t totalGapSize, int toAllocate, const vector<int> &initialBigStepLocations, const vector<double> &volSize, vector<int> &bigStepBools)
{
    if (totalGapSize == 0)
        return toAllocate;

    int toAlloc0 = toAllocate;
    //Compute gap sizes between current big steps. Recall 0 is a big step.
    vector<size_t> gapSizes;
    vector<double> gapWeights;
    vector<size_t> gapAllocations;

    double totalWeight = 0.0;
   
    for (size_t i = 1;i<initialBigStepLocations.size();i++)
    {
        size_t gapSize = initialBigStepLocations[i] - initialBigStepLocations[i-1] - 1;
        gapSizes.push_back(gapSize);

        double weight = getGapWeight(initialBigStepLocations[i-1], initialBigStepLocations[i], volSize);
        totalWeight += weight;
        gapWeights.push_back(weight);
    }
    if (gapWeights.size() == 0)
        return toAllocate;

    toAlloc0 = toAllocate;
    for (size_t i = 1;i<initialBigStepLocations.size();i++)
    {
        double weight = gapWeights[i-1];
        size_t gapSize = gapSizes[i-1];

        size_t gapAllocation = 0;
        gapAllocation = (int)((weight * (double)toAllocate) / totalWeight);
        toAlloc0 -= gapAllocation;
        if (toAlloc0 < 0)
        {
            gapAllocation -= (- toAlloc0);
            gapAllocation = Maths::max(gapAllocation, 0);
        }
        gapAllocations.push_back(gapAllocation);
    }

    if (toAlloc0 > 0)
    {
        vector<GapInfo> gis;
        for (size_t i = 0;i<gapSizes.size();i++)
        {
            GapInfo gi;
            gi.gapID = i;
            gi.weight = gapWeights[i];
            gis.push_back(gi);        
        }

        std::sort(gis.begin(), gis.end());
        for (size_t i = 0;i<gis.size() && toAlloc0 > 0;i++)
        {
            const GapInfo &gi = gis[i];
            gapAllocations[gi.gapID]++;
            toAlloc0--;
        }

    }

    size_t leftPos = 0;
    for (size_t i = 0;i<gapSizes.size() && toAllocate>0;i++)
    {
        int gapSize = (int)gapSizes[i];
        int gapAlloc = (int)gapAllocations[i];
        double space = (double)(gapSize - gapAlloc) / (double)(gapAlloc + 1);
        if (space < 1.0)
            space = 1.0;

        for (int j = 0;j<gapAlloc && toAllocate>0;j++)
        {
            size_t pos = (size_t)(leftPos + (space + 1.0) * (double)(j+1));
            if (!bigStepBools[pos])
            {
                bigStepBools[pos] = true;
                toAllocate--;
            }
        }
        leftPos += gapSize + 1;
    }

    return toAllocate;
}

size_t SRMEquityDiffuseStrictSecondOrder::addTimePointsDivideLargestGap(int toAllocate, const vector<double> &volSize, vector<int> &bigStepBools)
{
    while(toAllocate > 0)
    {
        //Find largest gap
        double maxGapWeight = -1;
        int maxGapCentre = -1;
        vector<size_t> gapSizes;
        vector<size_t> gapAllocations;

        int left = 0;
        int right = 0;
        for (size_t i = 0;i<bigStepBools.size();i++)
        {
            if (bigStepBools[i])
            {
                right = i;
                double gapWeight = getGapWeight(left, right, volSize);
                if (gapWeight > maxGapWeight)
                {
                    maxGapCentre = (right + left) / 2;
                    maxGapWeight = gapWeight;
                }
                left = right;            
            }

            /*
            gapWeight += volSize[i];
            if (bigStepBools[i])
            {
                right = i;
                double gapSize = gapWeight * (right - left - 1);
                if (gapSize > maxGapSize)
                {
                    maxGapCentre = (right + left) / 2;
                    maxGapSize = gapSize;
                }
                left = right;                
                gapWeight = volSize[i];
            }
            */
        }

        if (maxGapWeight <= 0)
            break;

        if (bigStepBools[maxGapCentre])
        {
            break;
        }


        bigStepBools[maxGapCentre] = true;
        toAllocate--;
    }

    return toAllocate;
}

size_t SRMEquityDiffuseStrictSecondOrder::addTimePointsInOrder(int toAllocate, vector<int> &bigStepBools)
{

    //only a few remainders left...
    for (size_t i = 0;i<bigStepBools.size() && toAllocate != 0;i++)
    {
        if (!bigStepBools[i])
        {
            bigStepBools[i] = true;
            toAllocate--;
        }
    }
    return toAllocate;
}

size_t SRMEquityDiffuseStrictSecondOrder::addTimePointsFillLargeGaps(double maxGapYearFrac, vector<int> &bigStepBools)
{
    size_t numAdded = 0;
    double gapPeriod = 0.0;
   
    int left = 0;
    int right = 0;
    for (size_t i = 0;i<bigStepBools.size();i++)
    {
        if (bigStepBools[i])
        {
            right = i;
            if (gapPeriod > maxGapYearFrac)
            {
                int gapCentre = (right + left) / 2;
                bigStepBools[gapCentre] = true;
            }      
            left = right;          
            gapPeriod = 0.0;
        }
        gapPeriod += yearFrac[i];
    }
    return numAdded;
}

DRLIB_END_NAMESPACE  


