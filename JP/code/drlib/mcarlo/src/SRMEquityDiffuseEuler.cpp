#include "edginc/config.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/SRMEquityDiffuseEuler.hpp"
#include "edginc/IQMCRNGManager.hpp"

DRLIB_BEGIN_NAMESPACE  

void SRMEquityDiffuseEuler::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    static const string method("SRMEquityDiffuseEuler::generatePath");

    int numDates = this->sigmaEQ.size();

    randoms = rngMgr->getCorrelatedRandoms(randomIndex); // save pointer to randoms
    LnE = origLnE; // reset
    spotEqPos = spotEqStart; // position in spotEQ/spotEQIndexes
    expSpotEqPos = 0;
    stopIdx = Maths::min(spotEqIndexes[spotEqPos],
        expSpotEqIndexes[expSpotEqPos]);
    divPos = 0; // position in exDivIndexes


    // loop over all dates
    for (int simDateIdx = 0; simDateIdx < numDates; /* increment in body */){

        // update LnE and go ex- any continuous amounts

        //sigmaEQ is the local vol, including the spotVol factor
        double sigmaEQ = this->sigmaEQ[simDateIdx];

        double sqyf = sqrtYearFrac[simDateIdx];
        double yf = yearFrac[simDateIdx];
        double x_random = randoms[simDateIdx];


        //We diffuse LnE. LnE is log spot equity minus the log of the money market account e^{\int r_t}.
        //The diffusion for LnE does not therefore include the r_t term.

        //We compute LnE at simDateIdx+1 in terms of parameter values at simDateIdx.
        LnE += sigmaEQ * (sqrtYearFrac[simDateIdx] * randoms[simDateIdx] -
            0.5 * sigmaEQ * yearFrac[simDateIdx]);
        LnE += CfactorArrayEQ[simDateIdx];

        /* CUPS adjust if necessary
        actually could use sigmaFX as a flag here */
        if (ccyTreatment != ccyVanilla) {
            LnE -= corrEqFX * (*sigmaFX)[simDateIdx] * sigmaEQ * yearFrac[simDateIdx];
        }

        // store LnE if needed
        simDateIdx++; // makes indexing easier for dates

        // Some things about dividends - note this is AFTER the "move to next step" (simDateIdx has been ++d)
        // so that a div on this date is treated as ex- and not cum-
        // eqdiffuse.c lines 445 - 492
        if (simDateIdx == exDivIndexes[divPos]) {
            // no dollar divs as such, no pv around div pmt date vs div ex date
            // so the formula is a bit simpler and no need to go via actual fwd eq
            LnE += lnDivYield[divPos];
            divPos++;
        }

        if (simDateIdx + todayIdx == stopIdx){
            // must save the current spot EQ (or actually log of it here)

            //As remarked above, we need to add on \int{r_t}.
            //We have already incremented simDateIdx, and *irLnMONEY[n]
            //is the integral from 0 to n+1 of r_t (note irLnMoney[0] != 0 in general).
            //Hence r_t is integrated to the same date as LnE,
            //the current value of simDateIdx.
            double lnSpotEQ = LnE + (*irLnMONEY)[simDateIdx-1];
            if (ccyTreatment == ccyStruck) {
                lnSpotEQ += (*spotFXFullPath)[simDateIdx];
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
        // compute new value for sigmaEQ
        if (simDateIdx < (int) SpotVol.size()){
            double A_Smile = EqSmile_a1[simDateIdx]; /* skew  */
            double B_Smile = EqSmile_a2[simDateIdx]; /* smile */
            double spotvol = SpotVol[simDateIdx];
            if (A_Smile == 0 && B_Smile == 0){ /* no smile, can save time */
                sigmaEQ = spotvol;
            } else {

                //See the remark above as to why we are adding (*irLnMONEY)[simDateIdx-1].
                //Similar reasoning leads to the addition of LnFwdEq[simDateIdx-1]:
                //LnFwdEq[n] is the log fwd at time n+1 (LnFwdEq[0] != 0 in general).

                double x = LnE + (*irLnMONEY)[simDateIdx-1] - LnFwdEq[simDateIdx-1];
                double C_Smile = EqSmile_a3[simDateIdx];
                double w = C_Smile * x;
                const double SMILE_CUTOFF = 100.0;
                if (w < SMILE_CUTOFF && w > - SMILE_CUTOFF) {
                    double y = exp(w);
                    double z = 1.0/y;
                    /* loc vol = spotvol * (1 + A * Tanh(C * x) +
                    *                   B * ( 1 - 1/Cosh(C * x)) ) */
                    sigmaEQ = spotvol * (1.0 + A_Smile * (y-z)/(y+z)
                        + B_Smile * (1.0 - 2.0/(y+z)));
                } else if (w <= - SMILE_CUTOFF) {
                    /* y = 0 */
                    sigmaEQ = spotvol * (1.0 - A_Smile + B_Smile);
                } else {
                    /* z = 0 */
                    sigmaEQ = spotvol * (1.0 + A_Smile + B_Smile);
                }
            }
            this->sigmaEQ[simDateIdx] = Maths::max(sigmaEQ, 0.01 * spotvol); // save it
        }
    }
    QLIB_VERIFY((size_t)spotEqPos == spotEq.size(), Format::toString("spotEqPos ended with the wrong value %l, %l expected", (size_t)spotEqPos, spotEq.size()));
    QLIB_VERIFY((size_t)expSpotEqPos == expSpotEq.size(), Format::toString("expSpotEqPos ended with the wrong value %l, %l expected", (size_t)expSpotEqPos, expSpotEq.size()));
}


DRLIB_END_NAMESPACE  

