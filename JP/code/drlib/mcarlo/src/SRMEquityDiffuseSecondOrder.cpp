#include "edginc/config.hpp"
#include "edginc/SRMEquityDiffuse.hpp"
#include "edginc/SRMEQVol.hpp"
#include "edginc/SRMEquityUtil.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/Maths.hpp"
#include "edginc/CriticalDateCollector.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/SRMEquityDiffuseSecondOrder.hpp"
#include "edginc/IQMCRNGManager.hpp"

DRLIB_BEGIN_NAMESPACE 



//JDee: I _think_ there is a bug in the maths of the 2nd order approx: BUGFIX_FACTOR should be 0.5
//Need to lookup in Kloeden and Platen.
#define BUGFIX_FACTOR 0.5

/** Simulate a new path = CalcStateQty_EQ + DiffuseEQ */
void SRMEquityDiffuseSecondOrder::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
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
    for (int simDateIdx = 0; simDateIdx < numDates; ++simDateIdx/* DO NOT increment in body ;)*/){

        double delta_t = yearFrac[simDateIdx];
        double sqrt_delta_t = sqrtYearFrac[simDateIdx];
        double x_random = randoms[simDateIdx];

        // compute new value for vol voldash, voldashdash
        if (simDateIdx > 0){
            double A_Smile = EqSmile_a1[simDateIdx]; /* skew  */
            double B_Smile = EqSmile_a2[simDateIdx]; /* smile */
            double spotvol = SpotVol[simDateIdx];
            if (A_Smile == 0 && B_Smile == 0){ /* no smile, can save time */
                sigmaEQ = spotvol;
                volDashEQ = 0.0;
                volDashDashEQ = 0.0;
            } else {
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
        double mu_2 = CfactorArrayEQ[simDateIdx];
        double mu_dash = - sigmaEQ * volDashEQ; //first derivative of 'drift' with respect to VOL
        double mu_dashdash = - volDashEQ * volDashEQ - sigmaEQ * volDashDashEQ;

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
        double first_order = sigmaEQ * x_random * sqrt_delta_t + mu_dt;

        //TODO: only compute if smile is nonzero
        double second_order = 0;
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
        second_order += sigmaEQ * volDashEQ * delta_t * 0.5 * (x_random * x_random - 1 - mu_dt);
        second_order += mu_dt * volDashEQ * x_random * sqrt_delta_t * 0.5       //mu already has a delta_t factor
            + sigmaEQSquare * (volDashDashEQ * 0.5
            - volDashEQ) * x_random * sqrt_delta_t * delta_t * 0.5;
        second_order -= (volDashEQ * volDashEQ + sigmaEQ * volDashDashEQ)  * sigmaEQSquare * delta_t * delta_t * BUGFIX_FACTOR * 0.5;

        LnE += first_order;
        LnE += second_order;

        /* End of main formula for 2nd order discretisation */

        // Store LnE if needed
        if (simDateIdx + 1 + todayIdx == stopIdx){
            // must save the current spot EQ (or actually log of it here)
            double lnSpotEQ = LnE + (*irLnMONEY)[simDateIdx];
            if (ccyTreatment == ccyStruck) {
                lnSpotEQ += (*spotFXFullPath)[simDateIdx+1];
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

void SRMEquityDiffuseSecondOrder::finalize(DateTimeArrayConstSP allDates)
{
    SRMEquityDiffuse::finalize(allDates);

    //Perform any extra initialization here...
 }

DRLIB_END_NAMESPACE  


