/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File:        grid1d.c        
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#include <math.h>
#include <stdio.h>

#include "crmultiq.h"

 
/*-----------------------------------------------------------------------------
 * Q3SimpsonPricer1D
 *
 * Integrate payoff times density using Simpson's rule.
 * 
 */
int Q3SimpsonPricer1D(
    PAYOFF  *optPayoff, /* (I) option payoff parameters */
    FPAYOFF *payFunc,   /* (I) option payoff function   */
    double  *grid,      /* (I) array of gaussian points */
    double  *yield,     /* (I) array of yield values    */
    double  *dens,      /* (I) array of density values  */
    double   densNorm,  /* (I) density normalization    */
    double  *premium    /* (O) fwd premium/rate         */
    )
{
    static char routine[] = "Q3SimpsonPricer1D";
    int         status    = FAILURE; 

    double paydens[Q3_1D_INT_PTS];
    double currPt[2];
    double prevStep, nextStep;
    double payoff;    
    long   i;
    
    /* fill (payoff * density) */
    prevStep = 0.;
    nextStep = 0.;
    for (i = 0; i < Q3_1D_INT_PTS; i++) 
    {
        prevStep = nextStep;
        if(i < Q3_1D_INT_PTS - 1) 
            nextStep = yield[i+1] - yield[i];

        currPt[0] = yield[i];
        currPt[1] = grid[i];

        if ((*payFunc)(
            optPayoff,
            currPt,
            MAX(prevStep, nextStep) * Q3_SMOOTH_FACTOR,
            &payoff) == FAILURE) goto RETURN;

        paydens[i] = dens[i] * payoff;

    }

    /* integrate */
    *premium = Q3SimpsIntegral(
        Q3_1D_STEP,
        Q3_1D_INT_PTS,
        paydens) / densNorm;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* Q3SimpsonPricer1D */


/*-----------------------------------------------------------------------------
 * Q3PayYield1D
 *
 */
int Q3Pay1D_Yield(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    /* avoid warnings about unused inputs */
    smooth=smooth; 
    prm=prm;
    *payoff = pt[0];
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3PayVnl1D
 *
 * Pays vanilla.
 *
 */
int Q3Pay1D_Vnl(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield = pt[0];

    *payoff = Q3SmoothMAX((yield - prm->strike) * prm->cop, smooth);
    
    return SUCCESS;
}



/*-----------------------------------------------------------------------------
 * Q3Pay1D_VnlNull
 *
 * Prices Vanilla on first variable, ignoring the second.
 *
 * Expects: params = (strike).
 */
int Q3Pay1D_VnlNull(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yInd      = pt[0];
    double strikeInd = prm->params[0];

    *payoff = Q3SmoothMAX((yInd - strikeInd) * prm->cop, smooth);
    
    return SUCCESS;
}
