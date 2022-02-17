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

#include "q3.h"

 
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
        Q3ErrMsg("%s: Failed\n", routine);
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
    smooth; prm;
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
 * Q3Pay1D_Pow
 *
 * Pays exponent or rate.
 *
 */
int Q3Pay1D_Pow(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield = pt[0];
    double expon = (prm->params)[0];
	smooth       = smooth;

    *payoff = pow(yield, expon);
    
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3PayBin1D
 *
 * Binary payoff.
 *
 */
int Q3Pay1D_Bin(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield = pt[0];

    *payoff = Q3SmoothStepFcn((yield - prm->strike) * prm->cop, smooth);    

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3PayAnn1D
 *
 * Annuity payoff.
 */
int Q3Pay1D_Ann(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double strike = prm->strike;
    double cop    = prm->cop;

    *payoff = Q3SmoothMAX((yield - strike) * cop, smooth);
    if (fabs((prm->params)[0]) < TINY) 
    {
        *payoff *= (1. - pow(1. + yield * (prm->params)[1], -1)) / yield ;
    } 
    else 
    {
        *payoff *= (1. - pow(1. + yield / (prm->params)[0], 
                             -(prm->params)[0] * (prm->params)[1])) / yield;
    }
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_Tec
 */
int Q3Pay1D_Tec(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double strike = prm->strike;
    double cop    = prm->cop;
    double base;

    base = MAX(1. + (yield + (prm->params)[1]) * (prm->params)[2], 0.);
    *payoff = Q3SmoothMAX(
        (pow(base, (prm->params)[0]) - 
         pow(1. + strike, (prm->params)[0])) * cop, smooth);

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_TecFwd
 */
int Q3Pay1D_TecFwd(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double base;

    smooth;

    base = MAX(1. + (yield + (prm->params)[1]) * (prm->params)[2], 0.);
    *payoff = pow(base, (prm->params)[0]) - 1.;

    return SUCCESS;
} /* Q3Pay1D_TecFwd */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FixedRibIn
 */
int Q3Pay1D_FixedRibIn(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double lbar   = prm->params[0]; /* low  barrier */
    double hbar   = prm->params[1]; /* high barrier */
    double leps   = prm->params[2]; /* low  epsilon */
    double heps   = prm->params[3]; /* high epsilon */

    *payoff = 0.;

    /* compute inside RIB payoff ___/   \___ */
    /* low barrier: distinguish between zero and non-zero epsilon */
    if (fabs(leps) < TINY) 
    {
        *payoff += Q3SmoothStepFcn(yield - lbar, smooth);
    } 
    else 
    {
        *payoff += (Q3SmoothMAX(yield - lbar + leps, smooth) -
                    Q3SmoothMAX(yield - lbar, smooth)) / leps;
    }
         
    /* high barrier: distinguish between zero and non-zero epsilon */
    if (fabs(heps) < TINY) 
    {
        *payoff -= Q3SmoothStepFcn(yield - hbar, smooth);
    } 
    else 
    {
        *payoff -= (Q3SmoothMAX(yield - hbar, smooth) -
                    Q3SmoothMAX(yield - hbar - heps, smooth)) / heps;
    }

    return SUCCESS;

} /* Q3Pay1D_FixedRibIn */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FixedRibOut
 */
int Q3Pay1D_FixedRibOut(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield  = pt[0];
    double lbar   = prm->params[0]; /* low  barrier */
    double hbar   = prm->params[1]; /* high barrier */
    double leps   = prm->params[2]; /* low  epsilon */
    double heps   = prm->params[3]; /* high epsilon */

    *payoff = 0.;

    /* compute inside RIB payoff ___/   \___ */
    /* low barrier: distinguish between zero and non-zero epsilon */
    if (fabs(leps) < TINY) 
    {
        *payoff += Q3SmoothStepFcn(yield - lbar, smooth);
    } 
    else 
    {
        *payoff += (Q3SmoothMAX(yield - lbar + leps, smooth) -
                    Q3SmoothMAX(yield - lbar, smooth)) / leps;
    }
         
    /* high barrier: distinguish between zero and non-zero epsilon */
    if (fabs(heps) < TINY) 
    {
        *payoff -= Q3SmoothStepFcn(yield - hbar, smooth);
    } 
    else 
    {
        *payoff -= (Q3SmoothMAX(yield - hbar, smooth) -
                    Q3SmoothMAX(yield - hbar - heps, smooth)) / heps;
    }

    *payoff = 1. - *payoff;

    return SUCCESS;

} /* Q3Pay1D_FixedRibOut */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FlrFlrOrCapCapEmbedFlt
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is FixedRate + Call/Put(strike1d).
 * expects params = (spread, strike, leverage)
 */
int Q3Pay1D_FlrFlrOrCapCapEmbedFlt(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO((prm->params[2]));
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double fixedRate, strikeCnd, priceCnd;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strike for 1d problem */
    strikeCnd  = (cop * (yInd - strikeInd) < 0 ? strikeInd : yInd);
    
    /* fixed rate for 1d problem */
    fixedRate = strikeCnd;

    /* adjust strike for leverage and spread */
    strikeCnd -= spread;

    /* price conditional option */
    if (Q3MQLevPricer(
        mqCnd,
        Q3_COP(cop),
        strikeCnd,       
        leverage,
        &priceCnd) == FAILURE) goto RETURN;

    /* option price adjusted for leverage */
    *payoff = fixedRate + cop * priceCnd;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FlrFlrOrCapCap
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is FixedRate + Call/Put(strike1d).
 * expects params = (spread, strike, leverage)
 */
int Q3Pay1D_FlrFlrOrCapCap(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];    
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double fixedRate, strikeCnd, priceCnd;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strike for 1d problem */
    strikeCnd  = (cop * (yInd - strikeInd) < 0 ? strikeInd : yInd);
    
    /* fixed rate for 1d problem */
    fixedRate = (cop * (yInd - strikeInd) < 0 ? cop * (strikeInd - yInd) : 0.);

    /* adjust strike for leverage and spread */
    strikeCnd -= spread;
        
    if (Q3MQLevPricer(
        mqCnd,
        Q3_COP(cop),
        strikeCnd,
        leverage,
        &priceCnd) == FAILURE) goto RETURN;

    /* option price adjusted for leverage */
    *payoff = fixedRate + priceCnd;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FlrCapOrCapFlrEmbedFlt
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is fixedRate + Put/Call Spread.
 * expects params = (spread, strike, leverage)
 */
int Q3Pay1D_FlrCapOrCapFlrEmbedFlt(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    /* independent variables */
    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double fixedRate, priceCnd1, priceCnd2;
    double strikeCnd1, strikeCnd2;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strikes for 1d problem */
    strikeCnd1 = yInd      - spread;
    strikeCnd2 = strikeInd - spread;
    
    if (cop * (yInd - strikeInd) < 0)
    {
        if (Q3MQLevPricer(
            mqCnd,
            Q3_COP(cop),
            strikeCnd1,       
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (Q3MQLevPricer(
            mqCnd,
            Q3_COP(cop),
            strikeCnd2,
            leverage,
            &priceCnd2) == FAILURE) goto RETURN;
        
        fixedRate = yInd;
    }
    else
    {
        priceCnd1  = 0.;
        priceCnd2  = 0.;
        fixedRate = strikeInd;
    }

    /* option prices adjusted for leverage */
    *payoff = fixedRate + cop * (priceCnd1 - priceCnd2);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_FlrCapOrCapFlr
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is fixedRate + Put/Call Spread.
 * expects params = (spread, strike, leverage)
 */
int Q3Pay1D_FlrCapOrCapFlr(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = (prm->params)[0];
    double strikeInd = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double fixedRate, priceCnd1, priceCnd2;
    double strikeCnd1, strikeCnd2;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strikes for 1d problem */
    strikeCnd1 = yInd      - spread;
    strikeCnd2 = strikeInd - spread;

    if (cop * (yInd - strikeInd) < 0)
    {
        if (Q3MQLevPricer(
            mqCnd,
            Q3_COP(cop),
            strikeCnd1,
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (Q3MQLevPricer(
            mqCnd,
            Q3_COP(cop),
            strikeCnd2,
            leverage,
            &priceCnd2) == FAILURE) goto RETURN;
        
        fixedRate = 0.;
    }
    else
    {
        priceCnd1  = 0.;
        priceCnd2  = 0.;
        fixedRate = cop * (strikeInd - yInd);
    }

    /* option prices adjusted for leverage */
    *payoff = fixedRate + priceCnd1 - priceCnd2;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_Sum
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is Max(CoP*(L*R1 + R2 - strikeInd), 0).
 * Expects ordering: params = (strike, leverage).
 */
int Q3Pay1D_Sum(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double leverage  = SHIFT_ZERO(prm->params[1]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double strikeCnd, priceCnd;

    smooth;    

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strike for 1d problem */
    strikeCnd = strikeInd - yInd;

    /* price option on conditioned variable */
    if (Q3MQLevPricer(
        mqCnd,
        Q3_COP(cop),
        strikeCnd,
        leverage,
        &priceCnd) == FAILURE) goto RETURN;

    /* option price adjusted for leverage */
    *payoff = priceCnd;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_Prod
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Expects ordering: params = (strike).
 */
int Q3Pay1D_Prod(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);


    double muMQSav, sigMQSav;
    double strikeCnd, priceCnd, leverage;

    smooth;    

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* leverage is conditioning rate */
    leverage = SHIFT_ZERO(yInd); 

    /* strike for 1d problem */
    strikeCnd = strikeInd;

    /* price option on conditioned variable */
    if (Q3MQLevPricer(
        mqCnd,
        Q3_COP(cop),
        strikeCnd,       
        leverage,
        &priceCnd) == FAILURE) goto RETURN;

    /* option price adjusted for leverage */
    *payoff = priceCnd;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_Perc
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Expects ordering: params = (strike).
 */
int Q3Pay1D_Perc(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double strikeCnd, priceCnd, leverage;

    smooth;    

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* leverage */
    leverage = yInd - strikeInd;

    /* conditional strike */
    strikeCnd = 0.;

    /* price option on conditioned variable */
    if (cop * leverage >= 0.)
    {
        if (Q3MQPricer(
            mqCnd,
            Q3_CALL,
            strikeCnd,       
            &priceCnd) == FAILURE) goto RETURN;
    }
    else
    {
        if (Q3MQPricer(
            mqCnd,
            Q3_PUT,
            strikeCnd,       
            &priceCnd) == FAILURE) goto RETURN;
    }

    /* option price adjusted for leverage */
    *payoff = fabs(leverage) * priceCnd;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_YldYld
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is yield2 * Integral(yield1).
 */
int Q3Pay1D_YldYld(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double call, put;

    smooth;    

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* price option on conditioned variable */
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        0.,       
        &call) == FAILURE) goto RETURN;

    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        0.,       
        &put) == FAILURE) goto RETURN;

    /* option price adjusted for leverage */
    *payoff = yInd * (call - put);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

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


/*-----------------------------------------------------------------------------
 * Q3Pay1D_YldNull
 *
 * Test function for bivariate product.  Prices yield on first variable,
 * ignoring the second.
 *
 */
int Q3Pay1D_YldNull(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yInd  = pt[0];
    smooth; prm;

    *payoff = yInd;

    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * 
 * Pays MIN(MAX(leverage*yCnd+spread,F),C) * W(yInd,LB,HB,LEPS,HEPS) with 
 * integration over yCnd.
 * W() is a weighting function smoothed by epsilon parameters LEPS, HEPS.
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 * 
 */
int Q3Pay1D_MinMaxIn(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_ZERO((prm->params)[6]);
    double heps      = SHIFT_ZERO((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double weight    = 0.;

    double muMQSav, sigMQSav;
    double fixedRate, priceCnd1, priceCnd2;
    double strikeCnd1, strikeCnd2;

    /* probability inside barriers */
    weight = Q3ProbInBarriers(
        yInd,
        lb,
        hb,
        leps,
        heps,
        smooth);

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strikes for 1d problem */
    strikeCnd1 = floor - spread;
    strikeCnd2 = cap   - spread;

    if (floor < cap)
    {
        if (Q3MQLevPricer(
            mqCnd,
            Q3_CALL,
            strikeCnd1,       
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (Q3MQLevPricer(
            mqCnd,
            Q3_CALL,
            strikeCnd2,
            leverage,
            &priceCnd2) == FAILURE) goto RETURN;
        
        fixedRate = floor;
    }
    else
    {
        priceCnd1  = 0.;
        priceCnd2  = 0.;
        fixedRate  = cap;
    }

    *payoff  = fixedRate + priceCnd1 - priceCnd2;
    *payoff *= weight;
    
    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxIn */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yCnd+spread,F),C) * W(yInd,LB,HB,LEPS,HEPS) with 
 * integration over yCnd.
 * W() is a weighting function smoothed by epsilon parameters LEPS, HEPS.
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 */
int Q3Pay1D_MinMaxOut(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_ZERO((prm->params)[6]);
    double heps      = SHIFT_ZERO((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double weight    = 0.;

    double muMQSav, sigMQSav;
    double fixedRate, priceCnd1, priceCnd2;
    double strikeCnd1, strikeCnd2;

    /* probability inside barriers */
    weight = Q3ProbOutBarriers(
        yInd,
        lb,
        hb,
        leps,
        heps,
        smooth);

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* strikes for 1d problem */
    strikeCnd1 = floor - spread;
    strikeCnd2 = cap   - spread;

    if (floor < cap)
    {
        if (Q3MQLevPricer(
            mqCnd,
            Q3_CALL,
            strikeCnd1,       
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (Q3MQLevPricer(
            mqCnd,
            Q3_CALL,
            strikeCnd2,       
            leverage,
            &priceCnd2) == FAILURE) goto RETURN;
        
        fixedRate = floor;
    }
    else
    {
        priceCnd1  = 0.;
        priceCnd2  = 0.;
        fixedRate  = cap;
    }

    /* option prices adjusted for leverage */
    *payoff  = fixedRate + priceCnd1 - priceCnd2;
    *payoff *= weight;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxOut */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxSpdIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<(yCnd-yInd)<HB, else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int Q3Pay1D_MinMaxSpdIn(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double loBnd, hiBnd;
    double price, priceBS;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    hiBnd = hb + yInd;
    loBnd = lb + yInd;

    /* calculate binary put spread */
    if (Q3MQBinPricer(
        mqCnd,
        Q3_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBS = price;
    
    if (Q3MQBinPricer(
        mqCnd,
        Q3_PUT,
        loBnd,
        &price) == FAILURE) goto RETURN;
    priceBS -= price;

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * priceBS;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxSpdIn */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxSpdOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if (yCnd-yInd)<LB || HB<(yCnd-yInd), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int Q3Pay1D_MinMaxSpdOut(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double loBnd, hiBnd;
    double price, priceBS;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    hiBnd = hb + yInd;
    loBnd = lb + yInd;

    /* calculate binary put spread */
    if (Q3MQBinPricer(
        mqCnd,
        Q3_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBS = price;
    
    if (Q3MQBinPricer(
        mqCnd,
        Q3_PUT,
        loBnd,
        &price) == FAILURE) goto RETURN;
    priceBS -= price;

    *payoff = MIN(MAX(leverage*yInd + spread, floor), cap) * (1. - priceBS);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxSpdOut */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxSpdIn_EPS
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<(yCnd-yInd)<HB, else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 */
int Q3Pay1D_MinMaxSpdIn_EPS(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_ZERO((prm->params)[6]);
    double heps      = SHIFT_ZERO((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav   = 0.;
    double sigMQSav  = 0.;

    double loBnd, hiBnd;
    double price, priceBPS, priceBCS;

    smooth;

    /* validate inputs */
    if ((fabs(leps) < TINY) || (fabs(heps) < TINY)) goto RETURN;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    hiBnd = hb + yInd;
    loBnd = lb + yInd;

    /* calculate binary call spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        loBnd - leps,
        &price) == FAILURE) goto RETURN;
    priceBCS = price;
    
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        loBnd,
        &price) == FAILURE) goto RETURN;
    priceBCS -= price;
    priceBCS /= leps;

    /* calculate binary put spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hiBnd + heps,
        &price) == FAILURE) goto RETURN;
    priceBPS = price;
    
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBPS -= price;
    priceBPS /= heps;


    *payoff = MIN(MAX(leverage*yInd + spread, floor), cap) 
        * (priceBCS + priceBPS - 1.);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    if (status == SUCCESS)
    {
        mqCnd->muMQ  = muMQSav;
        mqCnd->sigMQ = sigMQSav;
    }

    return status;

} /* Q3Pay1D_MinMaxSpdIn_EPS */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxSpdOut_EPS
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if (yCnd-yInd)<LB || HB<(yCnd-yInd), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 */
int Q3Pay1D_MinMaxSpdOut_EPS(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_ZERO((prm->params)[6]);
    double heps      = SHIFT_ZERO((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav   = 0.;
    double sigMQSav  = 0.;
    double loBnd, hiBnd;
    double price, priceBPS, priceBCS;

    smooth;

    /* validate inputs */
    if ((fabs(leps) < TINY) || (fabs(heps) < TINY)) goto RETURN;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    hiBnd = hb + yInd;
    loBnd = lb + yInd;

    /* calculate binary call spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        loBnd - leps,
        &price) == FAILURE) goto RETURN;
    priceBCS = price;
    
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        loBnd,
        &price) == FAILURE) goto RETURN;
    priceBCS -= price;
    priceBCS /= leps;

    /* calculate binary put spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hiBnd + heps,
        &price) == FAILURE) goto RETURN;
    priceBPS = price;
    
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBPS -= price;
    priceBPS /= heps;


    *payoff = MIN(MAX(leverage*yInd + spread, floor), cap) 
        * (2. - priceBCS - priceBPS);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    if (status == SUCCESS)
    {
        mqCnd->muMQ  = muMQSav;
        mqCnd->sigMQ = sigMQSav;
    }

    return status;

} /* Q3Pay1D_MinMaxSpdOut_EPS */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxIn_PayInd
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (lB,HB).
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<yCnd<HB (w/ epsilons), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, loEps, hiEps)
 */
int Q3Pay1D_MinMaxIn_PayInd(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_EPS((prm->params)[6]);
    double heps      = SHIFT_EPS((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, p1, p2;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate put spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hb + heps,
        &p1) == FAILURE) goto RETURN;

    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hb,
        &p2) == FAILURE) goto RETURN;

    putSpd = p1 - p2;
    
    /* compute call spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        lb - leps,
        &p1) == FAILURE) goto RETURN;

    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        lb,
        &p2) == FAILURE) goto RETURN;

    callSpd = p1 - p2;

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * 
	(callSpd/leps + putSpd/heps - 1.);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxIn_PayInd */


/*-----------------------------------------------------------------------------
 * Q3Pay1D_MinMaxOut_PayInd
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (lB,HB).
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if yCnd<LB||yCnd>HB (w/ epsilons), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, loEps, hiEps)
 */
int Q3Pay1D_MinMaxOut_PayInd(
    PAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    MQDATA *mqInd    = (prm->mq)[0];
    MQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = SHIFT_EPS((prm->params)[6]);
    double heps      = SHIFT_EPS((prm->params)[7]);
    double corr      = COLLAR(Q3_MIN_CORR, Q3_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, p1, p2;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate put spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hb + heps,
        &p1) == FAILURE) goto RETURN;

    if (Q3MQPricer(
        mqCnd,
        Q3_PUT,
        hb,
        &p2) == FAILURE) goto RETURN;

    putSpd = p1 - p2;
    
    /* compute call spread */
    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        lb - leps,
        &p1) == FAILURE) goto RETURN;

    if (Q3MQPricer(
        mqCnd,
        Q3_CALL,
        lb,
        &p2) == FAILURE) goto RETURN;

    callSpd = p1 - p2;

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * 
	(2. - callSpd/leps - putSpd/heps);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* Q3Pay1D_MinMaxOut_PayInd */
