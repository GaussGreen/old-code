/******************************************************************************
 * Module:      CRXQ
 * File:        crxqgrid1d.c        
 * Author:      Credit QRD Charles Morcom
 *****************************************************************************/
#include <math.h>
#include <stdio.h>

#include "crxq.h"
#include "q3.h"
#include "crxerror.h"

 
/*-----------------------------------------------------------------------------
 * CRXQSimpsonPricer1D
 *
 * Integrate payoff times density using Simpson's rule.
 * 
 */
int CRXQSimpsonPricer1D(
    CRXQPAYOFF  *optPayoff, /* (I) option payoff parameters */
    CRXQFPAYOFF *payFunc,   /* (I) option payoff function   */
    double  *grid,      /* (I) array of gaussian points */
    double  *yield,     /* (I) array of yield values    */
    double  *dens,      /* (I) array of density values  */
    double   densNorm,  /* (I) density normalization    */
    double  *premium    /* (O) fwd premium/rate         */
    )
{
    static char routine[] = "CRXQSimpsonPricer1D";
    int         status    = FAILURE; 

    double paydens[CRXQ_1D_INT_PTS];
    double currPt[2];
    double prevStep, nextStep;
    double payoff;    
    long   i;
    
    /* fill (payoff * density) */
    prevStep = 0.;
    nextStep = 0.;
    for (i = 0; i < CRXQ_1D_INT_PTS; i++) 
    {
        prevStep = nextStep;
        if(i < CRXQ_1D_INT_PTS - 1) 
            nextStep = yield[i+1] - yield[i];

        currPt[0] = yield[i];
        currPt[1] = grid[i];

        if ((*payFunc)(
            optPayoff,
            currPt,
            MAX(prevStep, nextStep) * CRXQ_SMOOTH_FACTOR,
            &payoff) == FAILURE) goto RETURN;

        paydens[i] = dens[i] * payoff;

    }
    /* integrate (from q3.h) and adjust for distribution normalization */
    *premium = Q3SimpsIntegral(
        CRXQ_1D_STEP,
        CRXQ_1D_INT_PTS,
        paydens) / densNorm;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQSimpsonPricer1D */

/******************************************************************************
 * FROM HERE TO THE END OF THE FILE FOLLOW ALL THE DIFFERENT PAYOFF FUNCTIONS
 * FOR ALL THE VARIOUS OPTION TYPES
 *****************************************************************************/

/*-----------------------------------------------------------------------------
 * CRXQPayYield1D
 * Pays just the yield itself
 */
int CRXQPay1D_Yield(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPayVnl1D
 * Pays vanilla call or put
 */
int CRXQPay1D_Vnl(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPay1D_Pow
 * Pays exponent of rate.
 */
int CRXQPay1D_Pow(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield = pt[0];
    double expon = (prm->params)[0];

    *payoff = pow(yield, expon);
    
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPay1D_Quad
 * Pays y*MAX(CoP *(y - K), 0) - call/put payoff multiplied by rate
 */
int CRXQPay1D_Quad(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    double yield = pt[0];

    *payoff = Q3SmoothMAX((yield - prm->strike) * prm->cop, smooth) * yield;
    
    return SUCCESS;
}


/*-----------------------------------------------------------------------------
 * CRXQPayBin1D
 * Binary payoff.
 */
int CRXQPay1D_Bin(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPayAnn1D
 * Annuity payoff. This is the simple swap annuity - would need to be 
 * changed for credit
 */
int CRXQPay1D_Ann(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPay1D_Tec
 */
int CRXQPay1D_Tec(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPay1D_TecFwd
 */
int CRXQPay1D_TecFwd(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
} /* CRXQPay1D_TecFwd */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_FixedRibIn
 */
int CRXQPay1D_FixedRibIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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

} /* CRXQPay1D_FixedRibIn */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_FixedRibOut
 */
int CRXQPay1D_FixedRibOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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

} /* CRXQPay1D_FixedRibOut */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_FlrFlrOrCapCapEmbedFlt
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is FixedRate + Call/Put(strike1d).
 * expects params = (spread, strike, leverage)
 */
int CRXQPay1D_FlrFlrOrCapCapEmbedFlt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO((prm->params[2]));
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQLevPricer(
        mqCnd,
        CRXQ_COP(cop),
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
 * CRXQPay1D_FlrFlrOrCapCap
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is FixedRate + Call/Put(strike1d).
 * expects params = (spread, strike, leverage)
 */
int CRXQPay1D_FlrFlrOrCapCap(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];    
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
        
    if (CRXQLevPricer(
        mqCnd,
        CRXQ_COP(cop),
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
 * CRXQPay1D_FlrCapOrCapFlrEmbedFlt
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is fixedRate + Put/Call Spread.
 * expects params = (spread, strike, leverage)
 */
int CRXQPay1D_FlrCapOrCapFlrEmbedFlt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    /* independent variables */
    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = prm->params[0];
    double strikeInd = prm->params[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
        if (CRXQLevPricer(
            mqCnd,
            CRXQ_COP(cop),
            strikeCnd1,       
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (CRXQLevPricer(
            mqCnd,
            CRXQ_COP(cop),
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
 * CRXQPay1D_FlrCapOrCapFlr
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is fixedRate + Put/Call Spread.
 * expects params = (spread, strike, leverage)
 */
int CRXQPay1D_FlrCapOrCapFlr(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double spread    = (prm->params)[0];
    double strikeInd = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
        if (CRXQLevPricer(
            mqCnd,
            CRXQ_COP(cop),
            strikeCnd1,
            leverage,
            &priceCnd1) == FAILURE) goto RETURN;

        if (CRXQLevPricer(
            mqCnd,
            CRXQ_COP(cop),
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
 * CRXQPay1D_Sum
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is Max(CoP*(L*R1 + R2 - strikeInd), 0).
 * Expects ordering: params = (strike, leverage).
 */
int CRXQPay1D_Sum(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double leverage  = SHIFT_ZERO(prm->params[1]);
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQLevPricer(
        mqCnd,
        CRXQ_COP(cop),
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
 * CRXQPay1D_Prod
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Expects ordering: params = (strike).
 */
int CRXQPay1D_Prod(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQLevPricer(
        mqCnd,
        CRXQ_COP(cop),
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
 * CRXQPay1D_Perc
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Expects ordering: params = (strike).
 */
int CRXQPay1D_Perc(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double strikeInd = prm->params[0];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            strikeCnd,       
            &priceCnd) == FAILURE) goto RETURN;
    }
    else
    {
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
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
 * CRXQPay1D_YldYld
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Payoff is yield2 * Integral(yield1).
 */
int CRXQPay1D_YldYld(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQPricer(
        mqCnd,
        CRXQ_CALL,
        0.,       
        &call) == FAILURE) goto RETURN;

    if (CRXQPricer(
        mqCnd,
        CRXQ_PUT,
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
 * CRXQPay1D_VnlNull
 *
 * Prices Vanilla on first variable, ignoring the second.
 *
 * Expects: params = (strike).
 */
int CRXQPay1D_VnlNull(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPay1D_YldNull
 *
 * Test function for bivariate product.  Prices yield on first variable,
 * ignoring the second.
 *
 */
int CRXQPay1D_YldNull(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
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
 * CRXQPay1D_InSpdRIB
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<(yCnd-yInd)<HB, else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int CRXQPay1D_InSpdRIB(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBS = price;
    
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
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

} /* CRXQPay1D_InSpdRIB */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutSpdRIB
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if (yCnd-yInd)<LB || HB<(yCnd-yInd), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int CRXQPay1D_OutSpdRIB(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        hiBnd,
        &price) == FAILURE) goto RETURN;
    priceBS = price;
    
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
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

} /* CRXQPay1D_OutSpdRIB */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_InSpdRIB_EPS
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<(yCnd-yInd)<HB, else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 */
int CRXQPay1D_InSpdRIB_EPS(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = (prm->params)[6];
    double heps      = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav   = 0.;
    double sigMQSav  = 0.;

    double loBnd, hiBnd;
    double price, priceBPS, priceBCS;

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

    if (fabs(heps) < CRXQ_MIN_EPS)
    { /* calculate binary put */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd,
	    &priceBPS) == FAILURE) goto RETURN;
    }
    else
    { /* calculate put spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd + heps,
	    &price) == FAILURE) goto RETURN;
	priceBPS = price;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd,
	    &price) == FAILURE) goto RETURN;
	priceBPS -= price;
	priceBPS /= heps;
    }

    if (fabs(leps) < CRXQ_MIN_EPS)
    { /* calculate binary call */ 
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd,
	    &priceBCS) == FAILURE) goto RETURN;
    }
    else
    { /* calculate call spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd - leps,
	    &price) == FAILURE) goto RETURN;
	priceBCS = price;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd,
	    &price) == FAILURE) goto RETURN;
	priceBCS -= price;
	priceBCS /= leps;
    }

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

} /* CRXQPay1D_InSpdRIB_EPS */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutSpdRIB_EPS
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if (yCnd-yInd)<LB || HB<(yCnd-yInd), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C, LEPS, HEPS)
 */
int CRXQPay1D_OutSpdRIB_EPS(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = (prm->params)[6];
    double heps      = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav   = 0.;
    double sigMQSav  = 0.;
    double loBnd, hiBnd;
    double price, priceBPS, priceBCS;

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

    if (fabs(heps) < CRXQ_MIN_EPS)
    { /* calculate binary put */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd,
	    &priceBPS) == FAILURE) goto RETURN;
    }
    else
    {
	/* calculate binary put spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd + heps,
	    &price) == FAILURE) goto RETURN;
	priceBPS = price;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hiBnd,
	    &price) == FAILURE) goto RETURN;
	priceBPS -= price;
	priceBPS /= heps;
    }

    if (fabs(leps) < CRXQ_MIN_EPS)
    { /* calculate binary call */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd,
	    &priceBCS) == FAILURE) goto RETURN;
    }
    else
    { /* calculate call spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd - leps,
	    &price) == FAILURE) goto RETURN;
	priceBCS = price;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    loBnd,
	    &price) == FAILURE) goto RETURN;
	priceBCS -= price;
	priceBCS /= leps;
    }

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

} /* CRXQPay1D_MinMaxSpdOut_EPS */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_InRIB_EPS
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
int CRXQPay1D_InRIB_EPS(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = (prm->params)[6];
    double heps      = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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

    if (fabs(heps) <  CRXQ_MIN_EPS)
    { /* calculate binary put */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb,
	    &putSpd) == FAILURE) goto RETURN;
    }
    else
    { /* calculate put spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb + heps,
	    &p1) == FAILURE) goto RETURN;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb,
	    &p2) == FAILURE) goto RETURN;

	putSpd = (p1 - p2) / heps;
    }
    
    if (fabs(leps) < CRXQ_MIN_EPS)
    { /* calculate binary call */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb,
	    &callSpd) == FAILURE) goto RETURN;
    }
    else
    { /* compute call spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb - leps,
	    &p1) == FAILURE) goto RETURN;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb,
	    &p2) == FAILURE) goto RETURN;
	
	callSpd = (p1 - p2) / leps;
    }

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * 
	(callSpd + putSpd - 1.);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InRIB_EPS */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutRIB_EPS
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
int CRXQPay1D_OutRIB_EPS(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double leps      = (prm->params)[6];
    double heps      = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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

    if (fabs(heps) < CRXQ_MIN_EPS)
    { /* calculate binary put */
	if( CRXQBinPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb,
	    &putSpd) == FAILURE) goto RETURN;
    }
    else
    { /* calculate put spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb + heps,
	    &p1) == FAILURE) goto RETURN;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_PUT,
	    hb,
	    &p2) == FAILURE) goto RETURN;
	
	putSpd = (p1 - p2) / heps;
    }

    if (fabs(leps) < CRXQ_MIN_EPS) 
    { /* calculate binary call */
	if (CRXQBinPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb,
	    &callSpd) == FAILURE) goto RETURN;
    }
    else
    { /* compute call spread */
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb - leps,
	    &p1) == FAILURE) goto RETURN;
	
	if (CRXQPricer(
	    mqCnd,
	    CRXQ_CALL,
	    lb,
	    &p2) == FAILURE) goto RETURN;
	
	callSpd = (p1 - p2) / leps;
    }

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * 
	(2. - callSpd - putSpd);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutRIB_EPS */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_InRIB
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (lB,HB).
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if LB<yCnd<HB (no epsilons).
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int CRXQPay1D_InRIB(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double priceBS, price;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        hb,
        &price) == FAILURE) goto RETURN;
    priceBS = price;

    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        lb,
        &price) == FAILURE) goto RETURN;
    priceBS -= price;

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * priceBS;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InRIB */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutRIB
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (lB,HB).
 *
 * Pays MIN(MAX(leverage*yInd+spread,F),C) if yCnd<LB||yCnd>HB (no epsilons), 
 * else pays 0. 
 *
 * expects params = (LB, HB, leverage, spread, F, C)
 */
int CRXQPay1D_OutRIB(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb        = (prm->params)[0];
    double hb        = (prm->params)[1];
    double leverage  = SHIFT_ZERO(prm->params[2]);
    double spread    = (prm->params)[3];
    double floor     = (prm->params)[4];
    double cap       = (prm->params)[5];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double priceBS, price;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        hb,
        &price) == FAILURE) goto RETURN;
    priceBS = price;

    if (CRXQBinPricer(
        mqCnd,
        CRXQ_PUT,
        lb,
        &price) == FAILURE) goto RETURN;
    priceBS -= price;

    *payoff = MIN(MAX(leverage*yInd+spread,floor),cap) * (1. - priceBS);

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutRIB */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_PercWgt
 *
 * Evaluates value of a 1d slice of 2d payoff: MAX(0,R1*(w1*R1+w2*R2-K))
 * 
 * Calculate as MAX(0,Y_ind*(w1*Y_ind+w2*Y_cnd-K))
 *
 * Expects ordering: params = (weight1, weight2, strike).
 */
int CRXQPay1D_PercWgt(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double cop       = prm->cop;
    double weight1   = prm->params[0];
    double weight2   = prm->params[1];
    double strikeInd = prm->params[2];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
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

    /* Handle weight2 = 0 case */
    if (fabs(weight2) < TINY)
    { 
	*payoff = yInd * (weight1 * yInd - strikeInd);
    }
    else
    {
	/* leverage */
	leverage = weight2 * yInd;
	
	/* conditional strike */
	strikeCnd = (strikeInd - weight1 * yInd) / weight2;
	
	/* price option on conditioned variable */
	if (cop * leverage >= 0.)
	{
	    if (CRXQPricer(
		mqCnd,
		CRXQ_CALL,
		strikeCnd,       
		&priceCnd) == FAILURE) goto RETURN;
	}
	else
	{
	    if (CRXQPricer(
		mqCnd,
		CRXQ_PUT,
		strikeCnd,       
		&priceCnd) == FAILURE) goto RETURN;
	}
	
	/* option price adjusted for leverage */
	*payoff = fabs(leverage) * priceCnd;
    }

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

}


/*-----------------------------------------------------------------------------
 * CRXQPay1D_InAndIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 and LB2<yCnd<HB2 (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_InAndIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbInBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    if(weight < TINY)
    {
        *payoff = 0.0;
    }
    else
    {
        /* apply conditioning */
        mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
        mqCnd->sigMQ *= corrCmp;

        /* compute conditional forward rate */
        mqCnd->calcFwd = TRUE;

        /* calculate binary put spread */
        if(fabs(heps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary put */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &putSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate put spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2+heps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &price1) == FAILURE) goto RETURN;
            putSpd = (price2 - price1) / heps2;
        }

        if (fabs(leps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary call */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &callSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate call spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2-leps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &price1) == FAILURE) goto RETURN;
            callSpd = (price2 - price1) / leps2;
        }

        *payoff = weight * (callSpd + putSpd - 1.0);
    }

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InAndIn */

/*-----------------------------------------------------------------------------
 * CRXQPay1D_InOrIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 Or LB2<yCnd<HB2 (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_InOrIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2, priceBS;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbInBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if(fabs(heps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary put */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &putSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate put spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2+heps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &price1) == FAILURE) goto RETURN;
        putSpd = (price2 - price1) / heps2;
    }

    if (fabs(leps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary call */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &callSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate call spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2-leps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &price1) == FAILURE) goto RETURN;
        callSpd = (price2 - price1) / leps2;
    }

    priceBS = callSpd + putSpd - 1.0;
    *payoff = weight + priceBS - weight * priceBS;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InOrIn */

/*-----------------------------------------------------------------------------
 * CRXQPay1D_InAndOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) and LB2<yCnd<HB2 (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_InAndOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbOutBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    if(weight < TINY)
    {
        *payoff = 0.0;
    }
    else
    {
        /* apply conditioning */
        mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
        mqCnd->sigMQ *= corrCmp;

        /* compute conditional forward rate */
        mqCnd->calcFwd = TRUE;

        /* calculate binary put spread */
        if(fabs(heps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary put */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &putSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate put spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2+heps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &price1) == FAILURE) goto RETURN;
            putSpd = (price2 - price1) / heps2;
        }

        if (fabs(leps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary call */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &callSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate call spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2-leps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &price1) == FAILURE) goto RETURN;
            callSpd = (price2 - price1) / leps2;
        }

        *payoff = weight * (callSpd + putSpd - 1.0);
    }

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InAndOut */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_InOrOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) or LB2<yCnd<HB2 (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_InOrOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2, priceBS;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbOutBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if(fabs(heps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary put */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &putSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate put spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2+heps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &price1) == FAILURE) goto RETURN;
        putSpd = (price2 - price1) / heps2;
    }

    if (fabs(leps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary call */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &callSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate call spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2-leps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &price1) == FAILURE) goto RETURN;
        callSpd = (price2 - price1) / leps2;
    }

    priceBS = callSpd + putSpd - 1.0;
    *payoff = weight + priceBS - weight * priceBS;
    
    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_InOrOut */

/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutAndIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 and (yCnd<LB2 or yCnd>HB2) (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_OutAndIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbInBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    if(weight < TINY)
    {
        *payoff = 0.0;
    }
    else
    {
        /* apply conditioning */
        mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
        mqCnd->sigMQ *= corrCmp;

        /* compute conditional forward rate */
        mqCnd->calcFwd = TRUE;

        /* calculate binary put spread */
        if(fabs(heps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary put */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &putSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate put spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2+heps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &price1) == FAILURE) goto RETURN;
            putSpd = (price2 - price1) / heps2;
        }

        if (fabs(leps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary call */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &callSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate call spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2-leps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &price1) == FAILURE) goto RETURN;
            callSpd = (price2 - price1) / leps2;
        }

        *payoff = weight * (2.0 - callSpd - putSpd);
    }

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutAndIn */

/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutOrIn
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if LB1<yInd<HB1 or (yCnd<LB2 or yCnd>HB2) (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_OutOrIn(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2, priceBS;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbInBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if(fabs(heps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary put */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &putSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate put spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2+heps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &price1) == FAILURE) goto RETURN;
        putSpd = (price2 - price1) / heps2;
    }

    if (fabs(leps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary call */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &callSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate call spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2-leps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &price1) == FAILURE) goto RETURN;
        callSpd = (price2 - price1) / leps2;
    }

    priceBS = 2.0 - callSpd - putSpd;
    *payoff = weight + priceBS - weight * priceBS;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutOrIn */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutAndOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) and (yCnd<LB2 or yCnd>HB2) (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_OutAndOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbOutBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    if(weight < TINY)
    {
        *payoff = 0.0;
    }
    else
    {
        /* apply conditioning */
        mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
        mqCnd->sigMQ *= corrCmp;

        /* compute conditional forward rate */
        mqCnd->calcFwd = TRUE;

        /* calculate binary put spread */
        if(fabs(heps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary put */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &putSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate put spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2+heps2,
                &price2) == FAILURE) goto RETURN;
        
            if (CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                hb2,
                &price1) == FAILURE) goto RETURN;
            putSpd = (price2 - price1) / heps2;
        }

        if (fabs(leps2) < CRXQ_MIN_EPS)
        {
            /* calculate binary call */
            if (CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &callSpd) == FAILURE) goto RETURN;
        }
        else
        {
            /* calculate call spread*/
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2-leps2,
                &price2) == FAILURE) goto RETURN;
            
            if (CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                lb2,
                &price1) == FAILURE) goto RETURN;
            callSpd = (price2 - price1) / leps2;
        }

        *payoff = weight * (2.0 - callSpd - putSpd);
    }

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutAndOut */


/*-----------------------------------------------------------------------------
 * CRXQPay1D_OutOrOut
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * Integrates over conditioned variable to find probability that it falls
 * inside range (LB2,HB2).
 *
 * Pays 1 if (yInd<LB1 or yInd>HB1) or (yCnd<LB2 or yCnd>HB2) (w epsilons).
 * else pays 0. 
 *
 * expects params = (LB1, HB1, LB2, HB2, Leps1, Heps1, Leps2, Heps2)
 */
int CRXQPay1D_OutOrOut(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double lb1       = (prm->params)[0];
    double hb1       = (prm->params)[1];
    double lb2       = (prm->params)[2];
    double hb2       = (prm->params)[3];
    double leps1     = (prm->params)[4];
    double heps1     = (prm->params)[5];
    double leps2     = (prm->params)[6];
    double heps2     = (prm->params)[7];
    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double muMQSav, sigMQSav;
    double putSpd, callSpd, price1, price2, priceBS;

    double weight   = 0.;

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    weight = Q3ProbOutBarriersSmooth(
        yInd,
        lb1,
        hb1,
        leps1,
        heps1,
        smooth);

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    /* calculate binary put spread */
    if(fabs(heps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary put */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &putSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate put spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2+heps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_PUT,
            hb2,
            &price1) == FAILURE) goto RETURN;
        putSpd = (price2 - price1) / heps2;
    }

    if (fabs(leps2) < CRXQ_MIN_EPS)
    {
        /* calculate binary call */
        if (CRXQBinPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &callSpd) == FAILURE) goto RETURN;
    }
    else
    {
        /* calculate call spread*/
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2-leps2,
            &price2) == FAILURE) goto RETURN;
        
        if (CRXQPricer(
            mqCnd,
            CRXQ_CALL,
            lb2,
            &price1) == FAILURE) goto RETURN;
        callSpd = (price2 - price1) / leps2;
    }

    priceBS = 2.0 - callSpd - putSpd;

    *payoff = weight + priceBS - weight * priceBS;

    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_OutOrOut */

/*-----------------------------------------------------------------------------
 * CRXQPay1D_ComplexSPD
 *
 * Evaluates value of a 1d slice of a 2d payoff.
 * payoff = max(K0, a00*R1 + a10*R2 + 
 *                  b0 * max(K1, a01*R1 + a11*R2) +
 *                  b1 * max(K2, a02*R1 + a12*R2) +
 *                  b2 * max(K3, a03*R1 + a13*R2) +
 *                  b3 * max(K4, a04*R1 + a14*R2))
 * Expect parameters:
 * (K0, K1, K2, K3, K4, a00, a10, a01, a11, a02, a12, a03, a13,
 *  a04, a14, b0, b1, b2, b3)
 */
int CRXQPay1D_ComplexSPD(
    CRXQPAYOFF *prm,    /* (I) payoff parameters    */
    double *pt,     /* (I) yield,gaussian point */
    double  smooth, /* (I) smoothing window     */
    double *payoff  /* (O) payoff               */
    )
{
    int status = FAILURE;

    double yInd      = pt[0];
    double xInd      = pt[1];
    CRXQDATA *mqInd    = (prm->mq)[0];
    CRXQDATA *mqCnd    = (prm->mq)[1];
    double muMQInd   = mqInd->muMQ;
    double sigMQInd  = SHIFT_ZERO(mqInd->sigMQ);

    double corr      = COLLAR(CRXQ_MIN_CORR, CRXQ_MAX_CORR, prm->corr);
    double corrCmp   = sqrt(1. - corr * corr);

    double K[5], b[4], a[2][5];
    double Strike[5], ZeroFlag[5];
    double muMQSav, sigMQSav;
    double price1, price2;
    int    i,j,idx;
    int    COP[5];
    int    NumOfIntervals;

    typedef struct{
        double Slope;
        double LB;
        double RB;
        double Intercep;
        double Root;
    }PieceLinear;

    PieceLinear LinFunc[5];

    smooth;

    /* save current values */
    muMQSav  = mqCnd->muMQ;
    sigMQSav = mqCnd->sigMQ;

    /* apply conditioning */
    mqCnd->muMQ  += corr * (xInd - muMQInd) * mqCnd->sigMQ / sigMQInd;
    mqCnd->sigMQ *= corrCmp;

    /* compute conditional forward rate */
    mqCnd->calcFwd = TRUE;

    idx = 0;
    for(i=0;i<5;i++)
    {
        K[i] = prm->params[idx++];
    }

    for(i=0;i<5;i++)
    {
        for(j=0;j<2;j++)
        {
            a[j][i] = prm->params[idx++];
        }
    }

    for(i=0;i<4;i++)
    {
        b[i] = prm->params[idx++];
        if(fabs(b[i]) < TINY) /* b[i] == 0*/
        {
            K[i+1] = 0.0;
            a[0][i+1] = 0.0;
            a[1][i+1] = 0.0;
        }
    }
   
    for(i=0;i<5;i++)
    {
        if(fabs(a[0][i]) < TINY)
        {
            ZeroFlag[i] = TRUE;
            COP[i] = CRXQ_CALL;
        }
        else
        {
            ZeroFlag[i] = FALSE; 
            if(a[0][i] > 0)
            {
                COP[i] = CRXQ_CALL; /*Call*/
            }
            else
            {
                COP[i] = CRXQ_PUT; /*Put*/
            }
        }
    }

    Strike[0] = a[1][0] * yInd - K[0];
    for(i=1;i<5;i++)
    {
        if(ZeroFlag[i] == FALSE)
        {
            Strike[i]  = (K[i] - a[1][i] * yInd) / a[0][i];
            Strike[0] += b[i-1] * K[i];
        }
        else
        {
            Strike[0] += b[i-1] * MAX(K[i], a[1][i]*yInd);
            Strike[i]  = -INF;
        }
    }

    /*sort Strikes[1:4] in ascending order*/

    for(i=1;i<5;i++)
    {
        idx = i;;
        for(j=i+1;j<5;j++)
        {
            if (Strike[j] < Strike[idx]) idx = j;
        }
        if( Q3SwapDouble( 
            &(Strike[idx]), 
            &(Strike[i])) == FAILURE) goto RETURN;
        if( Q3SwapInt( 
            &(COP[idx]), 
            &(COP[i])) == FAILURE) goto RETURN;
        if (Q3SwapDouble( 
            &(b[idx-1]), 
            &(b[i-1])) == FAILURE) goto RETURN;
        if (Q3SwapDouble( 
            &(a[0][idx]), 
            &(a[0][i])) == FAILURE) goto RETURN;
        if (Q3SwapDouble( 
            &(a[1][idx]), 
            &(a[1][i])) == FAILURE) goto RETURN;
    }
    
    /*search for first nonzero a[0][i]*/
    idx = 0;
    for(i=1;i<5;i++)
    {
        if(fabs(Strike[i] + INF)< TINY) 
        {
            idx = i;
        }
    }

    NumOfIntervals = 5 - idx;
    
    /*initialize piece-wise linear function*/
    idx++;
    /*First interval [-infty, Strike[1]].*/
    LinFunc[0].LB       = -INF;
    LinFunc[0].RB       = Strike[idx];
    LinFunc[0].Slope    = a[0][0];
    LinFunc[0].Intercep = Strike[0];
    for(i=idx;i<5;i++)
    {
        if (COP[i] == CRXQ_PUT)
        {
            LinFunc[0].Slope    += ( b[i-1] * a[0][i] );
            LinFunc[0].Intercep -= ( b[i-1] * a[0][i] * Strike[i] );
        }
    }
    if(fabs(LinFunc[0].Slope) > TINY)
    {
        LinFunc[0].Root = -LinFunc[0].Intercep / LinFunc[0].Slope;
    }
    else
    {
        LinFunc[0].Root = -INF;
    }


    /*intevals in the middle*/
    for(i=1;i<NumOfIntervals-1; i++)
    {
        LinFunc[i].LB       = Strike[idx];
        LinFunc[i].RB       = Strike[idx+1];
        LinFunc[i].Slope    = LinFunc[i-1].Slope;
        LinFunc[i].Intercep = LinFunc[i-1].Intercep;
        if(COP[idx] == CRXQ_CALL)
        {
            LinFunc[i].Slope    += b[idx-1] * a[0][idx];
            LinFunc[i].Intercep -= b[idx-1] * a[0][idx] * Strike[idx];
        }
        else
        {
            LinFunc[i].Slope    -= b[idx-1] * a[0][idx];
            LinFunc[i].Intercep += b[idx-1] * a[0][idx] * Strike[idx];
        }
        if(fabs(LinFunc[i].Slope)> TINY)
        {
            LinFunc[i].Root = -LinFunc[i].Intercep / LinFunc[i].Slope;
        }
        else
        {
             LinFunc[i].Root = -INF;
        }
        idx++;
    }


    /*last interval [strike, +infty]*/
    i = NumOfIntervals - 1;
    if (i==0) /*One interval only*/
    {
        LinFunc[0].LB = -INF;
        LinFunc[0].RB = INF;
    }
    else
    {
        LinFunc[i].LB       = Strike[idx];
        LinFunc[i].RB       = INF;
        LinFunc[i].Slope    = LinFunc[i-1].Slope;
        LinFunc[i].Intercep = LinFunc[i-1].Intercep;
        if(COP[idx] == CRXQ_CALL)
        {
            LinFunc[i].Slope    += b[idx-1] * a[0][idx];
            LinFunc[i].Intercep -= b[idx-1] * a[0][idx] * Strike[idx];
        }
        else
        {
            LinFunc[i].Slope    -= b[idx-1] * a[0][idx];
            LinFunc[i].Intercep += b[idx-1] * a[0][idx] * Strike[idx];
        }
        if(fabs(LinFunc[i].Slope)> TINY)
        {
            LinFunc[i].Root = -LinFunc[i].Intercep / LinFunc[i].Slope;
        }
        else
        {
            LinFunc[i].Root = -INF;
        }
    }
        
    /*Decompse the payoff function accroding to strike intervals */
    *payoff = 0.0;
    for(i=0;i<NumOfIntervals;i++)
    {
        if ((LinFunc[i].Root > LinFunc[i].LB) &&
            (LinFunc[i].Root < LinFunc[i].RB))
        {
            if(LinFunc[i].Slope > 0) /*Call*/
            {
                if(CRXQPricer(
                    mqCnd,
                    CRXQ_CALL,
                    LinFunc[i].Root,
                    &price1) == FAILURE) goto RETURN;

                if(CRXQPricer(
                    mqCnd,
                    CRXQ_CALL,
                    LinFunc[i].RB,
                    &price2) == FAILURE) goto RETURN;
                *payoff += LinFunc[i].Slope * (price1 - price2);

                if(CRXQBinPricer(
                    mqCnd,
                    CRXQ_CALL,
                    LinFunc[i].RB,
                    &price1) == FAILURE) goto RETURN;
                *payoff -= (LinFunc[i].Slope * LinFunc[i].RB + 
                            LinFunc[i].Intercep) * price1;
            }
            else
            {
                if(CRXQPricer(
                    mqCnd,
                    CRXQ_PUT,
                    LinFunc[i].Root,
                    &price1) == FAILURE) goto RETURN;

                if(CRXQPricer(
                    mqCnd,
                    CRXQ_PUT,
                    LinFunc[i].LB,
                    &price2) == FAILURE) goto RETURN;
                *payoff += -LinFunc[i].Slope * (price1 - price2);

                if(CRXQBinPricer(
                    mqCnd,
                    CRXQ_PUT,
                    LinFunc[i].LB,
                    &price1) == FAILURE) goto RETURN;
                *payoff -= (LinFunc[i].Slope * LinFunc[i].LB + 
                            LinFunc[i].Intercep) * price1;
            }
        }
        else if (LinFunc[i].Root == -INF && 
                 LinFunc[i].Intercep > TINY)
        {
            if(CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].RB,
                &price1) == FAILURE) goto RETURN;

            if(CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].LB,
                &price2) == FAILURE) goto RETURN;
            *payoff += LinFunc[i].Intercep * (price1 - price2);
        }
        else if((LinFunc[i].Slope > 0) &&
                (LinFunc[i].Root <= LinFunc[i].LB))
        {
            if(CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].Root,
                &price1) == FAILURE) goto RETURN;

            if(CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].RB,
                &price2) == FAILURE) goto RETURN;
            *payoff += LinFunc[i].Slope * (price1 - price2);

            if(CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].RB,
                &price1) == FAILURE) goto RETURN;
            *payoff -= (LinFunc[i].Slope * LinFunc[i].RB + 
                        LinFunc[i].Intercep) * price1;
                
            if(CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].Root,
                &price1) == FAILURE) goto RETURN;

            if(CRXQPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].LB,
                &price2) == FAILURE) goto RETURN;
            *payoff -= LinFunc[i].Slope * (price1 - price2);

            if(CRXQBinPricer(
                mqCnd,
                CRXQ_CALL,
                LinFunc[i].LB,
                &price1) == FAILURE) goto RETURN;
            *payoff += (LinFunc[i].Slope * LinFunc[i].LB + 
                        LinFunc[i].Intercep) * price1;
        }
        else if((LinFunc[i].Slope < 0) &&
                (LinFunc[i].Root >= LinFunc[i].RB))
        {
            if(CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].Root,
                &price1) == FAILURE) goto RETURN;

            if(CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].LB,
                &price2) == FAILURE) goto RETURN;
            *payoff += LinFunc[i].Slope * (price1 - price2);

            if(CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].LB,
                &price1) == FAILURE) goto RETURN;
            *payoff -= (LinFunc[i].Slope * LinFunc[i].LB + 
                        LinFunc[i].Intercep) * price1;

            if(CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].Root,
                &price1) == FAILURE) goto RETURN;

            if(CRXQPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].RB,
                &price2) == FAILURE) goto RETURN;
            *payoff -= LinFunc[i].Slope * (price1 - price2);

            if(CRXQBinPricer(
                mqCnd,
                CRXQ_PUT,
                LinFunc[i].RB,
                &price1) == FAILURE) goto RETURN;
            *payoff += (LinFunc[i].Slope * LinFunc[i].RB + 
                        LinFunc[i].Intercep) * price1;
        }
    }

    *payoff += K[0];
    
    status = SUCCESS;

 RETURN:
    
    /* restore values */
    mqCnd->muMQ  = muMQSav;
    mqCnd->sigMQ = sigMQSav;

    return status;

} /* CRXQPay1D_ComplexSPD */

