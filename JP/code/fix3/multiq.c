/******************************************************************************
 * Module:      Q3TMX
 * Submodule:
 * File: multiq.c       
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <float.h>

#include "q3_tmx.h"

/* Q3TMXBSQ macro: Q Black-Scholes on an interval */
#define  Q3TMXBSQ_INT(m,s,q,x,y) \
    (Q3TMXBSQInt((q)*(s),(q)*(m),(((x)-(m))/(s)),(((y)-(m))/(s))))

/* MQ pricer for ___/|___ or ___|\___ payoff -> see end of file */
static int MQSharkFin(double,double,double,double,double *,double *);


/*f----------------------------------------------------------------------------
 * MultiQ Black pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 */
int Q3TMXMQPricer(
    MQDATA            *mq,           /* (I) MQ data  */
    long              optType,       /* (I) option type */
    double            strike,        /* (I) option strike */
    double            *price)        /* (O) option price & vega */
{
    static char      routine[]="Q3TMXMQPricer";
    int              status = FAILURE;

    double fwdRate = mq->fwdRate;
    double muMQ    = mq->muMQ;
    double sigMQ   = mq->sigMQ;

    long   nq;          /* number of intervals                    */
    double strkMult;    /* strike/fwdRate                         */
    double strkMultMQ;  /* strike/fwdRate after affine adjustment */
    double *q;          /* pointer to q parameters                */
    double *k;          /* pointer to strike ratios               */

    /* distribution data */
    double x[Q3TMX_NBQ+1];  /* Gaussian nodes                          */
    double d[Q3TMX_NBQ];    /* derivative of mapping function at nodes */

    /* affine adjustment constants */
    double fwdC = mq->fwdC;
    double volC = mq->volC;

    double temp, xtemp, d2q, parity;
    long   i, idx;
    int    sign;

    /* check option type */
    if (optType != Q3TMX_CALL && optType != Q3TMX_PUT) goto RETURN;

    /* check forward */
    if (fwdRate < Q3TMX_MIN_FWD) goto RETURN;

    /* for zero volatility, price intrinsic */
    if (sigMQ < Q3TMX_MIN_VOL) 
    {
        if (mq->calcFwd == TRUE)
        {
            if (Q3TMXMQMap(
                mq,
                mq->muMQ,
                &fwdRate) == FAILURE) goto RETURN;
        }

        *price = MAX(Q3TMX_COP_TO_COP(optType) * (fwdRate - strike), 0);
        status = SUCCESS;
        goto RETURN;
    }

    /* distribution data at forward point */
    x[0] = 0.;
    d[0] = 1.;

    /* branch into low/high strikes (left/right); subtle      */
    /* choice needed in calibration: for ATM option, use left */
    /* intervals for put and high intervals for call          */
    strkMult = strike/fwdRate;

    /* apply affine adjustment */
    strkMultMQ = 1 + (strkMult - fwdC) / volC;

    if (strkMultMQ < 1.) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        sign = 1;
    }

    /* shift all q away from zero */
    for (i=0; i<nq; i++) if (fabs(q[i])<Q3TMX_Q_SHIFT) q[i] = Q3TMX_Q_SHIFT;

    /* compute MQ distribution data from MQDATA */
    for (i=0; i<nq-1; i++) 
    {
        temp = 1 + q[i] * (k[i+1] - k[i])/d[i];
        /* if (k[i+1] * fwdRate) is not actually reached stop */
        if (temp <= 0) 
        {
            nq = i+1;
            break;
        }
        d[i+1] = d[i] * temp;
        x[i+1] = x[i] + log(temp)/q[i];
    }

    /* Gaussian integration cutoff */
    x[nq] = muMQ + sign * (Q3TMX_CUM_NUM_STDEV) * sigMQ;
      
    /* find strike location */
    idx = 0;
    while (idx < nq - 1 && sign * (strkMultMQ - k[idx+1]) >= 0.) idx++;

    /* find normal point corresponding to strike */
    d2q = d[idx]/q[idx];
    temp = 1 + (strkMultMQ - k[idx])/d2q;
    xtemp = (temp > 0 ? x[idx] + log(temp)/q[idx] : x[nq]);

    /* key step(!): price put for low strikes, call for high strikes */

    /* last interval */
    temp = exp(-q[idx] * x[idx]);
    *price = 
        strkMultMQ         * Q3TMXBSQ_INT(muMQ, sigMQ, 0,      x[nq],    xtemp) -
        (k[idx] - d2q)     * Q3TMXBSQ_INT(muMQ, sigMQ, 0,      x[idx+1], xtemp) -
        d2q * temp         * Q3TMXBSQ_INT(muMQ, sigMQ, q[idx], x[idx+1], xtemp);

    /* other intervals */
    for (i=idx+1; i<nq; i++) 
    {
        d2q = d[i]/q[i];
        temp = exp(-q[i] * x[i]);
        *price -= 
            ((k[i] - d2q)  * Q3TMXBSQ_INT(muMQ, sigMQ, 0,      x[i+1],   x[i]) +
             d2q * temp     * Q3TMXBSQ_INT(muMQ, sigMQ, q[i],   x[i+1],   x[i]));
    }

    /* apply affine adjustment */
    *price *= volC;

    /* calculate put/call parity forward */
    if (mq->calcFwd == FALSE)
    {
        parity = 1.;
    }
    else
    {
        parity = 0.;

        /* do current side */
        for (i=0; i<nq; i++) 
        {
            d2q = d[i]/q[i];
            temp = exp(-q[i] * x[i]);
            parity 
                -= sign * (
                    (k[i] - d2q) * Q3TMXBSQ_INT(muMQ, sigMQ, 0, x[i+1], x[i]) +
                    d2q * temp * Q3TMXBSQ_INT(muMQ, sigMQ, q[i], x[i+1], x[i])
                    );
        }

        /* do other side */
        if (sign > 0)
        {
            nq   = mq->nbQL;
            q    = mq->qL;
            k    = mq->kL;
            sign = -1;
        }
        else
        {
            nq   = mq->nbQR;
            q    = mq->qR;
            k    = mq->kR;
            sign = 1;
        }

        /* shift all q away from zero */
        for (i=0; i<nq; i++) if (fabs(q[i])<Q3TMX_Q_SHIFT) q[i] = Q3TMX_Q_SHIFT;

        /* compute MQ distribution data from MQDATA */
        x[0] = 0.;
        d[0] = 1.;
        for (i=0; i<nq-1; i++) 
        {
            temp = 1 + q[i] * (k[i+1] - k[i])/d[i];
            /* if (k[i+1] * fwdRate) is not actually reached stop */
            if (temp <= 0) 
            {
                nq = i+1;
                break;
            }
            d[i+1] = d[i] * temp;
            x[i+1] = x[i] + log(temp)/q[i];
        }

        /* Gaussian integration cutoff */
        x[nq] = muMQ + sign * (Q3TMX_CUM_NUM_STDEV) * mq->sigMQ;

        for (i=0; i<nq; i++) 
        {
            d2q = d[i]/q[i];
            temp = exp(-q[i] * x[i]);
            parity 
                -= sign * (
                    (k[i] - d2q) * Q3TMXBSQ_INT(muMQ, sigMQ, 0, x[i+1], x[i]) +
                    d2q * temp * Q3TMXBSQ_INT(muMQ, sigMQ, q[i], x[i+1], x[i])
                    );
        }

        /* affine adjustment */
        parity = volC * parity + fwdRate * (fwdC - volC);

    } /* calculate put/call parity */

    /* put/call parity: this has to be done in the original variable */
    if (strkMult < fwdC) 
    {
        /* put-call parity for low strikes */
        if (optType == Q3TMX_CALL) *price += (parity - strkMult);
    } 
    else 
    {
        /* put-call parity for high strikes */
        if (optType == Q3TMX_PUT)  *price -= (parity - strkMult);
    }

    /* rescale by fwd */
    *price *= fwdRate;
    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3TMXMQPricer */


/*f----------------------------------------------------------------------------
 * MultiQ: mapping function
 *
 * Evaluates mapping function at normal point x given MQDATA data
 */
int Q3TMXMQMap(
    MQDATA            *mq,      /* (I) MQ data       */
    double             xval,    /* (I) normal point  */
    double            *yield)   /* (O) mapped yield  */
{
    static char      routine[]="Q3TMXMQMap";
    int              status = FAILURE;

    /* distribution data */
    double x[Q3TMX_NBQ];    /* Gaussian nodes                          */
    double d[Q3TMX_NBQ];    /* derivative of mapping function at nodes */

    /* affine adjustment constants */
    double fwdC = mq->fwdC;
    double volC = mq->volC;

    double fwdRate = mq->fwdRate;

    long   nq; /* number of intervals      */
    double *q; /* pointer to q parameters  */
    double *k; /* pointer to strike ratios */
    
    double temp;
    long   i, idx;
    int    sign;

    /* check forward */
    if (fwdRate < Q3TMX_MIN_FWD) goto RETURN;

    /* for zero volatility, return forward */
    if (mq->sigMQ < Q3TMX_MIN_VOL) 
    {
        *yield = fwdRate;
        return (SUCCESS);
    }    

    /* branch into low/high strikes (left/right) */
    if (xval <= 0.)
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        sign = 1;
    }

    /* distribution data at forward point */
    x[0] = 0.;
    d[0] = 1.;

    /* shift all q away from zero */
    for (i=0; i<nq; i++) if (fabs(q[i])<Q3TMX_Q_SHIFT) q[i] = Q3TMX_Q_SHIFT;

    /* compute MQ distribution data from MQDATA */
    for (i=0; i<nq-1;i++) 
    {
        temp = 1 + q[i] * (k[i+1] - k[i])/d[i];
        /* if (k[i+1] * fwdRate) is not actually reached stop */
        if (temp <= 0) 
        {
            nq = i+1;
            break;
        }
        d[i+1] = d[i] * temp;
        x[i+1] = x[i] + log(temp)/q[i];
    }

    /* find x location */
    idx = 0;
    while (idx < nq - 1 && sign * (xval - x[idx+1]) >= 0.) idx++;

    /* evaluate mapping function */
    *yield = k[idx] + d[idx] * (exp(q[idx] * (xval - x[idx])) - 1)/q[idx];

    /* apply affine adjustment */
    *yield = fwdC + volC * (*yield - 1.);

    /* rescale by forward */
    *yield *= fwdRate;

    status = SUCCESS;

 RETURN:
    
    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3TMXMQMap */


/*f----------------------------------------------------------------------------
 * MultiQ: initialize data for calibration
 *
 * Determine 3 target option prices and expected fwd.
 *           target[0] = expected rate
 *           target[1] = ATM_CALL
 *           target[2]  = L1_CALL
 *           target[3]  = R1_PUT
 *
 */
int Q3TMXMQTargetSV(
    SVDATA            *sv,             /* (I) SV data given */    
    double            *mqGuess,        /* (I) MQ initial guess */
    MQDATA            *mq)             /* (O) MQ data */    
{
    static char      routine[]="Q3TMXMQTargetSV";
    int              status = FAILURE;

    double totVol, totVovSq, vovSq2Vol, q;
    double strikeATM, strikeL1, strikeR1, strikeR2;
    double priceRatio, yieldBdry;

    /* check if calibration needed */
    totVol = mq->sigATM * sqrt(mq->expiry);
    if (totVol < Q3TMX_MIN_VOL) 
    {
        mq->muMQ = 0.;
        mq->sigMQ = 0.;
        mq->putR[0] = mq->callL[0] = mq->putR[1] = mq->callL[1] = 0;
        status = SUCCESS;
        goto RETURN;
    } 

    /* initial guess for internal MQ parameters: use analytic or given guess */
    /* check positivity of sigma guess                                       */
    if (mqGuess == NULL) 
    {
        /* analytic */
        totVol    = sv->sigATM * sqrt(sv->expiry);
        totVovSq  = sv->volVolSV * sv->volVolSV * sv->expiry;
        vovSq2Vol = (totVol < Q3TMX_MIN_VOL)?0:(totVovSq/totVol);
        q         = 1. - sv->q;
    
        mq->sigMQ = totVol / (1+totVovSq/sqrt(3.));
        mq->muMQ  = -0.5 * q * totVol * mq->sigMQ;
        mq->qL[0] = q + vovSq2Vol * (q/3. -1.);    
        mq->qR[0] = mq->qL[0] + 2. * vovSq2Vol;
    } 
    else 
    {
        /* check positivity of sigma guess */
        if (mqGuess[0] < Q3TMX_TINY) 
        {
            Q3TMXErrMsg("%s: First guess must be > 0.\n", routine);
        }
        /* use given guess */
        mq->sigMQ = mqGuess[0];
        mq->muMQ  = mqGuess[1];
        mq->qL[0] = mqGuess[2];
        mq->qR[0] = mqGuess[3];
    }
       
    /* update tails */
    mq->qL[1] = (1.- mq->tauL) * (1. + mq->qL[0] * (mq->kL[1]-1.))/mq->kL[1];
    mq->qR[1] = (1.- mq->tauR) * (1. + mq->qR[0] * (mq->kR[1]-1.))/mq->kR[1];

    /* update negative tail */
    mq->qL[2] = (1. + mq->qL[0] * (mq->kL[1]-1.)
                    + mq->qL[1] * (mq->kL[2]-mq->kL[1])) / (mq->tauN+mq->kL[2]);

    /* crossover points */
    strikeATM  = 1 * mq->fwdRate;
    strikeL1   = mq->kL[1] * mq->fwdRate;
    strikeR1   = mq->kR[1] * mq->fwdRate;
    strikeR2   = mq->kR[2] * mq->fwdRate;

    /* check positivity and ordering of strikes */
    if (strikeL1 <= Q3TMX_MIN_FWD ||
        strikeL1 >= mq->fwdRate   ||
        strikeR1 <= mq->fwdRate   ||
        strikeR1 >= strikeR2) 
    {
        Q3TMXErrMsg("%s: SV strikes not correctly positioned.\n", routine);
        goto RETURN;
    }

    /* compute stoch vol option prices as targets */
    if (Q3TMXSVPricer(sv, Q3TMX_CALL, strikeATM, &(mq->callL[0]))== FAILURE ||
        Q3TMXSVPricer(sv, Q3TMX_CALL, strikeL1,  &(mq->callL[1]))== FAILURE ||
        Q3TMXSVPricer(sv, Q3TMX_PUT,  strikeR1,  &(mq->putR[1])) == FAILURE) 
        goto RETURN;
    mq->optATM = mq->putR[0] = mq->callL[0];

    /* bounds for positivity of call spreads */
    mq->muOverSigMax =  Q3TMXNormCumInv ((mq->callL[1] - mq->callL[0]) / 
                   (strikeATM - strikeL1));
    mq->muOverSigMin = -Q3TMXNormCumInv ((mq->putR[1] - mq->putR[0]) /
                                    (strikeR1 - strikeATM));
    if (mq->muOverSigMax < mq->muOverSigMin + Q3TMX_TINY)
    {
        Q3TMXErrMsg ("Implied Stoch Vol OTM prices are not arbitrage-free.\n");
        goto RETURN;
    }     

    /* check arbitrage-free restrictions on tau parameters */
    if (mq->tauL < 1. - Q3TMX_TINY) 
    {
        priceRatio = (mq->callL[1] - mq->fwdRate + strikeL1) / mq->callL[0];

        /* determine lower bound of yield distribution */
        switch (mq->nbQL)
        {
        case 1: yieldBdry = -999999; break;
        case 2: yieldBdry = strikeL1 * mq->tauL / (mq->tauL - 1.); break;
        case 3: yieldBdry = -mq->tauN * mq->fwdRate; break;
        default: Q3TMXErrMsg ("Boundary conditions not implemented for "
                           ">=4 left Q intervals\n"); goto RETURN;
        }
        
        if ((strikeL1 - yieldBdry) / (mq->fwdRate - yieldBdry) < priceRatio) 
        {
            Q3TMXErrMsg("%s: MultiQ is not arbitrage-free for TauL = %f.\n",
                     routine, mq->tauL);
            goto RETURN;
        }
    }
 
    if (mq->tauR > 1. + Q3TMX_TINY) 
    {
        priceRatio = (mq->putR[1] + mq->fwdRate - strikeR1) / mq->putR[0];

        /* determine upper bound of yield distribution */
        switch (mq->nbQL)
        {
        case 1: yieldBdry = 999999; break;
        case 2: yieldBdry = strikeR1 * mq->tauR / (mq->tauR - 1.); break;
        case 3: yieldBdry = ((mq->tauR > mq->kR[2] / (mq->kR[2] - mq->kR[1])) ?
                            strikeR1 * mq->tauR / (mq->tauR - 1.) : 999999); break;
        default: Q3TMXErrMsg ("Boundary conditions not implemented for "
                           ">=4 right Q intervals\n"); goto RETURN;
        }

        if ((yieldBdry - strikeR1) / (yieldBdry - mq->fwdRate) < priceRatio) 
        {
            Q3TMXErrMsg("%s: MultiQ is not arbitrage-free for TauR = %f.\n",
                     routine, mq->tauR);
            goto RETURN;
        }
    }

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;

} /* Q3TMXMQTargetSV */


/*f----------------------------------------------------------------------------
 * Q3TMXSVToMQ
 *
 * Multiq: calibrate MultiQ parameters to Stoch Vol interface 
 */
int Q3TMXSVToMQ(
    double             fwdRate,         /* (I) fwd rate */
    double             sigATM,          /* (I) atm vol */
    double             expiry,          /* (I) in years */
    double             *smile,          /* (I) smile: MQ smile */
    MQDATA             *mq)             /* (0) MQ data  */
{
    static char      routine[]="Q3TMXSVToMQ";
    int              status = FAILURE;

    SVDATA sv;

    /* smile: initialize SV */
    if (Q3TMXSVSmileInit(
        fwdRate,
        sigATM,
        expiry,
        smile,
        &sv) == FAILURE) goto RETURN;

    /* calibrate SV */
    if (Q3TMXSVCalib(&sv) == FAILURE) goto RETURN;

    /* smile: initialize MQ structure */
    if (Q3TMXMQSmileInit(
        fwdRate,
        sigATM,
        expiry,
        smile,
        mq) == FAILURE) goto RETURN;

    /* set initial MQ params & define targets for calibration */
    if (Q3TMXMQTargetSV(
        &sv,
        NULL,
        mq) == FAILURE) goto RETURN;

    /* calibrate MQ */
    if (Q3TMXMQCalib(mq) == FAILURE) goto RETURN;

    status = SUCCESS;

 RETURN:

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /* Q3TMXSVToMQ */


/* MQ pricer for ___/|___ or ___|\___ payoff: light version for speed.
 * To be used in calibration only.
 * Only works on middle intervals, based on k[0] = 1, x0 = 0.
 * Does not work for q = 0. 
 * Returns premium divided by forward and the derivative wrt q.
 */
static int MQSharkFin(
    double muMQ,          /* (I) MQ internal mean */
    double sigMQ,         /* (I) MQ internal standard deviation */
    double q,             /* (I) current value of q */
    double k1,            /* (I) first strike ratio right/left of forward */
    double *price,        /* (O) premium */
    double *deriv)        /* (O) derivative wrt q */
{
    /* middle interval assumptions */
    double x0 = 0., k0 = 1.;

    /* q almost zero cutoff */
    double qCutoff = Q3TMX_Q_SHIFT / sigMQ;

    /* calculation variables */
    double temp, x1, bsInt0, bsIntQ, xInf;
    int    sign = (k1 < 1. ? -1 : 1) ;
    
    /* pick point at infinity */
    xInf         = muMQ + sign * Q3TMX_CUM_NUM_STDEV * sigMQ;

    if (fabs(q) < qCutoff) 
    {
        double q0, qf1, qf2;

        /* evaluate at qCutoff */
        q0 = qCutoff;

        /* update inverse mapping of strike */ 
        temp     = 1. + q0 * (k1 - k0);                 
        x1       = (temp > 0. ? x0 + log(temp) / q0 : xInf );  
                                                         
        /* partial Q3TMXBSQ_INT results */  
        bsInt0   = Q3TMXBSQ_INT(muMQ, sigMQ, 0,  x1, x0); 
        bsIntQ   = Q3TMXBSQ_INT(muMQ, sigMQ, q0, x1, x0);   
       
        /* price target option  */ 
        qf1      = (1. - k1) * bsInt0 + (bsIntQ - bsInt0) / q0;   

        /* evaluate at -qCutoff */
        q0 = -qCutoff;

        /* update inverse mapping of strike */ 
        temp     = 1. + q0 * (k1 - k0);                 
        x1       = (temp > 0. ? x0 + log(temp) / q0 : xInf );  
                                                         
        /* partial Q3TMXBSQ_INT results */  
        bsInt0   = Q3TMXBSQ_INT(muMQ, sigMQ, 0,  x1, x0); 
        bsIntQ   = Q3TMXBSQ_INT(muMQ, sigMQ, q0, x1, x0);   
       
        /* price target option  */ 
        qf2      = (1. - k1) * bsInt0 + (bsIntQ - bsInt0) / q0;   

        /* price target option by linear interpolation */
        *price = 0.5 * ((1. - q/qCutoff) * qf2 + (1. + q/qCutoff) * qf1);
        
        /* numerical derivative wrt q */
        *deriv = 0.5 * (qf1 - qf2) / qCutoff;
    } 
    else 
    {
        /* calculation variables used in derivative */
        double sigMQSq      = sigMQ * sigMQ;
        double muPlusSigSqQ = muMQ + q * sigMQSq;
        double invQ         = 1. /q;

        /* update inverse mapping of strike */ 
        temp     = 1. + q * (k1 - k0);                 
        x1       = (temp > 0. ? x0 + log(temp) / q: xInf );  
                                                         
        /* partial Q3TMXBSQ_INT results */  
        bsInt0   = Q3TMXBSQ_INT(muMQ, sigMQ, 0, x1, x0); 
        bsIntQ   = Q3TMXBSQ_INT(muMQ, sigMQ, q, x1, x0);   
       
        /* price target option  */ 
        *price   = (1. - k1) * bsInt0 + invQ * (bsIntQ - bsInt0);   

        /* analytical derivative wrt q */
        *deriv  = (invQ / q)             * (bsInt0 - bsIntQ) +
            invQ * (muPlusSigSqQ  *  bsIntQ +
                    sigMQ  * Q3TMXNormDens(muMQ/sigMQ) * 
                    (exp(x1 * (muPlusSigSqQ - 0.5 * x1) / sigMQSq) -
                     exp(x0 * (muPlusSigSqQ - 0.5 * x0) / sigMQSq)));
    } /* if fabs(q) */

    return SUCCESS;
}   /* MQSharkFin */


/*f----------------------------------------------------------------------------
 * Q3TMXMQCalibMidQs
 *
 * Uses Newton-Raphson to match call/put spreads assuming the mean and std 
 * are known and that the q's have been assigned initial guesses
 */
int Q3TMXMQCalibMidQs(
    MQDATA          *mq)       /* (I/O) smeasure with MQ data        */
{
    static char      routine[]="Q3TMXMQCalibMidQs";
    int              status = FAILURE;

    /* local read-only MQ data */
    double fwdRate   = mq->fwdRate;
    double muMQ      = mq->muMQ;
    double sigMQ     = mq->sigMQ;
    int    calibType = mq->calibType;

    /* middle q's (modified by procedure) */
    double *q0 = NULL;

    /* used in Newton-Raphson */
    double q, qf, qdf;
    int    qcount, qfound[2];

    double k, target;
    int    i;

    /* fixed number of steps: unconditional success */
    if (calibType == Q3TMX_MQ_FIX) 
    {
        qfound[0] = TRUE;
        qfound[1] = TRUE;
    } 
    else 
    {
        qfound[0] = FALSE;
        qfound[1] = FALSE;
    }

    /* successive Newton-Raphson for each middle q: 
     * calibrate qL0 for i = 0 and qR0 for i = 1   */
    for (i=0; i<=1; i++) 
    {
        /* select strike, q parameter, and set up target option price 
         * note: target is divided by fwdRate */
        if (i == 0) 
        {
            k       = mq->kL[1];
            target  = (mq->callL[1] - mq->callL[0]) / fwdRate - 
                      (1. - k) * Q3TMXNormCum(muMQ/sigMQ);
            q0      = &(mq->qL[0]);
        } 
        else 
        {
            k       = mq->kR[1];
            target  = (mq->putR[1] -  mq->putR[0]) / fwdRate +
                      (1. - k) * Q3TMXNormCum(-muMQ/sigMQ);
            q0      = &(mq->qR[0]);
        }

        /* Is this feasible with the current muMQ and sigMQ? */
        if (target < 0) 
        {
            Q3TMXErrMsg("%s: Algorithm has exited feasible region.\n", routine);
            goto RETURN;
        }

        /* Newton-Raphson loop */
        for (qcount = 0; qcount < mq->qSteps; qcount++) 
        {
            /* compute target function and its derivative wrt q */
            MQSharkFin(muMQ, sigMQ, *q0, k, &qf, &qdf);  
                     
            /* Q3TMX_MQ_TOL: check error and exit if within tolerance */
            if ((calibType == Q3TMX_MQ_TOL) &&
                (fabs(qf - target) < mq->otmTol * mq->callL[0])) 
            { 
                qfound[i] = TRUE;
                break;
            }
                
            if (fabs(qdf) < Q3TMX_MQ_RESN) 
            {
                Q3TMXErrMsg("%s: df/dq = 0 in Newton-Raphson.\n", routine);
                goto RETURN;
            }
            *q0 -= (qf - target)/qdf; 
        } /* for qcount */
        
        /* update tails */
        q = *q0;
        if (i==0) 
        {
            mq->qL[1] = (1. - mq->tauL) * (1. + q * (k - 1.))/k;

            /* update negative tail as well */
            mq->qL[2] = (1. + q * (k - 1.) + mq->qL[1] * (mq->kL[2] - k)) 
                      / (mq->tauN + mq->kL[2]);
        } 
        else 
        {
            mq->qR[1] = (1. - mq->tauR) * (1. + q * (k - 1.))/k;
        }
    }  /* for i loop */

    if (qfound[0] && qfound[1]) status = SUCCESS; 

  RETURN:  

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3TMXMQCalibMidQs */


/*f----------------------------------------------------------------------------
 * Q3TMXMQCalibMean
 *
 * Uses Newton-Raphson to match forward assuming the std is known and the mean
 * has been assigned an initial guess. Recalibrates middle q's at each step.
 */
int Q3TMXMQCalibMean(
    MQDATA  *mq)  /* (I/O) MQ data */
{
    static char routine[] = "Q3TMXMQCalibMean";
    int         status    = FAILURE;

    /* local MQ data */
    double fwdRate   = mq->fwdRate;
    int    calibType = mq->calibType;
 
    /* Newton-Raphson variable */
    double *muMQ = &(mq->muMQ);

    /* used in Newton-Raphson */
    double mf, mff, mdf;
    int    mcount, mfound;
    double opt1, opt2;
    double qLBak, qRBak;
    double muStepMax, muStepMin;   

    /* fail is sigMQ <= 0 */
    if (mq->sigMQ < Q3TMX_TINY) goto RETURN;

    /* fixed number of steps: unconditional success */
    if (calibType == Q3TMX_MQ_FIX) 
    {
        mfound = TRUE;
    } 
    else 
    {
        mfound = FALSE;
    }   

    /* NR loop over normal mean */
    for (mcount = 0; mcount < mq->mSteps; mcount++) 
    { 
        /* store q's */
        qLBak = mq->qL[0];
        qRBak = mq->qR[0];

        /* vary mean and (if required) calibrate middle q's */
        *muMQ += mq->mDelta;
        if ((mq->calibQs) && (Q3TMXMQCalibMidQs(mq) == FAILURE)) goto RETURN;

        /* compute forward with these q's use left intervals for put and
         * right intervals for call */
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, fwdRate + Q3TMX_MQ_RESN, &opt1) == FAILURE ||
            Q3TMXMQPricer(mq, Q3TMX_PUT,  fwdRate - Q3TMX_MQ_RESN, &opt2) == FAILURE ) 
            goto RETURN;
        mff = opt1 - opt2;

        /* restore mean and q's, and (if required) recalibrate the q's */
        *muMQ -= mq->mDelta;
        mq->qL[0] = qLBak;
        mq->qR[0] = qRBak;
        if ((mq->calibQs) && (Q3TMXMQCalibMidQs(mq) == FAILURE)) goto RETURN;

        /* recompute forward with these q's; use left intervals for put and
         * right intervals for call */
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, fwdRate+Q3TMX_MQ_RESN, &opt1) == FAILURE ||
            Q3TMXMQPricer(mq, Q3TMX_PUT,  fwdRate-Q3TMX_MQ_RESN, &opt2) == FAILURE ) 
            goto RETURN;  
        mf = opt1 - opt2;

        /* Q3TMX_MQ_TOL: check error and exit if within tolerance */
        if ((calibType == Q3TMX_MQ_TOL) &&
            (fabs(mf) < mq->fwdTol * fwdRate)) { 
            mfound = TRUE;
            break;
        }

        /* NR step */
        mdf = (mf - mff)/mq->mDelta;
        if (fabs(mdf) < Q3TMX_MQ_RESN)  
        {
             Q3TMXErrMsg("%s: df/dm = 0 in Newton-Raphson.\n", routine);
             goto RETURN;
        }

        /* impose call spread constraint to NR step for muMQ */
        if (mq->calibSafe)
    {
        muStepMin = mq->sigMQ * mq->muOverSigMin - (*muMQ);
            muStepMax = mq->sigMQ * mq->muOverSigMax - (*muMQ);

            *muMQ += MAX (0.5 * muStepMin, MIN (0.5 * muStepMax, mf / mdf));
    }
        else
    {
            *muMQ += mf / mdf;
    }
    } /* for mcount */
    
    if (mfound) status = SUCCESS;

  RETURN:
    
    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3TMXCalibMean */


/*f----------------------------------------------------------------------------
 * Q3TMXMQCalib
 *
 * Given strikeL1, strikeR1, (external OTM strikes),
 * uses iterative Newton-Raphson to match 3 external option prices and expected.
 * NOTE: tail q dependence is implicit! 
 */
int Q3TMXMQCalib(
    MQDATA *mq)    /* (I/O) MQ data            */    
{
    static char routine[] = "Q3TMXMQCalib";
    int         status    = FAILURE;

    /* local MQ data */
    double fwdRate   = mq->fwdRate; 
    double sigATM    = mq->sigATM;
    double expiry    = mq->expiry;
    int    calibType = mq->calibType;

    double strikeL1, strikeR1, strikeR2;
    double stVolAtmPrice;  

    /* Newton-Raphson variable */
    double *sigMQ = &(mq->sigMQ);

    /* used in Newton-Raphson */
    double sf, sff, sdf;
    int    scount, sfound;
    double muBak, qLBak, qRBak;
    double sigStepMax, sigStepMin;

    double totVol; 

    /* compute total volatility; if too low, skip calibration and */
    /* set vol to zero (pricer then returns intrinsic value)      */
    totVol = sigATM * sqrt(expiry);
    if (totVol < Q3TMX_MIN_VOL_CALIB) 
    {
        mq->calibQs = FALSE;
        /* *sigMQ = 0.;
         status = SUCCESS;
         goto RETURN; */
    } 

    /* R-normalization */
    mq->callL[0] /= fwdRate; mq->callL[1] /= fwdRate;
    mq->putR[0]  /= fwdRate; mq->putR[1]  /= fwdRate;
    mq->optATM   /= fwdRate;
    mq->fwdRate = 1; 

    /* crossover points */
    strikeL1 = mq->kL[1] * mq->fwdRate;
    strikeR1 = mq->kR[1] * mq->fwdRate;
    strikeR2 = mq->kR[2] * mq->fwdRate;

    /* if not calibrating mid Qs, use different target */
    if (mq->calibQs == TRUE)
    {
        stVolAtmPrice = mq->callL[0];  
    }
    else
    {
        stVolAtmPrice = mq->optATM;
    }

    /* fixed number of steps: unconditional success in NR loop */
    if (calibType == Q3TMX_MQ_FIX) 
    {
        sfound = TRUE;
    } 
    else 
    {
        sfound = FALSE;
    }   

    /* iterative Newton-Raphson std > mean > q's */  
    for (scount = 0; scount < mq->sSteps; scount++) 
    { 
        /* NR loop over normal stdev; store other variables first */
        muBak = mq->muMQ;
        qLBak = mq->qL[0];
        qRBak = mq->qR[0];

        /* vary std and calibrate normal mean */
        *sigMQ += mq->sDelta;
        if (Q3TMXMQCalibMean(mq) == FAILURE) goto RETURN;

        /* compute ATM price with current MQ parameters */
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, mq->fwdRate, &sff) == FAILURE) goto RETURN;

        /* reset stdev, restore others, and recalibrate mean and q's */
        *sigMQ -= mq->sDelta;
        mq->muMQ  = muBak;
        mq->qL[0] = qLBak;
        mq->qR[0] = qRBak;
        if (Q3TMXMQCalibMean(mq) == FAILURE) goto RETURN;

        /* compute ATM price with current MQ parameters */
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, mq->fwdRate, &sf) == FAILURE) goto RETURN;

        /* Q3TMX_MQ_TOL: check error and exit if within tolerance */
        if ((calibType == Q3TMX_MQ_TOL) &&
            (fabs(sf - stVolAtmPrice) < mq->atmTol * stVolAtmPrice)) 
        { 
            sfound = TRUE;
            break;
        }
        
        /* Newton-Raphson step for stdev */
        sdf = (sf - sff)/mq->sDelta;
        if (fabs(sdf) < Q3TMX_MQ_RESN) 
        {
            Q3TMXErrMsg("%s: df/ds = 0 in Newton-Raphson.\n", routine);
            goto RETURN;
        }

        /* impose call spread constraint to NR step for sigMQ */
        if (mq->calibSafe)
    {
        if (mq->muOverSigMin > 0)
            {
            sigStepMax = mq->muMQ / (*sigMQ * mq->muOverSigMin) - 1.;
                sigStepMin = mq->muMQ / (*sigMQ * mq->muOverSigMax) - 1.;
        }
            else if (mq->muOverSigMax > 0)
        {
                sigStepMax = 2 * Q3TMX_S_MAX_STEP;
                sigStepMin = MAX (mq->muMQ / (*sigMQ * mq->muOverSigMin),
                                  mq->muMQ / (*sigMQ * mq->muOverSigMax)) - 1.;
        }
            else 
        {
                sigStepMax = mq->muMQ / (*sigMQ * mq->muOverSigMax) - 1.;
                sigStepMin = mq->muMQ / (*sigMQ * mq->muOverSigMin) - 1.;
        }

        /* add a priori bounds */
            sigStepMax = MIN (0.5 * sigStepMax,  Q3TMX_S_MAX_STEP);
            sigStepMin = MAX (0.5 * sigStepMin, -Q3TMX_S_MAX_STEP);
    }
        else
    {
            sigStepMax =  Q3TMX_S_MAX_STEP;
            sigStepMin = -Q3TMX_S_MAX_STEP;
    }

        /* NR step */
        *sigMQ += MAX(sigStepMin * (*sigMQ), MIN(sigStepMax * (*sigMQ), 
                      (sf - stVolAtmPrice)/sdf));

    } /* for scount */ 

    /* R-normalization */
    mq->fwdRate    = fwdRate; 
    mq->callL[0]  *= fwdRate; mq->callL[1] *= fwdRate;
    mq->putR[0]   *= fwdRate; mq->putR[1]  *= fwdRate;
    strikeL1      *= fwdRate; strikeR1     *= fwdRate;
    mq->optATM    *= fwdRate;    
    stVolAtmPrice *= fwdRate;
    
    /* apply affine adjustment to match forward and ATM exactly */
    if (mq->affineAdj == TRUE)
    {
        double opt1, opt2, fwdErr, K, fwdC, volC;
        
        /* compute forward calibration error */
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, fwdRate+Q3TMX_MQ_RESN, &opt1) == FAILURE ||
            Q3TMXMQPricer(mq, Q3TMX_PUT,  fwdRate+Q3TMX_MQ_RESN, &opt2) == FAILURE)
            goto RETURN;

        /* compute error on the forward */
        fwdErr = (opt1 - opt2) / fwdRate;

        /* option price at calibrated forward strike */
        K = (1 + fwdErr);
        if (Q3TMXMQPricer(mq, Q3TMX_CALL, K * fwdRate, &opt1) == FAILURE)
            goto RETURN; 

        /* compute affine coefficients */
        volC = stVolAtmPrice / opt1;
        fwdC = 1. + volC * (1. - K);

        /* update MQDATA */
        mq->fwdC  = fwdC;
        mq->volC  = volC;
    } 
    else
    {
        mq->fwdC = 1.;
        mq->volC = 1.;
    } /* affine adjustment */

    /* If calibration uses a fixed number of steps and mid Qs are 
     * being calibrated, check deltaL and deltaR vol;  with affine 
     * adjustment, no need to check ATM price.
     */
    if ((calibType == Q3TMX_MQ_FIX) && (mq->calibQs == TRUE)) 
    {
        double strike, priceSV, priceMQ;
        double sigOTMCalib, sigOTM, volErr;
        int    optType, i;

        for (i=0; i<=1; i++) 
        {
            if (i==0) 
            {
                strike = strikeL1;
                priceSV = mq->callL[1];
                optType = Q3TMX_CALL;
            } 
            else 
            {
                strike = strikeR1;
                priceSV = mq->putR[1];
                optType = Q3TMX_PUT;
            }

            if (Q3TMXBSQImpVol (
                fwdRate, 
                strike, 
                expiry, 
                1., 
                priceSV,
                optType, 
                sigATM, 
                &sigOTM) == FAILURE) goto RETURN;
            if (Q3TMXMQPricer(
                mq, 
                optType, 
                strike, 
                &priceMQ) == FAILURE) goto RETURN;
            if (Q3TMXBSQImpVol (
                fwdRate, 
                strike, 
                expiry, 
                1., 
                priceMQ,
                optType, 
                sigATM, 
                &sigOTMCalib) == FAILURE) goto RETURN;
   
            volErr = fabs (sigOTMCalib - sigOTM);
            if (volErr > mq->otmTol) 
            {        
                Q3TMXErrMsg ("%s: Delta%c Vol error (%f) is larger than %f.\n", 
                    routine, (i==0 ? 'L' : 'R'), volErr, mq->otmTol);
                goto RETURN;
            } 
        } /* for i */
    } /* end OTM vol check */

    if (sfound) status = SUCCESS;   

  RETURN:

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
    
    return status;
} /* Q3TMXMQCalib */


/*f----------------------------------------------------------------------------
 * Q3TMXMQSmileInit
 *
 * Multiq: initialize MQ smile.
 */
int Q3TMXMQSmileInit(
    double             fwdRate,         /* (I) fwd rate */
    double             sigATM,          /* (I) atm vol */
    double             expiry,          /* (I) in years */
    double             *smile,          /* (I) smile: MQ smile */
    MQDATA             *mq)             /* (0) MQ data  */
{
    static char      routine[]="Q3TMXMQSmileInit";
    int              status = FAILURE;

    double deltaL   = smile[4];
    double tauL     = smile[5];
    double deltaR   = smile[6];
    double tauR     = smile[7];
    double normTail = smile[8];
    double deltaN   = smile[10];
    double tauN     = smile[11];
    double totVol, nt;

    /* check MultiQ delta parameters */
    if (deltaL < Q3TMX_MIN_DELTA_L || deltaL > Q3TMX_MAX_DELTA_L) 
    {
        Q3TMXErrMsg("%s: DeltaL (%f) must be between %1.2f and %1.2f.\n",
                 routine, deltaL, Q3TMX_MIN_DELTA_L, Q3TMX_MAX_DELTA_L);
        goto RETURN;
    }
    if (deltaR < Q3TMX_MIN_DELTA_R || deltaR > Q3TMX_MAX_DELTA_R) 
    {
        Q3TMXErrMsg("%s: DeltaR (%f) must be between %1.2f and %1.2f.\n",
                 routine, deltaR, Q3TMX_MIN_DELTA_R, Q3TMX_MAX_DELTA_R);
        goto RETURN;
    }

    /* check that normal tail starts beyond deltaR */
    nt = Q3TMXNormCumInv (deltaR);
    if (normTail < nt + Q3TMX_TINY) 
    {
        Q3TMXErrMsg ("%s: Normal cutoff is too low (must be above %f.\n",
                  routine, nt);
        goto RETURN;
    }

    /* check MultiQ tau parameters */
    if (tauL < Q3TMX_MIN_TAU_L) 
    {
        Q3TMXErrMsg("%s: TauL (%f) must be >= %1.2f.\n", 
                 routine, tauL, Q3TMX_MIN_TAU_L);
        goto RETURN;
    }

    /* initialize MQ */
    mq->fwdRate   = fwdRate;
    mq->expiry    = expiry;
    mq->sigATM    = sigATM;

    mq->tauL      = tauL;
    mq->tauR      = tauR;
    mq->nbQL      = 3;
    mq->nbQR      = 3;
    mq->kL[0]     = 1.;
    mq->kR[0]     = 1.;

    /* convert delta points to multiples of forward */
    totVol = sigATM * sqrt(expiry);
    mq->kL[1]     = exp(Q3TMXNormCumInv(deltaL) * totVol);
    mq->kR[1]     = exp(Q3TMXNormCumInv(deltaR) * totVol);

    /* set normal tail */
    mq->kR[2]     = exp(normTail           * totVol); 
    mq->qR[2]     = 0.;

    /* set left tail */
    mq->kL[2]     = (deltaN<1E-14 ? 0. : exp(Q3TMXNormCumInv(deltaN) * totVol));
    mq->tauN      = tauN;

    /* no affine adjustment during calibration */
    mq->fwdC = 1.;
    mq->volC = 1.;
    mq->affineAdj = TRUE;

    /* calibrate mid Qs */
    mq->calibQs = TRUE;
    
    /* do not recompute forward */
    mq->calcFwd = FALSE;

    /* decode convergence criteria for calibration */
    if (Q3TMXDecodeNCK (smile[9], mq) == FAILURE) goto RETURN;

    status = SUCCESS;

 RETURN:
    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /* Q3TMXMQSmileInit */


/* Switch between 5Q and 6Q -- default is 5Q */
#define Q3TMX_6Q 6

/*f----------------------------------------------------------------------------
 * MultiQ: decode numerical configuration key (NCK) 
 */
int Q3TMXDecodeNCK(
    double  nck,  /* (I) Nck: calibration type and details */
    MQDATA  *mq)  /* (O) MQ data                           */
{
    static char      routine[]="Q3TMXDecodeNCK";
    int              status = FAILURE;

    unsigned long key = (unsigned long) nck;
    int  i0, i1, i2, i3, i4, i5, i6, i7;

    /* extract numerical information */
    i7 = key % 10; key /= 10;
    i6 = key % 10; key /= 10;
    i5 = key % 10; key /= 10;
    i4 = key % 10; key /= 10;
    i3 = key % 10; key /= 10;
    i2 = key % 100; key /= 100;
    i1 = key % 100; key /= 100;

    i0 = key;

    /* number of q's -- 5Q by default, unless user enters '6' */
    mq->nbQL = (i3 == Q3TMX_6Q ? 3 : 2);

    /* use safe method for mu and sig? */
    mq->calibSafe = (i0==Q3TMX_MQ_FIX_SAFE || i0==Q3TMX_MQ_TOL_SAFE);

    /* process numerical details */
    switch(i0) 
    {
    case Q3TMX_MQ_FIX:      
    case Q3TMX_MQ_FIX_SAFE:
        mq->calibType = Q3TMX_MQ_FIX;

        /* read number of steps */
        mq->sSteps = i1;
        mq->mSteps = i2;
        mq->qSteps = i4;

        /* Step size for numerical derivatives */
        mq->mDelta = Q3TMX_M_DELTA;
        mq->sDelta = Q3TMX_S_DELTA;

        /* option price tolerances */
        mq->fwdTol = 0.0001  * pow(2, -i5);
        mq->atmTol = 0.01    * pow(2, -i6);
        mq->otmTol = 0.01    * pow(2, -i7);

        break;
    case Q3TMX_MQ_TOL:
    case Q3TMX_MQ_TOL_SAFE:
        mq->calibType = Q3TMX_MQ_TOL;

        /* use internal limits for number of steps */
        mq->sSteps = Q3TMX_S_STEPS;
        mq->mSteps = Q3TMX_M_STEPS;
        mq->qSteps = Q3TMX_Q_STEPS;

        /* Step size for numerical derivatives */
        mq->mDelta = pow(10, -i5);
        mq->sDelta = pow(10, -i6);

     
        /* option price tolerances */
        mq->fwdTol = pow(10, -i5);
        mq->atmTol = pow(10, -i6);
        mq->otmTol = pow(10, -i7);

        break;
    default:
      Q3TMXErrMsg("%s: Invalid calibration type.\n", routine);
      goto RETURN;
    }

    status = SUCCESS;

  RETURN:

    if (status == FAILURE) 
    {
        Q3TMXErrMsg("%s: Failed\n", routine);
    }
  
    return status;
} /* Q3TMXDecode NCK */


