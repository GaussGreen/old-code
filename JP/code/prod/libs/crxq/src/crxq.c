/******************************************************************************
 * Module:      CRXQ
 * Submodule:
 * File:        crxq.c  
 * Function:
 * Author:      Credit QRD
 * Revision:    $Header: $
 *****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <float.h>

#include <crxmultiq/include/crmultiq.h>
#include <crxflow/include/crxerror.h>
#include <crxflow/include/optbsq.h>

#include "crxq.h"

/*f----------------------------------------------------------------------------
 *
 * simple function to calculate (exp(qi(x-xi))-1)/qi (internal function)
 */
static double qxq(
    double q,                   /* (I) q smile parameter         */
    double x,                   /* (I) x                         */
    double x0)                  /* (I) x0                        */
{
    if(fabs(q) < TINY)
    {
        return x-x0;
    }
    else
    {
        return (exp(q*(x-x0))-1.0)/q;
    }
}

/*f---------------------------------------------------------------------------
 * Q3MQSolveQ
 *                                   exp(q * A) - 1
 * Newton-Raphson solution of q |->  -------------- + B = 0. (taken from Q3 lib)
 *                                         q
 */
static int Q3MQSolveMap4Q(
    double  A,         /* (I) */
    double  B,         /* (I) */
    double *q          /* (I/O) initial guess/result */
    )
{
    int         status = FAILURE;
    int         count  = 0;

    double      qq, f, dfdq;

    /* guess */
    qq = -2 * ((B / A) + 1) / A;
    qq = COLLAR(-CRXQ_QGUESS_LIMIT, CRXQ_QGUESS_LIMIT, qq);

    /* use expansion in vicinity of q = 0 */
    if (fabs(qq) < CRXQ_Q_SHIFT)
        f = A * (1. +.5 * A * qq) + B;
    else
        f = ((exp(qq * A) - 1.) / qq) + B;

    while(fabs(f) > CRXQ_QSOLV_TOL)
    {
        if (++count > CRXQ_QSOLV_ITER ) goto RETURN;

        if (fabs(qq) < CRXQ_Q_SHIFT)
            dfdq = A * A / 2;
        else
            dfdq = ((A * exp(qq * A)) / qq) - ((exp(qq * A) - 1.) / (qq * qq));
        
        if (fabs(dfdq) > TINY)
            qq -= f / dfdq;
        else
            goto RETURN;

        if (fabs(qq) < CRXQ_Q_SHIFT)
            f = A * (1. +.5 * A * qq) + B;
        else
            f = ((exp(qq * A) - 1.) / qq) + B;
    }

    *q = qq;

    status = SUCCESS;

 RETURN:

    return(status);

} /* Q3MQSolveMap4Q */
/*f----------------------------------------------------------------------------
 *BSQ macro: Q Black-Scholes on an interval BSQInt from optbsq.h */
#define  BSQ_INT(m,s,q,x,y) \
    (BSQInt((q)*(s),(q)*(m),(((x)-(m))/(s)),(((y)-(m))/(s))))

/*f----------------------------------------------------------------------------
 * Single Q Black Pricer
 *
 * Price of a Call/Put with Vega using Q version of Black&Scholes.
 */
int CRXQBSQPricer(
    double  Y,        /* (I) Fwd yield                     */
    double  K,        /* (I) Strike                        */
    double  T,        /* (I) Option expiration in years    */
    double  s,        /* (I) Annualized volatility         */
    double  Q,        /* (I) Q weight                      */
    long    I,        /* (I) Instrument                    */
    double *P)        /* (O) Price & Vega                  */
{
    static char routine[] = "CRXQBSQPricer";

    double C;          /* Call price                        */
    double V = 0.;     /* Vega                              */
    double d;          /* d in N(d) in Black & Scholes      */
    double st;         /* Sigma * sqrt(T) * Q               */
    double Y1, K1;     /* Modified fwd and strike           */
    double r;          /* Adjusted strike                   */
    long   vOn;        /* Vega calc on or off               */

    if (Y < CRXQ_MIN_FWD) 
    {
        DR_Error("%s failed: forward %lf less than minimum permitted %lf.\n", 
            routine, Y, CRXQ_MIN_FWD);
        return FAILURE;
    }

    /* compute vega? */
    vOn = (CRXQ_PAY_TYPE(I) == CRXQ_VEGA);

    /* for zero volatility return intrinsic */
    if (s * sqrt(T) < CRXQ_MIN_VOL) 
    {
        C = MAX(Y-K, 0.);
        V = 0.;

        switch(CRXQ_COP_TYPE(I)) 
        {
        case CRXQ_CALL: 
            P[0] = C; 
            if (vOn) P[1] = V; 
            break;
        case CRXQ_PUT:            
            P[0] = C -(Y-K); 
            if (vOn) P[1] = V; 
            break;
        default:
            DR_Error("%s failed: illegal option type %d.\n", 
                routine, CRXQ_COP_TYPE(I));
            return FAILURE;
        }

        return SUCCESS;
    }

    if (fabs(Q) > CRXQ_Q_SHIFT) 
    {
        st  = s * sqrt(T) * Q;
        r   = (K/Y-1.) * Q + 1.;
        /* if K is not in the range of the distribution */
        if (r < TINY) 
        {
            C = MAX(Y-K, 0.);
            V = 0.;
        } 
        else 
        {
            d   = - log(r) / st;
            st /= 2.;

            Y1  = Y / Q;
            K1  = K - Y + Y1;

            C   = Y1 * NormCum(d+st) - K1 * NormCum(d-st);        
            if (vOn) V = Y * exp(-(d*d+st*st)/2.) * sqrt(r*T) * INVSQRT2PI;
        }
    } 
    else 
    {
        st  = s * sqrt(T);
        d   = - (K/Y-1.)/st;

        C   = Y * st * exp(-d*d/2) * INVSQRT2PI + (Y-K) * NormCum(d);
        if (vOn) V = Y * exp(-(d*d)/2.) * sqrt(T) * INVSQRT2PI;
    }

    switch(CRXQ_COP_TYPE(I)) 
    {
    case CRXQ_CALL: 
        P[0] = C; 
        if (vOn) P[1] = V; 
        break;
    case CRXQ_PUT:            
        P[0] = C- (Y-K); 
        if (vOn) P[1] = V; 
        break;
    default:
        DR_Error("%s failed: X illegal option type %d.\n", 
                routine, CRXQ_COP_TYPE(I));
        return FAILURE;
    }

    return SUCCESS;

}  /* CRXQBSQPricer */


/*f----------------------------------------------------------------------------
 * MultiQ Black pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 */
int CRXQPricer(
    const CRXQDATA  *mq,           /* (I) MQ data                            */
    long             optType,      /* (I) option type                        */
    double           strike,       /* (I) option strike                      */
    double          *price)        /* (O) option price & vega                */
{
    static char      routine[]="CRXQPricer";
    int              status = FAILURE;

    double           fwdRate = mq->fwdRate;
    double           muMQ    = mq->muMQ;
    double           sigMQ   = mq->sigMQ;

    long   nq;                     /* number of intervals                    */
    double strkMult;               /* strike/fwdRate                         */
    double strkMultMQ;             /* strike/fwdRate after affine adjustment */
    const double *q = NULL;        /* pointer to q parameters                */
    const double *k = NULL;        /* pointer to strike ratios               */
    const double *d = NULL;        /* pointer to strike slopes               */
    double *x = NULL;              /* pointer to gaussian space boundaries   */

    /* affine adjustment constants */
    double K = mq->K;
    double C = mq->C;
    double KC = K*C;

    double temp, xtemp, d2q, parity;
    long   i, idx;
    int    sign;

    /* check option type */
    if (optType != CRXQ_CALL && optType != CRXQ_PUT) goto RETURN;

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) goto RETURN;

    /* for zero volatility, price intrinsic */
    if (sigMQ < CRXQ_MIN_VOL) 
    {
       if( mq->calcFwd == TRUE)
        {
            if(CRXQMap(
                mq,
                mq->muMQ,
                &fwdRate) == FAILURE) goto RETURN;
        }
        *price = MAX(CRXQ_COP_TO_COP(optType) * (fwdRate - strike), 0);
        status = SUCCESS;
        goto RETURN;
    }

    /* branch into low/high strikes (left/right); subtle      */
    /* choice needed in calibration: for ATM option, use left */
    /* intervals for put and high intervals for call          */
    strkMult = strike/fwdRate;

    /* apply affine adjustment */
    strkMultMQ = (K - 1)/(KC) + strkMult/ (KC);

    if (strkMultMQ < 1.) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        d = mq->dL;
        x = (double*)mq->xL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        d = mq->dR;
        x = (double*)mq->xR;
        sign = 1;
    }

    /* Gaussian integration cutoff for the last value at "infinity" */
    x[nq] = muMQ + sign * (CRXQ_CUM_NUM_STDEV) * sigMQ;
      
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
        strkMultMQ         * BSQ_INT(muMQ, sigMQ, 0,      x[nq],    xtemp) -
        (k[idx] - d2q)     * BSQ_INT(muMQ, sigMQ, 0,      x[idx+1], xtemp) -
        d2q * temp         * BSQ_INT(muMQ, sigMQ, q[idx], x[idx+1], xtemp);

    /* other intervals */
    for (i=idx+1; i<nq; i++) 
    {
        d2q = d[i]/q[i];
        temp = exp(-q[i] * x[i]);
        *price -= 
            ((k[i] - d2q)  * BSQ_INT(muMQ, sigMQ, 0,      x[i+1],   x[i]) +
             d2q * temp     * BSQ_INT(muMQ, sigMQ, q[i],   x[i+1],   x[i]));
    }

   /* apply affine adjustment */
    *price *= KC;  

    /* calculate put/call parity forward */
    if (mq->calcFwd == FALSE)
    {
        parity = 1.;
    }    
    else
    {
        parity = 0.;
        if ((strkMult<(1. + K*(C-1.)) && optType==CRXQ_CALL)
            || (strkMult>=(1. + K*(C-1.)) && optType==CRXQ_PUT))
        {

            /* Compute integral of H(x)dm(x)= E[H(x)] */

            /* do current side */
            for (i=0; i<nq; i++) 
            {
                d2q = d[i]/q[i];
                temp = exp(-q[i] * x[i]);
                parity 
                    -= sign * (
                        (k[i] - d2q) * BSQ_INT(muMQ, sigMQ, 0, x[i+1], x[i]) +
                        d2q * temp * BSQ_INT(muMQ, sigMQ, q[i], x[i+1], x[i])
                        );
            }

            /* do other side */
            if (sign > 0)
            {
                nq   = mq->nbQL;
                q    = mq->qL;
                k    = mq->kL;
                d = mq->dL;
                x = (double*)mq->xL;
                sign = -1;
            }
            else
            {
                nq   = mq->nbQR;
                q    = mq->qR;
                k    = mq->kR;
                d = mq->dR;
                x = (double*)mq->xR;
                sign = 1;
            }

            /* Gaussian integration cutoff */
            x[nq] = muMQ + sign * (CRXQ_CUM_NUM_STDEV) * mq->sigMQ;

            for (i=0; i<nq; i++) 
            {
                d2q = d[i]/q[i];
                temp = exp(-q[i] * x[i]);
                parity 
                    -= sign * (
                        (k[i] - d2q) * BSQ_INT(muMQ, sigMQ, 0, x[i+1], x[i]) +
                        d2q * temp * BSQ_INT(muMQ, sigMQ, q[i], x[i+1], x[i])
                        );
            }

            /* affine adjustment */

            parity = K*(C*parity - 1.0) + 1.0;

        } /* calculate put/call parity */
    }

    /* put/call parity: this has to be done in the original variable */
    if (strkMult < (1. + K * (C-1.)))
    {
        /* put-call parity for low strikes */
        if (optType == CRXQ_CALL) *price += (parity - strkMult);
    } 
    else 
    {
        /* put-call parity for high strikes */
        if (optType == CRXQ_PUT)  *price -= (parity - strkMult);
    }

    /* rescale by fwd */
    *price *= fwdRate;

    status = SUCCESS;
 RETURN:

    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* CRXQPricer */

/*f----------------------------------------------------------------------------
 * CRXQLevPricer
 *
 * MultiQ Black pricer with leverage     
 *
 * Uses Q parameters and the normal mean/stdev to compute call/put price
 */
int CRXQLevPricer(
    const CRXQDATA    *mq,                  /* (I) MQ data             */
    long               optType,             /* (I) option type         */
    double             strike,              /* (I) option strike       */
    double             leverage,            /* (I) leverage            */
    double            *price)               /* (O) option price & vega */
{
    static char routine[] = "CRXQLevPricer";
    int         status    = FAILURE;

    double      lPrice, lStrike;
    long        lOptType;

    /* handle -ve leverage */
    if (leverage < 0.)
    {
        if (optType == CRXQ_CALL)
        {
            lOptType = CRXQ_PUT;
        }
        else if(optType == CRXQ_PUT)
        {
            lOptType = CRXQ_CALL;
        }
        else
        {
            goto RETURN;
        }
    }
    else
    {
        lOptType = optType;
    }

    /* adjust strike */
    lStrike = strike / leverage;

    if (CRXQPricer(
         mq,
         lOptType,
         lStrike,
         &lPrice) == FAILURE) goto RETURN;

    *price = fabs(leverage) * lPrice;
        
    status = SUCCESS;
    
 RETURN:

    if (status == FAILURE)
    {
        DR_Error("%s: Failed\n", routine);
    }

    return status;

} /* CRXQLevPricer */


/*f----------------------------------------------------------------------------
 * CRXQBinPricer
 *
 * MultiQ binary pricer     
 *
 * Uses Q parameters and the normal mean/stdev to compute binary price
 */
int CRXQBinPricer(
    const CRXQDATA *mq,                 /* (I) MQ data             */
    long            optType,            /* (I) call/put            */
    double          strike,             /* (I) option strike       */
    double         *price)              /* (O) option price & vega */
{
    static char     routine[] ="CRXQBinPricer";
    int             status    = FAILURE;

    double          fwdRate   = mq->fwdRate;
    double          muMQ      = mq->muMQ;
    double          sigMQ     = mq->sigMQ;

    double          strkMultMQ, strkMult, temp, xtemp, d2q;
    const double    *q = NULL;
    const double	*k = NULL;
    const double	*d = NULL;
    double          *x = NULL;
    long            nq, idx;         
    int             sign;

    /* affine adjustment constants */
    double C = mq->C;
    double K = mq->K;
    double KC = K*C;

    /* check option type */
    if (optType != CRXQ_CALL && optType != CRXQ_PUT) 
    {
        DR_Error("%s failed: illegal option type - must be CRXQ_CALL or \
                   CRXQ_PUT\n", routine);
        goto RETURN;
    }

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) 
    {
        DR_Error("%s failed: forward rate (%lf) too low\n", routine, fwdRate);
        goto RETURN;
    }

    /* branch into low/high strikes (left/right); subtle
     * choice needed in calibration: for ATM option, use left
     * intervals for put and high intervals for call */
    strkMult = strike/fwdRate;

    /* apply affine adjustment */
    strkMultMQ = (K - 1)/(KC) + strkMult/ (KC);

    if (strkMultMQ < 1.0) 
    {
        nq = mq->nbQL;
        q = mq->qL;
        k = mq->kL;
        x = (double*)mq->xL;
        d = mq->dL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        k = mq->kR;
        x = (double*)mq->xR;
        d = mq->dR;
        sign = 1;
    }
  
    /* Gaussian integration cutoff for the last value at "infinity" */
    x[nq] = muMQ + sign * (CRXQ_CUM_NUM_STDEV) * sigMQ;
      
    /* find strike location */
    idx = 0;
    while (idx < nq - 1 && sign * (strkMultMQ - k[idx+1]) > 0.) idx++;
    
    /* find normal point corresponding to strike */
    d2q = d[idx]/q[idx];
    temp = 1 + (strkMultMQ - k[idx])/d2q;
    xtemp = (temp > 0 ? x[idx] + log(temp)/q[idx] : x[nq]);

    /* for zero volatility, price intrinsic */
    if (sigMQ > CRXQ_MIN_VOL) 
    {
        *price = NormCum((xtemp - muMQ) / sigMQ);
    }
    else
    {
        *price = (strkMult >= 1.0 ? 1. : 0.);
    }

    if (optType == CRXQ_CALL) *price = 1. - *price;

    status = SUCCESS;

 RETURN:
    
    return status;

} /* CRXQBinPricer */


/*f----------------------------------------------------------------------------
 * CRXQCalibrateAffineAdjustments
 *
 * Calibration of affine adjustment parameters K and C to match forward and
 * ATM call price 
 * 
 * Used for the exact and approximate q distributions
 */
int CRXQCalibrateAffineAdjustments(
    CRXQDATA*     mq)        /* (I/O) MQ data */
{
    static char   routine[] = "CRXQCalibrateAffineAdjustments";
    int           status    = FAILURE;

    double        expecH = 0.0;     /* expected value of H(x) */
    double        atmCallPrice;     /* the calculated ATM call price */
    int           i;
    double        d2q;

    mq->K = 1.0;
    mq->C = 1.0;

    /* for zero vol: must make sure that H(muMQ) maps to the forward 
       Set K=1.0, since degenerate case. */
    if (mq->sigMQ < TINY)
    {
        mq->C = 1.0;
        mq->K = 1.0;
        mq->optATM = 0.0;
        if (CRXQMap(
             mq,
             mq->muMQ,
             &expecH) == FAILURE)
        {
            DR_Error("%s failed: unable to calculate mapping of muMQ for zero \
                          vol calibration.\n", routine);
            goto RETURN;
        }
        mq->C = 1/expecH;
        status = SUCCESS;
        goto RETURN;
    }

    /* Compute the expectation of the mapping function H(x) */
    /* Left side */
    for (i=0; i<mq->nbQL; i++)
    {
        d2q = mq->dL[i]/mq->qL[i];
        expecH += (mq->kL[i] - d2q) * fabs(BSQ_INT(mq->muMQ,
                                                   mq->sigMQ,
                                                   0.0,
                                                   mq->xL[i],
                                                   mq->xL[i+1]));

        expecH += d2q*exp(-mq->qL[i]*mq->xL[i]) * fabs(BSQ_INT(mq->muMQ,
                                                               mq->sigMQ,
                                                               mq->qL[i],
                                                               mq->xL[i],
                                                               mq->xL[i+1]));
    }
    /* Right side */
    for (i=0; i<mq->nbQR; i++)
    {
        d2q = mq->dR[i]/mq->qR[i];
        expecH += (mq->kR[i] - d2q) * fabs(BSQ_INT(mq->muMQ,
                                                   mq->sigMQ,
                                                   0.0,
                                                   mq->xR[i],
                                                   mq->xR[i+1]));

        expecH += d2q*exp(-mq->qR[i]*mq->xR[i]) * fabs(BSQ_INT(mq->muMQ,
                                                               mq->sigMQ,
                                                               mq->qR[i],
                                                               mq->xR[i],
                                                               mq->xR[i+1]));
    }
    /* Calculate C parameter */
    mq->C = 1/expecH;

    /* Since E[Y] = fwdRate, (K (C expecH-1) + 1) = 1. Also, if fwdX is the point
       in gaussian X-space corresonding to the forward rate,
       fwdRate = fwdRate (K (C H(fwdX)-1) + 1)
       fwdX = H_inv(expecH)
       Since the price of the call is
       E[(Y-fwdRate)+] = fwdRate K C E[(H(x) - H(fwdX)+] we can price the call
        assuming that K=1 and then rescale to get the correct call price
       */
    mq->K = 1.0;

    if (CRXQPricer(
         mq,
         CRXQ_CALL,
         mq->fwdRate,
         &atmCallPrice) == FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM option with forward rate \
                  %lf.\n", routine, mq->fwdRate);
        goto RETURN;
    }

    /* The target price is the simple BS ATM option price */
    if (CRXQBSQPricer(
         mq->fwdRate,
         mq->fwdRate,
         mq->optExpy,
         mq->sigATM,
         1.0,
         CRXQ_CALL,
         &(mq->optATM)) == FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM option with forward rate=%lf, \
         expiry=%lf, vol=%lf.\n", routine, mq->fwdRate, mq->optExpy, mq->sigATM);
        goto RETURN;
    }

    /* Now set K so that the ATM option price matches, and the expectation is the
       forward rate */
    mq->K = (mq->optATM)/atmCallPrice;

    status = SUCCESS;
    RETURN:
  if (status == FAILURE) 
  {
    DR_Error("%s: failed\n", routine);
  }

  return status;
}


/*f----------------------------------------------------------------------------
 * CRXQCalibrateCRXQDATA
 *
 *      
 *
 * Calibrate a CRXQ data structure from an input set of variables
 */

int CRXQCalibrateCRXQDATA(
    double           optionExpiry,  /* (I) option expiry                      */
    double           forwardRate,   /* (I) forward rate                       */
    double           volATM,        /* (I) ATM vol                            */
    int              numberOfQs,    /* (I) number of q's specified            */
    const double    *inputQ,        /* (I) q parameters                       */
    const double    *inputDelta,    /* (I) middle value to both L and R, same */
    CRXQDATA*        mq)            /* (I/O) MQ data                          */
{
    static char      routine[] = "CRXQCalibrateCRXQDATA";
    int              status    = FAILURE;

    double           sigTotal  = volATM * sqrt(optionExpiry); 
    int              i;

    if (optionExpiry<TINY)
    {   
        DR_Error("%s failed: non-positive time to expiry.\n", routine);
        goto RETURN;
    } else if (volATM<TINY) 
    {
        DR_Error("%s failed: non-positive ATM vol.\n", routine);
        goto RETURN;
    } else if (numberOfQs<=0 || numberOfQs%2!=0) 
    {
        DR_Error("%s failed: You must provide an even number of \
                    q-parameters.\n", routine);
        goto RETURN;
    }

    /* need to check that the order of deltas is correct and no deltas are beyond
       the cutoff */
    for (i=0; i < numberOfQs-2; i++)
    {
        if (inputDelta[i+1]-inputDelta[i] < TINY)
        {
            DR_Error("%s failed: Deltas %lf [%d] and  %lf [%d] are too close \
                together for no. qs %d.\n", 
                routine, inputDelta[i], i, inputDelta[i], (i+1), numberOfQs);
          goto RETURN;
        }
    }
    /* note: add checks for delta cutoff and more resonable spacing */

    /* set-up */
    /* misc set-up */
    mq->nbQL = numberOfQs/2;
    mq->nbQR = mq->nbQL;
    mq->fwdRate   = forwardRate;
    mq->sigATM    = volATM;
    mq->optExpy   = optionExpiry;

    /* do not recompute forward */
    mq->calcFwd = FALSE;

    /* these are chosen. The calibration comes from the K and C params */
    mq->sigMQ     = sigTotal;
    mq->muMQ      = 0.0; /* this is by choice */

    /* left side */
    for(i = 0;i < mq->nbQL; i++)
    /* Reset all the x values to match the strike boundaries, deltaL.   */
    {
        mq->qL[i] = inputQ[mq->nbQL - i - 1];
        /* If the q's are zero, set them to a manageable small number */
        if (fabs(mq->qL[i])<CRXQ_Q_SHIFT) mq->qL[i] = CRXQ_Q_SHIFT;
        mq->xL[i] = sigTotal * NormCumInv(inputDelta[mq->nbQL - i - 1]);
    }
    mq->xL[mq->nbQL] = - sigTotal * CRXQ_CUM_NUM_STDEV; /* cut off beyond this */

    /* Build k and d coefficients using freedom to set first values 
     * so that H(0)=1 and H'(0)=1 and then
     * by requiring continuity and differentiability */
    mq->kL[0] = 1.0;
    mq->dL[0] = 1.0;
    for(i = 1;i < mq->nbQL;i++)
    {
        mq->dL[i] = mq->dL[i-1]*exp(mq->qL[i-1] * (mq->xL[i]-mq->xL[i-1]));
        mq->kL[i] = mq->kL[i-1] + mq->dL[i-1] * qxq(mq->qL[i-1],
                                                    mq->xL[i],
                                                    mq->xL[i-1]);
    }  
   
    /* right side */
    for(i = 0;i < mq->nbQR;i++)
    {
        mq->qR[i] = inputQ[mq->nbQR + i];
        /* If the q's are zero, set them to a manageable small number */
        if (fabs(mq->qR[i])<CRXQ_Q_SHIFT) mq->qR[i] = CRXQ_Q_SHIFT; 
        mq->xR[i] = sigTotal * NormCumInv(inputDelta[mq->nbQR -1 + i]);
    }
    mq->xR[mq->nbQR] = sigTotal * CRXQ_CUM_NUM_STDEV;
    mq->kR[0] = 1.0;
    mq->dR[0] = 1.0;
    for(i = 1;i < mq->nbQR;i++)
    {
        mq->dR[i] = mq->dR[i-1]*exp(mq->qR[i-1] * (mq->xR[i]-mq->xR[i-1]));
        mq->kR[i] = mq->kR[i-1] + mq->dR[i-1] * qxq(mq->qR[i-1],
                                                    mq->xR[i],
                                                    mq->xR[i-1]);
    } 
  
    /* Calibrate the affine adjustment parameters so that the expected value  
     * is the forward, and the ATM option price matches BS */
    if (CRXQCalibrateAffineAdjustments(
          mq) == FAILURE)
    {
        DR_Error("%s failed: unable to calculate K and C parameters \
                    in Q-calibration.", 
            routine);
        goto RETURN;
    }
 
    status = SUCCESS;
    

    RETURN:
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }

  return status;
}



/*f----------------------------------------------------------------------------
 * MultiQ: mapping function
 *
 * Evaluates mapping function at normal point x given CRXQDATA data
 */
int CRXQMap(
    const CRXQDATA    *mq,      /* (I) MQ data       */
    double             xval,    /* (I) normal point  */
    double            *yield    /* (O) mapped yield  */
    )
{
    static char      routine[]="CRMQMap";
    int              status = FAILURE;

    /* distribution data */

    double fwdRate = mq->fwdRate;

    long         nq;        /* number of intervals      */
    const double *q = NULL; /* pointer to q parameters  */
    const double *d = NULL; /* pointer to d parameters  */
    const double *k = NULL; /* pointer to strike ratios */
    const double *x = NULL; /* pointer to gaussian space strike boundaries */
   
    long   idx;
    int    sign;

    /* check forward */
    if (fwdRate < CRXQ_MIN_FWD) 
    {
        DR_Error("%s: fwdRate = %lf is too small.\n", routine, fwdRate);
        goto RETURN;
    }

    /* for zero volatility, return forward no matter what the point is */
    if (mq->sigMQ < CRXQ_MIN_VOL &&
        mq->calcFwd == FALSE)
    {
        *yield = fwdRate;
        status = SUCCESS;
        goto RETURN;
    }    

    /* branch into low/high strikes (left/right) */
    if (xval <= mq->xL[0])
    {
        nq = mq->nbQL;
        q = mq->qL;
        d = mq->dL;
        k = mq->kL;
        x = mq->xL;
        sign = -1;
    } 
    else 
    {
        nq = mq->nbQR;
        q = mq->qR;
        d = mq->dR;
        k = mq->kR;
        x = mq->xR;
        sign = 1;
    }
	
    /* find x location */
    idx = 0;
    while (idx < nq - 1 && sign * (xval - x[idx+1]) >= 0.) idx++;

    /* evaluate mapping function */

    *yield = k[idx] + d[idx] * (exp(q[idx]*(xval-x[idx])) - 1)/q[idx];

    *yield = mq->fwdRate*(mq->K*(mq->C * *yield - 1.0) + 1.0);
	
    status = SUCCESS;
 RETURN:
    if (status == FAILURE) 
    {
        DR_Error("%s: Failed\n", routine);
    }
    
    return status;
} /* CRXQMap */

/*-----------------------------------------------------------------------------
 * CRXQGridPricer
 *
 * Uses numerical grid integration to price a payoff on a multi-q distribution.
 *
 */
int CRXQGridPricer(
    CRXQPAYOFF  *optPayoff, /* (I) option payoff structure */
    CRXQFPAYOFF *payFunc,   /* (I) option payoff function  */
    CRXQDATA    *mq,        /* (I) measure data            */
    double      *premium    /* (O) fwd premium/rate        */
    )
{
    double grid[CRXQ_1D_INT_PTS];
    double yield[CRXQ_1D_INT_PTS];
    double dens[CRXQ_1D_INT_PTS];
    double point[2];
    double densNorm;
    double fwdRate;

    int status = FAILURE;

    /* check for degenerate case */
    if (mq->sigMQ > CRXQ_MIN_VOL)
    {
        /* calculate MQ yield density */
        if (CRXQDens(            
            mq,
            -CRXQ_1D_NUM_STDEV,   /* start */
            CRXQ_1D_NUM_STDEV,    /* end   */
            CRXQ_1D_INT_PTS,
            grid,
            yield,
            dens,
            &densNorm) == FAILURE) goto RETURN;

        /* integrate density * payoff (from crxq.h) */
        if (CRXQSimpsonPricer1D(
            optPayoff,
            payFunc,
            grid,
            yield,
            dens,
            densNorm,
            premium) == FAILURE) goto RETURN;
    }
    else
    {   
        fwdRate = mq->fwdRate;
        if( mq->calcFwd == TRUE)
        {
            if(CRXQMap(
               mq,
               mq->muMQ,
               &fwdRate) == FAILURE) goto RETURN;
        }
        point[0] = fwdRate;
        point[1] = mq->muMQ;
        if ((*payFunc)(
            optPayoff,
            point,
            0., 
            premium) == FAILURE) goto RETURN;
    }
    
    status = SUCCESS;

 RETURN:

    return status;

} /* CRXQGridPricer */


/*f----------------------------------------------------------------------------
 * CRXQDens
 *
 * Generates density for multi-q distribution from underlying measure.
 *
 */
int CRXQDens(
    const CRXQDATA *mq, /* (I) MQ measure                */
    double  start,      /* (I) start in std normal space */
    double  end,        /* (I) end in std normal space   */
    long    numGridPts, /* (I) num pts to generate       */
    double *grid,       /* (O) normal grid points        */
    double *yields,     /* (O) yield points              */
    double *dens,       /* (O) density                   */
    double *normC       /* (O) density normalization     */
    )
{
    static char routine[] = "CRXQDens";
    int         status    = FAILURE;

    double step, xval;
    long   i;

    step  = (end - start)/(numGridPts-1);
    xval  = start - step;

    /* calculate mapping and density */
    for (i = 0; i < numGridPts; i++)
    {
        xval   += step;
        grid[i] = mq->sigMQ * xval + mq->muMQ;
        if (CRXQMap(
            mq,
            grid[i],
            &(yields[i])) == FAILURE) goto RETURN;
             
        dens[i] = NormDens(xval);
    }

    /* normalize distribution */
    *normC = Q3SimpsIntegral(
        step,
        numGridPts,
        dens);
    if (*normC < TINY) 
    {
        DR_Error("%s: Non-positive total mass.\n", routine);
        goto RETURN;
    }

    status = SUCCESS;

 RETURN:

    return status;
} /* CRXQDens */

  
/*f----------------------------------------------------------------------------
 * CRXQApproximateFAWithMQ
 *
 *      
 *
 * Calibrate a CRXQDATA multi-q measure to approximate a forward-adjusted 
 * measure (used for faster bivariate pricing)
 */
int CRXQApproximateFAWithMQ(
    const CRXQFADATA* fa,   /* (I) Given FA measure to be approximated      */
    int               numQ, /* (I) No. of qs on each side of fwd to fit     */
    CRXQDATA*         pa    /* (O) Resulting multi-q measure                */    
    )             
{

    static char routine[] = "CRXQApproximateFAWithMQ";
    int         status    = FAILURE;

    CRXQPAYOFF pf; 
    /* Numerical FA density */
    double grid[CRXQ_1D_INT_PTS];
    double yield[CRXQ_1D_INT_PTS];
    double dens[CRXQ_1D_INT_PTS];
    double densNorm;
    double cumProb, loProb, hiProb;
    double tailStrike[2];
    double tailProb[2] = {1. - CRXQ_TAIL_PROB, CRXQ_TAIL_PROB};
    long   iLo,iHi;
    double deltaX;
    double fQTarget;
    double strikeStepL;
    double strikeStepR;
    double* k = NULL;
    double* x = NULL;
    double* d = NULL;
    double* q = NULL;
    double* cdf = NULL;
    double cdfL[CRXQ_NBQ];
    double cdfR[CRXQ_NBQ];
    int    i, side, nbQ;
    int sign;

    if (numQ<0 || numQ>CRXQ_NBQ) 
    {
        DR_Error("%s failed: illegal numQ, %d, should be positive and <=%d.\n", 
            routine, numQ, CRXQ_NBQ);
        goto RETURN;
    }

    /* calculate the FA yield density                                */
    /* note: use CRXQFADens & CRXQSimpsonPricer1D separately instead */
    /* of CRXQFAPricer for efficiency from reuse of density          */
    if (CRXQFADens(            
        fa,
        -CRXQ_1D_NUM_STDEV,   /* start */
        CRXQ_1D_NUM_STDEV,    /* end   */
        CRXQ_1D_INT_PTS,
        &(grid[0]),
        &(yield[0]),
        &(dens[0]),
        &densNorm) == FAILURE) 
    {
        DR_Error("%s failed: unable to compute numerical FA density.\n", 
            routine);
        goto RETURN;
    }

    /* calculate fwd and set in new MQ */
    if (CRXQSimpsonPricer1D(
        &pf,
        CRXQPay1D_Yield,
        grid,
        yield,
        dens,
        densNorm,
        &(pa->fwdRate)) == FAILURE) goto RETURN;

    if (pa->fwdRate < TINY) 
    {
        DR_Error("%s failed: Adjusted forward: %f < TINY.\n", 
            routine, pa->fwdRate);        
        goto RETURN;
    }

    /* find strikes for tail probabilities */
    cumProb = 0.;
    hiProb  = tailProb[0] * densNorm / CRXQ_1D_STEP;
    loProb  = tailProb[1] * densNorm / CRXQ_1D_STEP;        
    iLo     = 0;
    iHi     = CRXQ_1D_INT_PTS - 1;

    for (i = 0; i < CRXQ_1D_INT_PTS; i++) 
    {
        /* tail probs */
        cumProb += dens[i];
        if (cumProb < loProb)
            iLo = i;
        else if (cumProb < hiProb)
            iHi = i;
    }

    tailStrike[0] = yield[iHi];
    tailStrike[1] = yield[iLo];

    /* determine location of the min and max strike points. 
       strike points are equally spaced over yield range.
       number of q intervals is CRXQ_NBQ on each side of the distribution */
    strikeStepL = (1.0-tailStrike[1]/pa->fwdRate)/CRXQ_NBQ;
    strikeStepR = (tailStrike[0]/pa->fwdRate - 1.0)/CRXQ_NBQ; 

    /* Calculate the ATM call price and implied BS vol for the FA measure */
    pf.cop = 1; /* call */
    pf.strike = pa->fwdRate;
    if (CRXQSimpsonPricer1D(
        &pf,
        CRXQPay1D_Vnl,
        grid,
        yield,
        dens,
        densNorm,
        &(pa->optATM)) == FAILURE) 
    {
        DR_Error("%s failed: unable to price ATM call under FA measure", routine);
        goto RETURN;
    }
    /* Compute the BS implied vol for the ATM option (BSQImpVol from optbsq.c) */
    if (BSQImpVol(
        &(pa->sigATM),
        pa->fwdRate,
        pa->fwdRate, 
        fa->mq->optExpy, 
        1.0,
        pa->optATM, 
        CRXQ_CALL, 
        fa->mq->sigATM) == FAILURE)
    {
        DR_Error("%s failed: unable to calculate implied vol for ATM option \
        with price %lf.\n", routine, pa->optATM);
        goto RETURN;
    }

    /* Set parameters that will stay unchanged in the approximated measure */
    /* Set q, k, x, d in pa measure to original q, k, x, d inputs (this will be
       used for zero vol case */
    
    pa->optExpy = fa->mq->optExpy;
    pa->sigMQ = pa->sigATM*sqrt(pa->optExpy);
    pa->muMQ = 0.0;
    pa->nbQL = CRXQ_NBQ;
    pa->nbQR = CRXQ_NBQ;

    for (i=0; i<fa->mq->nbQL; i++)
    {
        pa->qL[i] = fa->mq->qL[i];
        pa->kL[i] = fa->mq->kL[i];
        pa->dL[i] = fa->mq->dL[i];
    }

    for (i=0; i<fa->mq->nbQR; i++)
    {
        pa->qR[i] = fa->mq->qR[i];
        pa->kR[i] = fa->mq->kR[i];
        pa->dR[i] = fa->mq->dR[i];
    }

    for (i=0; i<fa->mq->nbQL+1; i++)
    {
        pa->xL[i] = fa->mq->xL[i];
    }

    for (i=0; i<fa->mq->nbQR+1; i++)
    {
        pa->xR[i] = fa->mq->xR[i];
    }

    pa->calcFwd = fa->mq->calcFwd;

    /* no affine adjustments before bootstrapping */
    pa->K = 1.0;
    pa->C = 1.0;

    /* if vol is too small match the forward and skip further calibration */
    if (pa->sigMQ < CRXQ_MIN_VOL_CALIB)
    {
        pa->sigMQ = 0.;
        pa->muMQ = 0.;
        status = SUCCESS;
        goto RETURN;
    } 

    /* Initialize strike points */
    for (i=1; i<CRXQ_NBQ; i++)
    {
        pa->kL[i] = pa->kL[i-1] - strikeStepL;
        pa->kR[i] = pa->kR[i-1] + strikeStepR;
    }

    /* compute binary puts on left side */
    for (i=0; i<CRXQ_NBQ; i++)
    {
        pf.cop = -1; /* binary put */
        pf.strike = pa->kL[i] * pa->fwdRate;
        if (CRXQSimpsonPricer1D(
            &pf,
            CRXQPay1D_Bin,
            grid,
            yield,
            dens,
            densNorm,
            &(cdfL[i])) == FAILURE);
    }

    /* compute binary puts on right side */

    for (i=0; i<CRXQ_NBQ; i++)
    {
        pf.cop = -1; /* binary put */
        pf.strike = pa->kR[i] * pa->fwdRate;
        if (CRXQSimpsonPricer1D(
            &pf,
            CRXQPay1D_Bin,
            grid,
            yield,
            dens,
            densNorm,
            &(cdfR[i])) == FAILURE);
    }

    /* calibrate the x, q, d values to match the binary put prices */
    /* initialize q values for NR iteration */
    for (side = 0; side<2; side++)
    {
        if (side==0)
        {
            /* LEFT SIDE */
            nbQ = pa->nbQL;
            x = pa->xL;
            q = pa->qL;
            d = pa->dL;
            k = pa->kL;
            cdf = cdfL;
            sign = -1;
        }
        else
        {
            /* RIGHT SIDE */
            nbQ = pa->nbQR;
            x = pa->xR;
            q = pa->qR;
            d = pa->dR;
            k = pa->kR;
            cdf = cdfR;
            sign = 1;
        }

        /*  set boundary conditions for implied q's  */

        pa->qR[CRXQ_NBQ - 1] = 1E-6; /* normal in upper tail */
        pa->qL[CRXQ_NBQ - 1] = 1E-6; /* normal in lower tail */

        x[0] = NormCumInv(cdf[0])*pa->sigMQ + pa->muMQ;
        for (i=1; i<nbQ; i++)
        {
            /* x[i] is determined by the CDF:
               H(x[i]) = k[i], and so
               CDF[k[i].fwdRate] = Phi((x[i]-muMQ)/sigMQ) */
            x[i] = NormCumInv(cdf[i])*pa->sigMQ + pa->muMQ;
            if (fabs(d[i-1])<TINY)
            {
                DR_Error("%s failed: zero d-value at i=%d on %c side.",
                            routine, i, (side==0 ? 'L' : 'R'));
                goto RETURN;
            }
            /* Note that d[0] is set already */
            /* calibrate d, q to match x and k: */
            deltaX = x[i] - x[i-1];
            fQTarget = (k[i-1] - k[i])/d[i-1];

            /* solve for q by Newton Raphson*/
            if (Q3MQSolveMap4Q(
                deltaX,
                fQTarget,
                &(q[i-1])) == FAILURE) goto RETURN;

            /* If the q's are zero, set them to a manageable small number */
            if (fabs(q[i-1])<CRXQ_Q_SHIFT) q[i-1] = CRXQ_Q_SHIFT;

            /* Differentiability of H(x) requires d[i] = d[i-1].exp(q[i-1].deltaX) */
            d[i] = d[i-1] * exp(q[i-1]*deltaX);
        }

    /* Gaussian integration cutoff for the last value at "infinity" */
        x[nbQ] = pa->muMQ + sign * (CRXQ_CUM_NUM_STDEV) * pa->sigMQ;
    }

    /* Finally, calibrate C and K to match the ATM call price and the forward. */
    if (CRXQCalibrateAffineAdjustments(
         pa) == FAILURE) 
    {
        DR_Error("%s failed: unable to calibrate affine parameters C, K.");
        goto RETURN;
    }
    
    status = SUCCESS;
 RETURN:
    return status;

} /* CRXQApproximateFAWithMQ */




