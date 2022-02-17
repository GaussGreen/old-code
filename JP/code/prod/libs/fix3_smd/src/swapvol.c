/****************************************************************************/
/*      Bootstrapping of swaption volatility into spot volatility.          */
/****************************************************************************/
/*      SWAPVOL.c                                                           */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fix123head.h"


#ifndef    MIN_SPOT_VOL_RATIO
#define    MIN_SPOT_VOL_RATIO    0.2
#endif

#ifndef    SPOT_VOL_FILTER_AMOUNT
#define    SPOT_VOL_FILTER_AMOUNT    0.95
#endif



/*****  Fix3_ExpDecay  ***********************************************************/
/*
*       Exponential decay function.
*/
double  Fix3_ExpDecay (  double  a,
                    double  t)
{
    double  x;


    /* If a=0 or t=0, Fix3_ExpDecay(a,t)=1 */
    if (fabs (a * t) < ERROR)
    {
        return (1.);
    }

    x = (1. - exp (- a * t)) / a / t;

    return (x);

}  /* Fix3_ExpDecay */




/*****  Fix3_ExpInt  ***********************************************************/
/*
*       Exponential decay function.
*/
double  Fix3_ExpInt (double  T1,
                     double  T2,
                     double  T3,
                     double  beta)
{
    double  x;

    if (fabs(beta) < TINY)
    {
        return(T2-T1);
    }
    else
    {
        x = (exp(-beta * (T3 - T2)) -
             exp(-beta * (T3 - T1)))/ beta;
    }

    return(x);
} /* Fix3_ExpInt */




/*****  Smd_MeanDrift  ***********************************************************/
/*
*       Smd Mean Drift
*/


double Smd_MeanDrift(double        X,            /* (I)   X     value  */
                     double        A,            /* (I)   Skew  level  */
                     double        B,            /* (I)   Smile level  */
                     double        C)            /* (I)   Intensity    */

{
    double mean_drift;


    mean_drift = 1 + A * (exp(C*X) - exp(-C*X))/(exp(C*X) + exp(-C*X)) 
                   + B * (1 - 2.0 / (exp(C*X) + exp(-C*X)));

    mean_drift *= X;

    return(mean_drift);
} 




/*****  Fix3_Tanh  ***********************************************************/
/*
*       Hyperbolic tangent function.
*/
double  Fix3_Tanh (  double  a )
{
    double  x;

    x = (exp(a) - exp(-a)) / 
        (exp(a) + exp(-a));


    return (x);

}  /* Fix3_Tanh */





/*****  Fix3_BFactor    *********************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*/
int     Fix3_BFactor (double  *B,            /* (O) B ceoff for each factor       */
                 long    SwapSt,        /* (I) Underlying swap start         */
                 long    SwapMat,       /* (I) Underlying swap maturity      */
                 char    DCC,           /* (I) Underlying day count conv.    */
                 char    Freq,          /* (I) Underlying frequency          */
                 int     NbFactor,      /* (I) Number of factors             */
                 double  *Beta,         /* (I) Mean reversions               */
                 double  Bbq,           /* (I) Backbone parameter            */
                 double  VolNorm,       /* (I) Normal volatility in backbone */
                 double  VolLogn,       /* (I) Lognorm volatility in backbone*/
                 int     NbZero,        /* (I) Number of zeros               */
                 double  *Zero,         /* (I) Zero rates                    */
                 long    *ZeroDate,     /* (I) Zero maturity dates           */
                 long    BaseDate)      /* (I) Zero curve base date          */
{

    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    long        CpnPmtDate;           /* Current cpn pmt date    */
    long        TempDates[2];         /* To construct cpn list   */

    double  S;                   /* Swap start in years from base date */
    double  TtoCpn;              /* Time to current coupon in years    */
    double  PrevT;               /* Same for previous coupon           */

    double  DCCFrac;             /* Current coupon day count fraction  */
    double  Annuity;             /* Annuity price                      */
    double  FwdYield;            /* Forward yield                      */

    double  ZerotoS;             /* Zero to swap start                 */
    double  ZerotoCpn=0;         /* Zero to current cpn                */
    double  PrevZero;

    double  BbqAdj;              /* Back bone adjustment in A coeff    */
    double  A[3];                /* A in Christian's memo              */
    int     j, k;
    int     status = FAILURE;    /* Error status = FAILURE initially   */




    if (SwapMat > ZeroDate[NbZero-1])
    {        
        DR_Error ("Not enough zeros to calculate B (Fix3_BFactor)!");
        goto FREE_MEM_AND_RETURN;
    }


    TempDates[0] = SwapSt;
    TempDates[1] = SwapMat;

    CpnEventList = DrNewEventListFromFreq (  2,
                                             TempDates,
                                             Freq,
                                             'F',   /* Always front stub */
                                             'N',   /* Dates in not required */
                                             NULL, NULL, NULL, NULL, NULL);

    if (CpnEventList == NULL)
    {
        goto FREE_MEM_AND_RETURN;
    }


    /* Zero to swap start */
    ZerotoS = ZeroPrice(SwapSt,         /* (I) Discount bond maturity     */
                        BaseDate,       /* (I) Value Date                 */
                        NbZero,         /* (I) Number of zeros            */
                        ZeroDate,       /* (I) Zero maturity dates        */
                        Zero);          /* (I) Zero rates                 */
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    S = Daysact (BaseDate, SwapSt) / 365.;

    PrevT = S;
    PrevZero = ZerotoS;

    Annuity = 0.;

    A[0] = A[1] = A[2] = 0.;
    B[0] = B[1] = B[2] = 0.;

    for (j = 1; j < CpnEventList->NbEntries; j++)
    {
        CpnPmtDate = CpnEventList->Dates[j];

        if (DrDayCountFraction (CpnEventList->Dates[j-1],
                                CpnEventList->Dates[j], 
                                DCC,
                                &DCCFrac) == FAILURE)
        {
            DR_Error ("Could not calculate day count fraction (Fix3_BFactor)!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Zero to current cpn */
        ZerotoCpn = ZeroPrice(CpnPmtDate,   /* (I) Discount bond maturity     */
                              BaseDate,     /* (I) Value Date                 */
                              NbZero,       /* (I) Number of zeros            */
                              ZeroDate,     /* (I) Zero maturity dates        */
                              Zero);        /* (I) Zero rates                 */
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;

        TtoCpn = Daysact (BaseDate, CpnPmtDate) / 365.;
        
        Annuity += DCCFrac * ZerotoCpn;
            
        for (k = 0; k < NbFactor; k++)
        {
            /* cf Christian's memo for Vladimir's approximation */
            BbqAdj =      (Bbq) * VolLogn * log (PrevZero/ZerotoCpn) 
                   + (1. - Bbq) * VolNorm * (TtoCpn - PrevT);

            A[k]  += exp (-Beta[k] * (PrevT-S)) * BbqAdj
                   * Fix3_ExpDecay (Beta[k], (TtoCpn-PrevT));

            B[k]  += A[k] * DCCFrac * ZerotoCpn;
        }

        PrevT    = TtoCpn;
        PrevZero = ZerotoCpn;

    }  /* for j */

    FwdYield = (ZerotoS - ZerotoCpn) / Annuity;

    for (k = 0; k < NbFactor; k++)
    {
        B[k] *= FwdYield;
        B[k] += A[k] * ZerotoCpn;
        B[k] /= (ZerotoS - ZerotoCpn);
    }


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* Fix3_BFactor */






/*****  SpotVol    *********************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     SpotVol    (
       double  Aweight[MAXNBDATE],    /* (O) Spot vol                        */
       long    VolBaseDate,           /* (I) Volatility curve base date      */
       int     NbVol,                 /* (I) Nb of points in vol curve       */
       long    *VolDate,              /* (I) Volatility dates                */
       double  *Vol,                  /* (I) Vol curve                       */
       int     *VolUsed,              /* (I) TRUE if vol used in calibration */
       char    Freq,                  /* (I) Frequency of underlying rate    */
       char    DCC,                   /* (I) Day count convention            */
       long    *SwapSt,               /* (I) Underlying swap start           */
       long    *SwapMat,              /* (I) Underlying swap maturity        */
       double  QLeft,                 /* (I) Left Q mapping coefficient      */
       double  QRight,                /* (I) Right Q mapping coefficient     */
       double  FwdSh,                 /* (I) Fwd shift mapping coefficient   */
       double  Bbq,                   /* (I) Backbone parameter              */
       double  VolNorm,               /* (I) Normal volatility in backbone   */
       double  VolLogn,               /* (I) Lognorm volatility in backbone  */
       int     NbFactor,              /* (I) Number of factors               */
       double  *Beta,                 /* (I) Mean reversions                 */   
       double  Alpha1[MAXNBDATE],     /* (I) Vol factor 1                    */
       double  Alpha2[MAXNBDATE],     /* (I) Vol factor 2                    */
       double  Rho[MAXNBDATE],        /* (I) Correlation                     */   
       int     SkipFlag,              /* (I) Skip calibration failure points */
       int     CalibFlag,             /* (I) Index calibration flag          */
       int     NbZero,                /* (I) Number of zeros                 */
       double  *Zero,                 /* (I) Zero rates                      */
       long    *ZeroDate,             /* (I) Zero maturity dates             */
       long    ZeroBaseDate)          /* (I) Zero curve base date            */
{

    double  VolT[MAXNBDATE];      /* Expiries in years                      */
    double  t;                    /* Time between two consecutive expiries  */
    double  T;                    /* Time to current expiry                 */
    double  pY[MAXNBDATE];        /* Fwd par yield                          */
    double  Annuity;              /* Fwd annuity                            */

    double  B[MAXNBDATE][2];      /* B in Christian memo                    */
    double  Anorm=0;              /* Norm of historical aweights            */ 
    double  M;                    /* Matrix for lambda system               */
    double  y;                    /* Vector for lambda system               */
    double  atmPr;                /* Market price of atm option             */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */
    

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, j, k;
    int     status = FAILURE;     /* Error status = FAILURE initially       */

    char    ErrorMsg[MAXBUFF];    /* Error message                          */
 

    /*
    *  Norm of aweight is known form volnorm,vollog 
    */
    if (IS_EQUAL(Bbq,1))
    {
        Anorm = VolLogn; 
    }
    else
    if (IS_EQUAL(Bbq,0))
    {
        Anorm = VolNorm; 
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        goto RETURN;
    }       




    if (CalibFlag == FALSE)
    {
        Aweight[0] = Alpha1[0];       
        
        return (SUCCESS);
    }



    for (i = 0; i < NbVol; i++)
    {
        /* Only used vol points */
        if (!VolUsed[i]) continue;

        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

        /* B factor */
        if (Fix3_BFactor (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        NbFactor,
                        Beta,
                        Bbq,
                        VolNorm,
                        VolLogn,
                        NbZero,
                        Zero,
                        ZeroDate,
                        ZeroBaseDate) == FAILURE)
        {
            goto RETURN;
        }

        if (Par_Yield_From_Dates (&(pY[i]),
                                  &Annuity,
                                  SwapSt[i],
                                  SwapMat[i],
                                  DCC,
                                  Freq,
                                  'F',
                                  NbZero,
                                  Zero,
                                  ZeroDate,
                                  ZeroBaseDate) == FAILURE)
        {
            goto RETURN;
        }
    }
    
    
    LastCalibIdx = -1;              /* no calibrated points so far */

    for (i = 0; i < NbVol; i++)                                                    
    {       
        if (!VolUsed[i]) continue;

        T = VolT[i];
        t = ((LastCalibIdx == -1) ? VolT[i] : (VolT[i] - VolT[LastCalibIdx]));



        /* Lognormal to X-space vol adjustment. The single and */
        /* two q cases are treated differently for consistency */
        /* with old model.                                     */
        if ((fabs(QLeft - QRight) < TINY) && (fabs(FwdSh) < TINY))
        {

            if (fabs(QLeft) > QCUTOFF)                                                 
                y = Normal_InvH (QLeft * (NormalH (.5 * sqrt(T) * Vol[i]) - .5) + .5) / (.5 * QLeft);
            else
                y = (2. * NormalH (.5 * sqrt(T) * Vol[i]) - 1.) * sqrt (2. * PI);
  
            y = SQUARE (y);
        }
        else
        {
            atmPr= Option_BS2Q (pY[i],pY[i],T,Vol[i],'C',1.,1.,0.);
            if (atmPr < 0.0)
            {
                sprintf (ErrorMsg, "Spot vol: problem in BS2Q price at %ld ",
                         VolDate[i]);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
            
            y    = ImpVol_BS2Q (pY[i],pY[i],T,atmPr,'C',QLeft,QRight,FwdSh,Vol[i]);
            if (y < 0.0)
            {
                sprintf (ErrorMsg, "Spot vol: problem in BS2Q implied vol at %ld ",
                         VolDate[i]);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
            
            y = T * SQUARE (y);
        }


        /* Upgrade integrals */
        for (j=0; j <= LastCalibIdx ; j++)
        {
            if (j == 0)
            {
                y -=       Alpha1[j] * Alpha1[j] * B[i][0] * B[i][0] * Aweight[j] * Aweight[j] 
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], 2.0 * Beta[0]);
                y -=       Alpha2[j] * Alpha2[j] * B[i][1] * B[i][1] * Aweight[j] * Aweight[j] 
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], 2.0 * Beta[1]);
                y -= 2.0 * Alpha1[j] * Alpha2[j] * B[i][0] * B[i][1] * Aweight[j] * Aweight[j] * Rho[j]
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], Beta[1] + Beta[0]);
            }
            else
            {
                y -=       Alpha1[j] * Alpha1[j] * B[i][0] * B[i][0] * Aweight[j] * Aweight[j] 
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], 2.0 * Beta[0]);
                y -=       Alpha2[j] * Alpha2[j] * B[i][1] * B[i][1] * Aweight[j] * Aweight[j] 
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], 2.0 * Beta[1]);
                y -= 2.0 * Alpha1[j] * Alpha2[j] * B[i][0] * B[i][1] * Aweight[j] * Aweight[j] * Rho[j]
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], Beta[1] + Beta[0]);
            }

        }

        M = 0.0;

        for (j=LastCalibIdx+1; j <= i ; j++)
        {
            if (j == 0)
            {
                M +=       Alpha1[j] * Alpha1[j] * B[i][0] * B[i][0]  
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], 2.0 * Beta[0]);
                M +=       Alpha2[j] * Alpha2[j] * B[i][1] * B[i][1]  
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], 2.0 * Beta[1]);
                M += 2.0 * Alpha1[j] * Alpha2[j] * B[i][0] * B[i][1] * Rho[j]  
                           * Fix3_ExpInt (0.,VolT[j],VolT[i], Beta[1] + Beta[0]);
            }
            else
            {
                M +=       Alpha1[j] * Alpha1[j] * B[i][0] * B[i][0]  
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], 2.0 * Beta[0]);
                M +=       Alpha2[j] * Alpha2[j] * B[i][1] * B[i][1]  
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], 2.0 * Beta[1]);
                M += 2.0 * Alpha1[j] * Alpha2[j] * B[i][0] * B[i][1] * Rho[j]  
                           * Fix3_ExpInt (VolT[j-1],VolT[j],VolT[i], Beta[1] + Beta[0]);
            }
        }
        
        
        /*
         *  Solve using Gauss-Jordan method
         */
        if (fabs(M) < TINY)
        {
            sprintf (ErrorMsg, "Spot vol: problem in bootstrapping %ld "
                                "volatility (negative variance)",
                     VolDate[i]);
            DR_Error (ErrorMsg);
            goto RETURN;            
        }

        y /= M;
        
        /*
         * Check if solution is positive 
         */
        if (y < TINY)
        {
            sprintf (ErrorMsg, "Spot vol: problem in bootstrapping %ld "
                               "volatility (negative value)", VolDate[i]);
            DR_Error (ErrorMsg);

            if (SkipFlag)
            {
                VolUsed[i] = FALSE;
                continue;
            }
            else
            {
                goto RETURN;
            }            
        }  /* if */


        /*
         *   Record first calibrated lambda as a initial lambda level 
         *   if on the first calibration. If initial level is set already, 
         *   test new lambda w/r to it. Note possibly small lambda.  
         */

        lambdaNew = sqrt(y);

        if (LastCalibIdx == -1)
        {
            lambdaI = lambdaNew;
            lambda  = lambdaNew;
        }
        else
        {   
            if ((lambdaNew/lambdaI < MIN_SPOT_VOL_RATIO) || 
                (lambdaNew/lambdaI > 1./MIN_SPOT_VOL_RATIO))
            {
                sprintf (ErrorMsg, "Spot vol: problem in bootstrapping %ld "
                                   "volatility (low vol)", VolDate[i]);
                DR_Error (ErrorMsg);

                if(SkipFlag)
                {
                    VolUsed[i] = FALSE;
                    continue;
                }
                else
                {
                    goto RETURN;
                }
            }
            lambda  = lambdaNew;
        }

        /*
         *  Relative weights of factors
         */
        for (k = LastCalibIdx+1; k <= i; k++)
        {
            Aweight[k] = lambda;
        }

        /* 
         * Reset LastCalibIdx to current index
         */
        LastCalibIdx = i;
       
    }  /* for i */


    for (k = LastCalibIdx+1; k < NbVol; k++)
    {
        Aweight[k] = lambda ;        
        
    }  /* for k */
    

    if (LastCalibIdx == -1)
    {
        sprintf (ErrorMsg, "Spot vol: none of points is calibrated!");
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return (status);

}  /* SpotVol */






/*****  Smd_SpotVol    *********************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Smd_SpotVol    (
       double  Aweight[6][MAXNBDATE], /* (O) Spot vol curve of each factor   */
       long    VolBaseDate,           /* (I) Volatility curve base date      */
       int     NbVol,                 /* (I) Nb of points in vol curve       */
       long    *VolDate,              /* (I) Volatility dates                */
       double  *Vol,                  /* (I) Vol curve                       */
       int     *VolUsed,              /* (I) TRUE if vol used in calibration */
       char    Freq,                  /* (I) Frequency of underlying rate    */
       char    DCC,                   /* (I) Day count convention            */
       long    *SwapSt,               /* (I) Underlying swap start           */
       long    *SwapMat,              /* (I) Underlying swap maturity        */
       double  QLeft,                 /* (I) Left Q mapping coefficient      */
       double  QRight,                /* (I) Right Q mapping coefficient     */
       double  FwdSh,                 /* (I) Fwd shift mapping coefficient   */
       double  Bbq,                   /* (I) Backbone parameter              */
       double  VolNorm,               /* (I) Normal volatility in backbone   */
       double  VolLogn,               /* (I) Lognorm volatility in backbone  */
       int     NbFactor,              /* (I) Number of factors               */
       double  *Alpha,                /* (I) Relative size factors           */
       double  *Beta,                 /* (I) Mean reversions                 */
       double  *Rho,                  /* (I) Correlation between factors     */
       double  Dfac,                  /* (I) Correlation term structure      */
       int     SkipFlag,              /* (I) Skip calibration failure points */
       int     CalibFlag,             /* (I) Index calibration flag          */
       int     NbZero,                /* (I) Number of zeros                 */
       double  *Zero,                 /* (I) Zero rates                      */
       long    *ZeroDate,             /* (I) Zero maturity dates             */
       long    ZeroBaseDate)          /* (I) Zero curve base date            */
{

    /* Local Aweight Array */
    double SmdSpotVol[MAXNBDATE];
    double Alpha1L[MAXNBDATE];
    double Alpha2L[MAXNBDATE];
    double SmdAlp1L[MAXNBDATE];
    double SmdAlp2L[MAXNBDATE];
    double SmdCorrel[MAXNBDATE];
    double RhoL[MAXNBDATE];

    double SmdVolfac;
    double expiry;
    
    int    i;
    int    status = FAILURE;


    for (i=0; i < NbVol; i++)
    {
        expiry       = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;
        SmdAlp1L[i]  = Alpha[0];
        SmdAlp2L[i]  = Alpha[1];

        /* Correlation term structure */
        SmdCorrel[i] = Rho[0] * (1.0 + Fix3_Tanh(Dfac * expiry));

        /* Correlation boundaries */
        SmdCorrel[i] = MIN(SmdCorrel[i],  0.95);
        SmdCorrel[i] = MAX(SmdCorrel[i], -0.95);
    }


    /* Get equivalent VNFM parameters */
    if (IS_EQUAL(Beta[0], Beta[1]))
    {
        SmdVolfac = Beta[1] / 0.001;
    }
    else
    {
        SmdVolfac = Beta[1] / (Beta[1] - Beta[0]);
    }

    if (SmdVolfac > 0.)
    {
        for (i=0; i < NbVol; i++)
        {
            expiry     = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;         

            Alpha1L[i] = SmdVolfac;
            Alpha2L[i] = sqrt(SQUARE(SmdAlp2L[i]/SmdAlp1L[i]) + SQUARE(SmdVolfac) - 2.0 * SmdCorrel[i]*SmdAlp2L[i]/SmdAlp1L[i]*SmdVolfac);
            RhoL[i]    = (SmdAlp2L[i]/SmdAlp1L[i] * SmdCorrel[i] - SmdVolfac) / Alpha2L[i];                     
            
            if (RhoL[i] >= 0.999 || RhoL[i] <= -0.999)
            {
                DR_Error("Smd_SpotVol: Can not booststrapp spot vol ! Correlation is out of range!");
                goto RETURN;
            }

        }
    }
    else
    {
        for (i=0; i < NbVol; i++)
        {
            expiry     = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;

            Alpha1L[i] = -SmdVolfac;
            Alpha2L[i] = sqrt(SQUARE(SmdAlp2L[i]/SmdAlp1L[i]) + SQUARE(SmdVolfac) - 2.0 * SmdCorrel[i]*SmdAlp2L[i]/SmdAlp1L[i]*SmdVolfac);
            RhoL[i]    = -(SmdAlp2L[i]/SmdAlp1L[i] * SmdCorrel[i] - SmdVolfac) / Alpha2L[i];
            
            if (RhoL[i] >= 0.999 || RhoL[i] <= -0.999)
            {
                DR_Error("Smd_SpotVol: Can not booststrapp spot vol ! Correlation is out of range!");
                goto RETURN;
            }
        }
    }


    /* Get Spot Vol */
    if (SpotVol    (SmdSpotVol,
                    VolBaseDate,  
                    NbVol,
                    VolDate, 
                    Vol,
                    VolUsed,  
                    Freq,    
                    DCC, 
                    SwapSt,    
                    SwapMat,   
                    QLeft,   
                    QRight,   
                    FwdSh,     
                    Bbq,     
                    VolNorm, 
                    VolLogn, 
                    NbFactor,  
                    Beta,    
                    Alpha1L,
                    Alpha2L,
                    RhoL,
                    SkipFlag,  
                    CalibFlag,   
                    NbZero, 
                    Zero, 
                    ZeroDate,  
                    ZeroBaseDate) == FAILURE)
    {
        goto RETURN;
    }


    /* Now populate SMD aweights */
    for (i=0; i < NbVol; i++)
    {
        expiry     = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;

        Aweight[0][i] = SmdSpotVol[i];
        Aweight[1][i] = SmdAlp2L[i]/SmdAlp1L[i] * SmdSpotVol[i] * SmdCorrel[i];
        Aweight[2][i] = SmdAlp2L[i]/SmdAlp1L[i] * SmdSpotVol[i] * sqrt(1.0 - SmdCorrel[i] * SmdCorrel[i]);
    }


    status = SUCCESS;

RETURN:

    return(status);

    
}  /* Smd_SpotVol */






/*****  Smd_Interp_SpotVol  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Smd_Interp_SpotVol (
       double  **AweightCurve,         /* (O) Factors aweight curves        */
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP)                   /* (I) Total number of time points   */
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  t, T;
    int     i, j, k;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */



    NbAweight = 3;
    


    /* Calculate expiries */
    for (j = 0; j < NbVol; j++)
    {       
        VolT[j] = Daysact (VolBaseDate, VolDate[j]) / 365.;
    }

    j = 0;                  
    t = T = 0.;        
    
    for (i = 0; i <= NbTP; i++)                                              
    {
        /* End of current time step in years from value date */
        T += Length[i];

        while ((j < NbVol - 1) && (T > VolT[j] + ERROR))
            j++;

        
        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        for (k = 0; k < NbAweight; k++)
        {
            AweightCurve[k][i] = Aweight[k][j];
        }
      
        
    }  /*for i */
    
    
    /* One extra value needed to calculate the jump size at period 0 */
    for (k = 0; k < NbAweight; k++)
    {
        AweightCurve[k][-1] = AweightCurve[k][0];
    }
 
        
    status = SUCCESS;

    return (status);

}  /* Fix3_Interp_SpotVol */






























/*****  IndexVol    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     IndexVol (
          double  *Vol,                  /* (O) Index vol curve              */
          long    SwapSt,                /* (I) Underlying swap start        */
          long    SwapMat,               /* (I) Underlying swap maturity     */
          char    Freq,                  /* (I) Frequency of underlying rate */
          char    DCC,                   /* (I) Day count convention         */
          int     NbVol,                 /* (I) Number of spot vol points    */
          long    VolBaseDate,           /* (I) Volatility curve base date   */
          long    *VolDate,              /* (I) Spot vol dates               */
          double  SpotVol[MAXNBDATE],    /* (I) Aweight curve                */
          double  QLeft,                 /* (I) Left Q mapping coefficient   */
          double  QRight,                /* (I) Right Q mapping coefficient  */
          double  FwdShift,              /* (I) Fwd shift mapping coefficient*/
          double  Bbq,                   /* (I) Backbone parameter           */
          double  VolNorm,               /* (I) Normal volatility in backbone*/
          double  VolLogn,               /* (I) Lognorm volatility in backbon*/
          int     NbFactor,              /* (I) Number of factors            */
          double  *Beta,                 /* (I) Mean reversions              */          
          double  *Alpha1,               /* (I) Factor 1 weight              */
          double  *Alpha2,               /* (I) Factor 1 weight              */
          double  *Rho,                  /* (I) Correlation                  */
          int     NbZero,                /* (I) Number of zeros              */
          double  *Zero,                 /* (I) Zero rates                   */
          long    *ZeroDate,             /* (I) Zero maturity dates          */
          long    ZeroBaseDate)          /* (I) Zero curve base date         */
{

    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  ParYield;
    double  Annuity;
    double  atmPr;
    double  B[3];
    double  T;

    double  M = 0.;                 /* Total variance                        */
    int     Bucket;
    int     i, j;

    int     status = FAILURE;  /* Error status = FAILURE initially      */
    char    ErrorMsg[MAXBUFF]; /* Error message                         */



    /* Avoid spot index expiration */
    if (SwapSt <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }
 

    /* Find B coefficients */
    if (Fix3_BFactor (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 NbFactor, 
                 Beta,     
                 Bbq,
                 VolNorm,
                 VolLogn,
                 NbZero,   
                 Zero,     
                 ZeroDate, 
                 ZeroBaseDate) == FAILURE)
    {
        goto RETURN;
    }

    /* Find par yield */


    if (Par_Yield_From_Dates (&ParYield,
                              &Annuity,
                              SwapSt,
                              SwapMat,
                              DCC,
                              Freq,
                              'F',
                              NbZero,
                              Zero,
                              ZeroDate,
                              ZeroBaseDate) == FAILURE)
    {
        goto RETURN;
    }


    /* Time to expiry in years */
    for (i=0; i < NbVol; i++)
    {
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;
    }

    T = Daysact (VolBaseDate, SwapSt) / 365.;


    /* Find lowest Bucket */
    for (i=0; i < NbVol; i++)
    {
        if (VolDate[i] >= SwapSt)
            break;
    }
    Bucket = i-1;

 


    /* Upgrade integrals */
    for (j=0; j <= Bucket; j++)
    {
        if (j == 0)
        {
            M +=       Alpha1[j] * Alpha1[j] * B[0] * B[0] * SpotVol[j] * SpotVol[j] 
                * Fix3_ExpInt (0.,VolT[j], T, 2.0 * Beta[0]);
            M +=       Alpha2[j] * Alpha2[j] * B[1] * B[1] * SpotVol[j] * SpotVol[j] 
                * Fix3_ExpInt (0.,VolT[j], T, 2.0 * Beta[1]);
            M += 2.0 * Alpha1[j] * Alpha2[j] * B[0] * B[1] * SpotVol[j] * SpotVol[j] * Rho[j]
                * Fix3_ExpInt (0.,VolT[j], T, Beta[1] + Beta[0]);
        }
        else
        {
            M +=       Alpha1[j] * Alpha1[j] * B[0] * B[0] * SpotVol[j] * SpotVol[j] 
                * Fix3_ExpInt (VolT[j-1],VolT[j],T, 2.0 * Beta[0]);
            M +=       Alpha2[j] * Alpha2[j] * B[1] * B[1] * SpotVol[j] * SpotVol[j] 
                * Fix3_ExpInt (VolT[j-1],VolT[j],T, 2.0 * Beta[1]);
            M += 2.0 * Alpha1[j] * Alpha2[j] * B[0] * B[1] * SpotVol[j] * SpotVol[j] * Rho[j]
                * Fix3_ExpInt (VolT[j-1],VolT[j],T, Beta[1] + Beta[0]);
        }
    }


    /* Account for extra period */

    if (SwapSt <= VolDate[0])
    {
        M +=       Alpha1[0] * Alpha1[0] * B[0] * B[0] * SpotVol[0] * SpotVol[0]
            * Fix3_ExpInt (0.,T, T, 2.0 * Beta[0]);
        M +=       Alpha2[0] * Alpha2[0] * B[1] * B[1] * SpotVol[0] * SpotVol[0] 
            * Fix3_ExpInt (0.,T, T, 2.0 * Beta[1]);
        M += 2.0 * Alpha1[0] * Alpha2[0] * B[0] * B[1] * SpotVol[0] * SpotVol[0] * Rho[0]
            * Fix3_ExpInt (0.,T, T, Beta[1] + Beta[0]);
    }
    else if (SwapSt > VolDate[NbVol-1])
    {
        M +=       Alpha1[NbVol-1] * Alpha1[NbVol-1] * B[0] * B[0] * SpotVol[NbVol-1] * SpotVol[NbVol-1] 
            * Fix3_ExpInt (VolT[NbVol-1], T, T, 2.0 * Beta[0]);
        M +=       Alpha2[NbVol-1] * Alpha2[NbVol-1] * B[1] * B[1] * SpotVol[NbVol-1] * SpotVol[NbVol-1]
            * Fix3_ExpInt (VolT[NbVol-1], T, T, 2.0 * Beta[1]);
        M += 2.0 * Alpha1[NbVol-1] * Alpha2[NbVol-1] * B[0] * B[1] * SpotVol[NbVol-1] * SpotVol[NbVol-1] * Rho[NbVol-1]
            * Fix3_ExpInt (VolT[NbVol-1], T, T, Beta[1] + Beta[0]);
    }
    else
    {
        Bucket ++ ;

        M +=       Alpha1[Bucket] * Alpha1[Bucket] * B[0] * B[0] * SpotVol[Bucket] * SpotVol[Bucket] 
            * Fix3_ExpInt (VolT[Bucket-1], T, T, 2.0 * Beta[0]);
        M +=       Alpha2[Bucket] * Alpha2[Bucket] * B[1] * B[1] * SpotVol[Bucket] * SpotVol[Bucket]
            * Fix3_ExpInt (VolT[Bucket-1], T, T, 2.0 * Beta[1]);
        M += 2.0 * Alpha1[Bucket] * Alpha2[Bucket] * B[0] * B[1] * SpotVol[Bucket] * SpotVol[Bucket] * Rho[Bucket]
            * Fix3_ExpInt (VolT[Bucket-1], T, T, Beta[1] + Beta[0]);
    }



    if (M < SQUARE(TINY))
    {
        sprintf (ErrorMsg, "Fix3_IndexVol: problem in integrating index volatility "
                 "at date %ld (Fix3_IndexVol)", SwapSt);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    M = sqrt (M / T);

    /* Convert from q measure back to log-normal */
    if ((fabs(QLeft - QRight) < TINY) && (fabs(FwdShift) < TINY))
    {
        M = M * sqrt (T);
        if (fabs(QLeft) > QCUTOFF)                                                 
            M = 2. * Normal_InvH ((NormalH (.5 * QLeft * M) - .5) / QLeft + .5);
        else
            M = 2. * Normal_InvH (.5 * (1. + M / sqrt(2.*PI)));

        *Vol = M / sqrt(T);
    }
    else
    {
        atmPr= Option_BS2Q (ParYield,ParYield,T,M,'C',QLeft,QRight,FwdShift);
        if (atmPr < 0.0)
        {
            sprintf (ErrorMsg, "Fix3_IndexVol: problem in BS2Q price at %ld ",
                     SwapSt);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        M    = ImpVol_BS2Q (ParYield,ParYield,T,atmPr,'C',1.,1.,0.,M);
        if (M < 0.0)
        {
            sprintf (ErrorMsg, "Fix3_IndexVol: problem in BS2Q implied vol at %ld ",
                     SwapSt);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        *Vol = M;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* IndexVol */









/*****  Fix3_IndexVol    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Smd_IndexVol (
          double  *Vol,                  /* (O) Index vol curve              */
          long    SwapSt,                /* (I) Underlying swap start        */
          long    SwapMat,               /* (I) Underlying swap maturity     */
          char    Freq,                  /* (I) Frequency of underlying rate */
          char    DCC,                   /* (I) Day count convention         */
          int     CalibFlag,             /* (I) Index calibration flag       */
          int     NbVol,                 /* (I) Number of spot vol points    */
          long    VolBaseDate,           /* (I) Volatility curve base date   */
          long    *VolDate,              /* (I) Spot vol dates               */
          double  Aweight[6][MAXNBDATE], /* (I) Aweight curve                */
          double  QLeft,                 /* (I) Left Q mapping coefficient   */
          double  QRight,                /* (I) Right Q mapping coefficient  */
          double  FwdShift,              /* (I) Fwd shift mapping coefficient*/
          double  Bbq,                   /* (I) Backbone parameter           */
          double  VolNorm,               /* (I) Normal volatility in backbone*/
          double  VolLogn,               /* (I) Lognorm volatility in backbon*/
          int     NbFactor,              /* (I) Number of factors            */
          double  *Alpha,                /* (I) Vol factors                  */
          double  *Beta,                 /* (I) Mean reversions              */
          double  *Rho,                  /* (I) Rho                          */
          double  Dfac,
          int     NbZero,                /* (I) Number of zeros              */
          double  *Zero,                 /* (I) Zero rates                   */
          long    *ZeroDate,             /* (I) Zero maturity dates          */
          long    ZeroBaseDate)          /* (I) Zero curve base date         */
{

    /* Local Aweight Array */
    double SmdSpotVol[MAXNBDATE];
    double Alpha1L[MAXNBDATE];
    double Alpha2L[MAXNBDATE];
    double SmdAlp1L[MAXNBDATE];
    double SmdAlp2L[MAXNBDATE];
    double RhoL   [MAXNBDATE];
    double SmdCorrel[MAXNBDATE];

    double SmdVolfac;

    double expiry;
    
    int    i;
    int    status = SUCCESS;


    
    for (i=0; i < NbVol; i++)
    {
        expiry       = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;
        SmdAlp1L[i]  = Alpha[0];
        SmdAlp2L[i]  = Alpha[1];

        /* Correlation term structure */
        SmdCorrel[i] = Rho[0] * (1.0 + Fix3_Tanh(Dfac * expiry));

        /* Correlation boundaries */
        SmdCorrel[i] = MIN(SmdCorrel[i],  0.95);
        SmdCorrel[i] = MAX(SmdCorrel[i], -0.95);


    }


    /* Get equivalent VNFM parameters */
    if (IS_EQUAL(Beta[0], Beta[1]))
    {
        SmdVolfac = Beta[1] / 0.001;
    }
    else
    {
        SmdVolfac = Beta[1] / (Beta[1] - Beta[0]);
    }

    if (SmdVolfac > 0.)
    {
        for (i=0; i < NbVol; i++)
        {
            expiry         = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;
            
            Alpha1L[i]     = SmdVolfac;
            Alpha2L[i]     = sqrt(SQUARE(SmdAlp2L[i]/SmdAlp1L[i]) + SQUARE(SmdVolfac) - 2.0 * SmdCorrel[i]*SmdAlp2L[i]/SmdAlp1L[i]*SmdVolfac);
            RhoL[i]        = (SmdAlp2L[i]/SmdAlp1L[i] * SmdCorrel[i] - SmdVolfac) / Alpha2L[i];
            SmdSpotVol[i]  = Aweight[0][i]; 
        }
    }
    else
    {
        for (i=0; i < NbVol; i++)
        {
            expiry     = Daysact(ZeroBaseDate, VolDate[i]) / 365.0;

            Alpha1L[i] = -SmdVolfac;
            Alpha2L[i] = sqrt(SQUARE(SmdAlp2L[i]/SmdAlp1L[i]) + SQUARE(SmdVolfac) - 2.0 * SmdCorrel[i]*SmdAlp2L[i]/SmdAlp1L[i]*SmdVolfac);
            RhoL[i]    = -(SmdAlp2L[i]/SmdAlp1L[i] * SmdCorrel[i] - SmdVolfac) / Alpha2L[i];

            SmdSpotVol[i]  = Aweight[0][i];

        }
    }

    if (IndexVol (Vol,
                  SwapSt,
                  SwapMat,               
                  Freq,   
                  DCC,                   
                  NbVol,
                  VolBaseDate,
                  VolDate, 
                  SmdSpotVol,
                  QLeft,   
                  QRight,  
                  FwdShift,
                  Bbq,  
                  VolNorm,  
                  VolLogn, 
                  NbFactor, 
                  Beta,  
                  Alpha1L,
                  Alpha2L,
                  RhoL, 
                  NbZero,   
                  Zero,   
                  ZeroDate, 
                  ZeroBaseDate) == FAILURE)
    {
        goto RETURN;
    }




    status = SUCCESS;

RETURN:

    return(status);
}






















