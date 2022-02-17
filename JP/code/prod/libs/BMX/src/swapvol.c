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
#include "bmx123head.h"


#ifndef    MIN_SPOT_VOL_RATIO
#define    MIN_SPOT_VOL_RATIO    0.2
#endif

#ifndef    SPOT_VOL_FILTER_AMOUNT
#define    SPOT_VOL_FILTER_AMOUNT    0.95
#endif



/*****  ExpDecay  ***********************************************************/
/*
*       Exponential decay function.
*/
double  ExpDecay (  double  a,
                    double  t)
{
    double  x;


    /* If a=0 or t=0, ExpDecay(a,t)=1 */
    if (fabs (a * t) < ERROR)
    {
        return (1.);
    }

    x = (1. - exp (- a * t)) / a / t;

    return (x);

}  /* ExpDecay */


/******* VFAC ******************************************************************
*
*   Utility function for integration of exp(beta * t) over [t1,t2]
*/
double  VFAC(double     t1,
             double     t2,
             double     beta)
{
    double x;

    if (fabs(beta) < TINY)
    {
        return (fabs(t2 - t1));
    }
    else
    {
        x = (exp(beta * t2) - exp(beta * t1)) / beta;
        return(x);
    }
}/*VFAC */





/*****  BFactor    *********************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*/
int     BFactor (double  *B,            /* (O) B ceoff for each factor       */
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
        DR_Error ("Not enough zeros to calculate B (BFactor)!");
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
            DR_Error ("Could not calculate day count fraction (BFactor)!");
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
                   * ExpDecay (Beta[k], (TtoCpn-PrevT));

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

}  /* BFactor */



/*****  SpotVol    *********************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     SpotVol    (
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
       double  Bbq,                   /* (I) Backbone parameter              */
       double  VolNorm,               /* (I) Normal volatility in backbone   */
       double  VolLogn,               /* (I) Lognorm volatility in backbone  */
       double  Smile[NBVOLPARS][MAXNBDATE], /* (I) Vol and smile params      */
       int     NbFactor,              /* (I) Number of factors               */
       double  *Alpha,                /* (I) Relative size factors           */
       double  *Beta,                 /* (I) Mean reversions                 */
       double  *Rho,                  /* (I) Correlation between factors     */
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

    double  B[MAXNBDATE][3];      /* B in Christian memo                    */
    double  L[3][3];              /* Integrals of lambda factors            */
    double  D[3][3];              /* aweight*B                              */
    double  aw[3][3];             /* Historical aweights                    */
    double  Anorm=0;              /* Norm of historical aweights            */ 
    double  M;                    /* Matrix for lambda system               */
    double  y;                    /* Vector for lambda system               */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, k, p, q;
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


     /*
     *  Constant Aweight numbers (determined historically).
     *  NOTE: alphas are NOT normalized, but aweights ARE
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++) 
            aw[p][q] = 0.;

    aw[0][0]  = Alpha[0] * 1.;
    aw[0][0] /= Anorm;

    if (NbFactor > 1) 
    {
        aw[1][0] = Alpha[1] * Rho[0];
        aw[1][1] = Alpha[1] * sqrt(1 - Rho[0] * Rho[0]);
        aw[1][0] /= Anorm;
        aw[1][1] /= Anorm;
    }

    if (NbFactor > 2)
    {
        aw[2][0] = Alpha[2] * Rho[1];
        aw[2][1] = Alpha[2] * (Rho[2] - Rho[0]*Rho[1])/sqrt(1 - Rho[0]*Rho[0]);
        aw[2][2] = Alpha[2] * sqrt(1.0 - Rho[0]*Rho[0] - Rho[1]*Rho[1] 
          - Rho[2]*Rho[2] + 2.*Rho[0]*Rho[1]*Rho[2]) / sqrt(1 - Rho[0]*Rho[0]);
        aw[2][0] /= Anorm;
        aw[2][1] /= Anorm;
        aw[2][2] /= Anorm;
    }


    if (CalibFlag == FALSE)
    {
        Aweight[0][0] = aw[0][0];
        
        if (NbFactor > 1)
        {
            Aweight[1][0] = aw[1][0];
            Aweight[2][0] = aw[1][1];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][0] = aw[2][0];
            Aweight[4][0] = aw[2][1];
            Aweight[5][0] = aw[2][2];
        }

        return (SUCCESS);
    }



    for (i = 0; i < NbVol; i++)
    {
        /* Only used vol points */
        if (!VolUsed[i]) continue;

        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

        /* B factor */
        if (BFactor (   B[i],
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


    /* 
     *  Bootstrapp swaption volatility.
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
        {
            L[p][q] = 0.;
            D[p][q] = 0.;
        }
    
    
    LastCalibIdx = -1;              /* no calibrated points so far */

    for (i = 0; i < NbVol; i++)                                                    
    {       
        if (!VolUsed[i]) continue;

        T = VolT[i];
        t = ((LastCalibIdx == -1) ? VolT[i] : (VolT[i] - VolT[LastCalibIdx]));
        
        D[0][0] = aw[0][0] * B[i][0];
        
        if (NbFactor > 1)
        {
            D[1][0] = aw[1][0] * B[i][1];
            D[1][1] = aw[1][1] * B[i][1];
        }
        
        if (NbFactor > 2)
        {
            D[2][0] = aw[2][0] * B[i][2];
            D[2][1] = aw[2][1] * B[i][2];
            D[2][2] = aw[2][2] * B[i][2];
        }

        /* TMX3: use Vnfm vol directly */
        y = T * SQUARE(Vol[i]);

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * ExpDecay (2. * Beta[0], t) * t;
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * Beta[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(Beta[0] + Beta[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * Beta[1] * t);
            
            M += D[1][0] * D[1][0] * ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * ExpDecay (Beta[0] + Beta[1], t) * t;
            
            M += D[1][1] * D[1][1] * ExpDecay (2. * Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(Beta[0] + Beta[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * Beta[2] * t);
            
            M += D[2][0] * D[2][0] * ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][1] * D[2][1] * ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][2] * D[2][2] * ExpDecay (2. * Beta[2], t) * t;
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
            Aweight[0][k] = lambda * aw[0][0];
            
            if (NbFactor > 1)
            {
                Aweight[1][k] = lambda * aw[1][0];
                Aweight[2][k] = lambda * aw[1][1];
            }
            
            if (NbFactor > 2)
            {
                Aweight[3][k] = lambda * aw[2][0];
                Aweight[4][k] = lambda * aw[2][1];
                Aweight[5][k] = lambda * aw[2][2];
            }
        }

        /* 
         * Reset LastCalibIdx to current index
         */
        LastCalibIdx = i;

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * Beta[0] * t);
        L[0][0] += lambda * lambda * ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += lambda * lambda * ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2.*lambda*lambda * ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += lambda * lambda * ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2.*lambda*lambda * ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2.*lambda*lambda * ExpDecay ((Beta[1]+Beta[2]), t) * t;
        }
    }  /* for i */


    for (k = LastCalibIdx+1; k < NbVol; k++)
    {
        Aweight[0][k] = lambda * aw[0][0];
        
        if (NbFactor > 1)
        {
            Aweight[1][k] = lambda * aw[1][0];
            Aweight[2][k] = lambda * aw[1][1];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][k] = lambda * aw[2][0];
            Aweight[4][k] = lambda * aw[2][1];
            Aweight[5][k] = lambda * aw[2][2];
        }
    }  /* for k */
    

    if (LastCalibIdx == -1)
    {
        sprintf (ErrorMsg, "Spot vol: none of points is calibrated!)");
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return (status);

}  /* SpotVol */



/*****  Interp_SpotVol  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Interp_SpotVol (
       double  **AweightCurve,         /* (O) Factors aweight curves        */
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       int     NbFactor,               /* (I) Number of factors             */
       double  *Beta,                  /* (I) Mean reversions               */
       double  Bbq,                    /* (I) Backbone parameter            */
       double  VolNorm,                /* (I) Normal volatility in backbone */
       double  VolLogn,                /* (I) Lognorm volatility in backbone*/
       int     CalibFlag,              /* (I) Index calibration flag        */
       double  *FwdRate,               /* (I) One period forward rate       */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP)                   /* (I) Total number of time points   */
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  PrevVolT;              /* Time to previous vol point */
    double  Cov[3][3];             /* Covariance matrix          */
    double  x, t, T;
    double  BbqAdj;               /* Backbone vol adjustment    */

    int     i, j, k, l, p, q;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

    char    ErrorMsg[MAXBUFF];     /* Error message              */



    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;


    if (CalibFlag == FALSE)
    {
        for (i = -1; i <= NbTP; i++)                                              
        {
            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][0];
            }
        }

        for (i = 0; i < NbTP; i++)                                           
        {
            for (k = 0; k < NbAweight; k++)
            {
                l = (k > 0) + (k > 2);

                BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                         (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
                AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*ExpDecay(Beta[l], Length[i+1]);
            }
        }
        
        return (SUCCESS);
    }


    for (j = 0; j < NbVol; j++)
    {       
        VolT[j] = Daysact (VolBaseDate, VolDate[j]) / 365.;
    }

    /* 
    *   Interpolate spot volatility curves at each node.
    */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
            Cov[p][q] = 0.;

    j = 0;                  
    t = T = 0.;        
    
    for (i = 0; i <= NbTP; i++)                                              
    {
        /* End of current time step in years from value date */
        T += Length[i];

        while ((j < NbVol - 1) && (T > VolT[j] + ERROR))
            j++;

        PrevVolT = ((j == 0)? 0. : VolT[j-1]);

        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        if (t >= PrevVolT)
        {
            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][j];
            }
        }
        /* The time step is straddling two buckets: */
        /* we recalculate the covariance matrix.    */
        else
        {
            Cov[0][0]  = Aweight[0][j-1] * Aweight[0][j-1] * ExpDecay (2. * Beta[0], VolT[j-1] - t) * (VolT[j-1] - t);
            Cov[0][0] += Aweight[0][j]   * Aweight[0][j]   * ExpDecay (2. * Beta[0], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[0] * (VolT[j-1] - t));
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = Aweight[0][j-1] * Aweight[1][j-1] * ExpDecay (Beta[0] + Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][0] += Aweight[0][j]   * Aweight[1][j]   * ExpDecay (Beta[0] + Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[1])*(VolT[j-1]-t));
            
                Cov[1][1]  = Aweight[1][j-1] * Aweight[1][j-1] * ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[1][j]   * Aweight[1][j]   * ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
                Cov[1][1] += Aweight[2][j-1] * Aweight[2][j-1] * ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[2][j]   * Aweight[2][j]   * ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = Aweight[0][j-1] * Aweight[3][j-1] * ExpDecay (Beta[0] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][0] += Aweight[0][j]   * Aweight[3][j]   * ExpDecay (Beta[0] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[2])*(VolT[j-1]-t));
            
                Cov[2][1]  = Aweight[1][j-1] * Aweight[3][j-1] * ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[1][j]   * Aweight[3][j]   * ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));
                Cov[2][1] += Aweight[2][j-1] * Aweight[4][j-1] * ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[2][j]   * Aweight[4][j]   * ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));

                Cov[2][2]  = Aweight[3][j-1] * Aweight[3][j-1] * ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[3][j]   * Aweight[3][j]   * ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[4][j-1] * Aweight[4][j-1] * ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[4][j]   * Aweight[4][j]   * ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[5][j-1] * Aweight[5][j-1] * ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[5][j]   * Aweight[5][j]   * ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
            }


            x = Cov[0][0] / ExpDecay (2. * Beta[0], T - t) / (T - t);

            if (x < TINY)
            {
                sprintf (ErrorMsg, "Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                DR_Error (ErrorMsg);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / ExpDecay (Beta[0] + Beta[1], T - t) / (T - t);

                x = Cov[1][1] / ExpDecay (2. * Beta[1], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

                if (x < TINY)
                {
                    sprintf (ErrorMsg, "Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                    DR_Error (ErrorMsg);
                    goto RETURN;
                }
                        
                AweightCurve[2][i] = sqrt (x);
            }

            if (NbFactor > 2)
            {
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / ExpDecay (Beta[0] + Beta[2], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / ExpDecay (Beta[1] + Beta[2], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / ExpDecay (2. * Beta[2], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

                if (x < TINY)
                {
                    sprintf (ErrorMsg, "Interp_SpotVol: problem in "
                              "interpolating spot volatility at node #%d!", i);
                    DR_Error (ErrorMsg);
                    goto RETURN;
                }
                        
                AweightCurve[5][i] = sqrt (x);
            }
        }  /* if then else */           
    
        t += Length[i];

    }  /*for i */
    
    
    /* One extra value needed to calculate the jump size at period 0 */
    for (k = 0; k < NbAweight; k++)
    {
        AweightCurve[k][-1] = AweightCurve[k][0];
    }


    /* 
    *   Adjustment to the spot volatility as FwdRate is not an instantaneous 
    *   rate but a Length[] discrete rate.
    */
    for (i = 0; i < NbTP; i++)                                           
    {
        /* Note period i volatility corresponds to period i+1 rate */
        for (k = 0; k < NbAweight; k++)
        {
            l = (k > 0) + (k > 2);

            BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                     (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*ExpDecay(Beta[l], Length[i+1]);
        }
    }

        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Interp_SpotVol */



/*****  IndexVol    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     IndexVol (
          double  *Vol,                  /* (I/O) Index vol curve            */
          long    SwapSt,                /* (I) Underlying swap start        */
          long    SwapMat,               /* (I) Underlying swap maturity     */
          char    Freq,                  /* (I) Frequency of underlying rate */
          char    DCC,                   /* (I) Day count convention         */
          int     CalibFlag,             /* (I) Index calibration flag       */
          int     NbVol,                 /* (I) Number of spot vol points    */
          long    VolBaseDate,           /* (I) Volatility curve base date   */
          long    *VolDate,              /* (I) Spot vol dates               */
          double  Aweight[6][MAXNBDATE], /* (I) Aweight curve                */
          double  Bbq,                   /* (I) Backbone parameter           */
          double  VolNorm,               /* (I) Normal volatility in backbone*/
          double  VolLogn,               /* (I) Lognorm volatility in backbon*/
          int     NbFactor,              /* (I) Number of factors            */
          double  *Beta,                 /* (I) Mean reversions              */
          int     NbZero,                /* (I) Number of zeros              */
          double  *Zero,                 /* (I) Zero rates                   */
          long    *ZeroDate,             /* (I) Zero maturity dates          */
          long    ZeroBaseDate)          /* (I) Zero curve base date         */
{

    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  t;                 /* Time between two consecutive expiries */
    double  T=0.;              /* Time to current expiry                */
    double  ParYield;
    double  Annuity;
    double  atmPr;

    double  L[3][3];           /* Integrals of factor spot vol          */
    double  B[3];              /* B in Christian memo                   */
    double  D[3][3];           /* aweight*B                             */
    double  M;                 /* Total variance                        */

    long    EndDate;           /* End of current integration bucket     */
    int     i, j, p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially      */
    char    ErrorMsg[MAXBUFF]; /* Error message                         */
 

    /* Find B coefficients */
    if (BFactor (B,
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

    /* Constant spot volatility case */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
            D[p][q] = 0.;

    if (CalibFlag == FALSE)
    {
        /* Avoid spot index expiration; note we use base date of zero */
        /* curve as there is no volatility information available.     */
        if (SwapSt <= ZeroBaseDate) 
        {
            *Vol = 0.00;
            return(SUCCESS);
        }
    
        t = Daysact (ZeroBaseDate, SwapSt) / 365.;

        /* Find D coefficients */
        D[0][0] = Aweight[0][0] * B[0];
    
        if (NbFactor > 1)
        {
            D[1][0] = Aweight[1][0] * B[1];
            D[1][1] = Aweight[2][0] * B[1];
        }
    
        if (NbFactor > 2)
        {
            D[2][0] = Aweight[3][0] * B[2];
            D[2][1] = Aweight[4][0] * B[2];
            D[2][2] = Aweight[5][0] * B[2];
        }
    
        M = D[0][0] * D[0][0] * ExpDecay (2. * Beta[0], t) * t;
    
        if (NbFactor > 1) 
        {
            M += D[1][0] * D[1][0] * ExpDecay (2. * Beta[1], t) * t;
            M += D[1][1] * D[1][1] * ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * ExpDecay (Beta[0] + Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            M += D[2][0] * D[2][0] * ExpDecay (2. * Beta[2], t) * t;
            M += D[2][1] * D[2][1] * ExpDecay (2. * Beta[2], t) * t;
            M += D[2][2] * D[2][2] * ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * ExpDecay (Beta[1] + Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * ExpDecay (Beta[1] + Beta[2], t) * t;
        }

        
        if (M < SQUARE(TINY))
        {
            sprintf (ErrorMsg, "IndexVol: problem in integrating index volatility "
                     "at date %ld (IndexVol)", SwapSt);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        
        /* TMX3: use Vnfm directly */
        *Vol = sqrt (M / t);

        return (SUCCESS);

    }  /* if */


    /* Avoid spot index expiration */
    if (SwapSt <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }
    

    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++) 
            L[i][j] = 0.;
    
    
    for (i = 0; i < NbVol; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate = (i < NbVol-1) ? MIN(SwapSt, VolDate[i]) : SwapSt;


        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;

        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));
        

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * Beta[0] * t);
        L[0][0] += Aweight[0][i] * Aweight[0][i] * ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += Aweight[1][i] * Aweight[1][i] * ExpDecay (2. * Beta[1], t) * t;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[4][i] * Aweight[4][i] * ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[5][i] * Aweight[5][i] * ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * ExpDecay ((Beta[1]+Beta[2]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == SwapSt)
        {
            break;
        }
    }  /* for i */


    
    M = B[0] * B[0] * L[0][0];
    
    if (NbFactor > 1) 
    {
        M += B[1] * B[1] * L[1][1];
        M += B[0] * B[1] * L[0][1];
    }
    
    if (NbFactor > 2)
    {
        M += B[2] * B[2] * L[2][2];
        M += B[0] * B[2] * L[0][2];
        M += B[1] * B[2] * L[1][2];
    }
    

    if (M < SQUARE(TINY))
    {
        sprintf (ErrorMsg, "IndexVol: problem in integrating index volatility "
                 "at date %ld (IndexVol)", SwapSt);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* TMX3: use Vnfm directly */
    *Vol = sqrt (M / T);

    status = SUCCESS;

    RETURN:

    return (status);

}  /* IndexVol */





/* NEW FUNCTIONS ADDED BY BM */


/*****  TreeVol    *********************************************************/
/*
*       Build the Aweights.
*/

int  TreeVol(MKTVOL_DATA     *mktvol_data0,         /* (I/0) mkt vol data structure  */
             MKTVOL_DATA     *mktvol_data1,         /* (I/0) mkt vol data structure  */             
             TREE_DATA       *tree_data,            /* (I)   tree data structure     */
             T_CURVE         *t_curve)              /* (I)   t_curve stucture data   */

{
    /* Variables */
    int          i;
    long         t;
    double       *SpotVol1, *SpotVol2;  
    int          status = FAILURE;


    /* Tree weights calculation variables */
    double       var1       = 0.,
                 var2       = 0.,
                 cov        = 0.,
                 corrL      = 0.,
                 biv_correl = 0.;

    /* Previous and current time for tree weights calculation */ 
    double       PrevT, CurrT;

    /* Mean-reversions */
    double Beta1            = mktvol_data0->Beta[0];
    double Beta2            = mktvol_data0->Beta[1];



    
    /* Initialize memory */
    SpotVol1 = (double*) DR_Array(DOUBLE, -1, tree_data->NbTP+1);
    SpotVol2 = (double*) DR_Array(DOUBLE, -1, tree_data->NbTP+1);

    if (SpotVol1 == NULL  ||
        SpotVol2 == NULL)
    {
        DR_Error("TreeVol : Can't allocate memory!");
        return(status);
    }


    /* 1 FACTOR CASE */
    if (tree_data->NbFactor == 1 )   
    {
        /* Bootstrap spot vol */
        if (SpotVol (   mktvol_data0->Aweight,
                        mktvol_data0->BaseDate,
                        mktvol_data0->NbVol,
                        mktvol_data0->VolDate,
                        mktvol_data0->Vol[0],
                        mktvol_data0->VolUsed,
                        mktvol_data0->Freq,
                        mktvol_data0->DCC,
                        mktvol_data0->SwapSt,
                        mktvol_data0->SwapMat,
                        mktvol_data0->Bbq,
                        mktvol_data0->VolNorm,
                        mktvol_data0->VolLogn,
                        mktvol_data0->Vol,
                        1,                              /* 1 Factor */
                        mktvol_data0->Alpha,
                        mktvol_data0->Beta,
                        mktvol_data0->Rho,
                        mktvol_data0->SkipFlag,
                        mktvol_data0->CalibFlag,
                        (t_curve[tree_data->CvDiff]).NbZero,
                        (t_curve[tree_data->CvDiff]).Zero,
                        (t_curve[tree_data->CvDiff]).ZeroDate,
                        (t_curve[tree_data->CvDiff]).ValueDate) == FAILURE)
        {
            goto RETURN;
        }

        /* Interpolate spot volatility curves */
        if (Interp_SpotVol (tree_data->Aweight,
                            mktvol_data0->BaseDate,
                            mktvol_data0->NbVol,
                            mktvol_data0->VolDate,
                            mktvol_data0->Aweight,
                            1,                             /* 1 factor */  
                            mktvol_data0->Beta,
                            mktvol_data0->Bbq,
                            mktvol_data0->VolNorm,
                            mktvol_data0->VolLogn,
                            mktvol_data0->CalibFlag,
                            tree_data->FwdRate[tree_data->CvDiff],
                            tree_data->Length,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }

    }

    /* 2 FACTOR CASE */
    else if (tree_data->NbFactor == 2)        /* Calibrate two sets of spot vols */
    {

        /* We now calibrate the two mkt vol data structures with aq one factor model 
           each. The vols are then orthogonalized to fit the correlation target       */


        /* Read bivariate correlation first */
        if (Correl_Input_W(tree_data, mktvol_data0) == FAILURE)
            goto RETURN;

        /* FIRST INDEX RATE BOOSTRAPPING */
        /* Calculate the vols for the first factor */
        if (SpotVol (   mktvol_data0->Aweight,
                        mktvol_data0->BaseDate,
                        mktvol_data0->NbVol,
                        mktvol_data0->VolDate,
                        mktvol_data0->Vol[0],
                        mktvol_data0->VolUsed,
                        mktvol_data0->Freq,
                        mktvol_data0->DCC,
                        mktvol_data0->SwapSt,
                        mktvol_data0->SwapMat,
                        mktvol_data0->Bbq,
                        mktvol_data0->VolNorm,
                        mktvol_data0->VolLogn,
                        mktvol_data0->Vol,
                        1,                      /* 1 Factor always */
                        mktvol_data0->Alpha,
                        mktvol_data0->Beta,
                        mktvol_data0->Rho,
                        mktvol_data0->SkipFlag,
                        mktvol_data0->CalibFlag,
                        (t_curve[tree_data->CvDiff]).NbZero,
                        (t_curve[tree_data->CvDiff]).Zero,
                        (t_curve[tree_data->CvDiff]).ZeroDate,
                        (t_curve[tree_data->CvDiff]).ValueDate) == FAILURE)
        {
            goto RETURN;
        }
        
        /* Interpolate spot volatility curves */
        if (Interp_SpotVol (&SpotVol1,
                            mktvol_data0->BaseDate,
                            mktvol_data0->NbVol,
                            mktvol_data0->VolDate,
                            mktvol_data0->Aweight,
                            1,                    /* 1 factor always */
                            mktvol_data0->Beta,
                            mktvol_data0->Bbq,
                            mktvol_data0->VolNorm,
                            mktvol_data0->VolLogn,
                            mktvol_data0->CalibFlag,
                            tree_data->FwdRate[tree_data->CvDiff],
                            tree_data->Length,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }


        /* SECOND INDEX RATE BOOSTRAPPING */
        mktvol_data1->Beta[0]   = mktvol_data1->Beta[1];
        mktvol_data1->Alpha[0]  = mktvol_data1->Alpha[1];
        


        /* Calculate the vols for the first factor */
        if (SpotVol (   mktvol_data1->Aweight,
                        mktvol_data1->BaseDate,
                        mktvol_data1->NbVol,
                        mktvol_data1->VolDate,
                        mktvol_data1->Vol[0],
                        mktvol_data1->VolUsed,
                        mktvol_data1->Freq,
                        mktvol_data1->DCC,
                        mktvol_data1->SwapSt,
                        mktvol_data1->SwapMat,
                        mktvol_data1->Bbq,
                        mktvol_data1->VolNorm,
                        mktvol_data1->VolLogn,
                        mktvol_data1->Vol,
                        1,                      /* 1 Factor always */
                        mktvol_data1->Alpha,
                        mktvol_data1->Beta,
                        mktvol_data1->Rho,
                        mktvol_data1->SkipFlag,
                        mktvol_data1->CalibFlag,
                        (t_curve[tree_data->CvDiff]).NbZero,
                        (t_curve[tree_data->CvDiff]).Zero,
                        (t_curve[tree_data->CvDiff]).ZeroDate,
                        (t_curve[tree_data->CvDiff]).ValueDate) == FAILURE)
        {
            goto RETURN;
        }
        
        /* Interpolate spot volatility curves */
        if (Interp_SpotVol (&SpotVol2,
                            mktvol_data1->BaseDate,
                            mktvol_data1->NbVol,
                            mktvol_data1->VolDate,
                            mktvol_data1->Aweight,
                            1,                    /* 1 factor always */
                            mktvol_data1->Beta,
                            mktvol_data1->Bbq,
                            mktvol_data1->VolNorm,
                            mktvol_data1->VolLogn,
                            mktvol_data1->CalibFlag,
                            tree_data->FwdRate[tree_data->CvDiff],
                            tree_data->Length,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }



        /* NOW ORTHOGONALIZE THE VOLS AND CALCULATE THE TREE
           WEIGHTS TO FIT THE TARGET CORRELATION BETWEEN DRIVERS */

        /* Initial bivariate correlation */
        biv_correl = tree_data->Biv_Correl1[0];

        /* Populate tree weights at time -1 */
        tree_data->Aweight[0][-1] = SpotVol1[-1];
        tree_data->Aweight[1][-1] = biv_correl * SpotVol2[-1];
        tree_data->Aweight[2][-1] = sqrt(1 - biv_correl * biv_correl) * SpotVol2[-1];

        for (t=1; t <= tree_data->NbTP; t++)
        {
            /* Previous and current time step */
            PrevT = Daysact(tree_data->TPDate[0], tree_data->TPDate[t-1]) / 365.0;
            CurrT = Daysact(tree_data->TPDate[0], tree_data->TPDate[t])   / 365.0;


            /* First gaussian driver variance at time t */
            var1 += SQUARE(SpotVol1[t-1]) * VFAC(PrevT, CurrT, 2.0 * Beta1);

            
            /* Second gaussian driver variance at time t */
            var2 += SQUARE(SpotVol2[t-1]) * VFAC(PrevT, CurrT, 2.0 * Beta2);

            /* Interpolated bivariate correlation */
            tableinterp (CurrT,&biv_correl, tree_data->Biv_Expiry, tree_data->Biv_Correl1, tree_data->NbCorrel); 

            /* Now calculate the tree weights correlation */
            corrL = (biv_correl * sqrt(var1 * var2) - cov) /
                    (VFAC(PrevT, CurrT, Beta1 + Beta2) * SpotVol1[t-1] * SpotVol2[t-1]);

            /* Bounds for corrL */
            corrL = MIN(0.999, corrL);
            corrL = MAX(0.001, corrL);

        
            /* Populate tree weights */
            tree_data->Aweight[0][t-1] = SpotVol1[t-1];
            tree_data->Aweight[1][t-1] = corrL * SpotVol2[t-1];
            tree_data->Aweight[2][t-1] = sqrt(1 - corrL * corrL) * SpotVol2[t-1];


            /* Update covariance */
            cov   += tree_data->Aweight[0][t-1] * tree_data->Aweight[1][t-1] * VFAC(PrevT, CurrT, Beta1 + Beta2);

            corrL = cov / sqrt(var1 * var2);
        }

        /* Populate tree weights at terminal points*/
        tree_data->Aweight[0][tree_data->NbTP] = SpotVol1[tree_data->NbTP];
        tree_data->Aweight[1][tree_data->NbTP] = corrL * SpotVol2[tree_data->NbTP];
        tree_data->Aweight[2][tree_data->NbTP] = sqrt(1 - corrL * corrL) * SpotVol2[tree_data->NbTP];
    }


    /* Now Initialize all the directions for the second yield mapping */
    tree_data->Map_dir1[0] = 0.0;
    tree_data->Map_dir1[1] = 1.0;

    tree_data->Map_dir2[0] = 0.0;
    tree_data->Map_dir2[1] = 1.0;

    /* Initialize the numeraire mapping direction */
    tree_data->NmrMap_dir[0] = NUM_MAP_DIR1;
    tree_data->NmrMap_dir[1] = NUM_MAP_DIR2;

                   

    status = SUCCESS;

RETURN:

    Free_DR_Array(SpotVol1, DOUBLE, -1, tree_data->NbTP + 1);
    Free_DR_Array(SpotVol2, DOUBLE, -1, tree_data->NbTP + 1);



    return(status);
}





/****** Correl_Input_W  *******************************************
 /*
  * Read the correlation inputs
  */

int     Correl_Input_W(TREE_DATA     *tree_data, 
                       MKTVOL_DATA   *mvd      )
{
    int     i; 
    int     readerror;          /* Reading error status */
    char    ErrorMsg[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;


    /* Open the correlation spec file */
    stream = fopen ("correl.dat", "r");
    if (stream == NULL)
    {
        tree_data->NbCorrel = NB_CORREL_MAX;
        
        for (i=0; i < NB_CORREL_MAX; i++)
        {
            tree_data->Biv_Expiry[i]   = i;
            tree_data->Biv_Correl1[i]  = mvd->Rho[0];
            tree_data->Biv_Correl2[i]  = mvd->Rho[0];
        }

        status = SUCCESS;

        goto RETURN;
    }

    /* Read the number of correlations */
    if (FindAndSkipComLine (stream, "Number of correlations", "(Correl_Input_W)", "correl.dat") == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%ld \n", &(tree_data->NbCorrel));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of correlations in %s! (Mode)", "correl.dat");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Check the number of correlations is not too high */
    if (tree_data->NbCorrel > NB_CORREL_MAX)
    {
        DR_Error("The number of correlation exceeds the limit!");
        goto RETURN;
    }


    /* Read the first and second maturity and the target correlation */
    if (FindAndSkipComLine (stream, "Correlations", "Correl_Input_W", "correl.dat") == FAILURE)
    {
        goto RETURN;
    }

    for (i=0; i < tree_data->NbCorrel; i++)
    {
        readerror = fscanf (stream, "%lf %lf %lf \n",
                           &(tree_data->Biv_Expiry[i]),
                           &(tree_data->Biv_Correl1[i]),
                           &(tree_data->Biv_Correl2[i]));
        if (readerror != 3)
        {
            sprintf (ErrorMsg, "Could not read correlation number %s! (Correl_Input_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    /* Check the correlations  */
    for (i=0; i < tree_data->NbCorrel; i++)
    {
        if (fabs(tree_data->Biv_Correl1[i]) > 1.0 || fabs(tree_data->Biv_Correl2[i]) > 1.0)
        {
            sprintf (ErrorMsg, " The correlation number %s has to be between -1.0 and 1.0! (Correl_Input_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }



    
    status = SUCCESS;
        
RETURN:
    
    if (stream != NULL)
    {
        fclose (stream);
    }
    return (status);

}  /*Correl_Input_W*/

        











