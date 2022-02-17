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
#include "cupslib.h"





/*****  Hyb3_ExpDecay  ***********************************************************/
/*
*       Exponential decay function.
*/
double  Hyb3_ExpDecay (  double  a,
                    double  t)
{
    double  x;


    /* If a=0 or t=0, Hyb3_ExpDecay(a,t)=1 */
    if (fabs (a * t) < ERROR)
    {
        return (1.);
    }

    x = (1. - exp (- a * t)) / a / t;

    return (x);

}  /* Hyb3_ExpDecay */



/*****  Hyb3_BFactor    *********************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*/
int     Hyb3_BFactor (double  *B,            /* (O) B ceoff for each factor       */
                 long    SwapSt,        /* (I) Underlying swap start         */
                 long    SwapMat,       /* (I) Underlying swap maturity      */
                 char    DCC,           /* (I) Underlying day count conv.    */
                 char    Freq,          /* (I) Underlying frequency          */
                 int     NbFactor,      /* (I) Number of factors             */
                 double  *Beta,         /* (I) Mean reversions               */
                 T_CURVE const* crv)
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

    double  A[3];                /* A in Christian's memo              */
    int     j, k;
    int     status = FAILURE;    /* Error status = FAILURE initially   */




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

    ZerotoS = PrevZero = GetZeroPrice(SwapSt, crv);
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    S = Daysact (crv->ValueDate, SwapSt) / 365.;

    PrevT = S;
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
            DR_Error ("Could not calculate day count fraction (Hyb3_BFactor)!");
            goto FREE_MEM_AND_RETURN;
        }

        /* Zero to current cpn */
        ZerotoCpn = GetZeroPrice(CpnPmtDate, crv);
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;

        TtoCpn = Daysact (crv->ValueDate, CpnPmtDate) / 365.;
        
        Annuity += DCCFrac * ZerotoCpn;
            
        for (k = 0; k < NbFactor; k++)
        {
            /* cf Christian's memo for Vladimir's approximation */
            A[k] += exp (-Beta[k] * (PrevT-S)) * log (PrevZero/ZerotoCpn)
                     * Hyb3_ExpDecay (Beta[k], (TtoCpn-PrevT));

            B[k] += A[k] * DCCFrac * ZerotoCpn;
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

}  /* Hyb3_BFactor */



/*****  Hyb3_SpotVol    *********************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Hyb3_SpotVol    (
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
       int     NbFactor,              /* (I) Number of factors               */
       double  *Alpha,                /* (I) Relative size factors           */
       double  *Beta,                 /* (I) Mean reversions                 */
       double  *Rho,                  /* (I) Correlation between factors     */
       int     SkipFlag,              /* (I) Skip calibration failure points */
       int     CalibFlag,             /* (I) Index calibration flag          */
       T_CURVE const* crv)            
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
    double  M;                    /* Matrix for lambda system               */
    double  y;                    /* Vector for lambda system               */
    double  atmPr;                /* Market price of atm option             */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, k, p, q;
    int     status = FAILURE;     /* Error status = FAILURE initially       */

     /*
     *  Constant Aweight numbers (determined historically).
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++) 
            aw[p][q] = 0.;

    aw[0][0] = Alpha[0] * 1.;

    if (NbFactor > 1) 
    {
        aw[1][0] = Alpha[1] * Rho[0];
        aw[1][1] = Alpha[1] * sqrt(1 - Rho[0] * Rho[0]);
    }

    if (NbFactor > 2)
    {
        aw[2][0] = Alpha[2] * Rho[1];
        aw[2][1] = Alpha[2] * (Rho[2] - Rho[0]*Rho[1])/sqrt(1 - Rho[0]*Rho[0]);
        aw[2][2] = Alpha[2] * sqrt(1.0 - Rho[0]*Rho[0] - Rho[1]*Rho[1] 
          - Rho[2]*Rho[2] + 2.*Rho[0]*Rho[1]*Rho[2]) / sqrt(1 - Rho[0]*Rho[0]);
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
        if (Hyb3_BFactor (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        NbFactor,
                        Beta,
                        crv) != SUCCESS)
        {
            goto RETURN;
        }

        if (ParYieldFromDates (&(pY[i]),
                               &Annuity,
                               SwapSt[i],
                               SwapMat[i],
                               DCC,
                               Freq,
                               'F',
                               crv) != SUCCESS)
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
            y = ImpVol_BS2Q (pY[i],pY[i],T,atmPr,'C',
                             QLeft,QRight,FwdSh,Vol[i]); 
	    if (y < 0.0)
            {
                DR_Error("Unable to convert vol to 2Q measure.\n");
                goto RETURN;
            }
            y = T * SQUARE (y);
        }

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * Hyb3_ExpDecay (2. * Beta[0], t) * t;
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * Beta[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(Beta[0] + Beta[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * Beta[1] * t);
            
            M += D[1][0] * D[1][0] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * Hyb3_ExpDecay (Beta[0] + Beta[1], t) * t;
            
            M += D[1][1] * D[1][1] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(Beta[0] + Beta[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * Beta[2] * t);
            
            M += D[2][0] * D[2][0] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * Hyb3_ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * Hyb3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][1] * D[2][1] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * Hyb3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][2] * D[2][2] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
        }
        
        /*
         *  Solve using Gauss-Jordan method
         */
        if (fabs(M) < TINY)
        {
            DR_Error ("Spot vol: problem in bootstrapping %ld "
                                "volatility (negative variance)",
                     YMDDateFromIRDate(VolDate[i]));
            goto RETURN;            
        }

        y /= M;
        
        /*
         * Check if solution is positive 
         */
        if (y < TINY)
        {
            DR_Error ("Spot vol: problem in bootstrapping %ld "
                               "volatility (negative value)", YMDDateFromIRDate(VolDate[i]));

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
            if ((lambdaNew/lambdaI < 0.2) || (lambdaNew/lambdaI > 5.0))
            {
                DR_Error ("Spot vol: problem in bootstrapping %ld "
                                   "volatility (low vol)", YMDDateFromIRDate(VolDate[i]));

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
        L[0][0] += lambda * lambda * Hyb3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += lambda * lambda * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2.*lambda*lambda * Hyb3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += lambda * lambda * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2.*lambda*lambda * Hyb3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2.*lambda*lambda * Hyb3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
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
        DR_Error ("Spot vol: none of points is calibrated!)");
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_SpotVol */



/*****  Hyb3_Interp_SpotVol  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Hyb3_Interp_SpotVol (
       double  **AweightCurve,         /* (O) Factors aweight curves        */
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       int     NbFactor,               /* (I) Number of factors             */
       double  *Beta,                  /* (I) Mean reversions               */
       int     CalibFlag,              /* (I) Index calibration flag        */
       double  *FwdRate,               /* (I) One period forward rate       */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP)                   /* (I) Total number of time points   */
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  PrevVolT;              /* Time to previous vol point */
    double  Cov[3][3];             /* Covariance matrix          */
    double  x, t, T;

    int     i, j, k, l, p, q;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

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

                AweightCurve[k][i] *= (1.+FwdRate[i+1]) 
                                    * log (1.+FwdRate[i+1]) / FwdRate[i+1]
                                     * Hyb3_ExpDecay (Beta[l], Length[i+1]);
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
            Cov[0][0]  = Aweight[0][j-1] * Aweight[0][j-1] * Hyb3_ExpDecay (2. * Beta[0], VolT[j-1] - t) * (VolT[j-1] - t);
            Cov[0][0] += Aweight[0][j]   * Aweight[0][j]   * Hyb3_ExpDecay (2. * Beta[0], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[0] * (VolT[j-1] - t));
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = Aweight[0][j-1] * Aweight[1][j-1] * Hyb3_ExpDecay (Beta[0] + Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][0] += Aweight[0][j]   * Aweight[1][j]   * Hyb3_ExpDecay (Beta[0] + Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[1])*(VolT[j-1]-t));
            
                Cov[1][1]  = Aweight[1][j-1] * Aweight[1][j-1] * Hyb3_ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[1][j]   * Aweight[1][j]   * Hyb3_ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
                Cov[1][1] += Aweight[2][j-1] * Aweight[2][j-1] * Hyb3_ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[2][j]   * Aweight[2][j]   * Hyb3_ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = Aweight[0][j-1] * Aweight[3][j-1] * Hyb3_ExpDecay (Beta[0] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][0] += Aweight[0][j]   * Aweight[3][j]   * Hyb3_ExpDecay (Beta[0] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[2])*(VolT[j-1]-t));
            
                Cov[2][1]  = Aweight[1][j-1] * Aweight[3][j-1] * Hyb3_ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[1][j]   * Aweight[3][j]   * Hyb3_ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));
                Cov[2][1] += Aweight[2][j-1] * Aweight[4][j-1] * Hyb3_ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[2][j]   * Aweight[4][j]   * Hyb3_ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));

                Cov[2][2]  = Aweight[3][j-1] * Aweight[3][j-1] * Hyb3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[3][j]   * Aweight[3][j]   * Hyb3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[4][j-1] * Aweight[4][j-1] * Hyb3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[4][j]   * Aweight[4][j]   * Hyb3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[5][j-1] * Aweight[5][j-1] * Hyb3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[5][j]   * Aweight[5][j]   * Hyb3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
            }


            x = Cov[0][0] / Hyb3_ExpDecay (2. * Beta[0], T - t) / (T - t);

            if (x < TINY)
            {
                DR_Error ("Hyb3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / Hyb3_ExpDecay (Beta[0] + Beta[1], T - t) / (T - t);

                x = Cov[1][1] / Hyb3_ExpDecay (2. * Beta[1], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

                if (x < TINY)
                {
                    DR_Error ("Hyb3_Interp_SpotVol: problem in interpolating spot volatility at node #%d(second factor)!", i);
                    goto RETURN;
                }
                        
                AweightCurve[2][i] = sqrt (x);
            }

            if (NbFactor > 2)
            {
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / Hyb3_ExpDecay (Beta[0] + Beta[2], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / Hyb3_ExpDecay (Beta[1] + Beta[2], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / Hyb3_ExpDecay (2. * Beta[2], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

                if (x < TINY)
                {
                    DR_Error ("Hyb3_Interp_SpotVol: problem in "
                              "interpolating spot volatility at node #%d (third factor)!", i);
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

            AweightCurve[k][i] *= (1.+FwdRate[i+1]) * log (1.+FwdRate[i+1])
                          / FwdRate[i+1] * Hyb3_ExpDecay (Beta[l], Length[i+1]);
        }
    }

        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Interp_SpotVol */



/*****  Hyb3_IndexVol    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Hyb3_IndexVol (
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
          int     NbFactor,              /* (I) Number of factors            */
          double  *Beta,                 /* (I) Mean reversions              */
          T_CURVE const* crv)
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

    /* Find B coefficients */
    if (Hyb3_BFactor (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 NbFactor, 
                 Beta,     
                 crv) != SUCCESS)
    {
        goto RETURN;
    }

    /* Find par yield */


    if (ParYieldFromDates (&ParYield,
                           &Annuity,
                           SwapSt,
                           SwapMat,
                           DCC,
                           Freq,
                           'N',
                           crv) != SUCCESS)
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
        if (SwapSt <= crv->ValueDate) 
        {
            *Vol = 0.00;
            return(SUCCESS);
        }
    
        t = Daysact (crv->ValueDate, SwapSt) / 365.;

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
    
        M = D[0][0] * D[0][0] * Hyb3_ExpDecay (2. * Beta[0], t) * t;
    
        if (NbFactor > 1) 
        {
            M += D[1][0] * D[1][0] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            M += D[1][1] * D[1][1] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * Hyb3_ExpDecay (Beta[0] + Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            M += D[2][0] * D[2][0] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            M += D[2][1] * D[2][1] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            M += D[2][2] * D[2][2] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * Hyb3_ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * Hyb3_ExpDecay (Beta[1] + Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * Hyb3_ExpDecay (Beta[1] + Beta[2], t) * t;
        }

        
        if (M < SQUARE(TINY))
        {
            DR_Error ("Hyb3_IndexVol: problem in integrating index volatility "
                     "at date %ld (Hyb3_IndexVol)", YMDDateFromIRDate(SwapSt));
            goto RETURN;
        }
        
        M = sqrt (M / t);

        /* Convert from q measure back to log-normal */
        if ((fabs(QLeft - QRight) < TINY) && (fabs(FwdShift) < TINY))
        {
            M = M * sqrt (t);
            if (fabs(QLeft) > QCUTOFF)                                                 
                M = 2. * Normal_InvH ((NormalH (.5 * QLeft * M) - .5) / QLeft + .5);
            else
                M = 2. * Normal_InvH (.5 * (1. + M / sqrt(2.*PI)));

            *Vol = M / sqrt(t);
        }
        else
        {
            atmPr= Option_BS2Q (ParYield,ParYield,t,M,'C',QLeft,QRight,FwdShift);
            M = ImpVol_BS2Q (ParYield,ParYield,t,atmPr,'C',
                             1.,1.,0.,M); 
	    if( M < 0.0)
            {
                DR_Error("Unable to convert 2Q vol to lognormal measure.\n");
                goto RETURN;
            }
            *Vol = M;
        }


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
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Hyb3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += Aweight[1][i] * Aweight[1][i] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * Hyb3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * Hyb3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[4][i] * Aweight[4][i] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[5][i] * Aweight[5][i] * Hyb3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * Hyb3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * Hyb3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * Hyb3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
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
        DR_Error ("Hyb3_IndexVol: problem in integrating index volatility "
                 "at date %ld (Hyb3_IndexVol)", YMDDateFromIRDate(SwapSt));
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
        M = ImpVol_BS2Q (ParYield,ParYield,T,atmPr,'C',
                        1.,1.,0.,M);
	if ( M <  0.0 )
        {
            DR_Error("Unable to convert vol to 2Q measure.\n");
            goto RETURN;
        }

        *Vol = M;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_IndexVol */



/*****  Hyb3_IndexSpotVol *********************************************************/
/*
*       Spot volatility of an index. 
*       If spot volatility of instantaneous rate is time dependent, this 
*       routine requires the relevant Aweight[6] (i.e. at time SwapSt).
*/
int  Hyb3_IndexSpotVol(double  *Vol,         /* (O) Index vol curve               */
                  char    Freq,         /* (I) Frequency of underlying rate  */
                  char    DCC,          /* (I) Day count convention          */
                  long    SwapSt,       /* (I) Underlying swap start         */
                  long    SwapMat,      /* (I) Underlying swap maturity      */
                  int     NbFactor,     /* (I) Number of factors             */
                  double  *Aweight,     /* (I) Relative size factors         */
                  double  *Beta,        /* (I) Mean reversions               */
                  T_CURVE const* crv)
{

    double  B[3];              /* B in Christian memo                   */
    double  D[3][3];           /* aweight*B                             */
    double  M;                 /* Matrix for lambda system              */

    int     p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially      */

    /* Avoid spot index expiration */
    if (SwapSt < crv->ValueDate) 
    {
        *Vol = 0.00;
        return(SUCCESS);
    }
    
    /* Find B coefficients */
    if (Hyb3_BFactor (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 NbFactor, 
                 Beta,     
                 crv) != SUCCESS)
    {
        goto RETURN;
    }
    
    /* Find D coefficients */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
            D[p][q] = 0.;

    D[0][0] = Aweight[0] * B[0];
    
    if (NbFactor > 1)
    {
        D[1][0] = Aweight[1] * B[1];
        D[1][1] = Aweight[2] * B[1];
    }
    
    if (NbFactor > 2)
    {
        D[2][0] = Aweight[3] * B[2];
        D[2][1] = Aweight[4] * B[2];
        D[2][2] = Aweight[5] * B[2];
    }
    
    M = D[0][0] * D[0][0];
    
    if (NbFactor > 1) 
    {
        M += D[1][0] * D[1][0];
        M += 2. * D[0][0] * D[1][0];
        M += D[1][1] * D[1][1];
    }
    
    if (NbFactor > 2)
    {
        M += D[2][0] * D[2][0];
        M += 2. * D[0][0] * D[2][0];
        M += 2. * D[1][0] * D[2][0];
        M += D[2][1] * D[2][1];
        M += 2. * D[1][1] * D[2][1];
        M += D[2][2] * D[2][2];
    }
    
    if (M < TINY)
    {
        DR_Error ("Hyb3_IndexSpotVol: problem in calculating spot index "
                 "volatility at date %ld ()", YMDDateFromIRDate(SwapSt));
        goto RETURN;
    }
    
    *Vol = sqrt(M);

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_IndexSpotVol */
