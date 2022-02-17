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



/*****  Fix3_BFactor_E2Q    *****************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*/
int     Fix3_BFactor_E2Q (
        double*         B,            /* (O) B ceoff for each factor            */
        long            SwapSt,       /* (I) Underlying swap start              */
        long            SwapMat,      /* (I) Underlying swap maturity           */
        char            DCC,          /* (I) Underlying day count conv          */
        char            Freq,         /* (I) Underlying frequency               */
        MKTVOL_DATA*	mktvol_data,  /* (I) Volatility data                    */
        T_CURVE const*  crv)          /* (I) zero curve                         */
{

    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmt                 */
    long        CpnPmtDate;           /* Current cpn pmt date                   */
    long        TempDates[2];         /* To construct cpn list                  */

    double      S;                   /* Swap start in years from base date      */
    double      TtoCpn;              /* Time to current coupon in years         */
    double      PrevT;               /* Same for previous coupon                */

    double      DCCFrac;             /* Current coupon day count fraction       */
    double      Annuity;             /* Annuity price                           */
    double      FwdYield;            /* Forward yield                           */

    double      ZerotoS;             /* Zero to swap start                      */
    double      ZerotoCpn=0;         /* Zero to current cpn                     */
    double      PrevZero;

    double      BbqAdj;              /* Back bone adjustment in A coeff         */
    double      A[3];                /* A in Christian's memo                   */
    int         j, k;
    int         status = FAILURE;    /* Error status = FAILURE initially        */

    /* temporary variables for volatility data parameters  */
    double const*   Beta;           /* (I) Mean reversions	               */
    double          Bbq;            /* (I) Backbone parameter              */
    double          Amap;           /* (I) A mapping parameter             */
    double          Bmap;           /* (I) B mapping parameter             */
    double          VolNorm;        /* (I) Normal volatility in backbone   */
    double          VolLogn;         /* (I) Lognorm volatility in backbone  */
    int             NbFactor;       /* (I) NbFactors                       */

    /* set temporary variables equal to mktvol_data parameters */
    Beta        = mktvol_data->Beta;
    Bbq         = mktvol_data->Bbq;
    Amap        = mktvol_data->Amap;
    Bmap        = mktvol_data->Bmap;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;

    TempDates[0] = SwapSt;
    TempDates[1] = SwapMat;

    CpnEventList = DrNewEventListFromFreq (  2,
                                             TempDates,
                                             Freq,
                                             'F',   /* Always front stub */
                                             'N',   /* Dates in not required */
                                             NULL, NULL, NULL, NULL, NULL);

    if (CpnEventList == NULL) goto FREE_MEM_AND_RETURN;


    /* Zero to swap start */
    ZerotoS = GetZeroPrice(SwapSt, crv);
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    S = Daysact (crv->ValueDate, SwapSt) / 365.;

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
            DR_Error ("Could not calculate day count fraction (Fix3_BFactor_E2Q)!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Zero to current cpn */
        ZerotoCpn = GetZeroPrice(CpnPmtDate, crv);
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;

        TtoCpn = Daysact (crv->ValueDate, CpnPmtDate) / 365.;
        
        Annuity += DCCFrac * ZerotoCpn;
            
        /* cf Christian's memo for Vladimir's approximation */
        if (mktvol_data->VolUnit == 1)
        {
            BbqAdj =      (Bbq) * VolLogn * log (PrevZero/ZerotoCpn) 
               + (1. - Bbq) * VolNorm * (TtoCpn - PrevT); 
        }
        else
        {
            BbqAdj =       Bbq * VolLogn * (TtoCpn-PrevT); 
        }

        for (k = 0; k < NbFactor; k++)
        {
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
        if (mktvol_data->VolUnit == 1)
        {
            B[k] *= FwdYield;
            B[k] += A[k] * ZerotoCpn;
            B[k] /= (ZerotoS - ZerotoCpn);
        }
        else
        {
            B[k] *= (ZerotoS - ZerotoCpn) / Annuity;
            B[k] += A[k] * ZerotoCpn;
            B[k] /= Annuity; 
        }
    }


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* Fix3_BFactor_E2Q */




/*****  Fix3_SpotVol_E2Q    *****************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Fix3_SpotVol_E2Q    (
       MKTVOL_DATA*     mktvol_data,      /* (I/O) Volatility data */
       T_CURVE          const* crv)       /* (I) Zero curve        */
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

   
    /* local variables of vol data structure */
    long            VolBaseDate;    /* (I) Volatility curve base date       */
    int             NbVol;          /* (I) Nb of points in vol curve        */
    long const*     VolDate;        /* (I) Volatility date                  */
    double const*   Vol;            /* (I) Vol curve                        */
    int             VolTypeFlag;    /* (I) Normal or Lognormal              */
    int*            VolUsed;        /* (I/O) TRUE if vol used in calibration*/
    char            Freq;           /* (I) Frequency of underlying rate     */
    char            DCC;            /* (I) Day count convention             */
    long const*     SwapSt;         /* (I) Underlying swap starr            */
    long const*     SwapMat;        /* (I) Underlying swap maturity         */
    double          FwdSh;          /* (I) Fwd shift mapping coefficient    */
    double          Bbq;            /* (I) Backbone parameter               */
    double          VolNorm;        /* (I) Normal volatility in backbone    */
    double          VolLogn;        /* (I) Lognorm volatility in backbone   */
    int             NbFactor;       /* (I) Number of factors	            */
    double const*   Alpha;          /* (I) Relative size factors            */
    double const*   Beta;           /* (I) Mean reversions                  */
    double const*   Rho;            /* (I) Correlation between factors      */
    int             SkipFlag;       /* (I) Skip calibration failure points  */
    int             CalibFlag;      /* (I) Index calibration flag           */
    double          QLeft, QRight;  /* (I) Q's                              */
    double          Amap, Bmap;     /* (I) mapping parameters               */
    double Aweight[6][MAXNBDATE];  /* (O) Spot vol curve of each factor   */



    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Vol         = mktvol_data->Vol;
    VolTypeFlag = mktvol_data->VolUnit;
    VolUsed     = mktvol_data->VolUsed;
    Freq        = mktvol_data->Freq;
    DCC         = mktvol_data->DCC;
    SwapSt      = mktvol_data->SwapSt;
    SwapMat     = mktvol_data->SwapMat;
    QLeft       = mktvol_data->QLeft;
    QRight      = mktvol_data->QRight;
    Amap        = mktvol_data->Amap;
    Bmap        = mktvol_data->Bmap;
    FwdSh       = mktvol_data->FwdShift;
    Bbq         = mktvol_data->Bbq;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    Alpha       = mktvol_data->Alpha;
    Beta        = mktvol_data->Beta;
    Rho         = mktvol_data->Rho;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;

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
        mktvol_data->Aweight[0][0] = aw[0][0];
        
        if (NbFactor > 1)
        {
            mktvol_data->Aweight[1][0] = aw[1][0];
            mktvol_data->Aweight[2][0] = aw[1][1];
        }
        
        if (NbFactor > 2)
        {
            mktvol_data->Aweight[3][0] = aw[2][0];
            mktvol_data->Aweight[4][0] = aw[2][1];
            mktvol_data->Aweight[5][0] = aw[2][2];
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
        if (Fix3_BFactor_E2Q (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        mktvol_data,
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

        if (Conv_MarketVol_To_XSpace(pY[i], 
                                     Vol[i], 
                                     VolDate[i],VolTypeFlag,
                                     &y,
                                     T,
                                     FwdSh,
                                     QLeft,
                                     QRight,
                                     Amap, 
                                     Bmap) == FAILURE) 
        {
            goto RETURN;
    }
         

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * Beta[0], t) * t;
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * Beta[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(Beta[0] + Beta[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * Beta[1] * t);
            
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (Beta[0] + Beta[1], t) * t;
            
            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(Beta[0] + Beta[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * Beta[2] * t);
            
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * Beta[2], t) * t;
        }
        
        /*
         *  Solve using Gauss-Jordan method
         */
        if (fabs(M) < TINY)
        {
            DR_Error ("Spot vol: problem in bootstrapping %ld "
                                "volatility (negative variance)",
                     VolDate[i]);
            goto RETURN;            
        }

        y /= M;
        
        /*
         * Check if solution is positive 
         */
        if (y < TINY)
        {
            DR_Error("Spot vol: problem in bootstrapping %ld "
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
            if ((lambdaNew/lambdaI < MIN_SPOT_VOL_RATIO) || 
                (lambdaNew/lambdaI > 1./MIN_SPOT_VOL_RATIO))
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
        L[0][0] += lambda * lambda * Fix3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += lambda * lambda * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += lambda * lambda * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
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
    
    /* store computed Aweights in the mkt vol structure CalibFlag = TRUE */   
    for (i = 0; i < NbVol; i++)
    {
        mktvol_data->Aweight[0][i] = Aweight[0][i];
        if (NbFactor > 1)
        {
            mktvol_data->Aweight[1][i] = Aweight[1][i];
            mktvol_data->Aweight[2][i] = Aweight[2][i];
        }
        if (NbFactor > 2)
        {
            mktvol_data->Aweight[3][i] = Aweight[3][i];
            mktvol_data->Aweight[4][i] = Aweight[4][i];
            mktvol_data->Aweight[5][i] = Aweight[5][i];
        }

    }
  
    if (LastCalibIdx == -1)
    {
        DR_Error ("Spot vol: none of points is calibrated!");
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_SpotVol_E2Q */



/*****  Fix3_Interp_SpotVol_E2Q  *************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Fix3_Interp_SpotVol_E2Q (FIX3_TREE_DATA*  tree_data,   /* (I/ O) Tree data   */
                                 MKTVOL_DATA*     mktvol_data) /* (I)Volatility data */
            
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  PrevVolT;              /* Time to previous vol point */
    double  Cov[3][3];             /* Covariance matrix          */
    double  x, t, T;

    int     i, j, k, l, p, q;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */
    
    /* local variables for tree and vol structure */
    double**         AweightCurve;         /* (O) Factors aweight curves        */
    long             VolBaseDate;           /* (I) Volatility curve base date    */
    int              NbVol;                 /* (I) Number of spot vol points     */
    long const*      VolDate;               /* (I) Spot vol dates                */
    double           Aweight[6][MAXNBDATE]; /* (I) Aweight curve of each factor  */
    int              NbFactor;              /* (I) Number of factors             */
    double const*    Beta;                  /* (I) Mean reversions               */
    double           Bbq;                   /* (I) Backbone parameter            */
    double           VolNorm;               /* (I) Normal volatility in backbone */
    double           VolLogn;               /* (I) Lognorm volatility in backbone*/
    int              CalibFlag;             /* (I) Index calibration flag        */
    double const*    FwdRate;               /* (I) One period forward rate       */
    double const*    Length;                /* (I) Length of each time step      */
    int              NbTP;                  /* (I) Total number of time points   */


    /* set temporary volatility param to mktvol_data param */
    VolBaseDate  = mktvol_data->BaseDate;
    NbVol        = mktvol_data->NbVol;
    VolDate      = mktvol_data->VolDate;
    NbFactor     = mktvol_data->NbFactor;
    Beta         = mktvol_data->Beta;
    Bbq          = mktvol_data->Bbq;
    VolNorm      = mktvol_data->VolNorm;
    VolLogn      = mktvol_data->VolLogn;
    CalibFlag    = mktvol_data->CalibFlag;
    FwdRate      = tree_data->FwdRate[tree_data->CvDiff];
    Length       = tree_data->Length;
    NbTP         = tree_data->NbTP;
    AweightCurve = tree_data->Aweight;
    
    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    /* set temporary Aweight array to mktvol_data Aweight */
    if (CalibFlag == TRUE)
    {
        for (i = 0; i < NbVol; i++)
        {
            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][i] = mktvol_data->Aweight[k][i]; 
            }
        }
    }
    else
    {
        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][0] = mktvol_data->Aweight[k][0]; 
        }
    }


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

                AweightCurve[k][i] *=(1.+ FwdRate[i+1]) * Fix3_ExpDecay(Beta[l], Length[i+1]);
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
            Cov[0][0]  = Aweight[0][j-1] * Aweight[0][j-1] * Fix3_ExpDecay (2. * Beta[0], VolT[j-1] - t) * (VolT[j-1] - t);
            Cov[0][0] += Aweight[0][j]   * Aweight[0][j]   * Fix3_ExpDecay (2. * Beta[0], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[0] * (VolT[j-1] - t));
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = Aweight[0][j-1] * Aweight[1][j-1] * Fix3_ExpDecay (Beta[0] + Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][0] += Aweight[0][j]   * Aweight[1][j]   * Fix3_ExpDecay (Beta[0] + Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[1])*(VolT[j-1]-t));
            
                Cov[1][1]  = Aweight[1][j-1] * Aweight[1][j-1] * Fix3_ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[1][j]   * Aweight[1][j]   * Fix3_ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
                Cov[1][1] += Aweight[2][j-1] * Aweight[2][j-1] * Fix3_ExpDecay (2. * Beta[1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[2][j]   * Aweight[2][j]   * Fix3_ExpDecay (2. * Beta[1], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1] * (VolT[j-1] - t));
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = Aweight[0][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (Beta[0] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][0] += Aweight[0][j]   * Aweight[3][j]   * Fix3_ExpDecay (Beta[0] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0]+Beta[2])*(VolT[j-1]-t));
            
                Cov[2][1]  = Aweight[1][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[1][j]   * Aweight[3][j]   * Fix3_ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));
                Cov[2][1] += Aweight[2][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (Beta[1] + Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[2][j]   * Aweight[4][j]   * Fix3_ExpDecay (Beta[1] + Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1]+Beta[2])*(VolT[j-1]-t));

                Cov[2][2]  = Aweight[3][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[3][j]   * Aweight[3][j]   * Fix3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[4][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[4][j]   * Aweight[4][j]   * Fix3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[5][j-1] * Aweight[5][j-1] * Fix3_ExpDecay (2. * Beta[2], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[5][j]   * Aweight[5][j]   * Fix3_ExpDecay (2. * Beta[2], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2] * (VolT[j-1] - t));
            }


            x = Cov[0][0] / Fix3_ExpDecay (2. * Beta[0], T - t) / (T - t);

            if (x < TINY)
            {
                DR_Error ("Fix3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / Fix3_ExpDecay (Beta[0] + Beta[1], T - t) / (T - t);

                x = Cov[1][1] / Fix3_ExpDecay (2. * Beta[1], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

                if (x < TINY)
                {
                    DR_Error ("Fix3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                    goto RETURN;
                }
                        
                AweightCurve[2][i] = sqrt (x);
            }

            if (NbFactor > 2)
            {
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / Fix3_ExpDecay (Beta[0] + Beta[2], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / Fix3_ExpDecay (Beta[1] + Beta[2], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / Fix3_ExpDecay (2. * Beta[2], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

                if (x < TINY)
                {
                    DR_Error ("Fix3_Interp_SpotVol: problem in "
                              "interpolating spot volatility at node #%d!", i);
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
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1]) * Fix3_ExpDecay(Beta[l], Length[i+1]);
        }
    }

        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Interp_SpotVol_E2Q */



/*****  Fix3_IndexVol_E2Q    ******************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Fix3_IndexVol_E2Q (
          double*       Vol,                  /* (O) Index vol curve              */
          long          SwapSt,               /* (I) Underlying swap start        */
          long          SwapMat,              /* (I) Underlying swap maturity     */
          char          Freq,                 /* (I) Frequency of underlying rate */
          char          DCC,                  /* (I) Day count convention         */        
          MKTVOL_DATA   *mktvol_data,         /* (I) Volatility data              */
          T_CURVE const* crv)                 /* (I) zero curve                   */
{

    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  t;                 /* Time between two consecutive expiries */
    double  T=0.;              /* Time to current expiry                */
    double  ParYield;
    double  Annuity;
   

    double  L[3][3];           /* Integrals of factor spot vol          */
    double  B[3];              /* B in Christian memo                   */
    double  D[3][3];           /* aweight*B                             */
    double  M;                 /* Total variance                        */

    long    EndDate;           /* End of current integration bucket     */
    int     i, j, p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially      */


    int             NbVol;          /* (I) Nb of points in vol curve        */
    long            VolBaseDate;    /* (I) Volatility base date             */
    long const*     VolDate;        /* (I) Volatility dates                 */
    int             VolTypeFlag;    /* (I) Normal or Lognormal              */
    double          FwdShift;       /* (I) Fwd shift mapping coefficient    */
    double          Bbq;            /* (I) Backbone parameter               */
    double          VolNorm;        /* (I) Normal volatility in backbone    */
    double          VolLogn;        /* (I) Lognorm volatility in backbone   */
    int             NbFactor;       /* (I) Number of factors                */
    double const*   Alpha;          /* (I) Relative size factors            */
    double const*   Beta;           /* (I) Mean reversions                  */
    double const*   Rho;            /* (I) Correlation between factors      */
    int             SkipFlag;       /* (I) Skip calibration failure points  */
    int             CalibFlag;      /* (I) Index calibration flag           */
    double          QLeft, QRight;  /* (I) Q's                              */
    double          Amap, Bmap;     /* (I) mapping parameters               */
    double   Aweight[6][MAXNBDATE]; /* (I) Spot vol curve of each factor    */
    


    VolBaseDate     = mktvol_data->BaseDate;
    NbVol           = mktvol_data->NbVol;
    VolDate         = mktvol_data->VolDate;
    VolTypeFlag     = mktvol_data->VolUnit;
    QLeft           = mktvol_data->QLeft;
    QRight          = mktvol_data->QRight;
    Amap            = mktvol_data->Amap;
    Bmap            = mktvol_data->Bmap;
    FwdShift        = mktvol_data->FwdShift;
    Bbq             = mktvol_data->Bbq;
    VolNorm         = mktvol_data->VolNorm;
    VolLogn         = mktvol_data->VolLogn;
    NbFactor        = mktvol_data->NbFactor;
    Alpha           = mktvol_data->Alpha;
    Beta            = mktvol_data->Beta;
    Rho             = mktvol_data->Rho;
    SkipFlag        = mktvol_data->SkipFlag;
    CalibFlag       = mktvol_data->CalibFlag;
    

    if (CalibFlag == TRUE)
    {
        for (i = 0; i < NbVol; i++)
        {
            Aweight[0][i] = mktvol_data->Aweight[0][i]; 
            if (NbFactor > 1)
            {
                Aweight[1][i] = mktvol_data->Aweight[1][i];
                Aweight[2][i] = mktvol_data->Aweight[2][i];
            }
            if (NbFactor > 2)
            {
                Aweight[3][i] = mktvol_data->Aweight[3][i];
                Aweight[4][i] = mktvol_data->Aweight[4][i];
                Aweight[5][i] = mktvol_data->Aweight[5][i];
            }

        }
    }
    else
    {
        Aweight[0][0] = mktvol_data->Aweight[0][0]; 
        if (NbFactor > 1)
        {
            Aweight[1][0] = mktvol_data->Aweight[1][0];
            Aweight[2][0] = mktvol_data->Aweight[2][0];
        }
        if (NbFactor > 2)
        {
            Aweight[3][0] = mktvol_data->Aweight[3][0];
            Aweight[4][0] = mktvol_data->Aweight[4][0];
            Aweight[5][0] = mktvol_data->Aweight[5][0];
        }
    }


    /* Find B coefficients */
    if (Fix3_BFactor_E2Q (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 mktvol_data,
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
                           'F',
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
    
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * Beta[0], t) * t;
    
        if (NbFactor > 1) 
        {
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (Beta[0] + Beta[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
        }

        
        if (M < SQUARE(TINY))
        {
            DR_Error ("Fix3_IndexVol: problem in integrating index volatility "
                     "at date %ld (Fix3_IndexVol)", YMDDateFromIRDate(SwapSt));
            goto RETURN;
        }
        
        M = sqrt (M / t);

        /* Convert from q measure back to log-normal */
        if (Conv_XSpaceVol_To_Market(ParYield,
                                     Vol,       
                                     VolTypeFlag,
                                     M,
                                     t,
                                     FwdShift,
                                     QLeft,
                                     QRight,
                                     Amap,
                                     Bmap) == FAILURE)
        {
            goto RETURN;
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
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
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
        DR_Error ("Fix3_IndexVol: problem in integrating index volatility "
                 "at date %ld (Fix3_IndexVol)", YMDDateFromIRDate(SwapSt));
        goto RETURN;
    }

    M = sqrt (M / T);

    /* Convert from q measure back to log-normal */
    if (Conv_XSpaceVol_To_Market(ParYield,
                                 Vol,        
                                 VolTypeFlag,
                                 M,
                                 T,
                                 FwdShift,
                                 QLeft,
                                 QRight,
                                 Amap,
                                 Bmap) == FAILURE) 
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_IndexVol_E2Q */


/*****  Fix3_Filtered_SpotVol_E2Q    **********************************************/
/*
*       Bootstrap input volatility curve to extract filtered spot vols.
*       Spot vols are not modified provided the ratio of the current spot vol
*       to the spot vol in the first bucket is greater than MIN_SPOT_VOL_RATIO
*       or less than 1/MIN_SPOT_VOL_RATIO. Otherwise the spot vols are smoothly
*       capped and floored, with the cap and floor levels parameterized by
*       SPOT_VOL_FILTER_AMOUNT.
*/
int    Fix3_Filtered_SpotVol_E2Q(
       MKTVOL_DATA			*mktvol_data,	/* (I/O) Volatility data		*/
       T_CURVE const*		crv)			/* (I) zero curve				*/
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
    double  lambdaI = 0.0001;     /* Relative weight, first value           */

    /* filter parameters */
    double  alpha = SPOT_VOL_FILTER_AMOUNT;
    double  lr    = MIN_SPOT_VOL_RATIO;
    double  ymin;                 /* Minimum variance                       */
    double  ymax;                 /* Maximum variance                       */
    double  lStep;                /* Low varaince step size                 */
    double  hStep;                /* High variance step size                */

    int     NbVolUsed;            /* Nb of vol points used in calibration   */
    int     i, k, p, q;
    int     status = FAILURE;     /* Error status = FAILURE initially       */



    /* local variables of vol data structure */
    long            VolBaseDate;    /* (I) Volatility curve base date       */
    int             NbVol;          /* (I) Nb of points in vol curve        */
    long const*     VolDate;        /* (I) Volatility date                  */
    double const*   Vol;            /* (I) Vol curve                        */
    int             VolTypeFlag;    /* (I) Normal or Lognormal              */
    int*            VolUsed;        /* (I/O) TRUE if vol used in calibration*/
    char            Freq;           /* (I) Frequency of underlying rate     */
    char            DCC;            /* (I) Day count convention             */
    long const*     SwapSt;         /* (I) Underlying swap starr            */
    long const*     SwapMat;        /* (I) Underlying swap maturity         */
    double          FwdSh;          /* (I) Fwd shift mapping coefficient    */
    double          Bbq;            /* (I) Backbone parameter               */
    double          VolNorm;        /* (I) Normal volatility in backbone    */
    double          VolLogn;        /* (I) Lognorm volatility in backbone   */
    int             NbFactor;       /* (I) Number of factors	            */
    double const*   Alpha;          /* (I) Relative size factors            */
    double const*   Beta;           /* (I) Mean reversions                  */
    double const*   Rho;            /* (I) Correlation between factors      */
    int             SkipFlag;       /* (I) Skip calibration failure points  */
    int             CalibFlag;      /* (I) Index calibration flag           */
    double          QLeft, QRight;  /* (I) Q's                              */
    double          Amap, Bmap;     /* (I) mapping parameters               */
    double Aweight[6][MAXNBDATE];  /* (O) Spot vol curve of each factor   */



    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Vol         = mktvol_data->Vol;
    VolTypeFlag = mktvol_data->VolUnit;
    VolUsed     = mktvol_data->VolUsed;
    Freq        = mktvol_data->Freq;
    DCC         = mktvol_data->DCC;
    SwapSt      = mktvol_data->SwapSt;
    SwapMat     = mktvol_data->SwapMat;
    QLeft       = mktvol_data->QLeft;
    QRight      = mktvol_data->QRight;
    Amap        = mktvol_data->Amap;
    Bmap        = mktvol_data->Bmap;
    FwdSh       = mktvol_data->FwdShift;
    Bbq         = mktvol_data->Bbq;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    Alpha       = mktvol_data->Alpha;
    Beta        = mktvol_data->Beta;
    Rho         = mktvol_data->Rho;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;

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
        mktvol_data->Aweight[0][0] = aw[0][0];
        
        if (NbFactor > 1)
        {
            mktvol_data->Aweight[1][0] = aw[1][0];
            mktvol_data->Aweight[2][0] = aw[1][1];
        }
        
        if (NbFactor > 2)
        {
            mktvol_data->Aweight[3][0] = aw[2][0];
            mktvol_data->Aweight[4][0] = aw[2][1];
            mktvol_data->Aweight[5][0] = aw[2][2];
        }

        return (SUCCESS);
    }

    /* set nb of vol points used in calibration */
    NbVolUsed = 0;
    while ( (NbVolUsed < NbVol) && (VolUsed[NbVolUsed] == TRUE) )
        NbVolUsed++;

    if (NbVolUsed == 0)
    {
        DR_Error("Filtered Spot vol: require at least one vol point "
                 "to calibrate!");
        goto RETURN;
    }

    /* check that all subsequent vol points are not used */ 
    for (i = NbVolUsed; i < NbVol; i++)
    {
        if (VolUsed[i] != FALSE)
        {
            DR_Error("Filtered Spot vol: cannot find last vol point "
                     "to calibrate!");
            goto RETURN;
        }
    }

    for (i = 0; i < NbVolUsed; i++)
    {
        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

        /* B factor */
        if (Fix3_BFactor_E2Q (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        mktvol_data,
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


    /* Set filter parameters */
    ymin  = alpha * lr;
    ymin *= ymin;
    ymax  = 1./ymin;

    lStep = lr*lr - ymin;
    hStep = ymax  - 1./(lr*lr);

    /* Bootstrap swaption volatility */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
        {
            L[p][q] = 0.;
            D[p][q] = 0.;
        }

    for (i = 0; i < NbVolUsed; i++)
    {
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));
        
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

        if (Conv_MarketVol_To_XSpace(pY[i],
                                     Vol[i],        
                                     VolDate[i], 
                                     VolTypeFlag,
                                     &y,
                                     T,
                                     FwdSh,
                                     QLeft,
                                     QRight,
                                     Amap,
                                     Bmap) == FAILURE) 
        {
            goto RETURN;
        }


        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * Beta[0], t) * t;
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * Beta[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(Beta[0] + Beta[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * Beta[1] * t);
            
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (Beta[0] + Beta[1], t) * t;
            
            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * Beta[1], t) * t;
        }

        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(Beta[0] + Beta[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * Beta[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(Beta[1] + Beta[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * Beta[2] * t);
            
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (Beta[0] + Beta[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (Beta[1] + Beta[2], t) * t;
            
            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * Beta[2], t) * t;
        }

        /* M must be stictly positive, and since it is not a function   */
        /* of spot vol, trap the error now rather than later in Fix3_SpotVol */
        if (M < TINY)
        {
            DR_Error ("Filtered spot vol: problem in bootstrapping %ld "
                                "volatility (negative variance)", YMDDateFromIRDate(VolDate[i]));
            goto RETURN;
        }

        y /= M;

        /* Filter variance ratio: note that initial bucket */
        /* is only filtered if the variance is too low     */
        y /= (lambdaI*lambdaI);

        /* Low variance filter */
        lambda = DrSmoothMax(y, ymin, lStep);

        /* set lambdaI for subsequent buckets */
        if (i == 0)
        {
            lambda  = lambdaI * sqrt(lambda);
            lambdaI = lambda;
        }
        else
        {
            /* High variance filter */
            lambda -= DrSmoothMax(y - ymax, 0, hStep);
            lambda  = lambdaI * sqrt(lambda);
        }

        /* Relative weights of factors */
        Aweight[0][i] = lambda * aw[0][0];
        
        if (NbFactor > 1)
        {
            Aweight[1][i] = lambda * aw[1][0];
            Aweight[2][i] = lambda * aw[1][1];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][i] = lambda * aw[2][0];
            Aweight[4][i] = lambda * aw[2][1];
            Aweight[5][i] = lambda * aw[2][2];
        }

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * Beta[0] * t);
        L[0][0] += lambda * lambda * Fix3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += lambda * lambda * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += lambda * lambda * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += 2.*lambda*lambda * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
        }
    }  /* for i */

    for (k = NbVolUsed; k < NbVol; k++)
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

     /* store computed Aweights in the mkt vol structure CalibFlag = TRUE */   
    for (i = 0; i < NbVol; i++)
    {
        mktvol_data->Aweight[0][i] = Aweight[0][i];

        if (NbFactor > 1)
        {
            mktvol_data->Aweight[1][i] = Aweight[1][i];
            mktvol_data->Aweight[2][i] = Aweight[2][i];
        }
        if (NbFactor > 2)
        {
            mktvol_data->Aweight[3][i] = Aweight[3][i];
            mktvol_data->Aweight[4][i] = Aweight[4][i];
            mktvol_data->Aweight[5][i] = Aweight[5][i];
        }
    }
 

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Filtered_SpotVol_E2Q */

 



/*****  GenIndexVol_E2Q    *****************************************************/
/*
*       Integrates spot volatility curve between two time points falling 
*       before or on the swap start date.
*/
int   Fix3_GenIndexVol_E2Q (
          double  *Vol,                  /* (O) Average vol bet 2 time points  */
          long    IStart,                /* (I) Integration start point        */
          long    IEnd,                  /* (I) Integration end point          */
          long    SwapStart,             /* (I) Underlying swap start          */
          long    SwapMat,               /* (I) Underlying swap maturity       */
          char    Freq,                  /* (I) Frequency of underlying rate   */
          char    DCC,                   /* (I) Day count convention           */
          MKTVOL_DATA *mktvol_data,      /* (I) Volatility data                */
          T_CURVE const* crv)
{

    double  VolT[MAXNBDATE];   /* Expiries in years							*/
    double  t;                 /* Time between two consecutive expiries     */
    double  T=0.;              /* Time to current expiry                    */
    double  Int01=0;           /* size of int from IStart to IEnd           */ 
    double  Int12=0;           /* size of int from IEnd to SwapStart        */ 
    double  Int02=0;           /* size of int from IStart to SwapStart      */ 
    double  ParYield;
    double  Annuity;

    double  L[3][3];           /* Integrals of factor spot vol              */
    double  B[3];              /* B in Christian memo                       */
    double  D[3][3];           /* aweight*B                                 */
    double  M;                 /* Total variance                            */

    long    EndDate;           /* End of current integration bucket         */
    int     StartIdx=0;
    int     i, j, p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially          */
    char    ErrorMsg[MAXBUFF]; /* Error message                             */
    
    /* temporary variables for volatility parameters                        */
    int             NbVol;          /* (I) Nb of points in vol curve        */
    long            VolBaseDate;    /* (I) Volatility base date             */
    long const*     VolDate;        /* (I) Volatility dates                 */
    int             VolTypeFlag;    /* (I) Normal or Lognormal              */
    double          FwdShift;       /* (I) Fwd shift mapping coefficient    */
    double          Bbq;            /* (I) Backbone parameter               */
    double          VolNorm;        /* (I) Normal volatility in backbone    */
    double          VolLogn;        /* (I) Lognorm volatility in backbone   */
    int             NbFactor;       /* (I) Number of factors                */
    double const*   Alpha;          /* (I) Relative size factors            */
    double const*   Beta;           /* (I) Mean reversions                  */
    double const*   Rho;            /* (I) Correlation between factors      */
    int             SkipFlag;       /* (I) Skip calibration failure points  */
    int             CalibFlag;      /* (I) Index calibration flag           */
    double          QLeft, QRight;  /* (I) Q's                              */
    double          Amap, Bmap;     /* (I) mapping parameters               */
    double   Aweight[6][MAXNBDATE];/* (O) Spot vol curve of each factor	*/
    


    /* Avoid spot index expiration; note we use base date of zero */
    /* curve as there is no volatility information available.     */
    if ((IEnd <= crv->ValueDate)|| ( IStart == IEnd))
    {
            *Vol = 0.00;
            return(SUCCESS);
    }
    
    /* set temporary variables to mktvol_data parameters */
    VolBaseDate     = mktvol_data->BaseDate;
    NbVol           = mktvol_data->NbVol;
    VolDate         = mktvol_data->VolDate;
    VolTypeFlag     = mktvol_data->VolUnit;
    QLeft           = mktvol_data->QLeft;
    QRight          = mktvol_data->QRight;
    Amap            = mktvol_data->Amap;
    Bmap            = mktvol_data->Bmap;
    FwdShift        = mktvol_data->FwdShift;
    Bbq             = mktvol_data->Bbq;
    VolNorm         = mktvol_data->VolNorm;
    VolLogn         = mktvol_data->VolLogn;
    NbFactor        = mktvol_data->NbFactor;
    Alpha           = mktvol_data->Alpha;
    Beta            = mktvol_data->Beta;
    Rho             = mktvol_data->Rho;
    SkipFlag        = mktvol_data->SkipFlag;
    CalibFlag       = mktvol_data->CalibFlag;

    /* set termporary Aweight to mktvol_data */
    if (CalibFlag == TRUE)
    {
        for (i = 0; i < NbVol; i++)
        {
            Aweight[0][i] = mktvol_data->Aweight[0][i]; 
            if (NbFactor > 1)
            {
                Aweight[1][i] = mktvol_data->Aweight[1][i];
                Aweight[2][i] = mktvol_data->Aweight[2][i];
            }
            if (NbFactor > 2)
            {
                Aweight[2][i] = mktvol_data->Aweight[2][i];
                Aweight[3][i] = mktvol_data->Aweight[3][i];
                Aweight[4][i] = mktvol_data->Aweight[4][i];
            }

        }
    }
    else
    {
        Aweight[0][0] = mktvol_data->Aweight[0][0]; 
        if (NbFactor > 1)
        {
            Aweight[1][0] = mktvol_data->Aweight[1][0];
            Aweight[2][0] = mktvol_data->Aweight[2][0];
        }
        if (NbFactor > 2)
        {
            Aweight[3][0] = mktvol_data->Aweight[3][0];
            Aweight[4][0] = mktvol_data->Aweight[4][0];
            Aweight[5][0] = mktvol_data->Aweight[5][0];
        }
    }


    Int01 = Daysact (IStart, IEnd)      / 365.;
    Int12 = Daysact (IEnd, SwapStart)   / 365.;
    Int02 = Daysact (IStart, SwapStart) / 365.;  
    
    /* Find B coefficients */
    if (Fix3_BFactor_E2Q (B,
                 SwapStart,
                 SwapMat,
                 DCC,
                 Freq,
                 mktvol_data,
                 crv) != SUCCESS)
    {
        goto RETURN;
    }

    /* Find par yield */
    if (ParYieldFromDates (&ParYield,
                           &Annuity,
                            SwapStart,
                            SwapMat,
                            DCC,
                            Freq,
                            'F',
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
    
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * Beta[0], Int01) * Int01 * 
            exp(-2 * Beta[0] * Int12);
    
        if (NbFactor > 1) 
        {
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * Beta[1], Int01) * Int01 * 
                 exp(-2 * Beta[1] * Int12);

            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * Beta[1], Int01) * Int01 * 
                exp(-2 * Beta[1] * Int12);

            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (Beta[0] + Beta[1], Int01) * 
                 Int01 * exp(-(Beta[0]+Beta[1]) * Int12);
        }
        
        if (NbFactor > 2)
        {
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * Beta[2], Int01) * Int01 * 
                exp(-2 * Beta[2] * Int12);

            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * Beta[2], Int01) * Int01 * 
                exp(-2 * Beta[2] * Int12);

            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * Beta[2], Int01) * Int01 * 
                exp(-2 * Beta[2] * Int12);

            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (Beta[0] + Beta[2], Int01) * 
                Int01 * exp(-(Beta[0]+Beta[2]) * Int12);

            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (Beta[1] + Beta[2], Int01) *
                Int01 * exp(-(Beta[1]+Beta[2]) * Int12);

            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (Beta[1] + Beta[2], Int01) * 
                Int01 * exp(-(Beta[1]+Beta[2]) * Int12);
        }

        
        if (M < SQUARE(TINY))
        {
            sprintf (ErrorMsg, "GenIndexVol: problem in integrating index volatility "
                     "bet dates %ld, and %ld", IStart, IEnd);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        
        M = sqrt (M / Int01);
        
        /* Convert from q measure back to log-normal */
        if (Conv_XSpaceVol_To_Market(ParYield,
                                    Vol,        
                                    VolTypeFlag,
                                    M,
                                    Int01,
                                    FwdShift,
                                    QLeft,
                                    QRight,
                                    Amap,
                                    Bmap) == FAILURE) 
        {
            goto RETURN;
}

        return (SUCCESS);

    }  /* if */
   
    /* Avoid spot index expiration */
    if (IEnd <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }

    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++) 
            L[i][j] = 0.;


    /*Look for start of integration bucket */ 
    while ((StartIdx < NbVol - 1) && (IStart > VolDate[StartIdx]))
        StartIdx ++;

    for (i = StartIdx; i < NbVol; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */

        EndDate = (i < NbVol-1) ? MIN(IEnd, VolDate[i]) : IEnd;

        /* Time to expiry in years */
        VolT[i] = Daysact (IStart, EndDate) / 365.;

        T = VolT[i];
        t = ((i == StartIdx) ? VolT[StartIdx] : (VolT[i]-VolT[i-1]));
        

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * Beta[0] * t) ;
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);

            L[1][1] += Aweight[1][i] * Aweight[1][i] * 
                   Fix3_ExpDecay (2. * Beta[1], t) * t;

            L[1][1] += Aweight[2][i] * Aweight[2][i] * 
                   Fix3_ExpDecay (2. * Beta[1], t) * t;

            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);

            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * 
                   Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * 
               Fix3_ExpDecay (2. * Beta[2], t) * t;

            L[2][2] += Aweight[4][i] * Aweight[4][i] * 
               Fix3_ExpDecay (2. * Beta[2], t) * t;

            L[2][2] += Aweight[5][i] * Aweight[5][i] * 
                Fix3_ExpDecay (2. * Beta[2], t) * t;

            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);

            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * 
                Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;

            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);

            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * 
                Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;

            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * 
                Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == IEnd)
        {
            break;
        }
    }  /* for i */

    M = B[0] * B[0] * L[0][0] * exp(-2. * Beta[0] * Int12);
    
    if (NbFactor > 1) 
    {
        M += B[1] * B[1] * L[1][1] * exp(-2. * Beta[1] * Int12);
        M += B[0] * B[1] * L[0][1] * exp(-(Beta[0]+Beta[1]) * Int12);
    }
    
    if (NbFactor > 2)
    {
        M += B[2] * B[2] * L[2][2] * exp(-2 * Beta[2] * Int12);
        M += B[0] * B[2] * L[0][2] * exp(-(Beta[0]+Beta[2]) * Int12);
        M += B[1] * B[2] * L[1][2] * exp(-(Beta[1]+Beta[2]) * Int12);
    }
    

    if (M < SQUARE(TINY))
    {
        sprintf (ErrorMsg, "GenIndexVol: problem in integrating  fwd index volatility "
                 "bet dates %ld  %ld",IStart, IEnd);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    M = sqrt (M / Int01);

    /* Convert from q measure back to log-normal */
    if (Conv_XSpaceVol_To_Market(ParYield,
                                Vol,        
                                VolTypeFlag,
                                M,
                                Int01,
                                FwdShift,
                                QLeft,
                                QRight,
                                Amap,
                                Bmap) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* GenIndexVol_E2Q */




