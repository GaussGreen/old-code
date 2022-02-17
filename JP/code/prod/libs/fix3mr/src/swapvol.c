/****************************************************************************/
/*      Bootstrapping of swaption volatility into spot volatility.          */
/****************************************************************************/
/*      SWAPVOL.c                                                           */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/swapvol.c,v 1.19 2003/09/12 17:31:57 dfung Exp $
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


    //printf("\n%ld", SwapSt);
    //fflush(stdout);

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
    ZerotoS = ZeroPrice(SwapSt,	        /* (I) Discount bond maturity     */
                        BaseDate,	    /* (I) Value Date                 */
                        NbZero,		    /* (I) Number of zeros            */
                        ZeroDate,       /* (I) Zero maturity dates        */
                        Zero);   		/* (I) Zero rates                 */
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
        ZerotoCpn = ZeroPrice(CpnPmtDate,	/* (I) Discount bond maturity     */
                              BaseDate,	    /* (I) Value Date                 */
                              NbZero,		/* (I) Number of zeros            */
                              ZeroDate,     /* (I) Zero maturity dates        */
                              Zero);   		/* (I) Zero rates                 */
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

    //printf(":  B = %15.10f", B[0]);
    //fflush(stdout);

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
       double  BetaBmk[3][MAXNBDATE], /* (O) Mr on all benchmark intvals     */
       long    *BmkDate,              /* (O) Bmk dates                       */
       int     *NbBmkMr,              /* (O) Nb of mr bmk dates              */
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
       int     NbMr,                  /* (I) Number of mr dates              */
       long    *MrDate,               /* (I) MR dates                        */
       double  Beta[3][MAXNBDATE],    /* (I) Mean reversions                 */
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
    double  atmPr;                /* Market price of atm option             */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, k, p, q;
    int     status = FAILURE;     /* Error status = FAILURE initially       */

    char    ErrorMsg[MAXBUFF];    /* Error message                          */

    double  BetaL[3];
    int     mrIdx, prevMrIdx, lastMrIdx;
    int     lastMrFlag;
    double  BetaAve[3];
    int     NbAweight;
    int     idx;
	int     startIdx;
	long    endInt;

    for (k = 0; k < NbFactor; k++)
    {
        BetaAve[k] = 0.;
        for (i = 0; i < NbMr; i++)
            BetaAve[k] += Beta[k][i];

        BetaAve[k] /= NbMr;
    }
    
    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

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

        for (i = 0; i < NbMr; i++)                                                    
        {       
            /* record mr on this bmark interval */
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][i] = Beta[k][i];
            }
        }

        return (SUCCESS);
    }

    
	startIdx = 0;

    for (i = 0; i < NbVol; i++)
    {
        /* Only used vol points */
        if (!VolUsed[i]) continue;

        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

		/* find average mean-reversion for the swap
		   in each factor */
		while ((startIdx < NbMr)&&(MrDate[startIdx] < SwapSt[i]))
			startIdx++;

		if (startIdx >= NbMr - 1)
		{
			for (k = 0; k < NbFactor; k++)
			    BetaAve[k] = Beta[k][NbMr - 1];
		}
		else
		{
			if (MrDate[startIdx] >= SwapMat[i])
			{
				for (k = 0; k < NbFactor; k++)
			        BetaAve[k] = Beta[k][startIdx];
			}
			else
			{
                for (k = 0; k < NbFactor; k++)
		            BetaAve[k] = Beta[k][startIdx] * Daysact (SwapSt[i], MrDate[startIdx]);
			
			    for (mrIdx = startIdx+1; mrIdx < NbMr; mrIdx++)
				{
				    endInt = MIN(MrDate[mrIdx], SwapMat[i]);
				    for (k = 0; k < NbFactor; k++)
				        BetaAve[k] += Beta[k][mrIdx] * Daysact (MrDate[mrIdx-1], endInt);

				    if (endInt == SwapMat[i])
					{
					    break;
					}
				}
                if (endInt < SwapMat[i])
                /* SwapMat[i] is beyond the last MrDate */
				{
				    for (k = 0; k < NbFactor; k++)
				        BetaAve[k] += Beta[k][NbMr - 1] * Daysact (endInt, SwapMat[i]);
				}

			    /* average: divide by swap tenor */
			    for (k = 0; k < NbFactor; k++)
			        BetaAve[k] /= Daysact(SwapSt[i], SwapMat[i]);
			}
		}

        /* B factor */
        /* if (BFactor (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        NbFactor,
                        BetaAve,
                        Bbq,
                        VolNorm,
                        VolLogn,
                        NbZero,
                        Zero,
                        ZeroDate,
                        ZeroBaseDate) == FAILURE)
        {
            goto RETURN;
        } */

		if (BFactor_New (B[i],
                         SwapSt[i],
                         SwapMat[i],
                         DCC,
                         Freq,
                         NbFactor,
                         NbMr,
                         MrDate,
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
    mrIdx = 0;
    prevMrIdx = -1;
    lastMrIdx = NbMr - 1;
    lastMrFlag = FALSE;

    for (i = 0; i < NbVol; i++)                                                    
    {       
        if (!VolUsed[i]) continue;

        T = VolT[i];
        t = ((LastCalibIdx == -1) ? VolT[i] : (VolT[i] - VolT[LastCalibIdx]));

        /* determine mr on the current interval.
           here we assume that all MrDates are VolDates! */
        if (lastMrFlag == FALSE)
        {
            while (MrDate[mrIdx] < VolDate[i])
            {
                mrIdx++;
                if (mrIdx == NbMr)
                {
                    mrIdx = lastMrIdx;
                    lastMrFlag = TRUE;
                    break;
                }
            }
        }

        if (mrIdx != prevMrIdx) /* this is always true for i = 0 */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaL[k] = Beta[k][mrIdx];
            }
        }

        /* record mr on this bmark interval */
        for (k = 0; k < NbFactor; k++)
        {
            BetaBmk[k][i] = BetaL[k];
        }
        
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

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * BetaL[0] * t);
        
        M = D[0][0] * D[0][0] * ExpDecay (2. * BetaL[0], t) * t;
        //printf("\ny = %15.10f, M = %15.10f", y, M);
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * BetaL[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(BetaL[0] + BetaL[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * BetaL[1] * t);
            
            M += D[1][0] * D[1][0] * ExpDecay (2. * BetaL[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * ExpDecay (BetaL[0] + BetaL[1], t) * t;
            
            M += D[1][1] * D[1][1] * ExpDecay (2. * BetaL[1], t) * t;
        }
        
        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * BetaL[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(BetaL[0] + BetaL[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * BetaL[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * BetaL[2] * t);
            
            M += D[2][0] * D[2][0] * ExpDecay (2. * BetaL[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * ExpDecay (BetaL[0] + BetaL[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * ExpDecay (BetaL[1] + BetaL[2], t) * t;
            
            M += D[2][1] * D[2][1] * ExpDecay (2. * BetaL[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * ExpDecay (BetaL[1] + BetaL[2], t) * t;
            
            M += D[2][2] * D[2][2] * ExpDecay (2. * BetaL[2], t) * t;
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
        prevMrIdx = mrIdx;

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaL[0] * t);
        L[0][0] += lambda * lambda * ExpDecay (2. * BetaL[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaL[1] * t);
            L[1][1] += lambda * lambda * ExpDecay (2. * BetaL[1], t) * t;
            L[0][1] *= exp (-(BetaL[0] + BetaL[1]) * t);
            L[0][1] += 2.*lambda*lambda * ExpDecay ((BetaL[0]+BetaL[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaL[2] * t);
            L[2][2] += lambda * lambda * ExpDecay (2. * BetaL[2], t) * t;
            L[0][2] *= exp (-(BetaL[0] + BetaL[2]) * t);
            L[0][2] += 2.*lambda*lambda * ExpDecay ((BetaL[0]+BetaL[2]), t) * t;
            L[1][2] *= exp (-(BetaL[1] + BetaL[2]) * t);
            L[1][2] += 2.*lambda*lambda * ExpDecay ((BetaL[1]+BetaL[2]), t) * t;
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

        /* determine mr on the current interval.
           here we assume that all MrDates are VolDates! */
        if (lastMrFlag == FALSE)
        {
            while (MrDate[mrIdx] < VolDate[i])
            {
                mrIdx++;
                if (mrIdx == NbMr)
                {
                    mrIdx = lastMrIdx;
                    lastMrFlag = TRUE;
                    break;
                }
            }
        }

        /* record mr on this bmark interval */
        for (i = 0; i < NbFactor; i++)
        {
            BetaBmk[i][k] = Beta[i][mrIdx];
        }
    }  /* for k */

    for (i = 0; i < NbVol; i++)
        BmkDate[i] = VolDate[i];

    /* record mr, Aweight, and BmkDate after the last vol date */
    if (MrDate[mrIdx] > VolDate[NbVol - 1])
    {
        for (idx = mrIdx; idx < NbMr; idx++)
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol + idx - mrIdx] = Beta[k][idx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol + idx - mrIdx] = Aweight[k][NbVol - 1];
            }

            BmkDate[NbVol + idx - mrIdx] = MrDate[idx];
        }

        for (k = 0; k < NbFactor; k++)
        {
            BetaBmk[k][NbVol + NbMr - mrIdx] = Beta[k][lastMrIdx];
        }

        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][NbVol + NbMr - mrIdx] = Aweight[k][NbVol - 1];
        }

        *NbBmkMr = NbVol + NbMr - mrIdx;
    }
    else
        if (MrDate[mrIdx] == VolDate[NbVol - 1])
        {
            for (idx = mrIdx + 1; idx < NbMr; idx++)
            {
                for (k = 0; k < NbFactor; k++)
                {
                    BetaBmk[k][NbVol + idx - mrIdx - 1] = Beta[k][idx];
                }

                for (k = 0; k < NbAweight; k++)
                {
                    Aweight[k][NbVol + idx - mrIdx - 1] = Aweight[k][NbVol - 1];
                }

                BmkDate[NbVol + idx - mrIdx - 1] = MrDate[idx];
            }

            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol + NbMr - mrIdx - 1] = Beta[k][lastMrIdx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol + NbMr - mrIdx - 1] = Aweight[k][NbVol - 1];
            }

            *NbBmkMr = NbVol + NbMr - mrIdx - 1;
        }
        else /* (MrDate[mrIdx] < VolDate[NbVol - 1]), so in particular  mrIdx == lastMrIdx */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol] = Beta[k][lastMrIdx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol] = Aweight[k][NbVol - 1];
            }

            *NbBmkMr = NbVol;
        }

    

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
       double  **BetaTD,               /* (O) Time-dependent mr in the tree */
       int     *NbMrInt,               /* (O) Number of MR intervals        */
       long    *MrSwitchT,             /* (O) Array of TP when MR changes   */
       double  **BetaInt,              /* (O) MR on each MR intval          */
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       int     NbFactor,               /* (I) Number of factors             */
       int     NbMr,                   /* (I) Number of mr dates           */
       long    *MrDate,                /* (I) MR dates                     */
       int     NbBmk,
       long    *BmkDate,
       double  Beta[3][MAXNBDATE],     /* (I) Mean reversions on bmark intvals */
       double  Bbq,                    /* (I) Backbone parameter            */
       double  VolNorm,                /* (I) Normal volatility in backbone */
       double  VolLogn,                /* (I) Lognorm volatility in backbone*/
       int     CalibFlag,              /* (I) Index calibration flag        */
       double  *FwdRate,               /* (I) One period forward rate       */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP,                   /* (I) Total number of time points   */
       long    ZeroBaseDate)           /* (I) Zero curve base date         */
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

    double  timeFrac;
    long    BaseDateL;
    long    BmkDateL[MAXNBDATE];
    int     NbBmkL;
    int     idx;
    double  prevBeta[3]; 
    int     sameMR;
    double  AweightL[6][MAXNBDATE];

    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    if (CalibFlag == FALSE)
    /* no vol info */
    {
        BaseDateL = ZeroBaseDate;
        for (i = 0; i < NbMr; i++)
        {
            BmkDateL[i] = MrDate[i];
            for (k = 0; k < NbAweight; k++)
                AweightL[k][i] = Aweight[k][0];
        }
        NbBmkL = NbMr;
    }
    else
    {
        BaseDateL = VolBaseDate;
        for (i = 0; i < NbBmk; i++)
        {
            BmkDateL[i] = BmkDate[i];
            for (k = 0; k < NbAweight; k++)
                AweightL[k][i] = Aweight[k][i];
        }
        
        NbBmkL = NbBmk;
    }


    for (j = 0; j < NbBmkL; j++)
    {       
        VolT[j] = Daysact (BaseDateL, BmkDateL[j]) / 365.;
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

        while ((j < NbBmkL - 1) && (T > VolT[j] + ERROR))
            j++;

        //printf("\ni =  %d:  j = %d;  ", i, j);

        PrevVolT = ((j == 0)? 0. : VolT[j-1]);

        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        if (t >= PrevVolT)
        {
            //printf("  No interp:");
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = Beta[k][j];
            }

            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][j];
            }
        }
        /* The time step is straddling two buckets: */
        /* we interp mean-reversion and             */
        /* recalculate the covariance matrix.       */
        else
        {
            //printf("  Interp:");
            /* interp mean-reversion */
            timeFrac = (T - PrevVolT)/(T - t); /* by construction, T-PrevVolT is > ERROR */
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = (1. - timeFrac) * Beta[k][j-1] + timeFrac * Beta[k][j];
            }

            //printf(" beta = %20.15lf ", BetaTD[0][i]);
            
            /* recalc covar matrix */
            Cov[0][0]  = AweightL[0][j-1] * AweightL[0][j-1] * ExpDecay (2. * Beta[0][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
            Cov[0][0] += AweightL[0][j]   * AweightL[0][j]   * ExpDecay (2. * Beta[0][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[0][j] * (VolT[j-1] - t));
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = AweightL[0][j-1] * AweightL[1][j-1] * ExpDecay (Beta[0][j-1] + Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][0] += AweightL[0][j]   * AweightL[1][j]   * ExpDecay (Beta[0][j] + Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0][j]+Beta[1][j])*(VolT[j-1] - t));
            
                Cov[1][1]  = AweightL[1][j-1] * AweightL[1][j-1] * ExpDecay (2. * Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += AweightL[1][j]   * AweightL[1][j]   * ExpDecay (2. * Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1][j] * (VolT[j-1] - t));
                Cov[1][1] += AweightL[2][j-1] * AweightL[2][j-1] * ExpDecay (2. * Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += AweightL[2][j]   * AweightL[2][j]   * ExpDecay (2. * Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[1][j] * (VolT[j-1] - t));
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = AweightL[0][j-1] * AweightL[3][j-1] * ExpDecay (Beta[0][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][0] += AweightL[0][j]   * AweightL[3][j]   * ExpDecay (Beta[0][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[0][j]+Beta[2][j])*(VolT[j-1] - t));
            
                Cov[2][1]  = AweightL[1][j-1] * AweightL[3][j-1] * ExpDecay (Beta[1][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += AweightL[1][j]   * AweightL[3][j]   * ExpDecay (Beta[1][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1][j]+Beta[2][j])*(VolT[j-1] - t));
                Cov[2][1] += AweightL[2][j-1] * AweightL[4][j-1] * ExpDecay (Beta[1][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += AweightL[2][j]   * AweightL[4][j]   * ExpDecay (Beta[1][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-(Beta[1][j]+Beta[2][j])*(VolT[j-1] - t));

                Cov[2][2]  = AweightL[3][j-1] * AweightL[3][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += AweightL[3][j]   * AweightL[3][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2][j] * (VolT[j-1] - t));
                Cov[2][2] += AweightL[4][j-1] * AweightL[4][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += AweightL[4][j]   * AweightL[4][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2][j] * (VolT[j-1] - t));
                Cov[2][2] += AweightL[5][j-1] * AweightL[5][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += AweightL[5][j]   * AweightL[5][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]) * exp (-2. * Beta[2][j] * (VolT[j-1] - t));
            }            

            x = Cov[0][0] / ExpDecay (2. * BetaTD[0][i], T - t) / (T - t);

            //printf(" x = %20.15lf; ", x); 

            if (x < TINY)
            {
                sprintf (ErrorMsg, "Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                DR_Error (ErrorMsg);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);
            //printf(" aweight = %20.15lf; ", AweightCurve[0][i]);
            //fflush(stdout);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / ExpDecay (BetaTD[0][i] + BetaTD[1][i], T - t) / (T - t);

                x = Cov[1][1] / ExpDecay (2. * BetaTD[1][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

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
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / ExpDecay (BetaTD[0][i] + BetaTD[2][i], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / ExpDecay (BetaTD[1][i] + BetaTD[2][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / ExpDecay (2. * BetaTD[2][i], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

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
    idx = 0;
    for (k = 0; k < NbFactor; k++)
    {
        prevBeta[k] = -999.;
    }
    for (i = 0; i < NbTP; i++)                                           
    {
        /* Note period i volatility corresponds to period i+1 rate */
        for (k = 0; k < NbAweight; k++)
        {
            l = (k > 0) + (k > 2);

            BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                     (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*ExpDecay(BetaTD[l][i], Length[i+1]);
            //printf("\nadj aweight[%d] = %20.15lf; ", i, AweightCurve[0][i]);
            //fflush(stdout);
        }
        
        sameMR = TRUE;
        for (k = 0; k < NbFactor; k++)
        {
            sameMR = sameMR &&(fabs(BetaTD[k][i] - prevBeta[k]) < TINY);   
        } 

        if (!sameMR)
        {
            MrSwitchT[idx] = i;
            for (k = 0; k < NbFactor; k++)
            {
                BetaInt[k][idx] = BetaTD[k][i];
                prevBeta[k] = BetaInt[k][idx];
            }
            idx++;
        } /* if */
    } /* for i */

    *NbMrInt = idx;

        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Interp_SpotVol */

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
          int     NbMr,                  /* (I) Number of mr dates           */
          long    *MrDate,               /* (I) MR dates                     */
          double  Beta[3][MAXNBDATE],    /* (I) Mean reversions betw mr dates */
          double  BetaBmk[3][MAXNBDATE],    /* (I) Mean reversions betw all bmk dates */
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

    double  AweightL[6][MAXNBDATE];
    long    BaseDateL;
    long    BmkDateL[MAXNBDATE];
    int     NbBmk;
    double  TtoSt, PrevT;
    int     k;
	int     startIdx, mrIdx;
	long    endInt;

    double  BetaAve[3];


    startIdx = 0;
	/* find average mean-reversion for the swap
	   in each factor */
	while ((startIdx < NbMr)&&(MrDate[startIdx] < SwapSt))
		startIdx++;

	if (startIdx >= NbMr - 1)
	{
		for (k = 0; k < NbFactor; k++)
		    BetaAve[k] = Beta[k][NbMr - 1];
	}
	else
	{
		if (MrDate[startIdx] >= SwapMat)
		{
			for (k = 0; k < NbFactor; k++)
		        BetaAve[k] = Beta[k][startIdx];
		}
		else
		{
               for (k = 0; k < NbFactor; k++)
	            BetaAve[k] = Beta[k][startIdx] * Daysact (SwapSt, MrDate[startIdx]);
		
		    for (mrIdx = startIdx+1; mrIdx < NbMr; mrIdx++)
			{
			    endInt = MIN(MrDate[mrIdx], SwapMat);
			    for (k = 0; k < NbFactor; k++)
			        BetaAve[k] += Beta[k][mrIdx] * Daysact (MrDate[mrIdx-1], endInt);
				
				if (endInt == SwapMat)
				{
				    break;
				}
			}
            
			if (endInt < SwapMat)
               /* SwapMat[i] is beyond the last MrDate */
			{
			    for (k = 0; k < NbFactor; k++)
			        BetaAve[k] += Beta[k][NbMr - 1] * Daysact (endInt, SwapMat);
			}
			    /* average: divide by swap tenor */
		    for (k = 0; k < NbFactor; k++)
		        BetaAve[k] /= Daysact(SwapSt, SwapMat);
		}
	}
 

    /* Find B coefficients */
    /* if (BFactor (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 NbFactor,
                 BetaAve,     
                 Bbq,
                 VolNorm,
                 VolLogn,
                 NbZero,   
                 Zero,     
                 ZeroDate, 
                 ZeroBaseDate) == FAILURE)
    {
        goto RETURN;
    } */

    if (BFactor_New (B,
                         SwapSt,
                         SwapMat,
                         DCC,
                         Freq,
                         NbFactor,
                         NbMr,
                         MrDate,
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

    /* printf("\n%ld, %ld, %20.15lf, ", SwapSt, SwapMat, B[0]);
        fflush(stdout); */

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
    /* no vol info */
    {
        BaseDateL = ZeroBaseDate;
        for (i = 0; i < NbMr; i++)
        {
            BmkDateL[i] = MrDate[i];
            for (k = 0; k < 6; k++)
                AweightL[k][i] = Aweight[k][0];
        }
        NbBmk = NbMr;
    }
    else
    {
        BaseDateL = VolBaseDate;
        for (i = 0; i < NbVol; i++)
        {
            BmkDateL[i] = VolDate[i];
            for (k = 0; k < 6; k++)
                AweightL[k][i] = Aweight[k][i];
        }
        
        if (MrDate[NbMr - 2] == VolDate[NbVol - 1])
        /* possible change of mr at the last vol date */
        {
            NbBmk = NbVol + 1;
            BmkDateL[NbVol] = MrDate[NbMr - 1];
            for (k = 0; k < 6; k++)
                AweightL[k][NbVol] = Aweight[k][NbVol - 1];
        }
        else
        /* flat extrapolation of mr over the last vol date */
        {
            NbBmk = NbVol;
        }
    }

    /* Avoid spot index expiration */
    if (SwapSt <= BaseDateL) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }

    TtoSt = Daysact (BaseDateL, SwapSt) / 365.;
    
    for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++) 
            L[i][j] = 0.;

    PrevT = 0;
    for (i = 0; i < NbBmk; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate = (i < NbBmk-1) ? MIN(SwapSt, BmkDateL[i]) : SwapSt;

        VolT[i] = Daysact (BaseDateL, EndDate) / 365.;
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));        

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaBmk[0][i] * t);
        L[0][0] += AweightL[0][i] * AweightL[0][i] * ExpDecay (2. * BetaBmk[0][i], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaBmk[1][i] * t);
            L[1][1] += AweightL[1][i] * AweightL[1][i] * ExpDecay (2. * BetaBmk[1][i], t) * t;
            L[1][1] += AweightL[2][i] * AweightL[2][i] * ExpDecay (2. * BetaBmk[1][i], t) * t;
            L[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);
            L[0][1] += 2. * AweightL[0][i] * AweightL[1][i] * ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L[2][2] += AweightL[3][i] * AweightL[3][i] * ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[2][2] += AweightL[4][i] * AweightL[4][i] * ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[2][2] += AweightL[5][i] * AweightL[5][i] * ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);
            L[0][2] += 2. * AweightL[0][i] * AweightL[3][i] * ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t;
            L[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);
            L[1][2] += 2. * AweightL[1][i] * AweightL[3][i] * ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
            L[1][2] += 2. * AweightL[2][i] * AweightL[4][i] * ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
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
            sprintf (ErrorMsg, "IndexVol: problem in BS2Q price at %ld ",
                     SwapSt);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        M    = ImpVol_BS2Q (ParYield,ParYield,T,atmPr,'C',1.,1.,0.,M);
        if (M < 0.0)
        {
            sprintf (ErrorMsg, "IndexVol: problem in BS2Q implied vol at %ld ",
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


/*****  Filtered_SpotVol    **************************************************/
/*
*       Bootstrap input volatility curve to extract filtered spot vols.
*       Spot vols are not modified provided the ratio of the current spot vol
*       to the spot vol in the first bucket is greater than MIN_SPOT_VOL_RATIO
*       or less than 1/MIN_SPOT_VOL_RATIO. Otherwise the spot vols are smoothly
*       capped and floored, with the cap and floor levels parameterized by
*       SPOT_VOL_FILTER_AMOUNT.
*/
int    Filtered_SpotVol(
       double  Aweight[6][MAXNBDATE], /* (O) Spot vol curve of each factor   */
       double  BetaBmk[3][MAXNBDATE], /* (O) Mr on all benchmark intvals     */
       long    *BmkDate,              /* (O) Bmk dates                       */
       int     *NbBmkMr,              /* (O) Nb of mr bmk dates              */
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
       int     NbMr,                  /* (I) Number of mr dates              */
       long    *MrDate,               /* (I) MR dates                        */
       double  Beta[3][MAXNBDATE],    /* (I) Mean reversions                 */
       double  *Rho,                  /* (I) Correlation between factors     */
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
	double  Bn[MAXNBDATE][3];      /* B in Christian memo                    */
    double  L[3][3];              /* Integrals of lambda factors            */
    double  D[3][3];              /* aweight*B                              */
    double  aw[3][3];             /* Historical aweights                    */
    double  Anorm=0;              /* Norm of historical aweights            */ 
    double  M;                    /* Matrix for lambda system               */
    double  y;                    /* Vector for lambda system               */
    double  atmPr;                /* Market price of atm option             */
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

    char    ErrorMsg[MAXBUFF];    /* Error message                          */

    double  BetaL[3];
    int     mrIdx, prevMrIdx, lastMrIdx;
    int     lastMrFlag;

    double  BetaAve[3];
    int     NbAweight;
    int     idx;
	int     startIdx;
	long    endInt;

    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

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

        for (i = 0; i < NbMr; i++)                                                    
        {       
            /* record mr on this bmark interval */
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][i] = Beta[k][i];
            }
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

    /* B-coefficients and par yields */
	startIdx = 0;

    for (i = 0; i < NbVolUsed; i++)
    {
        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

		/* find average mean-reversion for the swap
		   in each factor */
		while ((startIdx < NbMr)&&(MrDate[startIdx] < SwapSt[i]))
			startIdx++;

		if (startIdx >= NbMr - 1)
		{
			for (k = 0; k < NbFactor; k++)
			    BetaAve[k] = Beta[k][NbMr - 1];
		}
		else
		{
			if (MrDate[startIdx] >= SwapMat[i])
			{
				for (k = 0; k < NbFactor; k++)
			        BetaAve[k] = Beta[k][startIdx];
			}
			else
			{
                for (k = 0; k < NbFactor; k++)
		            BetaAve[k] = Beta[k][startIdx] * Daysact (SwapSt[i], MrDate[startIdx]);
			
			    for (mrIdx = startIdx+1; mrIdx < NbMr; mrIdx++)
				{
				    endInt = MIN(MrDate[mrIdx], SwapMat[i]);
				    for (k = 0; k < NbFactor; k++)
				        BetaAve[k] += Beta[k][mrIdx] * Daysact (MrDate[mrIdx-1], endInt);

				    if (endInt == SwapMat[i])
					{
					    break;
					}
				}
                if (endInt < SwapMat[i])
                /* SwapMat[i] is beyond the last MrDate */
				{
				    for (k = 0; k < NbFactor; k++)
				        BetaAve[k] += Beta[k][NbMr - 1] * Daysact (endInt, SwapMat[i]);
				}

			    /* average: divide by swap tenor */
			    for (k = 0; k < NbFactor; k++)
			        BetaAve[k] /= Daysact(SwapSt[i], SwapMat[i]);
			}
		}

        /* B factor */
       /* if (BFactor (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        NbFactor,
                        BetaAve,
                        Bbq,
                        VolNorm,
                        VolLogn,
                        NbZero,
                        Zero,
                        ZeroDate,
                        ZeroBaseDate) == FAILURE)
        {
            goto RETURN;
	   } */ 

		if (BFactor_New (B[i],
                         SwapSt[i],
                         SwapMat[i],
                         DCC,
                         Freq,
                         NbFactor,
                         NbMr,
                         MrDate,
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

    mrIdx = 0;
    prevMrIdx = -1;
    lastMrIdx = NbMr - 1;

    lastMrFlag = FALSE;

    for (i = 0; i < NbVolUsed; i++)
    {
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));

        /* determine mr on the current interval.
           here we assume that all MrDates are VolDates! */
        if (lastMrFlag == FALSE)
        {
            while (MrDate[mrIdx] < VolDate[i])
            {
                mrIdx++;
                if (mrIdx == NbMr)
                {
                    mrIdx = lastMrIdx;
                    lastMrFlag = TRUE;
                    break;
                }
            }
        }

        if (mrIdx != prevMrIdx) /* this is always true for i = 0 */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaL[k] = Beta[k][mrIdx];
            }
        }

        /* record mr on this bmark interval */
        for (k = 0; k < NbFactor; k++)
        {
            BetaBmk[k][i] = BetaL[k];
        }
        
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
            if (atmPr < 0.0)
            {
                sprintf (ErrorMsg, "Filter spot vol: problem in BS2Q price at %ld ",
                         VolDate[i]);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
                        
            y    = ImpVol_BS2Q (pY[i],pY[i],T,atmPr,'C',QLeft,QRight,FwdSh,Vol[i]);
            if (y < 0.0)
            {
                sprintf (ErrorMsg, "Filter spot vol: problem in BS2Q implied vol at %ld ",
                         VolDate[i]);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            y = T * SQUARE (y);
        }

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * BetaL[0] * t);
        
        M = D[0][0] * D[0][0] * ExpDecay (2. * BetaL[0], t) * t;
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][1] * exp (-2. * BetaL[1] * t);
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(BetaL[0] + BetaL[1]) * t);
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * BetaL[1] * t);
            
            M += D[1][0] * D[1][0] * ExpDecay (2. * BetaL[1], t) * t;
            M += 2. * D[0][0] * D[1][0] * ExpDecay (BetaL[0] + BetaL[1], t) * t;
            
            M += D[1][1] * D[1][1] * ExpDecay (2. * BetaL[1], t) * t;
        }

        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][2] * exp (-2. * BetaL[2] * t);
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(BetaL[0] + BetaL[2]) * t);
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t);
            y -= D[2][1] * D[2][1] * L[2][2] * exp (-2. * BetaL[2] * t);
            y -= D[1][1] * D[2][1] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t);
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * BetaL[2] * t);
            
            M += D[2][0] * D[2][0] * ExpDecay (2. * BetaL[2], t) * t;
            M += 2. * D[0][0] * D[2][0] * ExpDecay (BetaL[0] + BetaL[2], t) * t;
            M += 2. * D[1][0] * D[2][0] * ExpDecay (BetaL[1] + BetaL[2], t) * t;
            
            M += D[2][1] * D[2][1] * ExpDecay (2. * BetaL[2], t) * t;
            M += 2. * D[1][1] * D[2][1] * ExpDecay (BetaL[1] + BetaL[2], t) * t;
            
            M += D[2][2] * D[2][2] * ExpDecay (2. * BetaL[2], t) * t;
        }

        /* M must be stictly positive, and since it is not a function   */
        /* of spot vol, trap the error now rather than later in SpotVol */
        if (M < TINY)
        {
            sprintf (ErrorMsg, "Filtered spot vol: problem in bootstrapping %ld "
                                "volatility (negative variance)", VolDate[i]);
            DR_Error (ErrorMsg);
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

        /* reset mr index */
        prevMrIdx = mrIdx;
        
        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaL[0] * t);
        L[0][0] += lambda * lambda * ExpDecay (2. * BetaL[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaL[1] * t);
            L[1][1] += lambda * lambda * ExpDecay (2. * BetaL[1], t) * t;
            L[0][1] *= exp (-(BetaL[0] + BetaL[1]) * t);
            L[0][1] += 2.*lambda*lambda * ExpDecay ((BetaL[0]+BetaL[1]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaL[2] * t);
            L[2][2] += lambda * lambda * ExpDecay (2. * BetaL[2], t) * t;
            L[0][2] *= exp (-(BetaL[0] + BetaL[2]) * t);
            L[0][2] += 2.*lambda*lambda * ExpDecay ((BetaL[0]+BetaL[2]), t) * t;
            L[1][2] *= exp (-(BetaL[1] + BetaL[2]) * t);
            L[1][2] += 2.*lambda*lambda * ExpDecay ((BetaL[1]+BetaL[2]), t) * t;
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

        /* determine mr on the current interval.
           here we assume that all MrDates are VolDates! */
        if (lastMrFlag == FALSE)
        {
            while (MrDate[mrIdx] < VolDate[i])
            {
                mrIdx++;
                if (mrIdx == NbMr)
                {
                    mrIdx = lastMrIdx;
                    lastMrFlag = TRUE;
                    break;
                }
            }
        }

        /* record mr on this bmark interval */
        for (i = 0; i < NbFactor; i++)
        {
            BetaBmk[i][k] = Beta[i][mrIdx];
        }
    }  /* for k */

       for (i = 0; i < NbVol; i++)
        BmkDate[i] = VolDate[i];

    /* record mr, Aweight, and BmkDate after the last vol date */
    if (MrDate[mrIdx] > VolDate[NbVol - 1])
    {
        for (idx = mrIdx; idx < NbMr; idx++)
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol + idx - mrIdx] = Beta[k][idx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol + idx - mrIdx] = Aweight[k][NbVol - 1];
            }

            BmkDate[NbVol + idx - mrIdx] = MrDate[idx];
        }

        for (k = 0; k < NbFactor; k++)
        {
            BetaBmk[k][NbVol + NbMr - mrIdx] = Beta[k][lastMrIdx];
        }

        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][NbVol + NbMr - mrIdx] = Aweight[k][NbVol - 1];
        }

        *NbBmkMr = NbVol + NbMr - mrIdx;
    }
    else
        if (MrDate[mrIdx] == VolDate[NbVol - 1])
        {
            for (idx = mrIdx + 1; idx < NbMr; idx++)
            {
                for (k = 0; k < NbFactor; k++)
                {
                    BetaBmk[k][NbVol + idx - mrIdx - 1] = Beta[k][idx];
                }

                for (k = 0; k < NbAweight; k++)
                {
                    Aweight[k][NbVol + idx - mrIdx - 1] = Aweight[k][NbVol - 1];
                }

                BmkDate[NbVol + idx - mrIdx - 1] = MrDate[idx];
            }

            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol + NbMr - mrIdx - 1] = Beta[k][lastMrIdx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol + NbMr - mrIdx - 1] = Aweight[k][NbVol - 1];
            }

            *NbBmkMr = NbVol + NbMr - mrIdx - 1;
        }
        else /* (MrDate[mrIdx] < VolDate[NbVol - 1]), so in particular  mrIdx == lastMrIdx */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaBmk[k][NbVol] = Beta[k][lastMrIdx];
            }

            for (k = 0; k < NbAweight; k++)
            {
                Aweight[k][NbVol] = Aweight[k][NbVol - 1];
            }

            *NbBmkMr = NbVol;
        }


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Filtered_SpotVol */

/*****  Interp_SpotVol_New  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Interp_SpotVol_New (
       double  **AweightCurve,         /* (O) Factors aweight curves        */
       double  **BetaTD,               /* (O) Time-dependent mr in the tree */
       int     *NbMrInt,               /* (O) Number of MR intervals        */
       long    *MrSwitchT,             /* (O) Array of TP when MR changes   */
       double  **BetaInt,              /* (O) MR on each MR intval          */
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       int     NbFactor,               /* (I) Number of factors             */
       int     NbMr,                   /* (I) Number of mr dates           */
       long    *MrDate,                /* (I) MR dates                     */
       int     NbBmk,
       long    *BmkDate,
       double  Beta[3][MAXNBDATE],     /* (I) Mean reversions on bmark intvals */
       double  Bbq,                    /* (I) Backbone parameter            */
       double  VolNorm,                /* (I) Normal volatility in backbone */
       double  VolLogn,                /* (I) Lognorm volatility in backbone*/
       int     CalibFlag,              /* (I) Index calibration flag        */
       double  *FwdRate,               /* (I) One period forward rate       */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP,                   /* (I) Total number of time points   */
       long    ZeroBaseDate)           /* (I) Zero curve base date         */
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

    double  timeFrac;
    long    BaseDateL;
    long    BmkDateL[MAXNBDATE];
    int     NbBmkL;
    int     idx;
    double  prevBeta[3]; 
    int     sameMR;
    double  AweightL[6][MAXNBDATE];

    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    if (CalibFlag == FALSE)
    /* no vol info */
    {
        BaseDateL = ZeroBaseDate;
        for (i = 0; i < NbMr; i++)
        {
            BmkDateL[i] = MrDate[i];
            for (k = 0; k < NbAweight; k++)
                AweightL[k][i] = Aweight[k][0];
        }
        NbBmkL = NbMr;
    }
    else
    {
        BaseDateL = VolBaseDate;
        for (i = 0; i < NbBmk; i++)
        {
            BmkDateL[i] = BmkDate[i];
            for (k = 0; k < NbAweight; k++)
                AweightL[k][i] = Aweight[k][i];
        }
        
        NbBmkL = NbBmk;
    }


    for (j = 0; j < NbBmkL; j++)
    {       
        VolT[j] = Daysact (BaseDateL, BmkDateL[j]) / 365.;
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

        while ((j < NbBmkL - 1) && (T > VolT[j] + ERROR))
            j++;

        //printf("\ni =  %d:  j = %d;  ", i, j);

        PrevVolT = ((j == 0)? 0. : VolT[j-1]);

        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        if (t >= PrevVolT)
        {
            //printf("  No interp:");
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = Beta[k][j];
            }

            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][j];
            }
        }
        /* The time step is straddling two buckets: */
        /* we interp mean-reversion and             */
        /* recalculate the covariance matrix.       */
        else
        {
            //printf("  Interp:");
            /* interp mean-reversion */
            timeFrac = (T - PrevVolT)/(T - t); /* by construction, T-PrevVolT is > ERROR */
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = (1. - timeFrac) * Beta[k][j-1] + timeFrac * Beta[k][j];
            }

            //printf(" beta = %20.15lf ", BetaTD[0][i]);
            
            /* recalc covar matrix */
            Cov[0][0]  = AweightL[0][j-1] * AweightL[0][j-1] * ExpDecay (2. * Beta[0][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[0][j] * (T - VolT[j-1]));
            Cov[0][0] += AweightL[0][j]   * AweightL[0][j]   * ExpDecay (2. * Beta[0][j], T - VolT[j-1]) * (T - VolT[j-1]);
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = AweightL[0][j-1] * AweightL[1][j-1] * ExpDecay (Beta[0][j-1] + Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-(Beta[0][j]+Beta[1][j])*(T - VolT[j-1]));
                Cov[1][0] += AweightL[0][j]   * AweightL[1][j]   * ExpDecay (Beta[0][j] + Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]);
            
                Cov[1][1]  = AweightL[1][j-1] * AweightL[1][j-1] * ExpDecay (2. * Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[1][j] * (T - VolT[j-1]));
                Cov[1][1] += AweightL[1][j]   * AweightL[1][j]   * ExpDecay (2. * Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]);
                Cov[1][1] += AweightL[2][j-1] * AweightL[2][j-1] * ExpDecay (2. * Beta[1][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[1][j] * (T - VolT[j-1]));
                Cov[1][1] += AweightL[2][j]   * AweightL[2][j]   * ExpDecay (2. * Beta[1][j], T - VolT[j-1]) * (T - VolT[j-1]);
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = AweightL[0][j-1] * AweightL[3][j-1] * ExpDecay (Beta[0][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-(Beta[0][j]+Beta[2][j])*(T - VolT[j-1]));
                Cov[2][0] += AweightL[0][j]   * AweightL[3][j]   * ExpDecay (Beta[0][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);
            
                Cov[2][1]  = AweightL[1][j-1] * AweightL[3][j-1] * ExpDecay (Beta[1][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-(Beta[1][j]+Beta[2][j])*(T - VolT[j-1]));
                Cov[2][1] += AweightL[1][j]   * AweightL[3][j]   * ExpDecay (Beta[1][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);
                Cov[2][1] += AweightL[2][j-1] * AweightL[4][j-1] * ExpDecay (Beta[1][j-1] + Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-(Beta[1][j]+Beta[2][j])*(T - VolT[j-1]));
                Cov[2][1] += AweightL[2][j]   * AweightL[4][j]   * ExpDecay (Beta[1][j] + Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);

                Cov[2][2]  = AweightL[3][j-1] * AweightL[3][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[2][j] * (T - VolT[j-1]));
                Cov[2][2] += AweightL[3][j]   * AweightL[3][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);
                Cov[2][2] += AweightL[4][j-1] * AweightL[4][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[2][j] * (T - VolT[j-1]));
                Cov[2][2] += AweightL[4][j]   * AweightL[4][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);
                Cov[2][2] += AweightL[5][j-1] * AweightL[5][j-1] * ExpDecay (2. * Beta[2][j-1], VolT[j-1] - t) * (VolT[j-1] - t) * exp (-2. * Beta[2][j] * (T - VolT[j-1]));
                Cov[2][2] += AweightL[5][j]   * AweightL[5][j]   * ExpDecay (2. * Beta[2][j], T - VolT[j-1]) * (T - VolT[j-1]);
            }            

            x = Cov[0][0] / ExpDecay (2. * BetaTD[0][i], T - t) / (T - t);

            //printf(" x = %20.15lf; ", x); 

            if (x < TINY)
            {
                sprintf (ErrorMsg, "Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                DR_Error (ErrorMsg);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);
            //printf(" aweight = %20.15lf; ", AweightCurve[0][i]);
            //fflush(stdout);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / ExpDecay (BetaTD[0][i] + BetaTD[1][i], T - t) / (T - t);

                x = Cov[1][1] / ExpDecay (2. * BetaTD[1][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

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
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / ExpDecay (BetaTD[0][i] + BetaTD[2][i], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / ExpDecay (BetaTD[1][i] + BetaTD[2][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / ExpDecay (2. * BetaTD[2][i], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

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
    idx = 0;
    for (k = 0; k < NbFactor; k++)
    {
        prevBeta[k] = -999.;
    }
    for (i = 0; i < NbTP; i++)                                           
    {
        /* Note period i volatility corresponds to period i+1 rate */
        for (k = 0; k < NbAweight; k++)
        {
            l = (k > 0) + (k > 2);

            BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                     (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*ExpDecay(BetaTD[l][i], Length[i+1]);
            //printf("\nadj aweight[%d] = %20.15lf; ", i, AweightCurve[0][i]);
            //fflush(stdout);
        }
        
        sameMR = TRUE;
        for (k = 0; k < NbFactor; k++)
        {
            sameMR = sameMR &&(fabs(BetaTD[k][i] - prevBeta[k]) < TINY);   
        } 

        if (!sameMR)
        {
            MrSwitchT[idx] = i;
            for (k = 0; k < NbFactor; k++)
            {
                BetaInt[k][idx] = BetaTD[k][i];
                prevBeta[k] = BetaInt[k][idx];
            }
            idx++;
        } /* if */
    } /* for i */

    *NbMrInt = idx;

        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Interp_SpotVol_New */

/*****  BFactor_New    *********************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*/
int     BFactor_New (double  *B,            /* (O) B ceoff for each factor       */
                 long    SwapSt,        /* (I) Underlying swap start         */
                 long    SwapMat,       /* (I) Underlying swap maturity      */
                 char    DCC,           /* (I) Underlying day count conv.    */
                 char    Freq,          /* (I) Underlying frequency          */
                 int     NbFactor,      /* (I) Number of factors             */
                 int     NbMr,          /* (I) Number of mr dates              */
                 long    *MrDate,       /* (I) MR dates                      */
                 double  Beta[3][MAXNBDATE], /* (I) Mean reversions               */
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
	double  C[3];                
    int     j, k;
    int     status = FAILURE;    /* Error status = FAILURE initially   */

    int     mrIdx, firstMrIdx, lastMrIdx;
    double  TtoMrDate;
    long    PrevCpnDate;
    int     lastMrFlag;
    double  avgRate;


    /* printf("\n%ld", SwapSt);
    fflush(stdout); */

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
    ZerotoS = ZeroPrice(SwapSt,	        /* (I) Discount bond maturity     */
                        BaseDate,	    /* (I) Value Date                 */
                        NbZero,		    /* (I) Number of zeros            */
                        ZeroDate,       /* (I) Zero maturity dates        */
                        Zero);   		/* (I) Zero rates                 */
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    S = Daysact (BaseDate, SwapSt) / 365.;

    PrevT = S;
    PrevZero = ZerotoS;

    Annuity = 0.;

    A[0] = A[1] = A[2] = 0.;
    B[0] = B[1] = B[2] = 0.;
	C[0] = C[1] = C[2] = 1.;

    lastMrFlag = FALSE;
    firstMrIdx = 0;
    lastMrIdx = 0;
    PrevCpnDate = SwapSt;

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

        if (lastMrFlag == FALSE)
        {
            while(MrDate[firstMrIdx] <= PrevCpnDate)
            {
                firstMrIdx++;
                if (firstMrIdx == NbMr)
                {
                    lastMrFlag = TRUE;
                    firstMrIdx = NbMr - 1;
                    break;
                }
            } /* while */
            lastMrIdx = firstMrIdx;
        } /* if */
        
        if (lastMrFlag == FALSE)
        {
            while(MrDate[lastMrIdx] < CpnPmtDate)
            {
                lastMrIdx++;
                if (lastMrIdx == NbMr)
                {
                    lastMrIdx = NbMr - 1;
                    break;
                }
            }
        }
        
        /* Zero to current cpn */
        ZerotoCpn = ZeroPrice(CpnPmtDate,	/* (I) Discount bond maturity     */
                              BaseDate,	    /* (I) Value Date                 */
                              NbZero,		/* (I) Number of zeros            */
                              ZeroDate,     /* (I) Zero maturity dates        */
                              Zero);   		/* (I) Zero rates                 */
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;

        TtoCpn = Daysact (BaseDate, CpnPmtDate) / 365.;

        avgRate = log(PrevZero/ZerotoCpn)/(TtoCpn - PrevT);
        
        Annuity += DCCFrac * ZerotoCpn;
            
	    for (mrIdx = firstMrIdx; mrIdx < lastMrIdx; mrIdx++)
        {
			TtoMrDate = Daysact (BaseDate, MrDate[mrIdx]) / 365.;
            BbqAdj =      (Bbq) * VolLogn * avgRate * (TtoMrDate - PrevT)
                   + (1. - Bbq) * VolNorm * (TtoMrDate - PrevT);
			
            for (k = 0; k < NbFactor; k++)
            {
                A[k] += C[k] * BbqAdj
                   * ExpDecay (Beta[k][mrIdx], (TtoMrDate-PrevT));
				C[k] *= exp(-Beta[k][mrIdx]*(TtoMrDate-PrevT));
            } /* for k */

			PrevT = TtoMrDate;
		} /* for MrIdx */

        /* last mr intval: mrIdx = lastMrIdx; last date is CpnPmtDate */
        BbqAdj =      (Bbq) * VolLogn * avgRate * (TtoCpn - PrevT)
                   + (1. - Bbq) * VolNorm * (TtoCpn - PrevT);
		for (k = 0; k < NbFactor; k++)
        {
            A[k]  += C[k] * BbqAdj
                   * ExpDecay (Beta[k][mrIdx], (TtoCpn-PrevT));
			C[k] *= exp(-Beta[k][mrIdx]*(TtoCpn-PrevT));
		}
        

        for (k = 0; k < NbFactor; k++)
        {
		    B[k]  += A[k] * DCCFrac * ZerotoCpn;
		}

        PrevT    = TtoCpn;
        PrevZero = ZerotoCpn;
        PrevCpnDate = CpnPmtDate;

    }  /* for j */

    FwdYield = (ZerotoS - ZerotoCpn) / Annuity;

    for (k = 0; k < NbFactor; k++)
    {
        B[k] *= FwdYield;
        B[k] += A[k] * ZerotoCpn;
        B[k] /= (ZerotoS - ZerotoCpn);
    }

    /*     printf(":  B = %15.10f", B[0]);
    fflush(stdout); */

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* BFactor_New */




