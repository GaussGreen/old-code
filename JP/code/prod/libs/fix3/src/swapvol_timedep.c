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


/*****  QInterp   *********************************************************/
/*
*       Interpolation of Q's given
*      
*/
int QInterp ( double QDateFrac,/* (I) q date as frac of years from today  */
              double *Q,       /* (O) q at qdatefrac                      */
              long   Today,    /* (I) today date                          */   
              long   *QDates,  /* (I) q  dates                            */
              double *QValues, /* (I) q  values                           */
              int    nbValues) /* (I) nb input values                     */
{
    int i;
    double QInputDatesFrac[MAXNBDATE];
    int status = FAILURE;
    
    /* compute the years frac for input q dates */
    for (i = 0; i < nbValues; i++)
        QInputDatesFrac[i] = Daysact(Today,QDates[i]) / 365.0;

    /* linearly interpolate to find Q */
    tableinterp( QDateFrac,
                 Q,
                 QInputDatesFrac,
                 QValues,
                 nbValues);
          
          

    status = SUCCESS;
    return (status);
}




/*****  Fix3_BFactor_TimeDep    *********************************************************/
/*
*       Determine the value of the B coefficient (see Vladimir's memo)
*       Time dependent mean reversions, corr and weights
*/
int     Fix3_BFactor_TimeDep (
        double*         B,            /* (O) B ceoff for each factor       */
        long            SwapSt,       /* (I) Underlying swap start         */
        long            SwapMat,      /* (I) Underlying swap maturity      */
        char            DCC,          /* (I) Underlying day count conv.    */
        char            Freq,         /* (I) Underlying frequency          */
        MKTVOL_DATA* mktvol_data,     /* (I) Volatility data               */
        T_CURVE const* crv)           /* (I) zero curve                    */
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
    int     i, j, k;
    int     status = FAILURE;    /* Error status = FAILURE initially   */


    long     mrIdx, firstMrIdx, lastMrIdx;
    double  TtoMrDate;
    long    PrevCpnDate;
    int     lastMrFlag;
    double  avgRate;
    double  C[3];

    /* temporary variables for volatility data parameters */
    int             NbFactor;     /* (I) Number of factors             */
    int             NbMr;         /* (I) Number of mean reversion dates*/
    long            *MrDate;      /* (I) MR dates                      */
    double          Beta[3][MAXNBDATE];/* (I) Mean reversions          */
    double          Bbq;          /* (I) Backbone parameter            */
    double          VolNorm;      /* (I) Normal volatility in backbone */
    double          VolLogn;      /* (I) Lognorm volatility in backbone*/

    /* set temporary variables equal to mktvol_data parameters */
    NbFactor = mktvol_data->NbFactor;
    NbMr     = mktvol_data->NbTDInp;
    MrDate   = mktvol_data->TDInpDate;
    Bbq      = mktvol_data->Bbq;
    VolNorm  = mktvol_data->VolNorm;
    VolLogn  = mktvol_data->VolLogn;
    for (i = 0; i < NbMr; i++)
    {
        Beta[0][i] = mktvol_data->BetaTD[0][i];
        if (NbFactor > 1)
            Beta[1][i] = mktvol_data->BetaTD[1][i];
        if (NbFactor > 2)
            Beta[2][i] = mktvol_data->BetaTD[2][i];
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
    ZerotoS = GetZeroPrice(SwapSt, crv);
    if (ZerotoS < 0.0) goto FREE_MEM_AND_RETURN;

    S = Daysact (crv->ValueDate, SwapSt) / 365.;

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
            DR_Error ("Could not calculate day count fraction (Fix3_BFactor)!");
            goto FREE_MEM_AND_RETURN;
        }

        /* search for first index such that  MrDate[mrIdx-1] < PrevCpnDate <= MrDate[mrIdx]
        if not already past last mean reversion date */
        if(lastMrFlag == FALSE)
            firstMrIdx = GetDLOffset (NbMr, MrDate, PrevCpnDate, CbkHIGHER);
        if (firstMrIdx == NbMr - 1)
            lastMrFlag = TRUE;
        /* if all TD dates are smaller than target date use last value */
        if (firstMrIdx < 0)
            firstMrIdx = NbMr - 1;


        /*initialize last mean reversion index to first */
        lastMrIdx = firstMrIdx;

        /* search for first index such that MrDate[mrIdx-1] < CpnPmtDate <= MrDate[mrIdx]
        if not already past last mean reversion date */
        if(lastMrFlag == FALSE)
            lastMrIdx  = GetDLOffset (NbMr, MrDate, CpnPmtDate, CbkHIGHER);
      
        /* if all TD dates are smaller than target date use last value */
        if (lastMrIdx < 0)
            lastMrIdx = NbMr - 1 ;

        
        
       
        /* Zero to current cpn */
        ZerotoCpn = GetZeroPrice(CpnPmtDate, crv);
        if (ZerotoCpn < 0.0) goto FREE_MEM_AND_RETURN;

        TtoCpn = Daysact (crv->ValueDate, CpnPmtDate) / 365.;

        avgRate = log(PrevZero/ZerotoCpn)/(TtoCpn - PrevT);
        
        Annuity += DCCFrac * ZerotoCpn;

        for (mrIdx = firstMrIdx; mrIdx < lastMrIdx; mrIdx++)
        {
            TtoMrDate = Daysact (crv->ValueDate, MrDate[mrIdx]) / 365.;
            BbqAdj =      (Bbq) * VolLogn * avgRate * (TtoMrDate - PrevT)
                   + (1. - Bbq) * VolNorm * (TtoMrDate - PrevT);
			
            for (k = 0; k < NbFactor; k++)
            {
                A[k] += C[k] * BbqAdj
                   * Fix3_ExpDecay (Beta[k][mrIdx], (TtoMrDate-PrevT));
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
                   * Fix3_ExpDecay (Beta[k][mrIdx], (TtoCpn-PrevT));
            C[k]  *= exp(-Beta[k][mrIdx]*(TtoCpn-PrevT));
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


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* Fix3_BFactor_New */




/*****  Fix3_SpotVol    *********************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Fix3_SpotVol_TimeDep   (
       MKTVOL_DATA* mktvol_data,      /* (I/O) Volatility data */
       T_CURVE const* crv)
{

    double  VolT[MAXNBDATE];      /* Expiries in years                      */
    double  StartT[MAXNBDATE];    /* Starts in years                        */
    double  TDDateT[MAXNBDATE];   /* Time dep dates in years                */
    double  t;                    /* Time between two consecutive expiries  */
    double  t0;
    double  T;                    /* Time to current expiry                 */
    double  pY[MAXNBDATE];        /* Fwd par yield                          */
    double  Annuity;              /* Fwd annuity                            */

    double  B[MAXNBDATE][3];      /* B in Christian memo                    */
    double  L[3][3];              /* Integrals of lambda factors            */
    double  D[3][3];              /* aweight*B                              */
    double  E[3][3];              /* Integrals of exp                       */
    double  aw[3][3];             /* Historical aweights                    */
    double  M;                    /* Matrix for lambda system               */
    double  y;                    /* Vector for lambda system               */
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, j, k, p, q;
    int     status = FAILURE;     /* Error status = FAILURE initially       */


    double  BetaL[3];
	double  WeightL[3];
    long    mrIdx, mrIdx_st, prevMrIdx, lastMrIdx;
    int     lastMrFlag, lastMr_stFlag;


    /* temporary variables for volatility parameters */
    long             VolBaseDate;  /* (I) Volatility curve base date        */
    int              NbVol;        /* (I) Nb of points in vol curve         */
    long const*      VolDate;      /* (I) Volatility dates                  */
    double const*    Vol;          /* (I) Vol curve                         */    
    int*             VolUsed;      /* (I/O) TRUE if vol used in calibration */
    char             Freq;         /* (I) Frequency of underlying rate      */
    char             DCC;          /* (I) Day count convention              */
    long const*      SwapSt;       /* (I) Underlying swap start             */
    long const*      SwapMat;      /* (I) Underlying swap maturity          */
    double           QLeft;        /* (I) Left Q mapping coefficient        */
    double           QRight;       /* (I) Right Q mapping coefficient       */
    double           FwdSh;        /* (I) Fwd shift mapping coefficient     */
    double           Bbq;          /* (I) Backbone parameter                */
    int              VolTypeFlag;  /* (I) Volatility type                 */
    double           VolNorm;      /* (I) Normal volatility in backbone     */
    double           VolLogn;      /* (I) Lognorm volatility in backbone    */
    int              NbFactor;     /* (I) Number of factors                 */
    int              NbTDInp;      /* (I) Number of TD input dates          */
    long             *TDInpDate;   /* (I) TD Inp dates                      */
    int              SkipFlag;     /* (I) Skip calibration failure points   */
    int              CalibFlag;    /* (I) Index calibration flag            */
    long             *BmkDate;     /* (O) Bmk dates                         */
    int              *NbBmkMr;     /* (O) Nb of mr bmk dates                */
    long             EndDate;      /* (I) End of integration bucket         */
    double       BetaBmk[3][MAXNBDATE]; /* (O) Mr on all benchmark intvals  */
    double       FactorWeight[3][MAXNBDATE];/*(I) Time-dep factor weights   */
    double       Beta[3][MAXNBDATE];/* (I) Mean reversions                  */
    double       Rho[3][MAXNBDATE]; /* (I) Correlation between factors      */
    double   Aweight[6][MAXNBDATE]; /* (O) Spot vol curve of each factor    */
    
    
    /* set temporary variables to mktvol_data parameters */
    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
 
    VolDate     = mktvol_data->VolDate;
    Vol         = mktvol_data->Vol;
    VolUsed     = mktvol_data->VolUsed;
    Freq        = mktvol_data->Freq;
    DCC         = mktvol_data->DCC;
    SwapSt      = mktvol_data->SwapSt;
    SwapMat     = mktvol_data->SwapMat;
    Bbq         = mktvol_data->Bbq;
    VolTypeFlag = mktvol_data->VolUnit;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;
    TDInpDate   = mktvol_data->TDInpDate;
    NbBmkMr     = &(mktvol_data->NbBmkMr);
    NbTDInp     = mktvol_data->NbTDInp;
    BmkDate     = mktvol_data->BmkDate;
    

    for (i = 0; i < NbTDInp; i++)
    {
        FactorWeight[0][i] = mktvol_data->AlphaTD[0][i]; 
        Beta[0][i]         = mktvol_data->BetaTD[0][i];
        if (NbFactor > 1)
        {
            FactorWeight[1][i] = mktvol_data->AlphaTD[1][i];
            Beta[1][i]         = mktvol_data->BetaTD[1][i]; 
            Rho[0][i]          = mktvol_data->RhoTD[0][i];
        }
        if (NbFactor > 2)
        {
            FactorWeight[2][i] = mktvol_data->AlphaTD[2][i];
            Beta[2][i]         = mktvol_data->BetaTD[2][i];
            Rho[1][i]          = mktvol_data->RhoTD[1][i];
            Rho[2][i]          = mktvol_data->RhoTD[2][i];
        }

    }

     /*
     *  Constant Aweight numbers (determined historically).
     *  NOTE: alphas are NOT normalized, but aweights ARE
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++) 
            aw[p][q] = 0.;

    


    if (CalibFlag == FALSE)
    {
        *NbBmkMr = NbTDInp;
        for (i = 0; i < NbTDInp; i++)                                 
        {
            aw[0][0]  = 1.;

            if (NbFactor > 1) 
            {
                aw[1][0] =  Rho[0][i];
                aw[1][1] =  sqrt(1 - Rho[0][i] * Rho[0][i]);
            }

            if (NbFactor > 2)
            {
                aw[2][0] = Rho[1][i];
                aw[2][1] = (Rho[2][i] - Rho[0][i]*Rho[1][i])
                            /sqrt(1 - Rho[0][i]*Rho[0][i]);
                aw[2][2] = sqrt(1.0 - Rho[0][i]*Rho[0][i] - Rho[1][i]*Rho[1][i]              
                    - Rho[2][i]*Rho[2][i] + 2.*Rho[0][i]*Rho[1][i]*Rho[2][i])
                      / sqrt(1 - Rho[0][i]*Rho[0][i]);
            }


            /* record mr on this bmark interval */
            for (k = 0; k < NbFactor; k++)
            {
                mktvol_data->BmkDate[i] = TDInpDate[i];  
                mktvol_data->BetaBmk[k][i] = Beta[k][i];
            }

			/* record Aweight */
            mktvol_data->Aweight[0][i] = aw[0][0] * FactorWeight[0][i];
        
            if (NbFactor > 1)
            {
                mktvol_data->Aweight[1][i] = aw[1][0] * FactorWeight[1][i];
                mktvol_data->Aweight[2][i] = aw[1][1] * FactorWeight[1][i];
            }
        
            if (NbFactor > 2)
            {
                mktvol_data->Aweight[3][i] = aw[2][0] * FactorWeight[2][i];
                mktvol_data->Aweight[4][i] = aw[2][1] * FactorWeight[2][i];
                mktvol_data->Aweight[5][i] = aw[2][2] * FactorWeight[2][i];
            }
        }
      

        return (SUCCESS);
    } 


    for (i = 0; i < NbVol; i++)
    {
        /* Only used vol points */
        if (!VolUsed[i]) continue;

        /* Time to option expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;
        /* Time to swap start in years */
        StartT[i] = Daysact (VolBaseDate, SwapSt[i]) / 365.;

        /* B factor */
        if (Fix3_BFactor_TimeDep (   B[i],
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
            E[p][q] = 1.;
        }
    
    
    LastCalibIdx = -1;              /* no calibrated points so far */
    mrIdx = 0;
    mrIdx_st = 0;
    prevMrIdx = -1;
    lastMrIdx = NbTDInp - 1;
    lastMrFlag = FALSE;
    lastMr_stFlag = FALSE;

    for (i = 0; i < NbVol; i++)                                                    
    {       
        if (!VolUsed[i])
            continue;

        if ( QInterp ( VolT[i],
              &QLeft,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QLeftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }

        if ( QInterp ( VolT[i],
              &QRight,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QRightTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }


        if ( QInterp ( VolT[i],
              &FwdSh,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->FwdShiftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }


        T = VolT[i];
        t = ((LastCalibIdx == -1) ? VolT[i] : (VolT[i] - VolT[LastCalibIdx]));
       

        /* determine mr on the current interval*/

        /* search for first index such that
          TDInpDate[mrIdx-1] < VolDate[i] <= TDInpDate[mrIdx]*/
        if(lastMrFlag == FALSE)
           mrIdx = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, VolDate[i], CbkHIGHER);
        if( mrIdx == NbTDInp - 1)
            lastMrFlag = TRUE;

        /* if all time dependent input dates smaller than VolDate[i]
            use last value */
        if (mrIdx < 0)
             mrIdx = NbTDInp - 1;
       
        if (mrIdx != prevMrIdx) /* this is always true for i = 0 */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaL[k] = Beta[k][mrIdx];
				WeightL[k] = FactorWeight[k][mrIdx];
            }
        }


        /* search for first index such that
           TDInpDate[mrIdx_st-1] < SwapSt[i] <= TDInpDate[mrIdx_st]*/
        if(lastMr_stFlag == FALSE)
           mrIdx_st = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, SwapSt[i], CbkHIGHER);
        if( mrIdx_st == NbTDInp - 1)
            lastMr_stFlag = TRUE;

        /* if all time dependent input dates smaller than SwapSt[i]
            use last value */
        if (mrIdx_st < 0)
             mrIdx_st = NbTDInp - 1;

        /* update the integrals E=exp(-2\int_{VolDate[i]^SwapSt[i]} beta) */
        for (j = mrIdx; j <= mrIdx_st; j++)
        {  
            EndDate = (j < mrIdx_st) ? TDInpDate[j] : SwapSt[i];

            TDDateT[j] = Daysact (VolDate[i], EndDate) / 365.;

            t0 = ((j == mrIdx) ? TDDateT[mrIdx] : (TDDateT[j]-TDDateT[j-1]));
            E[0][0] *= exp (- 2 * Beta[0][j] * t0);

            if (NbFactor > 1)
            {
                E[1][1] *= exp(- 2 * Beta[1][j] * t0);
                E[0][1] *= exp(- (Beta[0][j] + Beta[1][j]) * t0);
            }
            if(NbFactor > 2)
            {
                E[0][2] *= exp(- (Beta[0][j] +  Beta[2][j]) * t0);
                E[1][2] *= exp(- (Beta[1][j] +  Beta[2][j]) * t0);
                E[2][2] *= exp(- 2 * Beta[2][j] * t0);
            }
        }
       
        aw[0][0]  = 1.;
      
        if (NbFactor > 1) 
        {
            aw[1][0] =  Rho[0][mrIdx];
            aw[1][1] =  sqrt(1 - Rho[0][mrIdx] * Rho[0][mrIdx]);
        }

        if (NbFactor > 2)
        {
            aw[2][0] = Rho[1][mrIdx];
            aw[2][1] = (Rho[2][mrIdx] - Rho[0][mrIdx]*Rho[1][mrIdx])
                /sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
            aw[2][2] = sqrt(1.0 - Rho[0][mrIdx]*Rho[0][mrIdx] - 
                Rho[1][mrIdx]*Rho[1][mrIdx] - Rho[2][mrIdx]*Rho[2][mrIdx]
                + 2.*Rho[0][mrIdx]*Rho[1][mrIdx]*Rho[2][mrIdx]) 
               / sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
        }

      
        /* record mt on this bmark interval*/
        for (k = 0; k < NbFactor; k++)
        {
            BetaBmk[k][i] = BetaL[k];
        }

        
        D[0][0] = B[i][0];
        
        if (NbFactor > 1)
        {
            D[1][0] = B[i][1];
            D[1][1] = B[i][1];
        }
        
        if (NbFactor > 2)
        {
            D[2][0] = B[i][2];
            D[2][1] = B[i][2];
            D[2][2] = B[i][2];
        }

        /* Lognormal to X-space vol adjustment. The single and */
        /* two q cases are treated differently for consistency */
        /* with old model.                                     */
        if (Conv_MarketVol_To_XSpace(pY[i],
                                    Vol[i],        
                                    VolDate[i],
                                    VolTypeFlag,
                                    &y,
                                    T,
                                    FwdSh,
                                    QLeft,
                                    QRight,
                                    0.0, 
                                    1.0) == FAILURE) 
        {
            goto RETURN;
        }
         
       
        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * BetaL[0] * t) * E[0][0];
        
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * BetaL[0], t) * t * WeightL[0] * WeightL[0] 
            * aw[0][0] *aw[0][0] * E[0][0];
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][0] * exp (-2. * BetaL[1] * t) * E[1][1];
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(BetaL[0] + BetaL[1]) * t) * E[0][1];
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * BetaL[1] * t) * E[1][1];
            
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][0] * aw[1][0] * E[1][1];
            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (BetaL[0] + BetaL[1], t) * t * 
                WeightL[0] * WeightL[1] * aw[0][0] * aw[1][0] * E[0][1];
            
            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][1] * aw[1][1] * E[1][1];
        }
        
        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][0] * exp (-2. * BetaL[2] * t) * E[2][2];
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(BetaL[0] + BetaL[2]) * t) * E[0][2];
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t) * E[1][2];
            y -= D[2][1] * D[2][1] * L[2][1] * exp (-2. * BetaL[2] * t) * E[2][2];
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * BetaL[2] * t) * E[2][2];
            
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * BetaL[2], t) * t * 
                WeightL[2] * WeightL[2] * aw[2][0] * aw[2][0] * E[2][2];
            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (BetaL[0] + BetaL[2], t) * t * 
                WeightL[0] * WeightL[2] * aw[0][0] * aw[2][0] * E[0][2];
            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (BetaL[1] + BetaL[2], t) * t *
                WeightL[1] * WeightL[2] * aw[1][0] * aw[2][0] * E[1][2];
            
            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * BetaL[2], t) * t *  
                WeightL[2] * WeightL[2] * aw[2][1] * aw[2][1] * E[2][2]; 
            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (BetaL[1] + BetaL[2], t) * t * 
                WeightL[1] * WeightL[2] * aw[1][1] * aw[2][1] * E[1][2];
            
            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * BetaL[2], t) * t * 
                WeightL[2] * WeightL[2] * aw[2][2] * aw[2][2] * E[2][2];
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
            Aweight[0][k] = lambda * aw[0][0] * WeightL[0];
            
            if (NbFactor > 1)
            {
                Aweight[1][k] = lambda * aw[1][0] * WeightL[1];
                Aweight[2][k] = lambda * aw[1][1] * WeightL[1];
            }
            
            if (NbFactor > 2)
            {
                Aweight[3][k] = lambda * aw[2][0] * WeightL[2];
                Aweight[4][k] = lambda * aw[2][1] * WeightL[2];
                Aweight[5][k] = lambda * aw[2][2] * WeightL[2];
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
        L[0][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[0], t) * t * 
            WeightL[0] * WeightL[0] * aw[0][0] * aw[0][0];

        if (NbFactor > 1)
        {
            L[1][0] *= exp (-2. * BetaL[1] * t);
            L[1][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][0] * aw[1][0];
            L[1][1] *= exp (-2. * BetaL[1] * t);
            L[1][1] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][1] * aw[1][1];
            L[0][1] *= exp (-(BetaL[0] + BetaL[1]) * t);
            L[0][1] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[0]+BetaL[1]), t) * t * 
                WeightL[0] * WeightL[1] * aw[0][0] * aw[1][0];
        }

        if (NbFactor > 2)
        {
            L[2][0] *= exp (-2. * BetaL[2] * t);
            L[2][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][0] * aw[2][0];
            L[2][1] *= exp (-2. * BetaL[2] * t);
            L[2][1] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][1] * aw[2][1];
            L[2][2] *= exp (-2. * BetaL[2] * t);
            L[2][2] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][2] * aw[2][2];
            L[0][2] *= exp (-(BetaL[0] + BetaL[2]) * t);
            L[0][2] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[0]+BetaL[2]), t) * t * 
                WeightL[0] * WeightL[2] * aw[0][0] * aw[2][0];
            L[1][2] *= exp (-(BetaL[1] + BetaL[2]) * t);
            L[1][2] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[1]+BetaL[2]), t) * t *
                WeightL[1] * WeightL[2] * (aw[1][0] * aw[2][0] + aw[1][1] * aw[2][1]);
        }
        
    }  /* for i */


    if (LastCalibIdx == -1)
    {
        DR_Error ("Spot vol: none of points is calibrated!");
        goto RETURN;
    }

    
    for (k = LastCalibIdx+1; k < NbVol; k++)
    {
        /* determine mr and weight on the current interval.
           here we assume that all MrDates are VolDates! */
         
        if(lastMrFlag == FALSE)
            mrIdx = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, VolDate[i], CbkHIGHER);
        
        if (mrIdx == NbTDInp - 1)
            lastMrFlag = TRUE;
       
        /* record mr on this bmark interval */
        for (i = 0; i < NbFactor; i++)
        {
            BetaBmk[i][k] = Beta[i][mrIdx];
        }

        aw[0][0]  = 1.;

        if (NbFactor > 1) 
        {
            aw[1][0] =  Rho[0][mrIdx];
            aw[1][1] =  sqrt(1 - Rho[0][mrIdx] * Rho[0][mrIdx]);
        }

        if (NbFactor > 2)
        {
            aw[2][0] = Rho[1][mrIdx];
            aw[2][1] = (Rho[2][mrIdx] - Rho[0][mrIdx]*Rho[1][mrIdx])/sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
            aw[2][2] = sqrt(1.0 - Rho[0][mrIdx]*Rho[0][mrIdx] - Rho[1][mrIdx]*Rho[1][mrIdx] 
              - Rho[2][mrIdx]*Rho[2][mrIdx] + 2.*Rho[0][mrIdx]*Rho[1][mrIdx]*Rho[2][mrIdx]) 
              / sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
        }


		/* record Aweight */
		Aweight[0][k] = lambda * aw[0][0] * FactorWeight[0][mrIdx];
        
        if (NbFactor > 1)
        {
            Aweight[1][k] = lambda * aw[1][0] * FactorWeight[1][mrIdx];
            Aweight[2][k] = lambda * aw[1][1] * FactorWeight[1][mrIdx];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][k] = lambda * aw[2][0] * FactorWeight[2][mrIdx];
            Aweight[4][k] = lambda * aw[2][1] * FactorWeight[2][mrIdx];
            Aweight[5][k] = lambda * aw[2][2] * FactorWeight[2][mrIdx];
        }

    }  /* for k */

    for (i = 0; i < NbVol; i++)
        BmkDate[i] = VolDate[i];

    /* record mr, Aweight, and BmkDate after the last vol date */
  
       
   
    *NbBmkMr = NbVol;
   
    /* store computed Aweights in the mkt vol structure */
    for (i = 0; i < *NbBmkMr; i++)
    {
        mktvol_data->Aweight[0][i] = Aweight[0][i];
        mktvol_data->BetaBmk[0][i] = BetaBmk[0][i];
        if (NbFactor > 1)
        {
            mktvol_data->BetaBmk[1][i] = BetaBmk[1][i];
            mktvol_data->Aweight[1][i] = Aweight[1][i];
            mktvol_data->Aweight[2][i] = Aweight[2][i];
        }
        if (NbFactor > 2)
        {
            mktvol_data->BetaBmk[2][i] = BetaBmk[2][i];
            mktvol_data->Aweight[3][i] = Aweight[3][i];
            mktvol_data->Aweight[4][i] = Aweight[4][i];
            mktvol_data->Aweight[5][i] = Aweight[5][i];
        }
        
    }


  
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_SpotVol */




/*****  Fix3_Interp_SpotVol_TimeDep  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point (time dependent factor mr, correlations and weights)
*/
int     Fix3_Interp_SpotVol_TimeDep (
                             FIX3_TREE_DATA*  tree_data,   /* (I/ O) Tree data   */
                             MKTVOL_DATA*     mktvol_data) /* (I)Volatility data */
           
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  PrevVolT;              /* Time to previous vol point */
    double  Cov[3][3];             /* Covariance matrix          */
    double  x, t, T;
    double  BbqAdj;               /* Backbone vol adjustment    */

    int     i, j, k, l, p, q;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

    double  timeFrac;
 
    /* temporary variables for volatility parameters */
    double**         AweightCurve;          /* (O) Factors aweight curves        */
    double**         BetaTD;                /* (O) Time-dependent mr in the tree */
    double*          QLeft;                 /* (O) Time-dep q left in tree       */
    double*          QRight;                /* (O) Time-dep q right in tree      */
    double*          FwdShift;              /* (O) Time-dep fwd shift in tree    */
    long             VolBaseDate;           /* (I) Volatility curve base date    */
    int              NbVol;                 /* (I) Number of spot vol points     */
    long const*      VolDate;               /* (I) Spot vol dates                */
    double           Aweight[6][MAXNBDATE]; /* (I) Aweight curve of each factor  */
    int              NbFactor;              /* (I) Number of factors             */
    int              NbBmk;                 /* (I) Number of bmk MR              */
    long             *BmkDate;              /* (I) Bmk MR dates                  */
    double           BetaBmk[3][MAXNBDATE]; /* (I) Mean reversions               */
    double           Bbq;                   /* (I) Backbone parameter            */
    double           VolNorm;               /* (I) Normal volatility in backbone */
    double           VolLogn;               /* (I) Lognorm volatility in backbone*/
    int              CalibFlag;             /* (I) Index calibration flag        */
    double const*    FwdRate;               /* (I) One period forward rate       */
    double const*    Length;                /* (I) Length of each time step      */
    int              NbTP;                  /* (I) Zero curve base date          */ 

    /* set temporary variables to mktvol_data param */
    VolBaseDate  = mktvol_data->BaseDate;
    NbVol        = mktvol_data->NbVol;
    VolDate      = mktvol_data->VolDate;
    NbFactor     = mktvol_data->NbFactor;
    NbBmk        = mktvol_data->NbBmkMr;
    BmkDate      = mktvol_data->BmkDate;
    Bbq          = mktvol_data->Bbq;
    VolNorm      = mktvol_data->VolNorm;
    VolLogn      = mktvol_data->VolLogn;
    CalibFlag    = mktvol_data->CalibFlag;
    FwdRate      = tree_data->FwdRate[tree_data->CvDiff];
    Length       = tree_data->Length;
    NbTP         = tree_data->NbTP;
    AweightCurve = tree_data->Aweight;
    BetaTD       = tree_data->BetaTD;
    QLeft        = tree_data->QLeft;
    QRight       = tree_data->QRight;
    FwdShift     = tree_data->FwdShift;
  
    


    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    for (i = 0; i < NbBmk; i++)
    {
        BetaBmk[0][i] = mktvol_data->BetaBmk[0][i];
        if (NbFactor > 1)
            BetaBmk[1][i] = mktvol_data->BetaBmk[1][i];
        if (NbFactor > 2)
            BetaBmk[2][i] = mktvol_data->BetaBmk[2][i];
    }

    /* set temporary Aweight array to mktvol_data Aweight */
    
    for (i = 0; i < NbBmk; i++)
    {
        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][i] = mktvol_data->Aweight[k][i]; 
        }
    }
   


    for (j = 0; j < NbBmk; j++)
    {       
        VolT[j] = Daysact (VolBaseDate, BmkDate[j]) / 365.;
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

        while ((j < NbBmk - 1) && (T > VolT[j] + ERROR))
            j++;

        /* find QLeft at the current time point */
        if ( QInterp ( T,
              &(QLeft[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QLeftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }
        

        /* find QRight at the current time point*/
        if ( QInterp ( T,
              &(QRight[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QRightTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }

        /* find Fwd shift at the current time point */
        if ( QInterp ( T,
              &(FwdShift[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->FwdShiftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }


      
        PrevVolT = ((j == 0)? 0. : VolT[j-1]);
         
        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        if (t >= PrevVolT)
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = BetaBmk[k][j];
            }
            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][j];
            }
        }
        /* The time step is straddling two buckets: */
        /* we recalculate the covariance matrix.    */
        else
        {
            timeFrac = (T- PrevVolT) / (T - t);
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = (1. - timeFrac) * BetaBmk[k][j-1] + timeFrac * BetaBmk[k][j];
            }

            Cov[0][0]  = Aweight[0][j-1] * Aweight[0][j-1] * 
                         Fix3_ExpDecay (2. * BetaBmk[0][j-1], VolT[j-1] - t) * (VolT[j-1] - t)
                         * exp (-2. * BetaBmk[0][j] * (T - VolT[j-1] ));
            Cov[0][0] += Aweight[0][j] * Aweight[0][j]   * Fix3_ExpDecay (2. * BetaBmk[0][j], T - VolT[j-1]) * 
                         (T - VolT[j-1]) ;
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = Aweight[0][j-1] * Aweight[1][j-1] * 
                             Fix3_ExpDecay (BetaBmk[0][j - 1] + BetaBmk[1][ j - 1], VolT[j-1] - t) 
                             * (VolT[j-1] - t)* exp (-(BetaBmk[0][j] + BetaBmk[1][j])*(T - VolT[j-1]));
                Cov[1][0] += Aweight[0][j]  * Aweight[1][j]  * Fix3_ExpDecay (BetaBmk[0][j] + BetaBmk[1][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) ;
            
                Cov[1][1]  = Aweight[1][j-1] * Aweight[1][j-1] * Fix3_ExpDecay (2. * BetaBmk[1][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-2. * BetaBmk[1][j] * (T - VolT[j-1] ));
                Cov[1][1] += Aweight[1][j]   * Aweight[1][j]   * Fix3_ExpDecay (2. * BetaBmk[1][j], T - VolT[j-1]) * (T - VolT[j-1]) ;
                Cov[1][1] += Aweight[2][j-1] * Aweight[2][j-1] * Fix3_ExpDecay (2. * BetaBmk[1][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-2. * BetaBmk[1][j] * (T - VolT[j-1] ));
                Cov[1][1] += Aweight[2][j]   * Aweight[2][j]   * Fix3_ExpDecay (2. * BetaBmk[1][j], T - VolT[j-1]) * 
                             (T - VolT[j-1]) ;
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = Aweight[0][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (BetaBmk[0][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-(BetaBmk[0][j]+BetaBmk[2][j])*(T - VolT[j-1]));
                Cov[2][0] += Aweight[0][j]   * Aweight[3][j]   * Fix3_ExpDecay (BetaBmk[0][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) ;
            
                Cov[2][1]  = Aweight[1][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (BetaBmk[1][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-(BetaBmk[1][j]+BetaBmk[2][j])*(T - VolT[j-1]));
                Cov[2][1] += Aweight[1][j]   * Aweight[3][j]   * Fix3_ExpDecay (BetaBmk[1][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) ;
                Cov[2][1] += Aweight[2][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (BetaBmk[1][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-(BetaBmk[1][j]+BetaBmk[2][j])*(T - VolT[j-1]));
                Cov[2][1] += Aweight[2][j]   * Aweight[4][j]   * Fix3_ExpDecay (BetaBmk[1][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) ;

                Cov[2][2]  = Aweight[3][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t)  * exp (-2. * BetaBmk[2][j] * (T - VolT[j-1]));
                Cov[2][2] += Aweight[3][j]   * Aweight[3][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]);
                Cov[2][2] += Aweight[4][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-2. * BetaBmk[2][j - 1] * (T - VolT[j-1] ));
                Cov[2][2] += Aweight[4][j]   * Aweight[4][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                             (T - VolT[j-1]) ;
                Cov[2][2] += Aweight[5][j-1] * Aweight[5][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * 
                    (VolT[j-1] - t) * exp (-2. * BetaBmk[2][j - 1] * (T - VolT[j-1] ));
                Cov[2][2] += Aweight[5][j]   * Aweight[5][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]);
            }


            x = Cov[0][0] / Fix3_ExpDecay (2. * BetaTD[0][i], T - t) / (T - t);

            if (x < TINY)
            {
                DR_Error ("Fix3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / Fix3_ExpDecay (BetaTD[0][i] + BetaTD[1][i], T - t) / (T - t);

                x = Cov[1][1] / Fix3_ExpDecay (2. * BetaTD[1][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

           

                if (x < TINY)
                {
                    DR_Error ("Fix3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                    goto RETURN;
                }
                        
                AweightCurve[2][i] = sqrt (x);
            }

            if (NbFactor > 2)
            {
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / Fix3_ExpDecay (BetaTD[0][i] + BetaTD[2][i], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / Fix3_ExpDecay (BetaTD[1][i] + BetaTD[2][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / Fix3_ExpDecay (2. * BetaTD[2][i], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

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

            BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                     (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*Fix3_ExpDecay(BetaTD[l][i], Length[i+1]);
        }
    }
  
    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Interp_SpotVol_TimeDep */






/*****  Fix3_IndexVol_TimeDep    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*       Time dependent mr, correlations and weights
*/
int     Fix3_IndexVol_TimeDep (
          double*       Vol,                  /* (O) Index vol curve              */
          long          OptExp,               /* (I) Option expiry                */
          long          SwapSt,               /* (I) Underlying swap start        */
          long          SwapMat,              /* (I) Underlying swap maturity     */
          char          Freq,                 /* (I) Frequency of underlying rate */
          char          DCC,                  /* (I) Day count convention         */        
          MKTVOL_DATA   *mktvol_data,         /* (I) Volatility data              */
          T_CURVE const* crv)                 /* (I) zero curve                   */
{

    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  BmkDateT[MAXNBDATE];/*Bmk dates in years                    */
    double  OptExpT;           /* Option expiry in years                */
    double  StartT;            /* Time to swap start in years           */
    double  t;                 /* Time between two consecutive expiries */
    double  T=0.;              /* Time to current expiry                */
    double  ParYield;
    double  Annuity;
    
    double  L[3][3];           /* Integrals of factor spot vol          */
    double  B[3];              /* B in Christian memo                   */
    double  D[3][3];           /* aweight*B                             */
    double  E[3][3];
    double  M;                 /* Total variance                        */

    long    EndDate;           /* End of current integration bucket     */
    int     i, j, p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially      */
    char    ErrorMsg[MAXBUFF]; /* Error message                         */

    long mrIdx , mrIdx_st, lastMrFlag, lastMr_stFlag;


    double  PrevT;

    /* temporary variables for volatility parameters */
    int              NbVol;        /* (I) Nb of points in vol curve       */
    long             VolBaseDate;  /* (I) Volatility base date            */
    long const*      VolDate;      /* (I) Volatility dates                */
    double           QLeft;        /* (I) Left Q mapping coefficient      */
    double           QRight;       /* (I) Right Q mapping coefficient     */
    double           FwdShift;     /* (I) Fwd shift mapping coefficient   */
    double           Bbq;          /* (I) Backbone parameter              */
    double           VolNorm;      /* (I) Normal volatility in backbone   */
    double           VolLogn;      /* (I) Lognorm volatility in backbone  */
    int              VolTypeFlag;  /* (I) Volatility type                 */
    int              NbFactor;     /* (I) Number of factors               */
    int              SkipFlag;     /* (I) Skip calibration failure points */
    int              CalibFlag;    /* (I) Index calibration flag          */
    double           Aweight[6][MAXNBDATE];/* (I) Spot vol curve          */ 
    long             *BmkDate;     /* (I) Bmk dates                       */
    int              NbBmkMr;     /* (I) Nb of mr bmk dates              */
    double       BetaBmk[3][MAXNBDATE]; /* (O) Mr on all benchmark intvals*/
   

    /* set termporary variables to mktvol_data parameters */
    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Bbq         = mktvol_data->Bbq;
    VolTypeFlag = mktvol_data->VolUnit;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;
    BmkDate     = mktvol_data->BmkDate;
    NbBmkMr    =  mktvol_data->NbBmkMr;
   
    
    
    for (i = 0; i < NbBmkMr; i++)
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
    

    for (i = 0; i < NbBmkMr; i++)
    {
        BmkDateT[i] = Daysact (VolBaseDate, BmkDate[i]) / 365.;

        BetaBmk[0][i] = mktvol_data->BetaBmk[0][i];
        if (NbFactor > 1)
            BetaBmk[1][i] = mktvol_data->BetaBmk[1][i];
        if (NbFactor > 2)
            BetaBmk[2][i] = mktvol_data->BetaBmk[2][i];
    }


    OptExpT = Daysact(VolBaseDate, OptExp) / 365.;

    /* Find B coefficients */
    if (Fix3_BFactor_TimeDep(B,
                        SwapSt,
                        SwapMat,
                        DCC,
                        Freq,
                        mktvol_data,
                        crv) == FAILURE)
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

    

    /* Avoid spot index expiration */
    if (SwapSt <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }

    StartT = Daysact (VolBaseDate, SwapSt) / 365.;
    
    for(i = 0; i < 3; i++)
    {
        for(j = 0; j < 3; j++) 
        {
            L[i][j] = 0.;
            E[i][j] = 1.;
        }
    }

    PrevT = 0;


    mrIdx = FALSE;
    mrIdx_st = FALSE;
    lastMrFlag = FALSE;
    lastMr_stFlag = FALSE;

    
    /* search for first index such that
      TDInpDate[mrIdx-1] < VolDate[i] <= TDInpDate[mrIdx]*/
    if(lastMrFlag == FALSE)
       mrIdx = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, OptExp, CbkHIGHER);
    if( mrIdx == NbBmkMr - 1)
        lastMrFlag = TRUE;

    /* if all time dependent input dates smaller than VolDate[i]
        use last value */
    if (mrIdx < 0)
         mrIdx = NbBmkMr - 1;

    /* search for first index such that
           TDInpDate[mrIdx_st-1] < SwapSt[i] <= TDInpDate[mrIdx_st]*/
    if(lastMr_stFlag == FALSE)
       mrIdx_st = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, SwapSt, CbkHIGHER);
    if( mrIdx_st == NbBmkMr - 1)
        lastMr_stFlag = TRUE;
    if (mrIdx_st < 0)
         mrIdx_st = NbBmkMr - 1;
  

    /* update the integrals E=exp(-2\int_{OptExp^SwapSt} beta) */
    for (i = mrIdx; i <= mrIdx_st; i++)
    {
        EndDate = (i < mrIdx_st) ?  BmkDate[i] : SwapSt;

        VolT[i] = Daysact (OptExp, EndDate) / 365.;

        t = ((i == mrIdx) ? VolT[mrIdx] : (VolT[i]-VolT[i-1]));
        E[0][0] *= exp (- 2 * BetaBmk[0][i] * t);

        if (NbFactor > 1)
        {
            E[1][1] *= exp(- 2 * BetaBmk[1][i] * t);
            E[0][1] *= exp(- (BetaBmk[0][i] + BetaBmk[1][i]) * t);
        }
        if(NbFactor > 2)
        {
            E[0][2] *= exp(- (BetaBmk[0][i] +  BetaBmk[2][i]) * t);
            E[1][2] *= exp(- (BetaBmk[1][i] +  BetaBmk[2][i]) * t);
            E[2][2] *= exp(- 2 * BetaBmk[2][i] * t);
        }
    }


    for (i = 0; i < NbBmkMr; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate = (i < NbBmkMr-1) ? MIN(OptExp, BmkDate[i]) : OptExp;

        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));        

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaBmk[0][i] * t);
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * BetaBmk[0][i], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaBmk[1][i] * t);
            L[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);
            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);
            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t;
            L[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);
            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == OptExp)
        {
            break;
        }
    }  /* for i */


    M = B[0] * B[0] * L[0][0] * E[0][0];
    
    if (NbFactor > 1) 
    {
        M += B[1] * B[1] * L[1][1] * E[1][1];
        M += B[0] * B[1] * L[0][1] * E[0][1];
    }
    
    if (NbFactor > 2)
    {
        M += B[2] * B[2] * L[2][2] * E[2][2];
        M += B[0] * B[2] * L[0][2] * E[0][2];
        M += B[1] * B[2] * L[1][2] * E[1][2];
    }
    

    if (M < SQUARE(TINY))
    {
        sprintf (ErrorMsg, "IndexVol: problem in integrating index volatility "
                 "at date %ld (IndexVol)", SwapSt);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    

    M = sqrt (M / T);
    
    
    /* find QLeft, QRight, FwdSh */
    if ( QInterp ( T,
          &QLeft,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->QLeftTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }

    if ( QInterp (T,
          &QRight,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->QRightTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }


    if ( QInterp (T,
          &FwdShift,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->FwdShiftTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }



    /* Convert from q measure back to market */
    if (Conv_XSpaceVol_To_Market(ParYield,
                                    Vol,        
                                    VolTypeFlag,
                                    M,
                                    T,
                                    FwdShift,
                                    QLeft,
                                    QRight,
                                    0.,
                                    1.) == FAILURE) 
    {
        goto RETURN;
    }

    
   
    status = SUCCESS;

    RETURN:

    return (status);

   
}  /* Fix3_IndexVol_TimeDep */


/*****  Fix3_IndexCorr_TimeDep    *********************************************************/
/*
*       Compute correlation between two rate indices (time dependent parameters)
*/
int     Fix3_IndexCorr_TimeDep (
          double*       Corr,                /* (O) Index vol curve              */
          long          OptExp1,             /* (I) Option expiry                */
          long          SwapSt1,             /* (I) Underlying swap start        */
          long          SwapMat1,            /* (I) Underlying swap maturity     */
          char          Freq1,               /* (I) Frequency of underlying rate */
          char          DCC1,                /* (I) Day count convention         */
          long          OptExp2,             /* (I) Option expiry                */        
          long          SwapSt2,             /* (I) Underlying swap start        */
          long          SwapMat2,            /* (I) Underlying swap maturity     */
          char          Freq2,               /* (I) Frequency of underlying rate */
          char          DCC2,                /* (I) Day count convention         */
          MKTVOL_DATA   *mktvol_data,        /* (I) Volatility data              */
          T_CURVE const* crv)           /* (I) zero curve */
{

    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  t;                 /* Time between two consecutive expiries */
    double  T=0.;              /* Time to current expiry                */
    double  ParYield1;
    double  Annuity1;
    double  ParYield2;
    double  Annuity2;
    double  Vol1, Vol2;

    double  L[3][3];           /* Integrals of factor spot vol          */
    double  L1[3][3];          /* Integrals of factor spot vol          */
    double  L2[3][3];          /* Integrals of factor spot vol          */
    double  B1[3];             /* B in Christian memo                   */
    double  B2[3];             /* B in Christian memo                   */
    double  D[3][3];           /* aweight*B                             */
    double  E1[3][3];          /* exponential integrals                 */
    double  E2[3][3];          /* exponential integrals                 */
    double  E[3];              /* exponential integrals                 */
    double  M;                 /* Total covariance                      */
    double  M1, M2;            /* Total variances                       */

    long    EndDate;           /* End of current integration bucket     */
    int     i, p, q;


    int     status = FAILURE;  /* Error status = FAILURE initially      */
  
   
    double  PrevT;
    long    minIdxS, maxIdxS;  /* index of bmk interval where Min, max between
                                  SwapSt1, SwapSt2 lie    */ 
    long    MinS, MaxS;        /* Min, max between SwapSt1, SwapSt2     */
    long    minIdxE;           /* index of bmk interval where Min between 
                                  Expiry1, Expiry2 lies                 */
    long    MinE;              /* Min expiry                            */


    int              NbVol;        /* (I) Nb of points in vol curve       */
    long             VolBaseDate;  /* (I) Volatility base date            */
    long const*      VolDate;      /* (I) Volatility dates                */
    double           Bbq;          /* (I) Backbone parameter              */
    double           VolNorm;      /* (I) Normal volatility in backbone   */
    double           VolLogn;      /* (I) Lognorm volatility in backbone  */
    int              NbFactor;     /* (I) Number of factors               */
    int              SkipFlag;     /* (I) Skip calibration failure points */
    int              CalibFlag;    /* (I) Index calibration flag          */
    double           Aweight[6][MAXNBDATE];/* (I) Spot vol curve          */ 
    long             *BmkDate;     /* (I) Bmk dates                       */
    int              NbBmkMr;      /* (I) Nb of mr bmk dates              */
    double       BetaBmk[3][MAXNBDATE]; /* (O) Mr on all benchmark intvals*/
  

    long mrIdx1 , mrIdx_st1, lastMrFlag1, lastMr_stFlag1; /*for 1st rate  */
    long mrIdx2 , mrIdx_st2, lastMrFlag2, lastMr_stFlag2; /*for 2nd rate  */

    /*initialize mr flags and indices */
    mrIdx1 = 0;
    mrIdx_st1 = 0;
    lastMrFlag1 = FALSE;
    lastMr_stFlag1 = FALSE;

    mrIdx2 = 0;
    mrIdx_st2 = 0;
    lastMrFlag2 = FALSE;
    lastMr_stFlag2 = FALSE;




    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Bbq         = mktvol_data->Bbq;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;
    BmkDate     = mktvol_data->BmkDate;
    NbBmkMr     = mktvol_data->NbBmkMr;
   

    for (i = 0; i < NbBmkMr; i++)
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
   

    for (i = 0; i < NbBmkMr; i++)
    {
        BetaBmk[0][i] = mktvol_data->BetaBmk[0][i];
        if (NbFactor > 1)
            BetaBmk[1][i] = mktvol_data->BetaBmk[1][i];
        if (NbFactor > 2)
            BetaBmk[2][i] = mktvol_data->BetaBmk[2][i];
    }

    /* Find B coefficients */
    if (Fix3_BFactor_TimeDep (B1,
                        SwapSt1,
                        SwapMat1,
                        DCC1,
                        Freq1,
                        mktvol_data,
                        crv) == FAILURE)
    {
         goto RETURN;
    }

    if (Fix3_BFactor_TimeDep (B2,
                        SwapSt2,
                        SwapMat2,
                        DCC2,
                        Freq2,
                        mktvol_data,
                        crv) == FAILURE)
    {
         goto RETURN;
    }

    /* Find par yield */


    if (ParYieldFromDates (&ParYield1,
                           &Annuity1,
                           SwapSt1,
                           SwapMat1,
                           DCC1,
                           Freq1,
                           'F',
                           crv) != SUCCESS)
    {
        goto RETURN;
    }

    
    if (ParYieldFromDates (&ParYield2,
                           &Annuity2,
                           SwapSt2,
                           SwapMat2,
                           DCC2,
                           Freq2,
                           'F',
                           crv) != SUCCESS)
    {
        goto RETURN;
    }


    /* Constant spot volatility case */
    for(p = 0; p < 3; p++)
    {
        for(q = 0; q < 3; q++)
        {
            D[p][q] = 1.;
            L[p][q] = 0.;
            L1[p][q] = 0.;
            L2[p][q] = 0.;
            E1[p][q] = 1.;
            E2[p][q] = 1.;
        }
        E[p] = 1.;
    }
        

    PrevT = 0;

     /* search for first index such that
      BmkDate[mrIdx-1] < OptExp1 <= BmkDate[mrIdx]*/
    if(lastMrFlag1 == FALSE)
       mrIdx1 = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, OptExp1, CbkHIGHER);
    if( mrIdx1 == NbBmkMr - 1)
        lastMrFlag1 = TRUE;

    /* if all time dependent input dates smaller than OptExp1 use last value */
    if (mrIdx1 < 0)
         mrIdx1 = NbBmkMr - 1;

    /* search for first index such that
           BmkDate[mrIdx_st-1] < SwapSt1 <= BmkDate[mrIdx_st]*/
    if(lastMr_stFlag1 == FALSE)
       mrIdx_st1 = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, SwapSt1, CbkHIGHER);
    if( mrIdx_st1 == NbBmkMr - 1)
        lastMr_stFlag1 = TRUE;
    /* if all time dependent input dates smaller than SwapSt1 use last value */
    if (mrIdx_st1 < 0)
         mrIdx_st1 = NbBmkMr - 1;
  

    /* update the integrals E=exp(-2\int_{OptExp1^SwapSt1} beta) */
    for (i = mrIdx1; i <= mrIdx_st1; i++)
    {
        EndDate = (i < mrIdx_st1) ?  BmkDate[i] : SwapSt1;

        VolT[i] = Daysact (OptExp1, EndDate) / 365.;

        t = ((i == mrIdx1) ? VolT[mrIdx1] : (VolT[i]-VolT[i-1]));
        E1[0][0] *= exp (- 2 * BetaBmk[0][i] * t);

        if (NbFactor > 1)
        {
            E1[1][1] *= exp(- 2 * BetaBmk[1][i] * t);
            E1[0][1] *= exp(- (BetaBmk[0][i] + BetaBmk[1][i]) * t);
        }
        if(NbFactor > 2)
        {
            E1[0][2] *= exp(- (BetaBmk[0][i] +  BetaBmk[2][i]) * t);
            E1[1][2] *= exp(- (BetaBmk[1][i] +  BetaBmk[2][i]) * t);
            E1[2][2] *= exp(- 2 * BetaBmk[2][i] * t);
        }
    }


    /* search for first index such that
     BmkDate[mrIdx-1] < OptExp2 <= BmkDate[mrIdx]*/
    if(lastMrFlag2 == FALSE)
       mrIdx2 = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, OptExp2, CbkHIGHER);
    if( mrIdx2 == NbBmkMr - 1)
        lastMrFlag2 = TRUE;
    /* if all benchmark dates smaller than OptExp2 use last value */
    if (mrIdx2 < 0)
         mrIdx2 = NbBmkMr - 1;

    /* search for first index such that
           BmkDate[mrIdx_st-1] < SwapSt2 <= BmkDate[mrIdx_st]*/
    if(lastMr_stFlag2 == FALSE)
        mrIdx_st2 = (NbBmkMr == 1) ? 0 : GetDLOffset (NbBmkMr, BmkDate, SwapSt2, CbkHIGHER);
    if( mrIdx_st2 == NbBmkMr - 1)
        lastMr_stFlag2 = TRUE;
    /* if all benchmark dates smaller than SwapSt2 use last value */
    if (mrIdx_st2 < 0)
        mrIdx_st2 = NbBmkMr - 1;
  

    /* update the integrals E=exp(-2\int_{OptExp2^SwapSt2} beta) */
    for (i = mrIdx2; i <= mrIdx_st2; i++)
    {
        EndDate = (i < mrIdx_st2) ?  BmkDate[i] : SwapSt2;

        VolT[i] = Daysact (OptExp2, EndDate) / 365.;

        t = ((i == mrIdx2) ? VolT[mrIdx2] : (VolT[i]-VolT[i-1]));
        E2[0][0] *= exp (- 2 * BetaBmk[0][i] * t);

        if (NbFactor > 1)
        {
            E2[1][1] *= exp(- 2 * BetaBmk[1][i] * t);
            E2[0][1] *= exp(- (BetaBmk[0][i] + BetaBmk[1][i]) * t);
        }
        if(NbFactor > 2)
        {
            E2[0][2] *= exp(- (BetaBmk[0][i] +  BetaBmk[2][i]) * t);
            E2[1][2] *= exp(- (BetaBmk[1][i] +  BetaBmk[2][i]) * t);
            E2[2][2] *= exp(- 2 * BetaBmk[2][i] * t);
        }
    }

    MinE    = MIN(OptExp1, OptExp2);
    MinS    = MIN(SwapSt1, SwapSt2);
    MaxS    = MAX(SwapSt1, SwapSt2);
    minIdxE = MIN(mrIdx1, mrIdx2);
    minIdxS = MIN(mrIdx_st1, mrIdx_st2);
    maxIdxS = MAX(mrIdx_st1, mrIdx_st2);
   

     /* need to compute exp(-\int_MIN(SwapSt1, SwapSt2)^max(SwapSt1, SwapSt2) (-\beta(u)du) */
    for (i = minIdxS ; i <= maxIdxS; i++)
    {
        EndDate = ((i < maxIdxS) ? BmkDate[i] : MaxS);
        VolT[i] = Daysact (MinS, EndDate) / 365.;
        t = ((i == minIdxS)  ? VolT[i] : (VolT[i] - VolT[i - 1]));

        E[0] *= exp(-BetaBmk[0][i] * t);
        if (NbFactor > 1)
            E[1] *= exp(-BetaBmk[1][i] * t);

        if (NbFactor > 2)
            E[2] *= exp(-BetaBmk[2][i] * t);
    }

      /* need to compute exp(-\int_MIN(OptExp1, OptExp2)^min(SwapSt1, SwapSt2) (-2 *\beta(u)du) */
    for (i = minIdxE ; i <= minIdxS; i++)
    {
        EndDate = ((i < minIdxS) ? BmkDate[i] : MinS);
        VolT[i] = Daysact (MinE, EndDate) / 365.;
        t = ((i == minIdxE)  ? VolT[i] : (VolT[i] - VolT[i - 1]));

        D[0][0] *= exp(-2 * BetaBmk[0][i] * t);
        if (NbFactor > 1)
        {
            D[1][1] *= exp(- 2 * BetaBmk[1][i] * t);
            D[0][1] *= exp(- (BetaBmk[0][i] + BetaBmk[1][i]) * t);
        }
        if(NbFactor > 2)
        {
            D[0][2] *= exp(- (BetaBmk[0][i] +  BetaBmk[2][i]) * t);
            D[1][2] *= exp(- (BetaBmk[1][i] +  BetaBmk[2][i]) * t);
            D[2][2] *= exp(- 2 * BetaBmk[2][i] * t);
        }
       
    }


    for (i = 0; i < NbBmkMr; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate = (i < NbBmkMr-1) ? MIN(MinE, BmkDate[i]) : MinE;

        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));        

        /* 
         *  Update integrals
         */


        L[0][0] *= exp (-2. * BetaBmk[0][i] * t) ;
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * BetaBmk[0][i], t) * t
            * E[0] * D[0][0];

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaBmk[1][i] * t) ;
            L[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t
                       * E[1] * D[1][1];;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t
                       * E[1] * D[1][1];

            L[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);
            L[0][1] += Aweight[0][i] *  Aweight[1][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t
                     * E[1] * D[0][1];
            L[0][1] += Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t
                     * E[0] * D[0][1];
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t
                       * E[2] * D[2][2];

            L[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t
                       * E[2] * D[2][2];

            L[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t
                        * E[2] * D[2][2];

            L[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);
            L[0][2] += Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t
                       * E[2] * D[0][2];
            L[0][2] += Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t
                       * E[0] * D[0][2];

            L[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);
            L[1][2] +=  Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t
                * E[2] * D[1][2];
            L[1][2] += Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t
                * E[1] * D[1][2] ;

            L[1][2] += Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t
                * E[2] * D[1][2];
            
            L[1][2] += Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t
                * E[1] * D[1][2];

        }

        /* End of integration loop */
        if (EndDate == MinE)
        {
            break;
        }
    }  /* for i */


   
    for (i = 0; i < NbBmkMr; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate = (i < NbBmkMr-1) ? MIN(OptExp1, BmkDate[i]) : OptExp1;

        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));        

        /* 
         *  Update integrals
         */
        L1[0][0] *= exp (-2. * BetaBmk[0][i] * t);
        L1[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * BetaBmk[0][i], t) * t;

        if (NbFactor > 1)
        {
            L1[1][1] *= exp (-2. * BetaBmk[1][i] * t);
            L1[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L1[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L1[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);
            L1[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t;
        }

        if (NbFactor > 2)
        {
            L1[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L1[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L1[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L1[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L1[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);
            L1[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t;
            L1[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);
            L1[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
            L1[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == OptExp1)
        {
            break;
        }
    }  /* for i */


    for (i = 0; i < NbBmkMr; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */
        EndDate= (i < NbBmkMr-1) ? MIN(SwapSt2, BmkDate[i]) : SwapSt2;

        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;
        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));        

        /* 
         *  Update integrals
         */
        L2[0][0] *= exp (-2. * BetaBmk[0][i] * t);
        L2[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * BetaBmk[0][i], t) * t;

        if (NbFactor > 1)
        {
            L2[1][1] *= exp (-2. * BetaBmk[1][i] * t);
            L2[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L2[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;
            L2[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);
            L2[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t;
        }

        if (NbFactor > 2)
        {
            L2[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L2[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L2[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L2[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;
            L2[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);
            L2[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t;
            L2[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);
            L2[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
            L2[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == OptExp2)
        {
            break;
        }
    }  /* for i */
   
    M = B1[0] * B2[0] * L[0][0];
    
    if (NbFactor > 1) 
    {
        M += B1[1] * B2[1] * L[1][1];
        M += B1[0] * B2[1] * L[0][1];
    }
    
    if (NbFactor > 2)
    {
        M += B1[2] * B2[2] * L[2][2];
        M += B1[0] * B2[2] * L[0][2];
        M += B1[1] * B2[2] * L[1][2];
    }
    
    M1 = B1[0] * B1[0] * L1[0][0] * E1[0][0];
    
    if (NbFactor > 1) 
    {
        M1 += B1[1] * B1[1] * L1[1][1] * E1[1][1];
        M1 += B1[0] * B1[1] * L1[0][1] * E1[0][1];
    }
    
    if (NbFactor > 2)
    {
        M1 += B1[2] * B1[2] * L1[2][2] * E1[2][2];
        M1 += B1[0] * B1[2] * L1[0][2] * E1[0][2];
        M1 += B1[1] * B1[2] * L1[1][2] * E1[1][2];
    }
    
    M2 = B2[0] * B2[0] * L2[0][0] * E2[0][0];
    
    if (NbFactor > 1) 
    {
        M2 += B2[1] * B2[1] * L2[1][1] * E2[1][1];
        M2 += B2[0] * B2[1] * L2[0][1] * E2[0][1];
    }
    
    if (NbFactor > 2)
    {
        M2 += B2[2] * B2[2] * L2[2][2] * E2[2][2];
        M2 += B2[0] * B2[2] * L2[0][2] * E2[0][2];
        M2 += B2[1] * B2[2] * L2[1][2] * E2[1][2];
    }


    if ((M1 < SQUARE(TINY)) || (M2 < SQUARE(TINY)))
    {
        DR_Error ( "IndexCorr: problem in integrating index volatility ");
        goto RETURN;
    }
    
    Vol1 = sqrt(M1);
    Vol2 = sqrt(M2);
    *Corr = M / (Vol1 * Vol2  );
    
    
    status = SUCCESS;

    RETURN:

    return (status);

   
}  /* Fix3_IndexCorr_TimeDep */

/*****  Fix3_Filtered_SpotVol_TimeDep    **************************************************/
/*
*       Bootstrap input volatility curve to extract filtered spot vols.
*       Spot vols are not modified provided the ratio of the current spot vol
*       to the spot vol in the first bucket is greater than MIN_SPOT_VOL_RATIO
*       or less than 1/MIN_SPOT_VOL_RATIO. Otherwise the spot vols are smoothly
*       capped and floored, with the cap and floor levels parameterized by
*       SPOT_VOL_FILTER_AMOUNT. (time dependent parameters )
*/
int    Fix3_Filtered_SpotVol_TimeDep(
       MKTVOL_DATA *mktvol_data,        /* (I/O) Volatility data            */
       T_CURVE const* crv)              /* (I) zero curve                   */
{
    double  VolT[MAXNBDATE];      /* Expiries in years                      */
    double  StartT[MAXNBDATE];    /* Swap start in years                    */
    double  TDDateT[MAXNBDATE];   /* Time dep dates in years                */
    double  t;                    /* Time between two consecutive expiries  */
    double  t0;
    double  T;                    /* Time to current expiry                 */
    double  pY[MAXNBDATE];        /* Fwd par yield                          */
    double  Annuity;              /* Fwd annuity                            */

    double  B[MAXNBDATE][3];      /* B in Christian memo                    */
    double  L[3][3];              /* Integrals of lambda factors            */
    double  D[3][3];              /* aweight*B                              */
    double  E[3][3];              /* Integrals of exp                       */
    double  aw[3][3];             /* Historical aweights                    */
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
    int     i, j, k, p, q;
    int     status = FAILURE;     /* Error status = FAILURE initially       */


    double  BetaL[3];
	double  WeightL[3];
    long    mrIdx, mrIdx_st, prevMrIdx, lastMrIdx;
    int     lastMrFlag, lastMr_stFlag;
    
	
    /* temporary variables for volatility data parameters*/
    double           Aweight[6][MAXNBDATE];/* (O) Spot vol curve of each factor   */
    double           BetaBmk[3][MAXNBDATE];/* (O) Mr on all benchmark intvals    */
    long             *BmkDate;             /* (O) Bmk dates                       */
    int              *NbBmkMr;             /* (O) Nb of mr bmk dates              */
    long             VolBaseDate;          /* (I) Volatility curve base date      */
    int              NbVol;                /* (I) Nb of points in vol curve       */
    long const*      VolDate;              /* (I) Volatility dates                */
    double const*    Vol;                  /* (I) Vol curve                       */
    int const*       VolUsed;              /* (I) TRUE if vol used in calibration */
    char             Freq;                 /* (I) Frequency of underlying rate    */
    char             DCC;                  /* (I) Day count convention            */
    long const*      SwapSt;               /* (I) Underlying swap start           */
    long const*      SwapMat;              /* (I) Underlying swap maturity        */
    double           QLeft;                /* (I) Left Q mapping coefficient      */
    double           QRight;               /* (I) Right Q mapping coefficient     */
    double           FwdSh;                /* (I) Fwd shift mapping coefficient   */
    double           Bbq;                  /* (I) Backbone parameter              */
    int              VolTypeFlag;          /* (I) Volatility type                 */
    double           VolNorm;              /* (I) Normal volatility in backbone   */
    double           VolLogn;              /* (I) Lognorm volatility in backbone  */
    int              NbFactor;             /* (I) Number of factors               */
    int              NbTDInp;              /* (I) Number of TD input dates        */
    long             *TDInpDate;           /* (I) TD input dates                  */
    long             EndDate;
    double           FactorWeight[3][MAXNBDATE]; /* (I) Time-dependent factor weights*/
    double           Beta[3][MAXNBDATE];    /* (I) Mean reversions                */
    double           Rho[3][MAXNBDATE];     /* (I) Correlation between factors    */
    int              CalibFlag;             /* (I) Index calibration flag         */
    int              SkipFlag;              /* (I) Skip calibration failure points*/



     /* initialize temporary variables to mktvol_data parameters */
    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Vol         = mktvol_data->Vol;
    VolUsed     = mktvol_data->VolUsed;
    Freq        = mktvol_data->Freq;
    DCC         = mktvol_data->DCC;
    SwapSt      = mktvol_data->SwapSt;
    SwapMat     = mktvol_data->SwapMat;
    Bbq         = mktvol_data->Bbq;
    VolTypeFlag = mktvol_data->VolUnit;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    NbFactor    = mktvol_data->NbFactor;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;

    NbTDInp     = mktvol_data->NbTDInp;
    TDInpDate   = mktvol_data->TDInpDate;
    BmkDate     = mktvol_data->BmkDate;
    NbBmkMr     = &(mktvol_data->NbBmkMr);
    

    for (i = 0; i < NbTDInp; i++)
    {
        FactorWeight[0][i] = mktvol_data->AlphaTD[0][i]; 
        Beta[0][i]         = mktvol_data->BetaTD[0][i];
        if (NbFactor > 1)
        {
            FactorWeight[1][i] = mktvol_data->AlphaTD[1][i];
            Beta[1][i]         = mktvol_data->BetaTD[1][i]; 
            Rho[0][i]          = mktvol_data->RhoTD[0][i];
        }
        if (NbFactor > 2)
        {
            FactorWeight[2][i] = mktvol_data->AlphaTD[2][i];
            Beta[2][i]         = mktvol_data->BetaTD[2][i];
            Rho[1][i]          = mktvol_data->RhoTD[1][i];
            Rho[2][i]          = mktvol_data->RhoTD[2][i];
        }

    }

     /*
     *  Constant Aweight numbers (determined historically).
     *  NOTE: alphas are NOT normalized, but aweights ARE
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++) 
            aw[p][q] = 0.;




    if (CalibFlag == FALSE)
    {
        *NbBmkMr = NbTDInp;
		for (i = 0; i < NbTDInp; i++)                                                    
        { 
            aw[0][0]  = 1.;

            if (NbFactor > 1) 
            {
                aw[1][0] =  Rho[0][i];
                aw[1][1] =  sqrt(1 - Rho[0][i] * Rho[0][i]);
            }

            if (NbFactor > 2)
            {
                aw[2][0] = Rho[1][i];
                aw[2][1] = (Rho[2][i] - Rho[0][i]*Rho[1][i])/sqrt(1 - Rho[0][i]*Rho[0][i]);
                aw[2][2] = sqrt(1.0 - Rho[0][i]*Rho[0][i] - Rho[1][i]*Rho[1][i] 
                  - Rho[2][i]*Rho[2][i] + 2.*Rho[0][i]*Rho[1][i]*Rho[2][i]) / sqrt(1 - Rho[0][i]*Rho[0][i]);
            }


            /* record mr on this bmark interval */
            for (k = 0; k < NbFactor; k++)
            {
                mktvol_data->BmkDate[i] = TDInpDate[i];  
                mktvol_data->BetaBmk[k][i] = Beta[k][i];
            }

			/* record Aweight */
            mktvol_data->Aweight[0][i] = aw[0][0] * FactorWeight[0][i];
        
            if (NbFactor > 1)
            {
                mktvol_data->Aweight[1][i] = aw[1][0] * FactorWeight[1][i];
                mktvol_data->Aweight[2][i] = aw[1][1] * FactorWeight[1][i];
            }
        
            if (NbFactor > 2)
            {
                mktvol_data->Aweight[3][i] = aw[2][0] * FactorWeight[2][i];
                mktvol_data->Aweight[4][i] = aw[2][1] * FactorWeight[2][i];
                mktvol_data->Aweight[5][i] = aw[2][2] * FactorWeight[2][i];
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

    for (i = 0; i < NbVolUsed; i++)
    {
        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;
        StartT[i] = Daysact (VolBaseDate, SwapSt[i]) / 365.;

        
         /* B factor */
        if (Fix3_BFactor_TimeDep (   B[i],
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
            E[p][q] = 1.;
        }

    mrIdx = 0;
    mrIdx_st = 0;
    prevMrIdx = -1;
    lastMrIdx = NbTDInp - 1;
    lastMrFlag = FALSE;
    lastMr_stFlag = FALSE;

    for (i = 0; i < NbVolUsed; i++)
    {
        /* find QLeft at the current VolDate */
        if ( QInterp ( VolT[i],
              &QLeft,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QLeftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }


        
        /* find QRight at the current VolDate */
        if ( QInterp ( VolT[i],
              &QRight,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QRightTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }

        /* find FwdSh at the current VolDate */
        if ( QInterp ( VolT[i],
              &FwdSh,
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->FwdShiftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }

        T = VolT[i];
        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));
        if(lastMrFlag == FALSE)
           mrIdx = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, VolDate[i], CbkHIGHER);
      
        if (mrIdx == NbTDInp -1)
            lastMrFlag = TRUE;

          /* if all time dependent input dates smaller than VolDate[i]
            use last value */
        if (mrIdx < 0)
             mrIdx = NbTDInp - 1;
       

        if (mrIdx != prevMrIdx) /* this is always true for i = 0 */
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaL[k] = Beta[k][mrIdx];
                WeightL[k] = FactorWeight[k][mrIdx];
            }
        }

         /* search for first index such that
           TDInpDate[mrIdx_st-1] < SwapSt[i] <= TDInpDate[mrIdx_st]*/
        if(lastMr_stFlag == FALSE)
           mrIdx_st = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, SwapSt[i], CbkHIGHER);
        if (mrIdx_st == NbTDInp -1)
            lastMr_stFlag = TRUE;


         /* if all time dependent input dates smaller than SwapSt[i]
            use last value */
        if (mrIdx_st < 0)
             mrIdx_st = NbTDInp - 1;


         /* update the integrals E=exp(-2\int_{OptExp^SwapSt} beta) */
        for (j = mrIdx; j <= mrIdx_st; j++)
        {  
            EndDate = (j < mrIdx_st ) ? TDInpDate[j] : SwapSt[i];

            TDDateT[j] = Daysact (VolDate[i], EndDate) / 365.;

            t0 = ((j == mrIdx) ? TDDateT[mrIdx] : (TDDateT[j]-TDDateT[j-1]));
            E[0][0] *= exp (- 2 * Beta[0][j] * t0);

            if (NbFactor > 1)
            {
                E[1][1] *= exp(- 2 * Beta[1][j] * t0);
                E[0][1] *= exp(- (Beta[0][j] + Beta[1][j]) * t0);
            }
            if(NbFactor > 2)
            {
                E[0][2] *= exp(- (Beta[0][j] +  Beta[2][j]) * t0);
                E[1][2] *= exp(- (Beta[1][j] +  Beta[2][j]) * t0);
                E[2][2] *= exp(- 2 * Beta[2][j] * t0);
            }
        }
     

        aw[0][0]  = 1.;
      
        if (NbFactor > 1) 
        {
            aw[1][0] =  Rho[0][mrIdx];
            aw[1][1] =  sqrt(1 - Rho[0][mrIdx] * Rho[0][mrIdx]);
        }

        if (NbFactor > 2)
        {
            aw[2][0] = Rho[1][mrIdx];
            aw[2][1] = (Rho[2][mrIdx] - Rho[0][mrIdx]*Rho[1][mrIdx])/sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
            aw[2][2] = sqrt(1.0 - Rho[0][mrIdx]*Rho[0][mrIdx] - Rho[1][mrIdx]*Rho[1][mrIdx] 
              - Rho[2][mrIdx]*Rho[2][mrIdx] + 2.*Rho[0][mrIdx]*Rho[1][mrIdx]*Rho[2][mrIdx]) 
              / sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
       
        }

        /* record mr on this bmark interval*/
        for (k = 0; k <NbFactor; k++)
        {
            BetaBmk[k][i] = BetaL[k];
        }

        
        D[0][0] = B[i][0];
        
        if (NbFactor > 1)
        {
            D[1][0] = B[i][1];
            D[1][1] = B[i][1];
        }
        
        if (NbFactor > 2)
        {
            D[2][0] = B[i][2];
            D[2][1] = B[i][2];
            D[2][2] = B[i][2];
        }
        
        

        /* Market to X-space vol adjustment. The single and */
        /* two q cases are treated differently for consistency */
        /* with old model.                                     */
        if (Conv_MarketVol_To_XSpace(pY[i],
                                    Vol[i],        
                                    VolDate[i], 
                                    VolTypeFlag,
                                    &y,
                                    T,
                                    FwdSh,
                                    QLeft,
                                    QRight,
                                    0.0,
                                    1.0) == FAILURE) 
        {
            goto RETURN;
        }
        
        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * BetaL[0] * t) * E[0][0];
        
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * BetaL[0], t) * t *
            WeightL[0] * WeightL[0] * aw[0][0] * aw[0][0] * E[0][0];
        
        if (NbFactor > 1) 
        {
            y -= D[1][0] * D[1][0] * L[1][0] * exp (-2. * BetaL[1] * t) * E[1][1];
            y -= D[0][0] * D[1][0] * L[0][1] * exp (-(BetaL[0] + BetaL[1]) * t) * E[0][1];
            y -= D[1][1] * D[1][1] * L[1][1] * exp (-2. * BetaL[1] * t) * E[1][1];
            
            M += D[1][0] * D[1][0] * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                 WeightL[1] * WeightL[1] * aw[1][0] * aw[1][0] * E[1][1] ;
            M += 2. * D[0][0] * D[1][0] * Fix3_ExpDecay (BetaL[0] + BetaL[1], t) * t * 
                 WeightL[0] * WeightL[1] * aw[0][0] *aw[1][0] * E[0][1];
            
            M += D[1][1] * D[1][1] * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                 WeightL[1] * WeightL[1] * aw[1][1] * aw[1][1] * E[1][1];
        }

        if (NbFactor > 2)
        {
            y -= D[2][0] * D[2][0] * L[2][0] * exp (-2. * BetaL[2] * t) * E[2][2];
            y -= D[0][0] * D[2][0] * L[0][2] * exp (-(BetaL[0] + BetaL[2]) * t) * E[0][2];
            y -= D[1][0] * D[2][0] * L[1][2] * exp (-(BetaL[1] + BetaL[2]) * t) * E[1][2];
            y -= D[2][1] * D[2][1] * L[2][1] * exp (-2. * BetaL[2] * t) * E[2][2];
            y -= D[2][2] * D[2][2] * L[2][2] * exp (-2. * BetaL[2] * t) * E[2][2];
            
            M += D[2][0] * D[2][0] * Fix3_ExpDecay (2. * BetaL[2], t) * t * 
                 WeightL[2] * WeightL[2] * aw[2][0] * aw[2][0] * E[2][2];
            M += 2. * D[0][0] * D[2][0] * Fix3_ExpDecay (BetaL[0] + BetaL[2], t) * t * 
                 WeightL[0] * WeightL[2] * aw[0][0] * aw[2][0] * E[0][2];
            M += 2. * D[1][0] * D[2][0] * Fix3_ExpDecay (BetaL[1] + BetaL[2], t) * t * 
                 WeightL[1] * WeightL[2] * aw[1][0] * aw[2][0] * E[1][2];
            
            M += D[2][1] * D[2][1] * Fix3_ExpDecay (2. * BetaL[2], t) * t * 
                 WeightL[2] * WeightL[2] * aw[2][1] * aw[2][1] *E[2][2];
            M += 2. * D[1][1] * D[2][1] * Fix3_ExpDecay (BetaL[1] + BetaL[2], t) * t * 
                 WeightL[1] * WeightL[2] * aw[1][1] * aw[2][1] * E[1][2];
            
            M += D[2][2] * D[2][2] * Fix3_ExpDecay (2. * BetaL[2], t) * t * 
                 WeightL[2] * WeightL[2] *  aw[2][2] * aw[2][2] * E[2][2];
        }

        /* M must be stictly positive, and since it is not a function   */
        /* of spot vol, trap the error now rather than later in Fix3_SpotVol */
        if (M < TINY)
        {
            DR_Error ("TimeDep Filtered spot vol: problem in bootstrapping %ld "
                      "volatility (negative or zero variance = %lf)", 
                      YMDDateFromIRDate(VolDate[i]), M);
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
        Aweight[0][i] = lambda * aw[0][0] * WeightL[0];
        
        if (NbFactor > 1)
        {
            Aweight[1][i] = lambda * aw[1][0] * WeightL[1];
            Aweight[2][i] = lambda * aw[1][1] * WeightL[1];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][i] = lambda * aw[2][0] * WeightL[2];
            Aweight[4][i] = lambda * aw[2][1] * WeightL[2];
            Aweight[5][i] = lambda * aw[2][2] * WeightL[2];
        }

        prevMrIdx = mrIdx;
        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaL[0] * t);
        L[0][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[0], t) * t * 
            WeightL[0] * WeightL[0] * aw[0][0] * aw[0][0];

        if (NbFactor > 1)
        {
            L[1][0] *= exp (-2. * BetaL[1] * t);
            L[1][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][0] * aw[1][0];
            L[1][1] *= exp (-2. * BetaL[1] * t);
            L[1][1] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[1], t) * t * 
                WeightL[1] * WeightL[1] * aw[1][1] * aw[1][1];
            L[0][1] *= exp (-(BetaL[0] + BetaL[1]) * t);
            L[0][1] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[0]+BetaL[1]), t) * t * 
                WeightL[0] * WeightL[1] * aw[0][0] * aw[1][0];
        }

        if (NbFactor > 2)
        {
            L[2][0] *= exp (-2. * BetaL[2] * t);
            L[2][0] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][0] * aw[2][0];
            L[2][1] *= exp (-2. * BetaL[2] * t);
            L[2][1] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][1] * aw[2][1];
            L[2][2] *= exp (-2. * BetaL[2] * t);
            L[2][2] += lambda * lambda * Fix3_ExpDecay (2. * BetaL[2], t) * t *
                WeightL[2] * WeightL[2] * aw[2][2] * aw[2][2];
            L[0][2] *= exp (-(BetaL[0] + BetaL[2]) * t);
            L[0][2] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[0]+BetaL[2]), t) * t * 
                WeightL[0] * WeightL[2] * aw[0][0] * aw[2][0];
            L[1][2] *= exp (-(BetaL[1] + BetaL[2]) * t);
            L[1][2] += 2.*lambda*lambda * Fix3_ExpDecay ((BetaL[1]+BetaL[2]), t) * t *
                WeightL[1] * WeightL[2] * (aw[1][0] * aw[2][0] + aw[1][1] * aw[2][1]);
        }
    }  /* for i */

    for (k = NbVolUsed; k < NbVol; k++)
    {
       /* determine mr on the current interval.
           here we assume that all MrDates are VolDates! */
        if(lastMrFlag == FALSE)
            mrIdx = (NbTDInp == 1) ? 0 : GetDLOffset (NbTDInp, TDInpDate, VolDate[i], CbkHIGHER);
      
        if (mrIdx == NbTDInp - 1)
            lastMrFlag = TRUE;
       
        /* record mr on this bmark interval*/
        for (i = 0; i < NbFactor; i++)
        {
            BetaBmk[i][k] = Beta[i][mrIdx];
        }

        aw[0][0]  = 1.;
      
        if (NbFactor > 1) 
        {
            aw[1][0] =  Rho[0][mrIdx];
            aw[1][1] =  sqrt(1 - Rho[0][mrIdx] * Rho[0][mrIdx]);
        }

        if (NbFactor > 2)
        {
            aw[2][0] = Rho[1][mrIdx];
            aw[2][1] = (Rho[2][mrIdx] - Rho[0][mrIdx]*Rho[1][mrIdx])/sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
            aw[2][2] = sqrt(1.0 - Rho[0][mrIdx]*Rho[0][mrIdx] - Rho[1][mrIdx]*Rho[1][mrIdx] 
              - Rho[2][mrIdx]*Rho[2][mrIdx] + 2.*Rho[0][mrIdx]*Rho[1][mrIdx]*Rho[2][mrIdx]) 
              / sqrt(1 - Rho[0][mrIdx]*Rho[0][mrIdx]);
       
        }


        Aweight[0][k] = lambda * aw[0][0] * FactorWeight[0][mrIdx];
        
        if (NbFactor > 1)
        {
            Aweight[1][k] = lambda * aw[1][0] * FactorWeight[1][mrIdx];
            Aweight[2][k] = lambda * aw[1][1] * FactorWeight[1][mrIdx];
        }
        
        if (NbFactor > 2)
        {
            Aweight[3][k] = lambda * aw[2][0] * FactorWeight[2][mrIdx];
            Aweight[4][k] = lambda * aw[2][1] * FactorWeight[2][mrIdx];
            Aweight[5][k] = lambda * aw[2][2] * FactorWeight[2][mrIdx];
        }
    }  /* for k */

    for (i = 0; i < NbVol; i++)
        BmkDate[i] = VolDate[i];

    *NbBmkMr = NbVol;
    for (i = 0; i < *NbBmkMr; i++)
    {
        mktvol_data->Aweight[0][i] = Aweight[0][i];
        mktvol_data->BetaBmk[0][i] = BetaBmk[0][i];
        if (NbFactor > 1)
        {
            mktvol_data->BetaBmk[1][i] = BetaBmk[1][i];
            mktvol_data->Aweight[1][i] = Aweight[1][i];
            mktvol_data->Aweight[2][i] = Aweight[2][i];
        }
        if (NbFactor > 2)
        {
            mktvol_data->BetaBmk[2][i] = BetaBmk[2][i];
            mktvol_data->Aweight[3][i] = Aweight[3][i];
            mktvol_data->Aweight[4][i] = Aweight[4][i];
            mktvol_data->Aweight[5][i] = Aweight[5][i];
        }
    }


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Filtered_SpotVol */




/*****  Fix3_Interp_SpotVol_TimeDepOld  *****************************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point. With the bug in classic fix3
*/
int     Fix3_Interp_SpotVol_TimeDepOld (
                             FIX3_TREE_DATA*  tree_data,   /* (I/ O) Tree data   */
                             MKTVOL_DATA*     mktvol_data) /* (I)Volatility data */
       
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  PrevVolT;              /* Time to previous vol point */
    double  Cov[3][3];             /* Covariance matrix          */
    double  x, t, T;
    double  BbqAdj;               /* Backbone vol adjustment    */

    int     i, j, k, l, p, q;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

    double  timeFrac;


    double**         AweightCurve;          /* (O) Factors aweight curves        */
    double**         BetaTD;                /* (O) Time-dependent mr in the tree */
    long             VolBaseDate;           /* (I) Volatility curve base date    */
    int              NbVol;                 /* (I) Number of spot vol points     */
    long const*      VolDate;               /* (I) Spot vol dates                */
    double           Aweight[6][MAXNBDATE]; /* (I) Aweight curve of each factor  */
    int              NbFactor;              /* (I) Number of factors             */
    int              NbBmk;                 /* (I) Number of bmk MR              */
    long             *BmkDate;              /* (I) Bmk MR dates                  */
    double           BetaBmk[3][MAXNBDATE]; /* (I) Mean reversions               */
    double*          QLeft;                 /* (O) Time-dep q left in tree       */
    double*          QRight;                /* (O) Time-dep q right in tree      */
    double*          FwdShift;              /* (O) Time-dep fwd shift in tree    */
    double           Bbq;                   /* (I) Backbone parameter            */
    double           VolNorm;               /* (I) Normal volatility in backbone */
    double           VolLogn;               /* (I) Lognorm volatility in backbone*/
    int              CalibFlag;             /* (I) Index calibration flag        */
    double const*    FwdRate;               /* (I) One period forward rate       */
    double const*    Length;                /* (I) Length of each time step      */
    int              NbTP;                  /* (I) Zero curve base date          */ 

    /* set temporary volatility param to mktvol_data param */
    VolBaseDate  = mktvol_data->BaseDate;
    NbVol        = mktvol_data->NbVol;
    VolDate      = mktvol_data->VolDate;
    NbFactor     = mktvol_data->NbFactor;
    NbBmk        = mktvol_data->NbBmkMr;
    BmkDate      = mktvol_data->BmkDate;
    Bbq          = mktvol_data->Bbq;
    VolNorm      = mktvol_data->VolNorm;
    VolLogn      = mktvol_data->VolLogn;
    CalibFlag    = mktvol_data->CalibFlag;
    FwdRate      = tree_data->FwdRate[tree_data->CvDiff];
    Length       = tree_data->Length;
    NbTP         = tree_data->NbTP;
    AweightCurve = tree_data->Aweight;
    BetaTD       = tree_data->BetaTD;
    QLeft        = tree_data->QLeft;
    QRight       = tree_data->QRight;
    FwdShift     = tree_data->FwdShift;
    
    
    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    for (i = 0; i < NbBmk; i++)
    {
        BetaBmk[0][i] = mktvol_data->BetaBmk[0][i];
        if (NbFactor > 1)
            BetaBmk[1][i] = mktvol_data->BetaBmk[1][i];
        if (NbFactor > 2)
            BetaBmk[2][i] = mktvol_data->BetaBmk[2][i];
    }




    /* set temporary Aweight array to mktvol_data Aweight */ 
    for (i = 0; i < NbBmk; i++)
    {
        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][i] = mktvol_data->Aweight[k][i]; 
        }
    }
    



    for (j = 0; j < NbBmk; j++)
    {       
        VolT[j] = Daysact (VolBaseDate, BmkDate[j]) / 365.;
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

         /* find QLeft at the current time point */
        if ( QInterp ( T,
              &(QLeft[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QLeftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }
        

        /* find QRight at the current time point*/
        if ( QInterp ( T,
              &(QRight[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->QRightTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }

        /* find Fwd shift at the current time point */
        if ( QInterp ( T,
              &(FwdShift[i]),
              mktvol_data->BaseDate,
              mktvol_data->SmileDate,
              mktvol_data->FwdShiftTD,
              mktvol_data->NbSmileDates) == FAILURE)
        {
            goto RETURN;
        }


        while ((j < NbBmk - 1) && (T > VolT[j] + ERROR))
            j++;

        PrevVolT = ((j == 0)? 0. : VolT[j-1]);

        /* The current time step is in between bucket j and j-1: */
        /* its Aweight is equal to the Aweight of bucket j.      */
        if (t >= PrevVolT)
        {
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = BetaBmk[k][j];
            }
            for (k = 0; k < NbAweight; k++)
            {
                AweightCurve[k][i] = Aweight[k][j];
            }
        }
        /* The time step is straddling two buckets: */
        /* we recalculate the covariance matrix.    */
        else
        {
            timeFrac = (T- PrevVolT) / (T - t);
            for (k = 0; k < NbFactor; k++)
            {
                BetaTD[k][i] = (1. - timeFrac) * BetaBmk[k][j-1] + timeFrac * BetaBmk[k][j];
            }

            Cov[0][0]  = Aweight[0][j-1] * Aweight[0][j-1] * 
                         Fix3_ExpDecay (2. * BetaBmk[0][j-1], VolT[j-1] - t) * (VolT[j-1] - t);
            Cov[0][0] += Aweight[0][j] * Aweight[0][j]   * Fix3_ExpDecay (2. * BetaBmk[0][j], T - VolT[j-1]) * 
                         (T - VolT[j-1]) * exp (-2. * BetaBmk[0][j] * (VolT[j-1] - t));
            
            if (NbFactor > 1)
            {
                Cov[1][0]  = Aweight[0][j-1] * Aweight[1][j-1] * 
                             Fix3_ExpDecay (BetaBmk[0][j - 1] + BetaBmk[1][ j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][0] += Aweight[0][j]  * Aweight[1][j]  * Fix3_ExpDecay (BetaBmk[0][j] + BetaBmk[1][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-(BetaBmk[0][j] + BetaBmk[1][j])*(VolT[j-1]-t));
            
                Cov[1][1]  = Aweight[1][j-1] * Aweight[1][j-1] * Fix3_ExpDecay (2. * BetaBmk[1][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[1][j]   * Aweight[1][j]   * Fix3_ExpDecay (2. * BetaBmk[1][j], T - VolT[j-1]) * (T - VolT[j-1]) 
                             * exp (-2. * BetaBmk[1][j] * (VolT[j-1] - t));
                Cov[1][1] += Aweight[2][j-1] * Aweight[2][j-1] * Fix3_ExpDecay (2. * BetaBmk[1][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[1][1] += Aweight[2][j]   * Aweight[2][j]   * Fix3_ExpDecay (2. * BetaBmk[1][j], T - VolT[j-1]) * 
                             (T - VolT[j-1]) * exp (-2. * BetaBmk[1][j] * (VolT[j-1] - t));
            }

            if (NbFactor > 2)
            {
                Cov[2][0]  = Aweight[0][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (BetaBmk[0][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][0] += Aweight[0][j]   * Aweight[3][j]   * Fix3_ExpDecay (BetaBmk[0][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-(BetaBmk[0][j]+BetaBmk[2][j])*(VolT[j-1]-t));
            
                Cov[2][1]  = Aweight[1][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (BetaBmk[1][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[1][j]   * Aweight[3][j]   * Fix3_ExpDecay (BetaBmk[1][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-(BetaBmk[1][j]+BetaBmk[2][j])*(VolT[j-1]-t));
                Cov[2][1] += Aweight[2][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (BetaBmk[1][j - 1] + BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][1] += Aweight[2][j]   * Aweight[4][j]   * Fix3_ExpDecay (BetaBmk[1][j] + BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-(BetaBmk[1][j]+BetaBmk[2][j])*(VolT[j-1]-t));

                Cov[2][2]  = Aweight[3][j-1] * Aweight[3][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[3][j]   * Aweight[3][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-2. * BetaBmk[2][j] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[4][j-1] * Aweight[4][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[4][j]   * Aweight[4][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                             (T - VolT[j-1]) * exp (-2. * BetaBmk[2][j - 1] * (VolT[j-1] - t));
                Cov[2][2] += Aweight[5][j-1] * Aweight[5][j-1] * Fix3_ExpDecay (2. * BetaBmk[2][j - 1], VolT[j-1] - t) * (VolT[j-1] - t);
                Cov[2][2] += Aweight[5][j]   * Aweight[5][j]   * Fix3_ExpDecay (2. * BetaBmk[2][j], T - VolT[j-1]) * 
                            (T - VolT[j-1]) * exp (-2. * BetaBmk[2][j - 1] * (VolT[j-1] - t));
            }


            x = Cov[0][0] / Fix3_ExpDecay (2. * BetaTD[0][i], T - t) / (T - t);

            if (x < TINY)
            {
                DR_Error ("Fix3_Interp_SpotVolTimeDepOld: problem in interpolating spot volatility at node #%d!", i);
                goto RETURN;
            }
                        
            AweightCurve[0][i] = sqrt (x);

            if (NbFactor > 1)
            {
                AweightCurve[1][i] = Cov[1][0] / AweightCurve[0][i] / Fix3_ExpDecay (BetaTD[0][i] + BetaTD[1][i], T - t) / (T - t);

                x = Cov[1][1] / Fix3_ExpDecay (2. * BetaTD[1][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[1][i];

                if (x < TINY)
                {
                    DR_Error ("Fix3_Interp_SpotVol: problem in interpolating spot volatility at node #%d!", i);
                    goto RETURN;
                }
                        
                AweightCurve[2][i] = sqrt (x);
            }

            if (NbFactor > 2)
            {
                AweightCurve[3][i] = Cov[2][0] / AweightCurve[0][i] / Fix3_ExpDecay (BetaTD[0][i] + BetaTD[2][i], T - t) / (T - t);

                AweightCurve[4][i] = (Cov[2][1] / Fix3_ExpDecay (BetaTD[1][i] + BetaTD[2][i], T - t) / (T - t) - AweightCurve[1][i] * AweightCurve[3][i]) / AweightCurve[2][i];

                x = Cov[2][2] / Fix3_ExpDecay (2. * BetaTD[2][i], T - t) / (T - t) - AweightCurve[3][i] * AweightCurve[3][i] - AweightCurve[4][i] * AweightCurve[4][i];

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

            BbqAdj = (Bbq*VolLogn*log(1.+FwdRate[i+1])+(1-Bbq)*VolNorm*Length[i+1])/
                     (Bbq*VolLogn*FwdRate[i+1]        +(1-Bbq)*VolNorm*Length[i+1]);
                
            AweightCurve[k][i]*=(1.+FwdRate[i+1])*BbqAdj*Fix3_ExpDecay(BetaTD[l][i], Length[i+1]);
        }
    }
  
        
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Interp_SpotVol_TimeDepOld */

/*****  GenIndexVol_TimeDep   *************************************************/
/*
*       Integrates spot volatility curve between two time points falling 
*       before or on the swap start date.
*/
int   Fix3_GenIndexVol_TimeDep (
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

    double  VolT[MAXNBDATE];   /* Expiries in years                          */
    double  t;                 /* Time between two consecutive expiries      */
    double  T=0.;              /* Time to current expiry                     */
    double  Int01=0;           /* size of int from IStart to IEnd            */ 
    double  Int12=0;           /* size of int from IEnd to SwapStart         */ 
    double  Int02=0;           /* size of int from IStart to SwapStart       */ 
    double  IntQ=0;            /* Time for Q interpolation                   */
    double  ParYield;
    double  Annuity;
  

    double  L[3][3];           /* Integrals of factor spot vol               */
    double  B[3];              /* B in Christian memo                        */
    double  D[3][3];           /* aweight*B                                  */
    double  M;                 /* Total variance                             */

    long    EndDate;           /* End of current integration bucket          */
    int     StartIdx=0;
    int     EndIdx = 0;
    int     i, j, p, q;

    int     status = FAILURE;  /* Error status = FAILURE initially           */
    char    ErrorMsg[MAXBUFF]; /* Error message                              */
    
    /* temporary variables for volatility parameters */
    int              NbVol;        /* (I) Nb of points in vol curve       */
    long             VolBaseDate;  /* (I) Volatility base date            */
    long const*      VolDate;      /* (I) Volatility dates                */
    double           QLeft;        /* (I) Left Q mapping coefficient      */
    double           QRight;       /* (I) Right Q mapping coefficient     */
    double           FwdShift;     /* (I) Fwd shift mapping coefficient   */
    double           Bbq;          /* (I) Backbone parameter              */
    double           VolNorm;      /* (I) Normal volatility in backbone   */
    double           VolLogn;      /* (I) Lognorm volatility in backbone  */
    int              VolTypeFlag;  /* (I) Normal or Lognormal             */
    int              NbFactor;     /* (I) Number of factors               */
    double const*    Alpha;        /* (I) Relative size factors           */
    double const*    Beta;         /* (I) Mean reversions                 */
    double const*    Rho;          /* (I) Correlation between factors     */
    int              SkipFlag;     /* (I) Skip calibration failure points */
    int              CalibFlag;    /* (I) Index calibration flag          */
    double   Aweight[6][MAXNBDATE];/* (O) Spot vol curve of each factor   */
    long             *BmkDate;     /* (I) Bmk dates                       */
    int              NbBmkMr;      /* (I) Nb of mr bmk dates              */
    double       BetaBmk[3][MAXNBDATE]; /* (O) Mr on all benchmark intvals*/
   
    


    /* Avoid spot index expiration; note we use base date of zero */
    /* curve as there is no volatility information available.     */
    if ((IEnd <= crv->ValueDate)|| ( IStart == IEnd))
    {
            *Vol = 0.00;
            return(SUCCESS);
    }
    
    /* set temporary variables to mktvol_data parameters */
    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    Bbq         = mktvol_data->Bbq;
    VolNorm     = mktvol_data->VolNorm;
    VolLogn     = mktvol_data->VolLogn;
    VolTypeFlag = mktvol_data->VolUnit;
    NbFactor    = mktvol_data->NbFactor;
    Alpha       = mktvol_data->Alpha;
    Beta        = mktvol_data->Beta;
    Rho         = mktvol_data->Rho;
    SkipFlag    = mktvol_data->SkipFlag;
    CalibFlag   = mktvol_data->CalibFlag;
    BmkDate     = mktvol_data->BmkDate;
    NbBmkMr     =  mktvol_data->NbBmkMr;

    
    for (i = 0; i < NbBmkMr; i++)
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
    
    

    for (i = 0; i < NbBmkMr; i++)
    {
        BetaBmk[0][i] = mktvol_data->BetaBmk[0][i];
        if (NbFactor > 1)
            BetaBmk[1][i] = mktvol_data->BetaBmk[1][i];
        if (NbFactor > 2)
            BetaBmk[2][i] = mktvol_data->BetaBmk[2][i];
    }


    Int01 = Daysact (IStart, IEnd)      / 365.;
    Int12 = Daysact (IEnd, SwapStart)   / 365.;
    Int02 = Daysact (IStart, SwapStart) / 365.;  
    
    /* Find B coefficients */
    if (Fix3_BFactor_TimeDep (B,
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
    StartIdx = GetDLOffset (NbBmkMr, BmkDate, IStart, CbkHIGHER);
    
    for (i = StartIdx; i < NbBmkMr; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */

        EndDate = (i < NbBmkMr-1) ? MIN(IEnd, BmkDate[i]) : IEnd;

        /* Time to expiry in years */
        VolT[i] = Daysact (IStart, EndDate) / 365.;

        T = VolT[i];
        t = ((i == StartIdx) ? VolT[StartIdx] : (VolT[i]-VolT[i-1]));
        

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * BetaBmk[0][i] * t) ;
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * BetaBmk[0][i], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * BetaBmk[1][i] * t);

            L[1][1] += Aweight[1][i] * Aweight[1][i] * 
                Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;

            L[1][1] += Aweight[2][i] * Aweight[2][i] *
                Fix3_ExpDecay (2. * BetaBmk[1][i], t) * t;

            L[0][1] *= exp (-(BetaBmk[0][i] + BetaBmk[1][i]) * t);

            L[0][1] += 2. * Aweight[0][i] * Aweight[1][i] * 
                Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[1][i]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * BetaBmk[2][i] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * 
                Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;

            L[2][2] += Aweight[4][i] * Aweight[4][i] * 
                Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;

            L[2][2] += Aweight[5][i] * Aweight[5][i] * 
                Fix3_ExpDecay (2. * BetaBmk[2][i], t) * t;

            L[0][2] *= exp (-(BetaBmk[0][i] + BetaBmk[2][i]) * t);

            L[0][2] += 2. * Aweight[0][i] * Aweight[3][i] * 
                Fix3_ExpDecay ((BetaBmk[0][i]+BetaBmk[2][i]), t) * t;

            L[1][2] *= exp (-(BetaBmk[1][i] + BetaBmk[2][i]) * t);

            L[1][2] += 2. * Aweight[1][i] * Aweight[3][i] * 
                Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;

            L[1][2] += 2. * Aweight[2][i] * Aweight[4][i] * 
                Fix3_ExpDecay ((BetaBmk[1][i]+BetaBmk[2][i]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == IEnd)
        {
            break;
        }
    }  /* for i */

    
   

    /* need to multiply with exp(-2\int_t {IEnd}^{SwapStart} beta(u)du)) */
    /*Look for start of integration bucket */ 
    EndIdx = GetDLOffset (NbBmkMr, BmkDate, IEnd, CbkHIGHER);

    for (i = EndIdx; i < NbBmkMr; i++)
    {
        EndDate = (i < NbBmkMr-1) ? MIN(SwapStart, BmkDate[i]) : SwapStart;

        /* Time to expiry in years */
        VolT[i] = Daysact (IEnd, EndDate) / 365.;

        t = ((i == EndIdx) ? VolT[EndIdx] : (VolT[i]-VolT[i-1]));
        L[0][0] *= exp (- 2 * BetaBmk[0][i] * t);

        if (NbFactor > 1)
        {
            L[1][1] *= exp(- 2 * BetaBmk[1][i] * t);
            L[0][1] *= exp(- (BetaBmk[0][i] + BetaBmk[1][i]) * t);
        }
        if(NbFactor > 2)
        {
            L[0][2] *= exp(- (BetaBmk[0][i] +  BetaBmk[2][i]) * t);
            L[1][2] *= exp(- (BetaBmk[1][i] +  BetaBmk[2][i]) * t);
            L[2][2] *= exp(- 2 * BetaBmk[2][i] * t);
        }

          /* End of integration loop */
        if (EndDate == SwapStart)
        {
            break;
        }
    }

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
        sprintf (ErrorMsg, "GenIndexVol: problem in integrating  fwd index volatility "
                 "bet dates %ld  %ld",IStart, IEnd);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    M = sqrt (M / Int01);

    IntQ = Daysact (VolBaseDate, IEnd) / 365.0;

     /* find QLeft, QRight, FwdSh */
    if ( QInterp ( IntQ,
          &QLeft,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->QLeftTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }

    if ( QInterp (IntQ,
          &QRight,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->QRightTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }


    if ( QInterp (IntQ,
          &FwdShift,
          mktvol_data->BaseDate,
          mktvol_data->SmileDate,
          mktvol_data->FwdShiftTD,
          mktvol_data->NbSmileDates) == FAILURE)
    {
        goto RETURN;
    }

    /* Convert from q measure back to market */
    if (Conv_XSpaceVol_To_Market(ParYield,
                                Vol,        
                                VolTypeFlag,
                                M,
                                Int01,
                                FwdShift,
                                QLeft,
                                QRight,
                                0.0,
                                1.0) == FAILURE) 
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_GenIndexVol_TimeDep */

