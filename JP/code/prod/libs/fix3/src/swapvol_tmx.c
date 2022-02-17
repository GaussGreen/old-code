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
#include "q3_tmx.h"


#ifndef    MIN_SPOT_VOL_RATIO
#define    MIN_SPOT_VOL_RATIO    0.2
#endif

#ifndef    SPOT_VOL_FILTER_AMOUNT
#define    SPOT_VOL_FILTER_AMOUNT    0.95
#endif



/*****   Fix3_BFactor_Tmx    *********************************************/
/*
*       Determine the Piterbarg weight for Libor Atm Vol bootstrapping
*/
int        Fix3_BFactor_Tmx (
                    double  *BFac,         /* (O) Piterbarg Weight          */                 
                    long    SwapSt,        /* (I) Underlying swap start     */
                    long    SwapMat,       /* (I) Underlying swap maturity  */
                    char    DCC,           /* (I) Underlying day count conv.*/
                    char    Freq,          /* (I) Underlying frequency      */
                    MKTVOL_DATA*    mktvol_data,  /* (I) Volatility data    */
                    T_CURVE const*  crv)          /* (I) zero curve         */
{
    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    long        TempDates[2];         /* To construct cpn list   */

    double  ParLibor0, AnnuityLibor0;  /* Annuity prices                    */

    double  fac, TtoCpn;

    int     i,j;
    int     status = FAILURE;    /* Error status = FAILURE initially   */

    /* temporary variables for volatility data parameters */
    double const*   Beta;         /* (I) Mean reversions               */
    double          Bbq;          /* (I) Backbone parameter            */
    double          VolNorm;      /* (I) Normal volatility in backbone */
    double          VolLogn;      /* (I) Lognorm volatility in backbone*/
    int             NbFactor;     /* (I) NbFactors                     */
    double  ParYield0;            /* (I) Fwd Yield                     */
    double  Annuity0;             /* (I) Fwd Annuity                   */

    /* set temporary variables equal to mktvol_data parameters */
    Beta     = mktvol_data->Beta;
    Bbq      = mktvol_data->Bbq;
    VolNorm  = mktvol_data->VolNorm;
    VolLogn  = mktvol_data->VolLogn;
    NbFactor = mktvol_data->NbFactor;


    if (ParYieldFromDates(&ParYield0,
                          &Annuity0,
                          SwapSt,
                          SwapMat,
                          DCC,
                          Freq,
                          'F',
                          crv) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Generate coupon dates for the underlying */
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

    /* Initialize the B Factor  */      
    for (i=0; i < NbFactor; i++)
        BFac[i] = 0.0;

    /* 1-Factor case */
    if (NbFactor == 1)
    {
        for (j = 1; j < CpnEventList->NbEntries; j++)
        {
            /* Compute par & annuity libor for the current interval */
            if (ParYieldFromDates (&(ParLibor0),
                                   &(AnnuityLibor0),
                                   CpnEventList->Dates[j-1],
                                   CpnEventList->Dates[j],
                                   DCC,
                                   Freq,
                                  'F',
                                   crv) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
            /* Decay fac */
            TtoCpn      = Daysact(SwapSt, CpnEventList->Dates[j-1]) / 365.;
            fac         = exp(-Beta[0] * TtoCpn);
            BFac[0]    += ParLibor0 * AnnuityLibor0 * fac;
        }  /* for j */
        
        BFac[0] /= (ParYield0 * Annuity0);
    }

    else if (NbFactor == 2)
    {
        ;
    }
    else if (NbFactor == 3)
    {
        ;
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    DrFreeEventList(CpnEventList);

    return (status);

}  /* Fix3_BFactor_Tmx */



/*****  Fix3_SpotVol_Tmx    **************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Fix3_SpotVol_Tmx (MKTVOL_DATA  *mktvol_data, /* (I/O) Volatility data*/
                          T_CURVE const* crv)        /* (I) Zero curve       */
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
    double  lambda = 1;           /* Relative weight, = 1 for no calib case */
    double  lambdaNew = 1;        /* Relative weight, current value         */
    double  lambdaI = 0;          /* Relative weight, first value           */
    long    NbSml;                /* Number of smile pts                    */
    double  SwapSmile[NBVOLPARS][MAXNBDATE];
    double  SmlT[MAXNBDATE];      /* Smile expiries                         */   
    double  NmrT;                 /* Numeraire date expiry                  */
    double  totVol;               /* totvol for numeraire date              */
    int     intvl;

    int     LastCalibIdx;         /* Index of last calibreted vol point     */
    int     i, k, p, q, s ,j;
    int     status = FAILURE;     /* Error status = FAILURE initially       */

    char    ErrorMsg[MAXBUFF];    /* Error message                          */   
    
    int     NbNmr;                 /* (I) Number of numeraire dates       */
    long    *NmrDate;              /* (I) Numeraire dates                 */
    int     *SmlIdx;               /* (I) Smile Index                     */
    long    VolBaseDate;           /* (I) Volatility curve base date      */
    int     NbVol;                 /* (I) Nb of points in vol curve       */
    long    *VolDate;              /* (I) Volatility dates                */
    int     *VolUsed;              /* (I) TRUE if vol used in calibration */
    char    Freq;                  /* (I) Frequency of underlying rate    */
    char    DCC;                   /* (I) Day count convention            */
    long    *SwapSt;               /* (I) Underlying swap start           */
    long    *SwapMat;              /* (I) Underlying swap maturity        */
    int     NbFactor;              /* (I) Number of factors               */
    double  *Alpha;                /* (I) Relative size factors           */
    double  *Beta;                 /* (I) Mean reversions                 */
    double  *Rho;                  /* (I) Correlation between factors     */
    int     SkipFlag;              /* (I) Skip calibration failure points */

    NbNmr       = mktvol_data->NbNmr;
    NmrDate     = mktvol_data->NmrDate;
    SmlIdx      = mktvol_data->SmlLiqDate;
    VolBaseDate = mktvol_data->BaseDate;
    NbVol       = mktvol_data->NbVol;
    VolDate     = mktvol_data->VolDate;
    VolUsed     = mktvol_data->VolUsed;
    Freq        = mktvol_data->Freq;
    DCC         = mktvol_data->DCC;
    SwapSt      = mktvol_data->SwapSt;
    SwapMat     = mktvol_data->SwapMat;
    NbFactor    = mktvol_data->NbFactor;
    Alpha       = mktvol_data->Alpha;
    Beta        = mktvol_data->Beta;
    Rho         = mktvol_data->Rho;
    SkipFlag    = mktvol_data->SkipFlag;

    /*
     *  Constant Aweight numbers (determined historically).
     *  NOTE: alphas are NOT normalized, but aweights ARE
     */
    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++) 
            aw[p][q] = 0.;

    aw[0][0]  = 1.0;


    for (i = 0; i < NbVol; i++)
    {
        /* Only used vol points */
        if (!VolUsed[i]) continue;

        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, VolDate[i]) / 365.;

        if (ParYieldFromDates (&(pY[i]),
                               &Annuity,
                               SwapSt[i],
                               SwapMat[i],
                               DCC,
                               Freq,
                               'F',
                               crv) == FAILURE)
        {
            goto RETURN;
        }


        /* B factor */
        if (Fix3_BFactor_Tmx (   
                        B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        mktvol_data,
                        crv) == FAILURE)
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
        
   
        y = T * SQUARE(mktvol_data->Vol[i]);

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * Fix3_ExpDecay (2. * Beta[0], t) * t;
  
        
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
            mktvol_data->Aweight[0][k] = lambda * aw[0][0];            
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


    }  /* for i */


    for (k = LastCalibIdx+1; k < NbVol; k++)
    {
        mktvol_data->Aweight[0][k] = lambda * aw[0][0];        
  
    }  /* for k */


    /* 
     *  DEDUCE NOW LIBOR SMILE.
     */

    if (mktvol_data->NmrLibVol == NULL)
    {
        status = SUCCESS;
        goto RETURN;
    }

    for(p = 0; p < 3; p++)
        for(q = 0; q < 3; q++)
        {
            L[p][q] = 0.;
            D[p][q] = 0.;
        }


    
    /* From Swap Smile to Libor Smile */
    NbSml = 0;
    for (i=0; i < NbVol; i++)
    {
        if (SmlIdx[i])
        {
            for (s=1; s < NBVOLPARS; s++)
            {
                SwapSmile[s][NbSml] = mktvol_data->Smile[s][i];
            }

            SmlT[NbSml] = Daysact (VolBaseDate, VolDate[i]) / 365.;

            NbSml ++;
        }
    }

    /* Curve frequency */
    intvl         = 12 / Conv_Freq(Freq);

    j = 0;
    for (i=1 ;i < NbNmr - 1; i++)
    {
        double RegParLibor, IrrParLibor;
        double RegZerotoMat, IrrZerotoMat;

        long RegLibMat = Nxtmth(NmrDate[i], intvl,1L);
        long IrrLibMat = NmrDate[i+1];


        NmrT = Daysact (NmrDate[0], NmrDate[i]) / 365.;

        while ((j < NbVol - 1) && (NmrT >= VolT[j] + ERROR))
            j++;

        /* Deduce OU standard deviation at NmrDate Ti */
        t = Daysact (NmrDate[i-1], NmrDate[i]) / 365.;


        if (NbFactor == 1)
        {
            L[0][0] *= exp (-2. * Beta[0] * t);
            L[0][0] += mktvol_data->Aweight[0][j] * mktvol_data->Aweight[0][j] * 
                Fix3_ExpDecay (2. * Beta[0], t) * t;       
        }
        else if (NbFactor == 2)
        {
            ;
        }
        else if (NbFactor == 3)
        {
            ;
        }



        /* Regular Libor */
        if (ParYieldFromDates (&(RegParLibor),
                               &Annuity,
                               NmrDate[i],
                               RegLibMat,
                               DCC,
                               Freq,
                               'F',
                               crv) == FAILURE)
        {
            goto RETURN;
        }


        /* Irregular Libor */
        if (ParYieldFromDates (&(IrrParLibor),
                               &Annuity,
                               NmrDate[i],
                               IrrLibMat,
                               DCC,
                               Freq,
                               'F',
                               crv) == FAILURE)
        {
            goto RETURN;
        }


        /* Zero to Regular libor maturity */
        RegZerotoMat = GetZeroPrice (RegLibMat, crv);
        
        /* Zero to Irregular libor maturity */
        IrrZerotoMat = GetZeroPrice (IrrLibMat, crv);
 
        /* Libor ATM Vol */
        mktvol_data->NmrLibVol[0][i]   = sqrt(L[0][0]/NmrT) * 
            RegParLibor/IrrParLibor * RegZerotoMat/IrrZerotoMat;
        

        /* Linearly interpolate the smile */
        for (s=1; s<NBVOLPARS; s++)
        {
            tableinterp(NmrT,
                        &(mktvol_data->NmrLibVol[s][i]),
                        SmlT,
                        SwapSmile[s],
                        NbSml);
        }
    }


    /* Apply Caps & Floors to smile to avoid MultiQ calibration errors */
    for (i=1 ;i < NbNmr - 1; i++)
    {
        NmrT   = Daysact (NmrDate[0], NmrDate[i]) / 365.;
        totVol = mktvol_data->NmrLibVol[0][i] * sqrt(NmrT);

        /* Skew limits */
        mktvol_data->NmrLibVol[1][i] = MAX(mktvol_data->NmrLibVol[1][i], 1. - Q3TMX_MAX_SKEW/totVol);
        mktvol_data->NmrLibVol[1][i] = MIN(mktvol_data->NmrLibVol[1][i], 1. + Q3TMX_MAX_SKEW/totVol);

        /* vvol limits */
        mktvol_data->NmrLibVol[2][i] = MAX(mktvol_data->NmrLibVol[2][i], TINY);
        mktvol_data->NmrLibVol[2][i] = MIN(mktvol_data->NmrLibVol[2][i], 
                                 MIN(Q3TMX_MAX_VOL_VOL/(sqrt(NmrT/3.)),TMX_MAX_VOL_VOL));
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

}  /* Fix3_SpotVol_Tmx */



/*****  Fix3_Interp_SpotVol_Tmx *********************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Fix3_Interp_SpotVol_Tmx (
                FIX3_TREE_DATA*  tree_data,   /* (I/ O) Tree data   */
                MKTVOL_DATA*     mktvol_data) /* (I)Volatility data */           
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  t, T;
    int     i, j, k;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

    double  **AweightCurve;         /* (O) Factors aweight curves        */
    long    VolBaseDate;            /* (I) Volatility curve base date    */
    int     NbVol;                  /* (I) Number of spot vol points     */
    long    *VolDate;               /* (I) Spot vol dates                */
    double  Aweight[6][MAXNBDATE];  /* (I) Aweight curve of each factor  */
    long    NbFactor;               /* (I) Number of factors             */
    double  *Length;                /* (I) Length of each time step      */
    int     NbTP;                   /* (I) Total number of time points   */

    AweightCurve = tree_data->Aweight;
    VolBaseDate  = mktvol_data->BaseDate;
    NbVol        = mktvol_data->NbVol;
    VolDate      = mktvol_data->VolDate;
    NbFactor     = mktvol_data->NbFactor;
    Length       = tree_data->Length;
    NbTP         = tree_data->NbTP;

    if (NbFactor == 1)
        NbAweight = 1;
    else if (NbFactor == 2)
        NbAweight = 3;
    else if (NbFactor == 3)
        NbAweight = 6;
    else
        NbAweight = 0;

    /* set temporary Aweight array to mktvol_data Aweight */
    for (i = 0; i < NbVol; i++)
    {
        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][i]=mktvol_data->Aweight[k][i];
        }
    }
    

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

}  /* Fix3_Interp_SpotVol_Tmx */


/*****  Fix3_IndexVol_Tmx    ******************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Fix3_IndexVol_Tmx (
            double          *Vol,            /* (O) Index vol curve              */
            long            OptExp,          /* (I) Option expiry                */
            long            SwapSt,          /* (I) Underlying swap start        */
            long            SwapMat,         /* (I) Underlying swap maturity     */
            char            Freq,            /* (I) Frequency of underlying rate */
            char            DCC,             /* (I) Day count convention         */
            MKTVOL_DATA*    mktvol_data,     /* (I) Volatility data              */
            T_CURVE const*  crv)             /* (I) Tree data                    */
{
    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  t;                 /* Time between two consecutive expiries */
    double  T=0.;              /* Time to current expiry                */
    double  ParYield;
    double  Annuity;
    double  Tenor, RegTenor;
    double  RegParLibor;
    double  RegZerotoMat, IrrZerotoMat;

    double  L[3][3];           /* Integrals of factor spot vol          */
    double  B[3];              /* B in Christian memo                   */
    double  M;                 /* Total variance                        */

    long    EndDate;           /* End of current integration bucket     */
    long    RegLibMat;  
    int     i, j;
    int     intvl;

    int     status = FAILURE;  /* Error status = FAILURE initially      */

    int           NbVol;                /* (I) Number of spot vol points    */
    long          VolBaseDate;          /* (I) Volatility curve base date   */
    long          *VolDate;             /* (I) Spot vol dates               */
    int           NbFactor;             /* (I) Number of factors            */
    double        *Beta;                /* (I) Mean reversions              */

    NbVol       = mktvol_data->NbVol;
    VolBaseDate = mktvol_data->BaseDate;
    VolDate     = mktvol_data->VolDate;
    NbFactor    = mktvol_data->NbFactor;
    Beta        = mktvol_data->Beta;


    /* Avoid spot index expiration */
    if (SwapSt <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
    }
    
    /* Find par yield */
    if (ParYieldFromDates (&ParYield,
                           &Annuity,
                           SwapSt,
                           SwapMat,
                           DCC,
                           Freq,
                           'F',
                           crv) == FAILURE)
    {
        goto RETURN;
    }
    
    /* Find B coefficient */
    if (Fix3_BFactor_Tmx (
                 B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 mktvol_data,
                 crv) == FAILURE)        
    {    
        goto RETURN;        
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
        L[0][0] += mktvol_data->Aweight[0][i] * mktvol_data->Aweight[0][i] * 
            Fix3_ExpDecay (2. * Beta[0], t) * t;

        /* End of integration loop */
        if (EndDate == SwapSt)
        {
            break;
        }
    }  /* for i */
    
    M = B[0] * B[0] * L[0][0];
        
    if (M < SQUARE(TINY))
    {
        DR_Error ("IndexVol: problem in integrating index volatility "
                 "at date %ld (IndexVol)", SwapSt);
        goto RETURN;
    }

    *Vol = sqrt (M / T);

    /* If the swap tenor is less than 1 freq, perform `reduced' libor correction */
    intvl     = 12 / Conv_Freq(Freq);
    RegLibMat = Nxtmth(SwapSt, intvl,1L);
    Tenor     = Daysact (SwapSt, SwapMat) / 365.;
    RegTenor  = Daysact (SwapSt, RegLibMat) / 365.;

    if (Tenor < RegTenor - TINY)
    {
        /* Regular Libor */
        if (ParYieldFromDates (&(RegParLibor),
                               &Annuity,
                               SwapSt,
                               RegLibMat,
                               DCC,
                               Freq,
                               'F',
                               crv) == FAILURE)
        {
            goto RETURN;
        }

        /* Zero to Regular libor maturity */
        RegZerotoMat = GetZeroPrice (RegLibMat, crv);
        
        /* Zero to Irregular libor maturity */
        IrrZerotoMat = GetZeroPrice (SwapMat, crv);

        *Vol *= RegParLibor/ParYield * RegZerotoMat/IrrZerotoMat;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_IndexVol_Tmx  */


/*********  Fix3_IndexLimits_Tmx ****************************************/ 
/*
 *  Estimates State Variable limits using a lognormal model.
 */
int     Fix3_IndexLimits_Tmx (
                 double      *MinYield,      /* (O) : Min yield (normalized)  */
                 double      *MaxYield,      /* (O) : Max yield (normalized)  */
                 long         NbStDev,       /* (I) : Nb of st dev            */         
                 long         RateReset,     /* (I) : Reset date              */           
                 long         SwapSt,        /* (I) : Swap start              */                 
                 int          IdxMat,        /* (I) : Rate Maturity           */                 
                 char         DCC,           /* (I) : Day count fraction      */                 
                 char         Freq,          /* (I) : Swap Frequency          */                 
                 MKTVOL_DATA *mvd,           /* (I) : mktvol data structure   */                 
                 T_CURVE const *t_curve)     /* (I) : t curve                 */            
{
    int     i, NbVol;     
    long    ValueDate;
    long    SwapMat;
    double  Expiry, ParYield, Annuity;
    double  *SwapStT;  
    double  Vol, Skew, VoV, Q, totVol, totVovSq, vovSq2Vol, sigMQ, muMQ, qL, qR; 
    int     status = FAILURE;

    /* Initializations */
    NbVol     = mvd->NbVol;
    ValueDate = t_curve->ValueDate;
    SwapMat   = Nxtmth (SwapSt,IdxMat,1L);
    Expiry    = Daysact(ValueDate, SwapSt) / 365.;    

    /* Memory Allocation */
    SwapStT = (double *) DR_Array(DOUBLE,0, NbVol);
    
    /* Market vol dates */
    for (i = 0; i < NbVol; i++)
    {
        SwapStT[i] = Daysact (ValueDate, mvd->SwapSt[i])/365.;
    }
   
    /* Compute forward from curve for index yield */  
    if (ParYieldFromDates (&ParYield,
                              &Annuity,
                              SwapSt, 
                              SwapMat,       
                              DCC,           
                              Freq,          
                              'F',  
                              t_curve) == FAILURE)
    {
        goto RETURN;
    }
    
    /* Annualized index volatility at reset date  for yield */  
    if (Fix3_IndexVol_Tmx (
                  &Vol,
                  SwapSt,/* option expiry assumed equal to swap start*/
                  SwapSt,                  
                  SwapMat,                  
                  Freq,                  
                  DCC,                                
                  mvd,
                  t_curve) == FAILURE)
    {
        goto RETURN;
    }

    /* Interpolate Smile info */
    tableinterp(Expiry, 
                &Skew, 
                SwapStT,
                mvd->Smile[1],
                NbVol);
    
    tableinterp(Expiry, 
                &VoV, 
                SwapStT,
                mvd->Smile[2],
                NbVol);
    
    /* Use MQ Initial Guess formulas for limit calculation */
    totVol    = Vol * sqrt(Expiry);
    totVovSq  = SQUARE(VoV) * Expiry;
    vovSq2Vol = (totVol < Q3TMX_MIN_VOL) ? 0 : (totVovSq / totVol);
    Q         = 1.0 - Skew;
    sigMQ     = totVol / (1.0 + totVovSq / sqrt(3.));
    muMQ      = -0.5 * Q * totVol * sigMQ;
    qL        = Q + vovSq2Vol * (Q / 3. - 1.);
    qR        = qL + 2. * vovSq2Vol;

    /* Left side */
    if (fabs(qL) > TINY)
    {
        *MinYield = ParYield * (1.0 + (exp(qL * (muMQ - NbStDev * sigMQ)) - 1.0) / qL);
    }
    else
    {
        *MinYield = ParYield * (1.0 + muMQ - NbStDev * sigMQ);
    }

    /* Right side */
    if (fabs(qR) > TINY)
    {
        *MaxYield = ParYield * (1.0 + (exp(qR * (muMQ + NbStDev * sigMQ)) - 1.0) / qR);
    }
    else
    {
        *MaxYield = ParYield * (1.0 + muMQ + NbStDev * sigMQ);
    }

    status = SUCCESS;

RETURN:

    Free_DR_Array(SwapStT, DOUBLE,  0, NbVol);
    
    return (status);

}  /* Fix3_IndexLimits_Tmx */


/*****  GenIndexVol_Tmx    *************************************************/
/*
*       Integrates spot volatility curve between two time points falling 
*       before or on the swap start date.
*/
int   Fix3_GenIndexVol_Tmx (
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
    double  ParYield;
    double  Annuity;
    double  Tenor, RegTenor;
    double  RegParLibor;
    double  RegZerotoMat, IrrZerotoMat;

    double  L[3][3];           /* Integrals of factor spot vol               */
    double  B[3];              /* B in Christian memo                        */
    double  M;                 /* Total variance                             */

    long    EndDate;           /* End of current integration bucket          */
    long    RegLibMat;  
    int     StartIdx=0;
    int     i, j;
    int     intvl;

    int     status = FAILURE;  /* Error status = FAILURE initially           */
    char    ErrorMsg[MAXBUFF]; /* Error message                              */

    /* local variables -- extract information from data structures */
    int     NbVol       = mktvol_data->NbVol;
    long    VolBaseDate = mktvol_data->BaseDate;
    long    *VolDate    = mktvol_data->VolDate;
    int     NbFactor    = mktvol_data->NbFactor;
    double  *Beta       = mktvol_data->Beta;
    

    /* Avoid spot index expiration; note we use base date of zero */
    /* curve as there is no volatility information available.     */
    if ((IEnd <= VolBaseDate) || (IStart == IEnd))
    {
        *Vol = 0.00;
        return(SUCCESS);
    }

    Int01 = Daysact (IStart, IEnd)      / 365.;
    Int12 = Daysact (IEnd, SwapStart)   / 365.;
    Int02 = Daysact (IStart, SwapStart) / 365.;  
    
    /* Find par yield */
    if (ParYieldFromDates (&ParYield,
                           &Annuity,
                           SwapStart,
                           SwapMat,
                           DCC,
                           Freq,
                           'F',
                           crv) == FAILURE)
    {
        goto RETURN;
    }
    
    /* Find B coefficient */
    if (Fix3_BFactor_Tmx (B,
                          SwapStart,
                          SwapMat,
                          DCC,
                          Freq,
                          mktvol_data,
                          crv) == FAILURE)        
    {    
        goto RETURN;        
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
        L[0][0] += mktvol_data->Aweight[0][i] * mktvol_data->Aweight[0][i] * Fix3_ExpDecay (2. * Beta[0], t) * t;

        /* End of integration loop */
        if (EndDate == IEnd)
        {
            break;
        }

    }  /* for i */

    M = B[0] * B[0] * L[0][0] * exp(-2. * Beta[0] * Int12);
    
    if (M < SQUARE(TINY))
    {
        sprintf (ErrorMsg, "GenIndexVol: problem in integrating  fwd index volatility "
                 "bet dates %ld  %ld",IStart, IEnd);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    *Vol = sqrt (M / Int01);

    /* If the swap tenor is less than 1 freq, perform `reduced' libor correction */
    intvl     = 12 / Conv_Freq(Freq);
    RegLibMat = Nxtmth(SwapStart, intvl,1L);
    Tenor     = Daysact (SwapStart, SwapMat) / 365.;
    RegTenor  = Daysact (SwapStart, RegLibMat) / 365.;

    if (Tenor < RegTenor - TINY)
    {
        /* Regular Libor */
        if (ParYieldFromDates (&(RegParLibor),
                               &Annuity,
                               SwapStart,
                               RegLibMat,
                               DCC,
                               Freq,
                               'F',
                               crv) == FAILURE)
        {
            goto RETURN;
        }

        /* Zero to Regular libor maturity */
        RegZerotoMat = GetZeroPrice (RegLibMat,
                                     crv);
        
        /* Zero to Irregular libor maturity */
        IrrZerotoMat = GetZeroPrice (SwapMat,
                                     crv);

        *Vol *= RegParLibor/ParYield * RegZerotoMat/IrrZerotoMat;
    }

    status = SUCCESS;

    RETURN:

    return (status);


} /* Fix3_GenIndexVol_Tmx */
