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
#include "tmx123head.h"
#include "q3.h"


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



/*****   BFactor    *****************************************************/
/*
*       Determine the Piterbarg weight for Libor Atm Vol bootstrapping
*/
int        BFactor (double  *BFac,         /* (O) Piterbarg Weight          */                 
                    long    SwapSt,        /* (I) Underlying swap start     */
                    long    SwapMat,       /* (I) Underlying swap maturity  */
                    char    DCC,           /* (I) Underlying day count conv.*/
                    char    Freq,          /* (I) Underlying frequency      */
                    double  ParYield0,     /* (I) Fwd Yield                 */
                    double  Annuity0,      /* (I) Fwd Annuity               */
                    double  *Beta,         /* (I) Mean reversion            */
                    long    NbFactor,      /* (I) Number of factor          */
                    int     NbZero,        /* (I) Number of zeros           */
                    double  *Zero,         /* (I) Zero rates                */
                    long    *ZeroDate,     /* (I) Zero maturity dates       */
                    long    ValueDate)     /* (I) Zero curve base date      */
{
    EVENT_LIST  *CpnEventList = NULL; /* Event list for cpn pmts */
    long        TempDates[2];         /* To construct cpn list   */

    double  ParLibor0, AnnuityLibor0;  /* Annuity prices                     */

    double  fac, TtoCpn;

    int     i,j;
    int     status = FAILURE;    /* Error status = FAILURE initially   */

    if (SwapMat > ZeroDate[NbZero-1])
    {        
        DR_Error ("Not enough zeros to calculate B (BFactor)!");
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
            if (Par_Yield_From_Dates (&(ParLibor0),
                                      &(AnnuityLibor0),
                                      CpnEventList->Dates[j-1],
                                      CpnEventList->Dates[j],
                                      DCC,
                                      Freq,
                                      'F',
                                      NbZero,
                                      Zero,
                                      ZeroDate,
                                      ValueDate) == FAILURE)
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
       double  NmrLibSmile[NBVOLPARS][MAXNBDATE],/* (O) Libor ATM volatility */
       int     NbNmr,                 /* (I) Number of numeraire dates       */
       long    *NmrDate,              /* (I) Numeraire dates                 */
       int     *SmlIdx,               /* (I) Smile Index                     */
       long    VolBaseDate,           /* (I) Volatility curve base date      */
       int     NbVol,                 /* (I) Nb of points in vol curve       */
       long    *VolDate,              /* (I) Volatility dates                */
       double  Vol[NBVOLPARS][MAXNBDATE],/* (I) Vol curve                    */
       int     *VolUsed,              /* (I) TRUE if vol used in calibration */
       char    Freq,                  /* (I) Frequency of underlying rate    */
       char    DCC,                   /* (I) Day count convention            */
       long    *SwapSt,               /* (I) Underlying swap start           */
       long    *SwapMat,              /* (I) Underlying swap maturity        */
       int     NbFactor,              /* (I) Number of factors               */
       double  *Alpha,                /* (I) Relative size factors           */
       double  *Beta,                 /* (I) Mean reversions                 */
       double  *Rho,                  /* (I) Correlation between factors     */
       int     SkipFlag,              /* (I) Skip calibration failure points */
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


        /* B factor */
        if (BFactor (   B[i],
                        SwapSt[i],
                        SwapMat[i],
                        DCC,
                        Freq,
                        pY[i],
                        Annuity,
                        Beta,
                        NbFactor,
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
        
   
        y = T * SQUARE(Vol[0][i]);

        y -= D[0][0] * D[0][0] * L[0][0] * exp(-2. * Beta[0] * t);
        
        M = D[0][0] * D[0][0] * ExpDecay (2. * Beta[0], t) * t;
  
        
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


    }  /* for i */


    for (k = LastCalibIdx+1; k < NbVol; k++)
    {
        Aweight[0][k] = lambda * aw[0][0];        
  
    }  /* for k */


    /* 
     *  DEDUCE NOW LIBOR SMILE.
     */

    if (NmrLibSmile == NULL)
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
                SwapSmile[s][NbSml] = Vol[s][i];
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
            L[0][0]            *= exp (-2. * Beta[0] * t);
            L[0][0]            += Aweight[0][j] * Aweight[0][j] * ExpDecay (2. * Beta[0], t) * t;       
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
        if (Par_Yield_From_Dates (&(RegParLibor),
                                  &Annuity,
                                  NmrDate[i],
                                  RegLibMat,
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


        /* Irregular Libor */
        if (Par_Yield_From_Dates (&(IrrParLibor),
                                  &Annuity,
                                  NmrDate[i],
                                  IrrLibMat,
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


        /* Zero to Regular libor maturity */
        RegZerotoMat = ZeroPrice (RegLibMat,
                                  ZeroBaseDate,
                                  NbZero,
                                  ZeroDate,
                                  Zero);


        
        /* Zero to Irregular libor maturity */
        IrrZerotoMat = ZeroPrice (IrrLibMat,
                                  ZeroBaseDate,
                                  NbZero,
                                  ZeroDate,
                                  Zero);


        /* Libor ATM Vol */
        NmrLibSmile[0][i]   = sqrt(L[0][0]/NmrT) * RegParLibor/IrrParLibor * RegZerotoMat/IrrZerotoMat;
        

        /* Linearly interpolate the smile */
        for (s=1; s<NBVOLPARS; s++)
        {
            tableinterp(NmrT,
                        &(NmrLibSmile[s][i]),
                        SmlT,
                        SwapSmile[s],
                        NbSml);
        }
    }


    /* Apply Caps & Floors to smile to avoid MultiQ calibration errors */
    for (i=1 ;i < NbNmr - 1; i++)
    {
        NmrT   = Daysact (NmrDate[0], NmrDate[i]) / 365.;
        totVol = NmrLibSmile[0][i] * sqrt(NmrT);

        /* Skew limits */
        NmrLibSmile[1][i] = MAX(NmrLibSmile[1][i], 1. - Q3_MAX_SKEW/totVol);
        NmrLibSmile[1][i] = MIN(NmrLibSmile[1][i], 1. + Q3_MAX_SKEW/totVol);

        /* vvol limits */
        NmrLibSmile[2][i] = MAX(NmrLibSmile[2][i], TINY);
        NmrLibSmile[2][i] = MIN(NmrLibSmile[2][i], 
                                 MIN(Q3_MAX_VOL_VOL/(sqrt(NmrT/3.)),TMX_MAX_VOL_VOL));
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
       long    VolBaseDate,            /* (I) Volatility curve base date    */
       int     NbVol,                  /* (I) Number of spot vol points     */
       long    *VolDate,               /* (I) Spot vol dates                */
       double  Aweight[6][MAXNBDATE],  /* (I) Aweight curve of each factor  */
       long    NbFactor,               /* (I) Number of factors             */
       double  *Length,                /* (I) Length of each time step      */
       int     NbTP)                   /* (I) Total number of time points   */
{

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  t, T;
    int     i, j, k;
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

}  /* Interp_SpotVol */


/*****  IndexVol    *********************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     IndexVol (
          double        *Vol,                 /* (O) Index vol curve              */
          long          SwapSt,               /* (I) Underlying swap start        */
          long          SwapMat,              /* (I) Underlying swap maturity     */
          char          Freq,                 /* (I) Frequency of underlying rate */
          char          DCC,                  /* (I) Day count convention         */
          int           NbVol,                /* (I) Number of spot vol points    */
          long          VolBaseDate,          /* (I) Volatility curve base date   */
          long          *VolDate,             /* (I) Spot vol dates               */
          double        Aweight[6][MAXNBDATE],/* (I) Aweight curve                */
          int           NbFactor,             /* (I) Number of factors            */
          double        *Beta,                /* (I) Mean reversions              */
          int           NbZero,               /* (I) Number of zeros              */       
          double        *Zero,                /* (I) Zero rates                   */       
          long          *ZeroDate,            /* (I) Zero maturity dates          */       
          long          ZeroBaseDate)         /* (I) Zero curve base date         */
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

    /* Avoid spot index expiration */
    if (SwapSt <= VolBaseDate) 
    {
        *Vol = 0.00;        
        return(SUCCESS);        
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
    
    /* Find B coefficient */
    if (BFactor (B,
                 SwapSt,
                 SwapMat,
                 DCC,
                 Freq,
                 ParYield,
                 Annuity,
                 Beta,
                 NbFactor,
                 NbZero,
                 Zero,
                 ZeroDate,
                 ZeroBaseDate) == FAILURE)        
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
        L[0][0] += Aweight[0][i] * Aweight[0][i] * ExpDecay (2. * Beta[0], t) * t;

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
        if (Par_Yield_From_Dates (&(RegParLibor),
                                  &Annuity,
                                  SwapSt,
                                  RegLibMat,
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

        /* Zero to Regular libor maturity */
        RegZerotoMat = ZeroPrice (RegLibMat,
                                  ZeroBaseDate,
                                  NbZero,
                                  ZeroDate,
                                  Zero);
        
        /* Zero to Irregular libor maturity */
        IrrZerotoMat = ZeroPrice (SwapMat,
                                  ZeroBaseDate,
                                  NbZero,
                                  ZeroDate,
                                  Zero);

        *Vol *= RegParLibor/ParYield * RegZerotoMat/IrrZerotoMat;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* IndexVol */


/*********  StateVar ****************************************************/ 
/*
 *  Estimates State Variable limits using a lognormal model.
 */
int     StateVar(double      *MinYield,      /* (O) : Min yield (normalised)  */
                 double      *MaxYield,      /* (O) : Max yield (normalised)  */
                 double       NbStDev,       /* (I) : Nb of st dev            */                 
                 long         SwapSt,        /* (I) : Swap start              */                 
                 long         SwapMat,       /* (I) : Swap Maturity           */                 
                 char         DCC,           /* (I) : Day count fraction      */                 
                 char         Freq,          /* (I) : Swap Frequency          */                 
                 int          NbFactor,      /* (I) : No factors in the tree  */
                 MKTVOL_DATA *mvd,           /* (I) : mktvol data struc       */                 
                 T_CURVE     *t_curve)       /* (I) : t curve                 */            
{
    int     i, NbVol;     
    long    ValueDate;
    double  Expiry, ParYield, Annuity;
    double  *SwapStT;  
    double  Vol, Skew, VoV, Q, totVol, totVovSq, vovSq2Vol, sigMQ, muMQ, qL, qR; 
    int     status = FAILURE;

    /* Initialisations */
    NbVol     = mvd->NbVol;
    ValueDate = t_curve->ValueDate;
    Expiry    = Daysact(ValueDate, SwapSt) / 365.;    

    /* Memory Allocation */
    SwapStT = (double *) DR_Array(DOUBLE,0, NbVol);
    
    /* Market vol dates */
    for (i = 0; i < NbVol; i++)
    {
        SwapStT[i] = Daysact (ValueDate, mvd->SwapSt[i])/365.;
    }
   
    /* Compute forward from curve for index yield */  
    if (Par_Yield_From_Dates (&ParYield,
                              &Annuity,
                              SwapSt, 
                              SwapMat,       
                              DCC,           
                              Freq,          
                              'F',  
                              t_curve->NbZero,        
                              t_curve->Zero,         
                              t_curve->ZeroDate,  
                              ValueDate) == FAILURE)
    {
        goto RETURN;
    }
    
    /* Annualized index volatility at reset date  for yield */  
    if (IndexVol (&Vol,
                  SwapSt,                  
                  SwapMat,                  
                  Freq,                  
                  DCC,                  
                  NbVol,                
                  mvd->BaseDate,                  
                  mvd->VolDate,                  
                  mvd->Aweight,             
                  NbFactor,                  
                  mvd->Beta,                  
                  t_curve->NbZero,                  
                  t_curve->Zero,                  
                  t_curve->ZeroDate,                
                  ValueDate) == FAILURE)
    {
        goto RETURN;
    }

    /* Interpolate Smile info */
    tableinterp(Expiry, 
                &Skew, 
                SwapStT,
                mvd->Vol[1],
                NbVol);
    
    tableinterp(Expiry, 
                &VoV, 
                SwapStT,
                mvd->Vol[2],
                NbVol);
    
    /* Use MQ Initial Guess formulas for limit calculation */
    totVol    = Vol * sqrt(Expiry);
    totVovSq  = SQUARE(VoV) * Expiry;
    vovSq2Vol = (totVol < Q3_MIN_VOL) ? 0 : (totVovSq / totVol);
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

}  /* StateVar */
