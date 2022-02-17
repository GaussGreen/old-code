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



/*****  Smd_SpotVol_2F_Timedep    *******************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Smd_SpotVol_2F_Timedep    (
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
       MKTVOL_DATA  *mktvol_data,     /* (I) Volatility data                 */
       T_CURVE const* t_curve)        /* (I)   Zero curve                    */
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
                        mktvol_data,
                        t_curve) == FAILURE)
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
                               t_curve) == FAILURE)
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

}  /* Smd_SpotVol_2F_TimeDep */




/*****  Fix3_SpotVol_Smd    ************************************************/
/*
*       Bootstrapp one volatility curve to extract spot volatilities.
*       Although we input both VolDate, the option expiry, and SwapSt, the
*       accrual start of the underlying, this routine only works if they are
*       identical.
*/
int     Fix3_SpotVol_Smd (
                MKTVOL_DATA*   mktvol_data,  /* (I/O) Volatility data     */
                T_CURVE const* t_curve)      /* (I)   Zero curve          */
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

    long    VolBaseDate;           /* (I) Volatility curve base date      */
    int     NbVol;                 /* (I) Nb of points in vol curve       */
    long    *VolDate;              /* (I) Volatility dates                */
    double  *Vol;                  /* (I) Vol curve                       */
    int     *VolUsed;              /* (I) TRUE if vol used in calibration */
    char    Freq;                  /* (I) Frequency of underlying rate    */
    char    DCC;                   /* (I) Day count convention            */
    long    *SwapSt;               /* (I) Underlying swap start           */
    long    *SwapMat;              /* (I) Underlying swap maturity        */
    double  QLeft;                 /* (I) Left Q mapping coefficient      */
    double  QRight;                /* (I) Right Q mapping coefficient     */
    double  FwdSh;                 /* (I) Fwd shift mapping coefficient   */
    double  Bbq;                   /* (I) Backbone parameter              */
    double  VolNorm;               /* (I) Normal volatility in backbone   */
    double  VolLogn;               /* (I) Lognorm volatility in backbone  */
    int     NbFactor;              /* (I) Number of factors               */
    double  *Alpha;                /* (I) Relative size factors           */
    double  *Beta;                 /* (I) Mean reversions                 */
    double  *Rho;                  /* (I) Correlation between factors     */
    double  Dfac;                  /* (I) Correlation term structure      */
    int     SkipFlag;              /* (I) Skip calibration failure points */
    int     CalibFlag;             /* (I) Index calibration flag          */

    VolBaseDate           = mktvol_data->BaseDate;
    NbVol                 = mktvol_data->NbVol;
    VolDate               = mktvol_data->VolDate;
    Vol                   = mktvol_data->Vol;
    VolUsed               = mktvol_data->VolUsed;
    Freq                  = mktvol_data->Freq;
    DCC                   = mktvol_data->DCC;
    SwapSt                = mktvol_data->SwapSt;
    SwapMat               = mktvol_data->SwapMat;
    QLeft                 = mktvol_data->QLeft;
    QRight                = mktvol_data->QRight;
    FwdSh                 = mktvol_data->FwdShift;
    Bbq                   = mktvol_data->Bbq;
    VolNorm               = mktvol_data->VolNorm;
    VolLogn               = mktvol_data->VolLogn;
    NbFactor              = mktvol_data->NbFactor;
    Alpha                 = mktvol_data->Alpha;
    Beta                  = mktvol_data->Beta;
    Rho                   = mktvol_data->Rho;
    Dfac                  = mktvol_data->Dfac;
    SkipFlag              = mktvol_data->SkipFlag;
    CalibFlag             = mktvol_data->CalibFlag;


    for (i=0; i < NbVol; i++)
    {
        expiry       = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;
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
            expiry     = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;         

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
            expiry     = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;

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
    if (Smd_SpotVol_2F_Timedep    (
                    SmdSpotVol,
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
                    mktvol_data,
                    t_curve) == FAILURE)
    {
        goto RETURN;
    }


    /* Now populate SMD aweights */
    for (i=0; i < NbVol; i++)
    {
        expiry     = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;

        mktvol_data->Aweight[0][i] = SmdSpotVol[i];
        mktvol_data->Aweight[1][i] = SmdAlp2L[i]/SmdAlp1L[i] * SmdSpotVol[i] * SmdCorrel[i];
        mktvol_data->Aweight[2][i] = SmdAlp2L[i]/SmdAlp1L[i] * SmdSpotVol[i] * sqrt(1.0 - SmdCorrel[i] * SmdCorrel[i]);
    }


    status = SUCCESS;

RETURN:

    return(status);

    
}  /* Fix3_SpotVol_Smd */






/*****  Fix3_Interp_SpotVol_Smd  *******************************************/
/*
*       Interpolate spot volatility and correlation curves at each time line
*       point.
*/
int     Fix3_Interp_SpotVol_Smd (
                FIX3_TREE_DATA*    tree_data,   /* (I/O) Tree data       */
                MKTVOL_DATA*       mktvol_data) /* (I) Volatility data   */
{
    double  **AweightCurve;         /* (O) Factors aweight curves        */
    long    VolBaseDate;            /* (I) Volatility curve base date    */
    int     NbVol;                  /* (I) Number of spot vol points     */
    long    *VolDate;               /* (I) Spot vol dates                */
    double  Aweight[6][MAXNBDATE];  /* (I) Aweight curve of each factor  */
    double  *Length;                /* (I) Length of each time step      */
    int     NbTP;                   /* (I) Total number of time points   */

    double  VolT[MAXNBDATE];       /* Time to vol point in years */
    double  t, T;
    int     i, j, k;
    int     NbAweight;             /* Nb of orthogonal weights   */
    int     status = FAILURE;      /* Error status               */

    AweightCurve = tree_data->Aweight;
    VolBaseDate  = mktvol_data->BaseDate;
    NbVol        = mktvol_data->NbVol;
    VolDate      = mktvol_data->VolDate;
    Length       = tree_data->Length;
    NbTP         = tree_data->NbTP;

    /* 2F only at this point */
    NbAweight = 3;

    /* set temporary Aweight array to mktvol_data Aweight */
    for (i = 0; i < NbVol; i++)
    {
        for (k = 0; k < NbAweight; k++)
        {
            Aweight[k][i] = mktvol_data->Aweight[k][i];
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

}  /* Fix3_Interp_SpotVol_Smd */



/*****  Smd_IndexVol_2F_Timedep    *******************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Smd_IndexVol_2F_Timedep (
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
          MKTVOL_DATA*   mktvol_data,    /* (I) Volatility data              */
          T_CURVE const* t_curve)        /* (I) Zero curve                   */
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
                 mktvol_data,
                 t_curve) == FAILURE)
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
                              t_curve) == FAILURE)
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

}  /* Smd_IndexVol_2F_Timedep */



/*****  Fix3_IndexVol_Smd    *************************************************/
/*
*       Integrate spot volatility curve to determine index volatility.
*/
int     Fix3_IndexVol_Smd (
          double  *Vol,                  /* (O) Index vol curve              */
          long    OptExp,                /* (I) Option expiry                */
          long    SwapSt,                /* (I) Underlying swap start        */
          long    SwapMat,               /* (I) Underlying swap maturity     */
          char    Freq,                  /* (I) Frequency of underlying rate */
          char    DCC,                   /* (I) Day count convention         */
          MKTVOL_DATA*   mktvol_data,    /* (I/O) Volatility data            */
          T_CURVE const* t_curve)        /* (I)   Zero curve                 */
{

    int     CalibFlag;             /* (I) Index calibration flag       */
    int     NbVol;                 /* (I) Number of spot vol points    */
    long    VolBaseDate;           /* (I) Volatility curve base date   */
    long    *VolDate;              /* (I) Spot vol dates               */
    double  Aweight[6][MAXNBDATE]; /* (I) Aweight curve                */
    double  QLeft;                 /* (I) Left Q mapping coefficient   */
    double  QRight;                /* (I) Right Q mapping coefficient  */
    double  FwdShift;              /* (I) Fwd shift mapping coefficient*/
    double  Bbq;                   /* (I) Backbone parameter           */
    double  VolNorm;               /* (I) Normal volatility in backbone*/
    double  VolLogn;               /* (I) Lognorm volatility in backbon*/
    int     NbFactor;              /* (I) Number of factors            */
    double  *Alpha;                /* (I) Vol factors                  */
    double  *Beta;                 /* (I) Mean reversions              */
    double  *Rho;                  /* (I) Rho                          */
    double  Dfac;


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

    int    NbAweight;
    
    int    i,k;
    int    status = SUCCESS;


    /* Extract local variables from mktvol_data */
    CalibFlag     = mktvol_data->CalibFlag;
    NbVol         = mktvol_data->NbVol;
    VolBaseDate   = mktvol_data->BaseDate;
    VolDate       = mktvol_data->VolDate;
    QLeft         = mktvol_data->QLeft;
    QRight        = mktvol_data->QRight;
    FwdShift      = mktvol_data->FwdShift;
    Bbq           = mktvol_data->Bbq;
    VolNorm       = mktvol_data->VolNorm;
    VolLogn       = mktvol_data->VolLogn;
    NbFactor      = mktvol_data->NbFactor;
    Alpha         = mktvol_data->Alpha;
    Beta          = mktvol_data->Beta;
    Rho           = mktvol_data->Rho;
    Dfac          = mktvol_data->Dfac;


    /* 2F only */
    NbAweight = 3;

    /* Set temporary Aweight array to mktvol_data Aweight */
    for (i = 0; i < NbVol; i++)
    {

        for (k = 0; k < NbAweight; k++)
        {

            Aweight[k][i] = mktvol_data->Aweight[k][i];
        }
    }
    
    for (i=0; i < NbVol; i++)
    {
        expiry       = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;
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
            expiry         = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;
            
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
            expiry     = Daysact(t_curve->ValueDate, VolDate[i]) / 365.0;

            Alpha1L[i] = -SmdVolfac;
            Alpha2L[i] = sqrt(SQUARE(SmdAlp2L[i]/SmdAlp1L[i]) + SQUARE(SmdVolfac) - 2.0 * SmdCorrel[i]*SmdAlp2L[i]/SmdAlp1L[i]*SmdVolfac);
            RhoL[i]    = -(SmdAlp2L[i]/SmdAlp1L[i] * SmdCorrel[i] - SmdVolfac) / Alpha2L[i];

            SmdSpotVol[i]  = Aweight[0][i];

        }
    }

    if (Smd_IndexVol_2F_Timedep (
                  Vol,
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
                  mktvol_data,
                  t_curve) == FAILURE)
    {
        goto RETURN;
    }




    status = SUCCESS;

RETURN:

    return(status);
} /* Fix3_IndexVol_Smd */






















