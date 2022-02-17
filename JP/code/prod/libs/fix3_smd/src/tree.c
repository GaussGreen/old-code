/****************************************************************************/
/*      Building of the tree: calculation of forward rates, spot volatility */
/*      interest rate drift at each time step and tree limits.              */
/****************************************************************************/
/*      TREE.c                                                              */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fix123head.h"

/* Width of the outer ellipse in std devs. Orignal code had outerStdDevs = 1
*  this caused centering of DFs to fail on long dated trades.
*  Modified JBL 30/1/03. */
#define OUTER_STD_DEVS 0


/*****  Fix3_Forward_Curve  ******************************************************/
/*
*       Interpolate zero and calculate forward at each node point i.e.
*       every Length[i] interval. When the value date is not the same as the 
*       zero curve value date we need to interpolate or even extrapolate 
*       backward.
*/
int     Fix3_Forward_Curve ( double  *ZeroCoupon, /* (O) Zero coupon price       */
                        double  *ZeroRate,   /* (O) Zero coupon rate        */
                        double  *FwdRate,    /* (O) Forward rate            */
                        long    ZValueDate,  /* (I) Zero curve value date   */
                        int     NbZero,      /* (I) Number of zeros         */
                        double  *Zero,       /* (I) Zero rates (SA ACT/365) */
                        long    *ZeroDate,   /* (I) Zero maturity dates     */
                        long    *TPDate,     /* (I) Date of each time point */
                        int     NbTP)        /* (I) Nb of time points       */
{
    int     i;
    double  inv_mat;            /* 1/Maturity of zero coupon bond   */
    int     status = FAILURE;   /* Error status = FAILURE initially */

    double  ZPrice_i,           /* ZeroPrice from ZValueDate to TPDate[i]   */
            ZPrice_i1;          /* ZeroPrice from ZValueDate to TPDate[i+1] */

    if (NbTP <= 0L) goto RETURN;

    /* compute the 1-period forwards */
    ZPrice_i = ZeroPrice(TPDate[0],ZValueDate, NbZero, ZeroDate, Zero);
    if (ZPrice_i < TINY) goto RETURN;

    for (i = 0; i <= NbTP; i++)
    {
        ZPrice_i1 = ZeroPrice(TPDate[i+1],ZValueDate, NbZero, ZeroDate, Zero);
        if (ZPrice_i1 < TINY) goto RETURN;

        FwdRate[i] = ZPrice_i/ZPrice_i1 - 1.0;

        ZPrice_i = ZPrice_i1;
    }

    /* compute ZeroCoupon price and ZeroRate that are based from TPDate[0] */
    ZeroCoupon[0] = 1.0;
    for (i = 0; i <= NbTP; i++)
    {
        ZeroCoupon[i+1] = ZeroCoupon[i] / (1.0 + FwdRate[i]);

        inv_mat  = 365.0 / Daysact (TPDate[0], TPDate[i+1]);
        ZeroRate[i+1] = pow(ZeroCoupon[i+1], -inv_mat) - 1.0;
    }
    ZeroRate[0] = ZeroRate[1];

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Fix3_Forward_Curve: failed.");
    }
    return (status);

}  /* Fix3_Forward_Curve */



/*****  Fix3_Tree_Limits  ********************************************************/
/*
*       Calculates tree limits: we calculate the volatility of each process
*       and cut the tree according to a parabola in each dimension.
*       Also allocate on the go memory for the tree limits.
*/
int     Fix3_Tree_Limits (int     *Width,         /* (O) Ellipsoid width         */
                        int     *HalfWidth,          /* (O) Ellipsoid width         */
                        int     *Top1,               /* (O) Upper limits of 1D tree */
                        int     *Bottom1,            /* (O) Lower limits            */
                        int     **Top2,              /* (O) Upper limits of 2D tree */
                        int     **Bottom2,           /* (O) Lower limits            */
                        int     ***Top3,             /* (O) Upper limits of 3D tree */
                        int     ***Bottom3,          /* (O) Lower limits            */
                        int     NbFactor,            /* (I) Number of factors       */
                        int     NbSigmaMax,          /* (I) Max nb of std devs      */
                        int     NbTP,                /* (I) Nb of time points       */
                        double  Afac,                /* (I) Skew level              */
                        double  Bfac,                /* (I) Smile level             */
                        double  Cfac,                /* (I) Intensity               */ 
                        double  Dfac,                /* (I) Decay                   */
                        long    *TPDate,             /* (I) Tree dates              */
                        double  *Beta,               /* (I) Mean reversions         */
                        double  **Aweight,           /* (I) Orthogonal weights      */
                        double  *Length,             /* (I) Length of time steps    */
                        double  *LengthJ)            /* (I) Time step for jump size */
{

   double  Var1,Var3;   /* Variance of factors               */
    double  Cov2,Cov3;   /* Covariance of factors             */
    double  Vol1,Vol3;   /* Correlation of factors            */
    double  Rho2,Rho3;   /* Volatility of factors             */
                                                                     
    double  Jump[6];            /* Jump sizes                        */
    double  Delta, x1;          /* Intermediate results              */
    double  a1 ;                /* Ellipsoid axis                    */
    double  du;                 /* Jump coefficient                  */
    double  gamma = 0;
    double  AlphaUp, AlphaDown;
    double  AdjFacUp, AdjFacDown;
    double  X1Min, X1Max;
    
    double  Var2Up, Var2Down;
    double  Vol2Up, Vol2Down;
    double  Var12Up, Var12Down;
    double  Rho1Up, Rho1Down;
    double  Cov1Up, Cov1Down;


                                                                     
    int     Center;             /* Center of the ellipse             */
    int     h;                  /* Local value of half ellipse width */
    int     h2    ;             /* Minimum width of ellipsoid        */

    int     t, i, j;            /* Time and node indices             */
    int     status = FAILURE;   /* Error status                      */ 


    /* Reset variables to 0 */
    Var1 = Var3        = 0.;
    Vol1 = Vol3        = 0.;
    Cov2 = Cov3        = 0.;
    Rho2 = Rho3        = 0.;

    Var2Up = Var2Down  = 0.;
    Vol2Up = Vol2Down  = 0.;
    Var12Up= Var12Down = 0.;
    Rho1Up = Rho1Down  = 0.;
    Cov1Up = Cov1Down  = 0.;


    /* 
     *  Start the tree with 5 nodes in each dimension.
     */

    for (i = 0; i < 3; i++)
    {
        Width[i]     = 5;
        HalfWidth[i] = 2;
    }

    Bottom1[0] = -2;
    Top1[0]    =  2;

    if (NbFactor >= 2)
    {  
        Bottom2[0] = (int *) DR_Array (INT, -2, 2);
        Top2[0]    = (int *) DR_Array (INT, -2, 2);

        if (  (Bottom2[0] == NULL)
           || (Top2[0]    == NULL))
        {
            goto RETURN;
        }

        if (NbFactor >= 3)
        {    
            Bottom3[0] = (int **) DR_Array (INT_PTR, -2, 2);
            Top3[0]    = (int **) DR_Array (INT_PTR, -2, 2);

            if (  (Bottom3[0] == NULL)
               || (Top3[0]    == NULL))
            {
                goto RETURN;
            }

        }

        for (i = -2; i <= 2; i++)
        {
            Bottom2[0][i] = -2;
            Top2[0][i]    =  2;
                
            if (NbFactor >= 3)
            {    
                Bottom3[0][i] = (int *) DR_Array (INT, -2, 2);
                Top3[0][i]    = (int *) DR_Array (INT, -2, 2);

                if (  (Bottom3[0][i] == NULL)
                   || (Top3[0][i]    == NULL))
                {
                    goto RETURN;
                }

                for (j = -2; j <= 2; j++)
                {
                    Bottom3[0][i][j] = -2;
                    Top3[0][i][j]    =  2;
                }
            }  /* if */
        }  /* for i */
    }  /* if */



    /* Initialize */
    AlphaUp     = 1.0;
    AlphaDown   = 1.0;

    if (IS_EQUAL(Beta[0], Beta[1]))
    {
        gamma = Beta[1] / 0.000001;
    }
    else
    {
        gamma = Beta[1] / (Beta[1] - Beta[0]);
    }
    
    AdjFacUp   = AlphaUp   * gamma;
    AdjFacDown = AlphaDown * gamma;

    for (t = 0; t <= NbTP; t++)                  
    {
        /* 
         *  Jump coefficients.
         */

        du = sqrt (JUMPCOEFF * LengthJ[t-1]);

        Jump[0] = Aweight[0][t-1] * du;
        Jump[1] = Aweight[1][t-1] * du;
        Jump[2] = Aweight[2][t-1] * du;
        Jump[3] = Aweight[3][t-1] * du;
        Jump[4] = Aweight[4][t-1] * du;
        Jump[5] = Aweight[5][t-1] * du;

        
        /*
         *  Limits of 1D tree.
         */

        if (t > 0)
        {
            /* fabs () is needed as Jumps can be negative */
            /* Also make sure there is at least 2 nodes   */
            h = (int) ceil (NbSigmaMax * Vol1 / fabs (Jump[0]));
            h = MAX (h,  2);

            /* No shift of axis in first dimension */
            Center = 0;

            Bottom1[t] = Center - h;
            Top1[t]    = Center + h;
            
            /* Record the maximum width over time */
            HalfWidth[0] = MAX (HalfWidth[0], h);
            Width[0]     = 2 * HalfWidth[0] + 1;
        }


        /* Approximate limits for the mean-drift */
        if ((NbFactor >= 2) && (t > 0))
        {
            X1Min        = -NbSigmaMax * Vol1;
            AlphaDown    = Smd_MeanDrift(X1Min, Afac,Bfac, Cfac);
            AlphaDown   /= X1Min;

            X1Max        = NbSigmaMax * Vol1;
            AlphaUp      = Smd_MeanDrift(X1Max, Afac,Bfac, Cfac);
            AlphaUp     /= X1Max;           
            
            AdjFacUp   = AlphaUp   * gamma;
            AdjFacDown = AlphaDown * gamma;
        }


        /* 
         *  Limits of 2D tree. They depend on first dimension.
         */

        if ((NbFactor >= 2) && (t > 0))
        {  
            /*
             *  Allocate memory for the 2D limits and the 3D matrices of limits.
             */

            Bottom2[t] = (int *) DR_Array (INT, Bottom1[t], Top1[t]);
            Top2[t]    = (int *) DR_Array (INT, Bottom1[t], Top1[t]);

            if (  (Bottom2[t] == NULL)
               || (Top2[t]    == NULL))
            {
                goto RETURN;
            }

            if (NbFactor >= 3)
            {    
                Bottom3[t] = (int **) DR_Array (INT_PTR, Bottom1[t], Top1[t]);
                Top3[t]    = (int **) DR_Array (INT_PTR, Bottom1[t], Top1[t]);

                if (  (Bottom3[t] == NULL)
                   || (Top3[t]    == NULL))
                {
                    goto RETURN;
                }
            }

            
            /*  Minimum half size of the ellipse: slope + 1 node. This is to
            *   avoid case(1) where there is no square of 9 nodes to branch to.
            *   (case(2) is correct):
            *                                       x
            *                                   x   x
            *               x               O   O   O
            *           x   x               O   O   O
            *       x   x   x   =======>    O   O   O
            *       x   x                   x   x
            *       x                       x
            *
            *           (1)                     (2)
            */
       
            for (i = Bottom1[t]; i <= Top1[t]; i++)
            {


                /* Slope of 2D ellipse axis */
                if (i <= 0)
                    a1 = (Rho1Up   * Vol2Up   / Vol1 * Jump[0] - Jump[1]) / fabs (Jump[2]);
                else
                    a1 = (Rho1Down * Vol2Down / Vol1 * Jump[0] - Jump[1]) / fabs (Jump[2]);     

                h2 = (int) ceil (fabs(a1)) + 2;     


                x1 = i * Jump[0] / Vol1;

                if (i <= 0)
                    Delta = MAX (TINY, (1. - Rho1Up  * Rho1Up)  * (NbSigmaMax*NbSigmaMax - x1*x1) * SQUARE (Vol2Up  / Jump[2]));
                else
                    Delta = MAX (TINY, (1. - Rho1Down* Rho1Down)* (NbSigmaMax*NbSigmaMax - x1*x1) * SQUARE (Vol2Down / Jump[2]));   


                Delta = sqrt (Delta);                   

 

                h = (int) ceil (Delta);
                h = MAX (h, h2);

                Center = NEAR_INT (i * a1);

                Bottom2[t][i] = Center - h;
                Top2[t][i]    = Center + h;
                
                HalfWidth[1] = MAX (HalfWidth[1], h);
                Width[1]     = 2 * HalfWidth[1] + 1;


                /* 
                 *  Limits of 3D tree. They depend on first two dimensions.
                */

                if (NbFactor >= 3)
                {
                }  /* if */
            }  /* for i */
        }  /* if */


        /*
         *  Compute correlation matrix at next time step.
         */

        Var1 *= exp (- 2. * Beta[0] * Length[t]);
        Var1 += Aweight[0][t] * Aweight[0][t] * Length[t] * Fix3_ExpDecay (2. * Beta[0], Length[t]);

        Vol1 = sqrt (Var1);

        if (NbFactor >= 2)
        {  
            Var2Up *= exp (- 2. * Beta[1] * Length[t]);
            Var2Up += SQUARE(Aweight[1][t] - AdjFacUp * Aweight[0][t]) * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);
            Var2Up += Aweight[2][t] * Aweight[2][t]                    * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);

            Var2Down *= exp (- 2. * Beta[1] * Length[t]);
            Var2Down += SQUARE(Aweight[1][t] - AdjFacDown * Aweight[0][t]) * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);
            Var2Down += Aweight[2][t] * Aweight[2][t]                      * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);
            

            Var12Up *= exp (- (Beta[0] + Beta[1]) * Length[t]);
            Var12Up += Aweight[0][t] * (Aweight[1][t] - AdjFacUp * Aweight[0][t]) * Length[t] * Fix3_ExpDecay (Beta[0] + Beta[1], Length[t]);

            Var12Down *= exp (- (Beta[0] + Beta[1]) * Length[t]);
            Var12Down += Aweight[0][t] * (Aweight[1][t] - AdjFacDown * Aweight[0][t]) * Length[t] * Fix3_ExpDecay (Beta[0] + Beta[1], Length[t]);

            Cov1Up   = Var12Up + AdjFacUp * Var1;
            Vol2Up   = sqrt (Var2Up + AdjFacUp * AdjFacUp * Var1 + 2.0 * AdjFacUp * Var12Up);

            Cov1Down = Var12Down + AdjFacDown * Var1;
            Vol2Down = sqrt (Var2Down + AdjFacDown * AdjFacDown * Var1 + 2.0 * AdjFacDown * Var12Down);
            

            Rho1Up   = Cov1Up   / Vol1 / Vol2Up;
            Rho1Down = Cov1Down / Vol1 / Vol2Down;      
        }
    

        if (NbFactor >= 3)
        {
            ;
        }

    }  /* for t */


    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Fix3_Tree_Limits: could not allocate memory for the tree!");
    }

    return (status);

}  /* Fix3_Tree_Limits */



/*****  Fix3_Build_Tree  *********************************************************/
/*
*       Main routine for the construction of the tree.
*/
int     Fix3_Build_Tree (T_CURVE         *t_curve,       /* (I) Zero curve      */
                    MKTVOL_DATA     *mktvol_data,   /* (I) Volatility data */
                    FIX3_TREE_DATA       *tree_data)     /* (O) Tree data       */
{

    int     i;                       /*                                  */
    int     status = FAILURE;        /* Error status = FAILURE initially */
                

    /* Calculation of index forward rate at each node in the tree */
    for (i = 0; i < 3; i++)
    {
        if (Fix3_Forward_Curve ( tree_data->ZeroCoupon[i],
                            tree_data->ZeroRate[i],
                            tree_data->FwdRate[i],
                            (t_curve[i]).ValueDate,
                            (t_curve[i]).NbZero,
                            (t_curve[i]).Zero,
                            (t_curve[i]).ZeroDate,
                            tree_data->TPDate,
                            tree_data->NbTP) == FAILURE)
        {
            goto RETURN;
        }
    }  /* for i */

    /* Bootstrap swaption volatility to get spot volatility curve   */
    /* Use diffused t_curve ! Perform filtering if necessary eg CET */
    if ( (mktvol_data->FilterSpotVolFlag == TRUE) &&
         (mktvol_data->SkipFlag == FALSE) )
    {
        if (Smd_SpotVol (   mktvol_data->Aweight,
                        mktvol_data->BaseDate,
                        mktvol_data->NbVol,
                        mktvol_data->VolDate,
                        mktvol_data->Vol,
                        mktvol_data->VolUsed,
                        mktvol_data->Freq,
                        mktvol_data->DCC,
                        mktvol_data->SwapSt,
                        mktvol_data->SwapMat,
                        mktvol_data->QLeft,
                        mktvol_data->QRight,
                        mktvol_data->FwdShift,
                        mktvol_data->Bbq,
                        mktvol_data->VolNorm,
                        mktvol_data->VolLogn,
                        tree_data->NbFactor,
                        mktvol_data->Alpha,
                        mktvol_data->Beta,
                        mktvol_data->Rho,
                        mktvol_data->Dfac,
                        mktvol_data->SkipFlag,
                        mktvol_data->CalibFlag,
                        (t_curve[tree_data->CvDiff]).NbZero,
                        (t_curve[tree_data->CvDiff]).Zero,
                        (t_curve[tree_data->CvDiff]).ZeroDate,
                        (t_curve[tree_data->CvDiff]).ValueDate) == FAILURE)
        {
            goto RETURN;
        }


    }
    else  /* no filtering */
    {
        if (Smd_SpotVol (   mktvol_data->Aweight,
                        mktvol_data->BaseDate,
                        mktvol_data->NbVol,
                        mktvol_data->VolDate,
                        mktvol_data->Vol,
                        mktvol_data->VolUsed,
                        mktvol_data->Freq,
                        mktvol_data->DCC,
                        mktvol_data->SwapSt,
                        mktvol_data->SwapMat,
                        mktvol_data->QLeft,
                        mktvol_data->QRight,
                        mktvol_data->FwdShift,
                        mktvol_data->Bbq,
                        mktvol_data->VolNorm,
                        mktvol_data->VolLogn,
                        tree_data->NbFactor,
                        mktvol_data->Alpha,
                        mktvol_data->Beta,
                        mktvol_data->Rho,
                        mktvol_data->Dfac,
                        mktvol_data->SkipFlag,
                        mktvol_data->CalibFlag,
                        (t_curve[tree_data->CvDiff]).NbZero,
                        (t_curve[tree_data->CvDiff]).Zero,
                        (t_curve[tree_data->CvDiff]).ZeroDate,
                        (t_curve[tree_data->CvDiff]).ValueDate) == FAILURE)
        {
            goto RETURN;
        }
    }

    /* Interpolate spot volatility curves */
    if (Smd_Interp_SpotVol (tree_data->Aweight,
                            mktvol_data->BaseDate,
                            mktvol_data->NbVol,
                            mktvol_data->VolDate,
                            mktvol_data->Aweight,
                            tree_data->Length,
                            tree_data->NbTP) == FAILURE)
    {
        goto RETURN;
    }
    


    /* Calculation of inner ellipse limits */
    if (Fix3_Tree_Limits (   tree_data->Width,
                        tree_data->HalfWidth,
                        tree_data->Top1,
                        tree_data->Bottom1,
                        tree_data->Top2,
                        tree_data->Bottom2,
                        tree_data->Top3,
                        tree_data->Bottom3,
                        tree_data->NbFactor,
                        (tree_data->NbSigmaMax+1),
                        tree_data->NbTP,
                        mktvol_data->Afac,
                        mktvol_data->Bfac,
                        mktvol_data->Cfac,
                        mktvol_data->Dfac,
                        tree_data->TPDate,
                        mktvol_data->Beta,
                        tree_data->Aweight,
                        tree_data->Length,
                        tree_data->LengthJ) == FAILURE)
    {
        goto RETURN;
    }

    /* Calculation of outer ellipse limits */
    /* Outer ellipse is one standard deviation away from inner ellipse */
    if (Fix3_Tree_Limits (   tree_data->Width,
                        tree_data->HalfWidth,
                        tree_data->OutTop1,
                        tree_data->OutBottom1,
                        tree_data->OutTop2,
                        tree_data->OutBottom2,
                        tree_data->OutTop3,
                        tree_data->OutBottom3,
                        tree_data->NbFactor,
                        (tree_data->NbSigmaMax+1+OUTER_STD_DEVS),
                        tree_data->NbTP,
                        mktvol_data->Afac,
                        mktvol_data->Bfac,
                        mktvol_data->Cfac,
                        mktvol_data->Dfac,
                        tree_data->TPDate,
                        mktvol_data->Beta,
                        tree_data->Aweight,
                        tree_data->Length,
                        tree_data->LengthJ) == FAILURE)
    {
        goto RETURN;
    }

    /* Calculation of the 1 period rate drift in the tree */
    /* Again use index zero and forward curves */
    if (tree_data->NbFactor == 1)                                               
    {
        if (Fix3_Drift_1D ( mktvol_data,
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Fix3_Drift_2D ( mktvol_data,         
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Fix3_Drift_3D ( mktvol_data,         
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Build_Tree */
