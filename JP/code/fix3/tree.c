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



/*****  Fix3_Forward_Curve  ******************************************************/
/*
*       Interpolate zero and calculate forward at each node point i.e.
*       every Length[i] interval. When the value date is not the same as the 
*       zero curve value date we need to interpolate or even extrapolate 
*       backward.
*/
int    	Fix3_Forward_Curve (
        double*         ZeroCoupon, /* (O) Zero coupon price       */
        double*         ZeroRate,   /* (O) Zero coupon rate        */
        double*         FwdRate,    /* (O) Forward rate            */
        double*         TermZero,   /* (O) Zero to terminal point  */
        T_CURVE const*  crv,        /* (I) Zero curve              */
        long const*     TPDate,     /* (I) Date of each time point */
        int             NbTP)       /* (I) Nb of time points       */
{
    int     i;
    double  inv_mat;            /* 1/Maturity of zero coupon bond   */
    int     status = FAILURE;   /* Error status = FAILURE initially */

    double  ZPrice_i,           /* ZeroPrice from ZValueDate to TPDate[i]   */
            ZPrice_i1;          /* ZeroPrice from ZValueDate to TPDate[i+1] */

    if (NbTP <= 0L) goto RETURN;

    /* compute the 1-period forwards */
    ZPrice_i = GetZeroPrice(TPDate[0], crv);
    if (ZPrice_i < TINY) goto RETURN;

    for (i = 0; i <= NbTP; i++)
    {
        ZPrice_i1 = GetZeroPrice(TPDate[i+1], crv);
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

    for (i = 0; i <=NbTP+1; i++)
    {
        TermZero[i] = ZeroCoupon[NbTP] / ZeroCoupon[i];
    }

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("Fix3_Forward_Curve: failed.");
    }
    return (status);

}  /* Fix3_Forward_Curve */



/*****  Fix3_Tree_Limits_Classic  ************************************************/
/*
*       Calculates tree limits: we calculate the volatility of each process
*       and cut the tree according to a parabola in each dimension.
*       Also allocate on the go memory for the tree limits.
*/
int     Fix3_Tree_Limits_Classic (int     *Width,/* (O) Ellipsoid width    */
                        int     *HalfWidth, /* (O) Ellipsoid width         */
                        int     *Top1,      /* (O) Upper limits of 1D tree */
                        int     *Bottom1,   /* (O) Lower limits            */
                        int     **Top2,     /* (O) Upper limits of 2D tree */
                        int     **Bottom2,  /* (O) Lower limits            */
                        int     ***Top3,    /* (O) Upper limits of 3D tree */
                        int     ***Bottom3, /* (O) Lower limits            */
                        int     NbSigmaMax, /* (I) Max nb of std devs      */
                        FIX3_TREE_DATA *tree_data,/* (I) Tree data         */
                        MKTVOL_DATA *mktvol_data) /* (I) Volatility data   */
{

    double  Var1, Var2, Var3;   /* Variance of factors               */
    double  Cov1, Cov2, Cov3;   /* Covariance of factors             */
    double  Vol1, Vol2, Vol3;   /* Correlation of factors            */
    double  Rho1, Rho2, Rho3;   /* Volatility of factors             */
                                                                     
    double  Jump[6];            /* Jump sizes                        */
    double  Delta, x1, x2, x3;  /* Intermediate results              */
    double  a1, a2, a3;         /* Ellipsoid axis                    */
    double  du;                 /* Jump coefficient                  */
                                                                     
    int     Center;             /* Center of the ellipse             */
    int     h;                  /* Local value of half ellipse width */
    int     h2, h3;             /* Minimum width of ellipsoid        */

    int     t, i, j;            /* Time and node indices             */
    int     status = FAILURE;   /* Error status                      */
   
    
    double *Beta;               /* Temporary mean reversion array    */
    int     NbFactor;           /* (I) Number of factors             */
    int     NbTP;               /* (I) Nb of time points             */
    double  **Aweight;          /* (I) Orthogonal weights            */
    double  *Length;            /* (I) Length of time steps          */
    double  *LengthJ;           /* (I) Time step for jump size       */

    /*initialize mean reversion array to mkvol_data one */
    Beta     = mktvol_data->Beta;
    NbFactor = mktvol_data->NbFactor;
    NbTP     = tree_data->NbTP;
    Aweight  = tree_data->Aweight;
    Length   = tree_data->Length;
    LengthJ  = tree_data->LengthJ; 

    Var1 = Var2 = Var3 = 0.;
    Vol1 = Vol2 = Vol3 = 0.;
    Cov1 = Cov2 = Cov3 = 0.;
    Rho1 = Rho2 = Rho3 = 0.;


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


            /* Slope of 2D ellipse axis */
            a1 = (Rho1 * Vol2 / Vol1 * Jump[0] - Jump[1]) / fabs (Jump[2]);

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

            h2 = (int) ceil (fabs(a1)) + 2;

            for (i = Bottom1[t]; i <= Top1[t]; i++)
            {
                x1 = i * Jump[0] / Vol1;

                Delta = MAX (TINY, (1. - Rho1*Rho1) * (NbSigmaMax*NbSigmaMax - x1*x1) * SQUARE (Vol2 / Jump[2]));
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
                    Bottom3[t][i] = (int *) DR_Array (INT, Bottom2[t][i], Top2[t][i]);
                    Top3[t][i]    = (int *) DR_Array (INT, Bottom2[t][i], Top2[t][i]);

                    if (  (Bottom3[t][i] == NULL)
                       || (Top3[t][i]    == NULL))
                    {
                        goto RETURN;
                    }

                    /* Slope of 3D ellipsoid axis */
                    a2 = (Vol3 * ((Rho2-Rho1*Rho3) / Vol1 * Jump[0]
                                + (Rho3-Rho1*Rho2) / Vol2 * Jump[1]) / (1.-Rho1*Rho1) - Jump[3]) / fabs(Jump[5]);
                    a3 = (Vol3 *  (Rho3-Rho1*Rho2) / Vol2 * Jump[2]  / (1.-Rho1*Rho1) - Jump[4]) / fabs(Jump[5]);

                    /*  Minimum half size of the ellipsoid */
                    h3 = (int) ceil (fabs(a2) + fabs(a3)) + 2;


                    for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
                    {
                        x2 = (i * Jump[1] + j * Jump[2]) / Vol2;
                        x3 = ((Rho2 - Rho1*Rho3) * x1 + (Rho3 - Rho1*Rho2) * x2) / (1. - Rho1*Rho1);

                        Delta =  NbSigmaMax*NbSigmaMax - (x1*x1 + x2*x2 - 2. * Rho1*x1*x2) / (1. - Rho1*Rho1);
                        Delta *= (1. - Rho1*Rho1 - Rho2*Rho2 - Rho3*Rho3 + 2.*Rho1*Rho2*Rho3) / (1. - Rho1*Rho1);
                        Delta = MAX (TINY, Delta * SQUARE (Vol3 / Jump[5]));
                        Delta = sqrt (Delta);

                        h = (int) ceil (Delta);
                        h = MAX (h, h3);

                        Center = NEAR_INT (a2 * i + a3 * j);

                        Bottom3[t][i][j] = Center - h;
                        Top3[t][i][j]    = Center + h;
                
                        HalfWidth[2] = MAX (HalfWidth[2], h);
                        Width[2]     = 2 * HalfWidth[2] + 1;
                    }
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
            Var2 *= exp (- 2. * Beta[1] * Length[t]);
            Var2 += Aweight[1][t] * Aweight[1][t] * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);
            Var2 += Aweight[2][t] * Aweight[2][t] * Length[t] * Fix3_ExpDecay (2. * Beta[1], Length[t]);

            Vol2 = sqrt (Var2);

            Cov1 *= exp (- (Beta[0] + Beta[1]) * Length[t]);
            Cov1 += Aweight[0][t] * Aweight[1][t] * Length[t] * Fix3_ExpDecay (Beta[0] + Beta[1], Length[t]);

            Rho1 = Cov1 / Vol1 / Vol2;
        }

        if (NbFactor >= 3)
        {    
            Var3 *= exp (- 2. * Beta[2] * Length[t]);
            Var3 += Aweight[3][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2], Length[t]);
            Var3 += Aweight[4][t] * Aweight[4][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2], Length[t]);
            Var3 += Aweight[5][t] * Aweight[5][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2], Length[t]);

            Vol3 = sqrt (Var3);

            Cov2 *= exp (- (Beta[0] + Beta[2]) * Length[t]);
            Cov2 += Aweight[0][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (Beta[0] + Beta[2], Length[t]);

            Rho2   = Cov2 / Vol1 / Vol3;

            Cov3 *= exp (- (Beta[1] + Beta[2]) * Length[t]);
            Cov3 += Aweight[1][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (Beta[1] + Beta[2], Length[t]);
            Cov3 += Aweight[2][t] * Aweight[4][t] * Length[t] * Fix3_ExpDecay (Beta[1] + Beta[2], Length[t]);

            Rho3 = Cov3 / Vol2 / Vol3;
        }
        
    }  /* for t */


    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Fix3_Tree_Limits: could not allocate memory for the tree!");
    }

    return (status);

}  /* Fix3_Tree_Limits_Classic */



/*****  Fix3_Build_Tree  *********************************************************/
/*
*       Main routine for the construction of the tree.
*/
int     Fix3_Build_Tree (
            T_CURVE const*  t_curve,       /* (I) Zero curve        */
            MKTVOL_DATA*    mktvol_data,   /* (I/O) Volatility data */
            FIX3_TREE_DATA* tree_data)     /* (O) Tree data         */
{

    int     i;                       /*                                  */
    int     status = FAILURE;        /* Error status = FAILURE initially */
                

    /* Calculation of index forward rate at each node in the tree */
    for (i = 0; i < 3; i++)
    {
        if (Fix3_Forward_Curve (tree_data->ZeroCoupon[i],
                                tree_data->ZeroRate[i],
                                tree_data->FwdRate[i],
                                tree_data->TermZero[i],
                                &t_curve[i],
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
        if (Fix3_Filtered_SpotVol (
                        mktvol_data,
                        &t_curve[tree_data->CvDiff]) != SUCCESS)
        {
            goto RETURN;
        }


    }
    else  /* no filtering */
    {
        if (Fix3_SpotVol (   
                        mktvol_data,
                        &t_curve[tree_data->CvDiff]) != SUCCESS)
        {
            goto RETURN;
        }
    }

    /* Interpolate spot volatility curves */
    if (Fix3_Interp_SpotVol (
                           tree_data,
                           mktvol_data) == FAILURE)
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
                        (tree_data->NbSigmaMax+1),
                        tree_data,
                        mktvol_data) == FAILURE)
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
                        (tree_data->NbSigmaMax+1+OUTER_STD_DEVS),
                        tree_data,
                        mktvol_data) == FAILURE)
    {
        goto RETURN;
    }

    /* Calculation of the 1 period rate drift in the tree */
    /* Again use index zero and forward curves */
    if (tree_data->NbFactor == 1)                                               
    {
        if (Fix3_Drift_1D (	mktvol_data,
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Fix3_Drift_2D (	mktvol_data,         
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Fix3_Drift_3D (	mktvol_data,         
                        tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }  /* if then else */

    /* Allocate and compute numeraire if necessary */
    if (mktvol_data->IsNmrModel)
    {
        /* Allocate memory for numeraire and mapped yield */
        if (Nmr_Alloc (mktvol_data,
                       tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Calculate and store numeraire values */
        if (Nmr_Calc (t_curve,
                      mktvol_data,
                      tree_data) == FAILURE)
        {
            goto RETURN;
        }
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Build_Tree */
