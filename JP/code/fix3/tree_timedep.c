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



/*****  Fix3_Tree_Limits  ********************************************************/
/*
*       Calculates tree limits: we calculate the volatility of each process
*       and cut the tree according to a parabola in each dimension.
*       Also allocate on the go memory for the tree limits.
*/
int     Fix3_Tree_Limits_TimeDep (   int     *Width,     /* (O) Ellipsoid width         */
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

    double **Beta;              /* Temporary mean reversion array    */
    int     NbFactor;           /* (I) Number of factors             */
    int     NbTP;               /* (I) Nb of time points             */
    double  **Aweight;          /* (I) Orthogonal weights            */
    double  *Length;            /* (I) Length of time steps          */
    double  *LengthJ;           /* (I) Time step for jump size       */
    /*initialize mean reversion array to mkvol_data one */
    NbFactor = tree_data->NbFactor;
    NbTP     = tree_data->NbTP;
    Aweight  = tree_data->Aweight;
    Beta     = tree_data->BetaTD;
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

        Var1 *= exp (- 2. * Beta[0][t] * Length[t]);
        Var1 += Aweight[0][t] * Aweight[0][t] * Length[t] * Fix3_ExpDecay (2. * Beta[0][t], Length[t]);

        Vol1 = sqrt (Var1);

        if (NbFactor >= 2)
        {  
            Var2 *= exp (- 2. * Beta[1][t] * Length[t]);
            Var2 += Aweight[1][t] * Aweight[1][t] * Length[t] * Fix3_ExpDecay (2. * Beta[1][t], Length[t]);
            Var2 += Aweight[2][t] * Aweight[2][t] * Length[t] * Fix3_ExpDecay (2. * Beta[1][t], Length[t]);

            Vol2 = sqrt (Var2);

            Cov1 *= exp (- (Beta[0][t] + Beta[1][t]) * Length[t]);
            Cov1 += Aweight[0][t] * Aweight[1][t] * Length[t] * Fix3_ExpDecay (Beta[0][t] + Beta[1][t], Length[t]);

            Rho1 = Cov1 / Vol1 / Vol2;
        }

        if (NbFactor >= 3)
        {    
            Var3 *= exp (- 2. * Beta[2][t] * Length[t]);
            Var3 += Aweight[3][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2][t], Length[t]);
            Var3 += Aweight[4][t] * Aweight[4][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2][t], Length[t]);
            Var3 += Aweight[5][t] * Aweight[5][t] * Length[t] * Fix3_ExpDecay (2. * Beta[2][t], Length[t]);

            Vol3 = sqrt (Var3);

            Cov2 *= exp (- (Beta[0][t] + Beta[2][t]) * Length[t]);
            Cov2 += Aweight[0][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (Beta[0][t] + Beta[2][t], Length[t]);

            Rho2   = Cov2 / Vol1 / Vol3;

            Cov3 *= exp (- (Beta[1][t] + Beta[2][t]) * Length[t]);
            Cov3 += Aweight[1][t] * Aweight[3][t] * Length[t] * Fix3_ExpDecay (Beta[1][t] + Beta[2][t], Length[t]);
            Cov3 += Aweight[2][t] * Aweight[4][t] * Length[t] * Fix3_ExpDecay (Beta[1][t] + Beta[2][t], Length[t]);

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

}  /* Fix3_Tree_Limits_Timedep */



