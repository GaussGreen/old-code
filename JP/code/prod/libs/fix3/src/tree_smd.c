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


/*****  Fix3_Tree_Limits_Smd  ******************************************************/
/*
*       Calculates tree limits: we calculate the volatility of each process
*       and cut the tree according to a parabola in each dimension.
*       Also allocate on the go memory for the tree limits.
*/
int     Fix3_Tree_Limits_Smd (
                        int     *Width,              /* (O) Ellipsoid width         */
                        int     *HalfWidth,          /* (O) Ellipsoid width         */
                        int     *Top1,               /* (O) Upper limits of 1D tree */
                        int     *Bottom1,            /* (O) Lower limits            */
                        int     **Top2,              /* (O) Upper limits of 2D tree */
                        int     **Bottom2,           /* (O) Lower limits            */
                        int     ***Top3,             /* (O) Upper limits of 3D tree */
                        int     ***Bottom3,          /* (O) Lower limits            */
                        int     NbSigmaMax,          /* (I) Max nb of std devs      */
                        FIX3_TREE_DATA *tree_data,   /* (I) Tree data               */
                        MKTVOL_DATA *mktvol_data)    /* (I) Volatility data         */         
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


    int     NbFactor  = tree_data->NbFactor;
    int     NbTP      = tree_data->NbTP;
    double  Afac      = mktvol_data->Afac;
    double  Bfac      = mktvol_data->Bfac;
    double  Cfac      = mktvol_data->Cfac;
    double  Dfac      = mktvol_data->Dfac;
    long    *TPDate   = tree_data->TPDate;
    double  *Beta     = mktvol_data->Beta;
    double  **Aweight = tree_data->Aweight;
    double  *Length   = tree_data->Length;
    double  *LengthJ  = tree_data->LengthJ;


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

}  /* Fix3_Tree_Limits_Smd */



