/****************************************************************************/
/*      Extrapolation of the correlations for several factor IR             */
/****************************************************************************/
/*      CORRELATIONS.C                                                      */
/****************************************************************************/
/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cupslib.h"


/***** Hyb3_AssetIrCorr *******************************************************
 *
 *      Produces the implied corr structure between
 *      Asset and the factors driving the IR curve.
 * 
 *               - 
 *      choice 1 = EXPONENTIAL: corr. between the exp. factors
 *      choice 2 = TRIANGULAR:. corr. between orthogonalised. factors
 *      choice 3 = PC        :  corr. between principal components
 *
 *      In case the YC dim is 1, Rho can be input as NULL.
 *      Rk: Asset is the equivalent of the FX rate in the doc
 *
 *      This function is Flat Forward ready
 */
int Hyb3_AssetIrCorr (
    
    double  *RhoAssetIr,    /*(O) corr between Asset and IR factors         */
    double  RhoAssetSwap,   /*(I) corr between Asset and swap rate          */
    double  *Rho,           /*(I) corr between domestic exponential factors */
    int     NbFactor,       /*(I) number of factors for domestic curve      */
    double  *Alpha,         /*(I) relative size factors for domestic curve  */
    double  *Beta,          /*(I) dom mean rev. coeffs                      */
    int     choice,         /*(I) specifies factors wrt which correl is 
                                  output: 1=expon.,2=triang;3=PC            */

    /* Swap rates specifications*/
    long    SwapSt,         /* (I) Swap start                               */
    long    SwapMat,        /* (I) Swap maturity                            */
    char    Dcc,            /* (I) underlying day count convention          */
    char    Freq,           /* (I) freq                                     */
    long    Today,          /* (I) reference date                           */
    /* Zero curves specifications*/
    T_CURVE *tcurve)
    
{
    static char routine[]= "Hyb3_AssetIrCorr";
    double  B[3];       /* B factors in Christian's paper           */
    double  *A = NULL;  /* Normalised weights of orthogonal factors */
    double  *DummyRhoAssetIr = NULL; /* Temporary Correlation array */

    double  sum;
    double  length; /* time to SwapSt in years */
    int     i,j,
            status = FAILURE;

    double  /* triangular transformation used to produce final */
            /* correlations in case of choice = 1 */

            Transform[3][3];
    
    double  
            /* variables used to test the validity of the claculated correls*/
            test_RhoAssetSwap,
            test_Vol;
            

    /* Some basic checks */

    if (RhoAssetIr == NULL ||
        Alpha      == NULL ||
        Beta       == NULL ||
        tcurve      == NULL)
    {
        DR_Error ("Invalid pointer inputs to %s\n", routine);
        return (status);
    }

    if (NbFactor != 1 && NbFactor != 2 && NbFactor !=3)
    {
        DR_Error ("%s: Invalid number of yield curve factors, " 
                    "must be 1,2,or 3\n", routine);
        return (status);
    }

    if ((NbFactor > 1) && (Rho == NULL))
    {
        DR_Error ("%s: Correlation input cannot be NULL if Nbfactor > 1\n", routine);
        return (status);
    }

    if (choice != 1 && choice != 2 && choice != 3)
    {
        DR_Error ("%s: Invalid choice number, should be 1 or 2, or 3\n", routine);
        return (status);
    }
    if (choice == 3)
    {
        DR_Error ("%s: choice 3 (Principal Component) not implemented!\n", routine);
        return (status);
    }
 
    if (NbFactor ==1) /*No need to extend correlation in this case */
    {
        RhoAssetIr[0] = RhoAssetSwap; /* +, as BFactor always >0*/
        status = SUCCESS;
        return (status);
    }
   
    A = (double*) DR_Array (DOUBLE,0,NbFactor - 1);
    DummyRhoAssetIr = (double*) DR_Array (DOUBLE,0,NbFactor - 1);
    
  
    if (DummyRhoAssetIr == NULL ||  A == NULL)
    {
        DR_Error("memory alloc.failure in %s\n", routine);
        goto RETURN;
    }

    if (Hyb3_BFactor (B,
                SwapSt, 
                SwapMat,
                Dcc,
                Freq,
                NbFactor,
                Beta,
                tcurve) == FAILURE)
    {

        goto RETURN;
    }
    
    
   
    length = Daysact (Today,SwapSt) / 365.0; 
    
    for (i = 0; i < NbFactor; i++)
    {          
        B[i] *= exp(- Beta[i] * length);
    }

    /* calculate now normalised weights */
    switch (choice)
    {
        case 1:
        case 2:

            if (Triangulation (Transform,NbFactor,Rho) == FAILURE)
            {
                DR_Error("%s: mis-specified YC correlations\n", routine);
                goto RETURN;
            }
            
            A[0] = Alpha[0] * B[0];
            
            if (NbFactor > 1)
            {                
               
                A[0] += Alpha[1] * B[1] * Transform[1][0];
                A[1]  = Alpha[1] * B[1] * Transform[1][1];

            }/* end of "NbFactor > 1"*/

            if (NbFactor > 2)
            {                
                A[0] += Alpha[2] * B[2] * Transform[2][0];

                A[1] += Alpha[2] * B[2] * Transform[2][1];

                A[2]  = Alpha[2] * B[2] * Transform[2][2];

            }/* end of NbFactor>2*/
            
            break;

        case 3: 
            /* Not Implemented*/
            goto RETURN;
/*
            if (Get3A(NbFactor,
                      Alpha,
                      Beta,
                      Rho,
                      Transform) == FAILURE)
            {
                DR_Error("Get3A failed in %s\n", routine);
                goto RETURN;
            }
            for (i = 0; i < NbFactor; i++)
            {
                for (j = 0; j < NbFactor; j++)
                {
                    A[i] += Alpha[j] * B[j] * Transform[j][i];
                } 
            } 
  */

    }/*end of switch on choice*/

    /* normalizing domestic weights*/

    sum = 0.;
    for (i = 0; i < NbFactor; i++)
        sum += A[i] * A[i];

    sum = sqrt(sum);

    if (sum < TINY)
    {
        DR_Error ("Problem in %s: sum of the normalising "
                    "weights is "
                    "zero, unable to normalise\n", routine);
        goto RETURN;
    }

    for (i = 0; i < NbFactor; i++)
        A[i] /= sum;
    
    /* DummyRhoAssetIr are the corr mentioned in doc*/ 
    for (i = 0; i < NbFactor; i++)
    {
        DummyRhoAssetIr[i] = A[i] * RhoAssetSwap;
    }
    

    if (choice == 1) /*Convert back to exponential factor*/
    {   
        for (i = 0; i < NbFactor; i++)
        {
            sum = 0.;
            for (j = 0; j < NbFactor; j++)
            {
                sum += Transform[i][j] * DummyRhoAssetIr[j];

            } /* for j */

            RhoAssetIr[i] = sum;

        }/* for i */
                
        /* now test whether we recover input correlations*/
        
        /* Corr Fx versus Dom*/
        test_RhoAssetSwap = 0.;
        test_Vol = 0.;

        for (i = 0; i < NbFactor; i++)
        {
            test_Vol += Alpha[i] * Alpha[i]* B[i]* B[i];
        }

        if (NbFactor  > 1)
        {
            test_Vol += 2.0 * Alpha[0] * Alpha[1] * B[0] * B[1]* Rho[0];
        }

        if (NbFactor > 2)
        {
            test_Vol += 2.0 * Alpha[0] * Alpha[2] * B[0] * B[2] * Rho[1];
            test_Vol += 2.0 * Alpha[1] * Alpha[2] * B[1] * B[2] * Rho[2];
        }

        if (test_Vol < TINY)
        {
            DR_Error ("Problem calculating test Vol for IR (%s)\n", routine);
            goto RETURN;
        }

        for (i = 0 ; i < NbFactor; i++)
        {
            test_RhoAssetSwap += Alpha[i] * B[i] * RhoAssetIr[i];
        }
        test_RhoAssetSwap /= sqrt(test_Vol);

        if ( fabs (test_RhoAssetSwap - RhoAssetSwap) > TINY )   
            
        {
            DR_Error ("%s failed to reproduce"
                    "input correlation coefficient\n", routine);
            goto RETURN;
        }

    }/* end of choice=1*/

    else if (choice == 2 || choice == 3)
    {
        /* the correlations in this case have already been obtained */

        for (i = 0; i < NbFactor; i++)
        {
            RhoAssetIr[i] = DummyRhoAssetIr[i];

        } /* for i*/

    }/* end of choice=2 or choice=3*/

    status = SUCCESS;

RETURN:

    Free_DR_Array (A,DOUBLE,0,NbFactor - 1);
    Free_DR_Array (DummyRhoAssetIr,DOUBLE,0,NbFactor - 1);
    return (status);

} /* Hyb3_AssetIrCorr */
/***** Hyb3_IR_IR_Corr *******************************************************
 *
 *      Produces the implied correlation structure between YC driving factors
 *
 *      choice for output :
 *
 *      1 = EXPONENTIAL FACTORS
 *      2 = TRIANGULAR
 *      3 = PC (Not implemented)
 * 
 *      Note: in case dimension of yield curve is 1, DomRho and ForRho
 *            can be both input as NULL.
 *
 *      This function is Flat Forward ready
 */
int Hyb3_IR_IR_Corr (
    
    double  *RhoDomFor,     /*(O) corr.between "dom" and "for" factors     */
    double  RhoDswapFswap,  /*(I) corr. between "dom" and "for" rates      */
    double  *DomRho,        /*(I) corr. between dom exponential factors    */
    double  *ForRho,        /*(I) corr. between for exponential factors    */
    int     DomNbFactor,    /*(I) number of factors for domestic curve     */
    int     ForNbFactor,    /*(I) number of factors for foreign curve      */
    double  *DomAlpha,      /*(I) relative size factors for domestic curve */
    double  *ForAlpha,      /*(I) relative size factors for foreign curve  */
    double  *DomBeta,       /*(I) dom mean rev. coeffs                     */
    double  *ForBeta,       /*(I) for mean rev.coeffs                      */
    int     choice,         /*(I) specifies factors wrt which correl is
                                  output: 1=expon.,2=triang;3=PC           */

    /* Swap rates specifications*/ 
    long    DomSwapSt,      /* (I) dom swap start                       */
    long    DomSwapMat,     /* (I) dom swap maturity                    */
    long    ForSwapSt,      /* (I) for swap start                       */
    long    ForSwapMat,     /* (I) for swap maturity                    */
    char    DomDcc,         /* (I) dom underlying day count convention  */
    char    ForDcc,         /* (I) for day count convention             */
    char    DomFreq,        /* (I) dom freq                             */
    char    ForFreq,        /* (I) for freq                             */
    long    DomRefDate,     /* (I) dom reference date                   */
    long    ForRefDate,     /* (I) for ref date                         */

    /* Zero curves specifications*/
    T_CURVE *DomTcurve,     /* (I) domestic zero curve                  */
    T_CURVE *ForTcurve)     /* (I) foreign zero curve                   */
{
    static char routine[] = "Hyb3_IR_IR_Corr";
    double  /* B factors in Christian's paper */
            DomB[3],    
            ForB[3],    

            /* Normalised weights of orthogonal factors*/
            *DomA = NULL,       /*normalised weights of dom.orthog.factors  */
            *ForA = NULL,       /*normalised weights of for.orthog.factors  */

            /* Correlation arrays outputted by call to Full_Correl_Mtx  */
            *DummyRhoDomFor = NULL;
            
    double  sum,cum,length;

    int     i,j,k,l,
            status = FAILURE;

    double  /* triangular transformations used to produce final */
            /* correlations in case of choice = 1 */

            ForTransform[3][3],
            DomTransform[3][3];

    double
            /* variables used to test the validity of the claculated correls*/
            test_RhoDswapFswap,
            test_DVol,
            test_FVol;

    /* Some basic checks */

    if (RhoDomFor   == NULL ||
        DomAlpha    == NULL ||
        ForAlpha    == NULL ||
        DomBeta     == NULL ||
        ForBeta     == NULL ||
        DomTcurve   == NULL ||
        ForTcurve   == NULL)
    {
        DR_Error ("Invalid pointer inputs to %s\n", routine);
        return (status);
    }
    
    if (DomNbFactor != 1 && DomNbFactor != 2 && DomNbFactor !=3)
    {
        DR_Error ("%s: Invalid number of domestic yield curve factors, " 
                    "must be 1,2,or 3\n", routine);
        return (status);
    }
    
    if ((DomNbFactor > 1) && (DomRho == NULL))
    {
        DR_Error ("%s: Cannot have DomNbFactor > 1 and DomRho = NULL\n", routine);
        return(status);
    }

    if (ForNbFactor != 1 && ForNbFactor != 2 && ForNbFactor !=3)
    {
        DR_Error ("%s: Invalid number of foreign yield curve factors, " 
                    "must be 1,2,or 3\n", routine);
        return (status);
    }

    

    if ((ForNbFactor > 1) && (ForRho == NULL))
    {
        DR_Error ("%s: Cannot have ForNbFactor > 1 and ForRho = NULL \n", routine);
        return(status);
    }
     

    if (choice != 1 && choice != 2 && choice != 3)
    {
        DR_Error ("%s: Invalid choice number, should be 1 or 2, or 3\n", routine);
        return (status);
    }
    if (choice == 3)
    {
        DR_Error ("%s: choice 3 (PCA) not supported\n", routine);
        return (status);
    }

    DomA = (double*) DR_Array (DOUBLE,0,DomNbFactor - 1);
    ForA = (double*) DR_Array (DOUBLE,0,ForNbFactor - 1);
    DummyRhoDomFor  = (double*) DR_Array (DOUBLE,
                                          0,DomNbFactor * ForNbFactor - 1);

    
    if (    DummyRhoDomFor == NULL 
        ||  DomA == NULL
        ||  ForA == NULL)
    {
        DR_Error("memory alloc.failure in %s\n", routine);
        goto RETURN;
    }

    if(DomNbFactor==1)
    {
        DomB[0]=1; /* No need to compute B factor in this case*/
    }
    else
    {
        if (Hyb3_BFactor (DomB,
                    DomSwapSt, 
                    DomSwapMat,
                    DomDcc,
                    DomFreq,
                    DomNbFactor,
                    DomBeta,
                    DomTcurve) == FAILURE)
        {
            goto RETURN;
        }
    

        length = Daysact(DomRefDate, DomSwapSt) / 365.0;
        for (i = 0; i < DomNbFactor; i++)
        {
            DomB[i] *= exp(- DomBeta[i] * length);
        }
    }

    if(ForNbFactor==1)
    {
        ForB[0]=1; /* No need to compute B factor in this case*/
    }
    else
    {
        if (Hyb3_BFactor (ForB,
                    ForSwapSt,
                    ForSwapMat,
                    ForDcc,
                    ForFreq,
                    ForNbFactor,
                    ForBeta,
                    ForTcurve) == FAILURE)
        {
            goto RETURN;
        }
    

        length = Daysact(ForRefDate, ForSwapSt) / 365.0;
        for (i = 0; i < ForNbFactor; i++)
        {
            ForB[i] *= exp(- ForBeta[i] * length);
        }
    }


    /* calculate now normalised weights of orthogonal factors*/
    switch (choice)
    {
        case 1:
        case 2:

            if (Triangulation (DomTransform,DomNbFactor,DomRho) == FAILURE)
            {
                DR_Error("%s: mis-specified domestic correlations\n", routine);
                goto RETURN;
            }

            if (Triangulation(ForTransform,ForNbFactor,ForRho) == FAILURE)
            {
                DR_Error("%s: mis_specified foreign correlations\n", routine);
                goto RETURN;
            }

            DomA[0] = DomAlpha[0] * DomB[0];
            ForA[0] = ForAlpha[0] * ForB[0];

            if (DomNbFactor > 1)
            {                
                /* domestic weights */

                DomA[0] += DomAlpha[1] * DomB[1] * DomTransform[1][0];
                DomA[1]  = DomAlpha[1] * DomB[1] * DomTransform[1][1];

            }/* end of "DomNbFactor > 1"*/

            if (DomNbFactor > 2)
            {                
                /* domestic weights */

                DomA[0] += DomAlpha[2] * DomB[2] * DomTransform[2][0];

                DomA[1] += DomAlpha[2] * DomB[2] * DomTransform[2][1];

                DomA[2]  = DomAlpha[2] * DomB[2] * DomTransform[2][2];

            }/* end of DomNbFactor>2*/

            if (ForNbFactor > 1)
            {                
                /* foreign weights */

                ForA[0] += ForAlpha[1] * ForB[1] * ForTransform[1][0];
                ForA[1]  = ForAlpha[1] * ForB[1] * ForTransform[1][1];

            }/* end of ForNbFactor>1*/

            if (ForNbFactor > 2)
            {
                /* foreign weights */

                ForA[0] += ForAlpha[2] * ForB[2] * ForTransform[2][0];

                ForA[1] += ForAlpha[2] * ForB[2] * ForTransform[2][1];

                ForA[2]  = ForAlpha[2] * ForB[2] * ForTransform[2][2];

            }/* end of ForNbFactor>2*/

            break;

        case 3:
            goto RETURN; /*Not implemented*/
      

    }/*end of switch on choice*/

    /* normalizing domestic weights*/

    sum = 0.;
    for (i = 0; i < DomNbFactor; i++)
        sum += DomA[i] * DomA[i];

    sum = sqrt(sum);

    if (sum < TINY)
    {
        DR_Error ("%s: Problem in Customised_Correl_Mtx: sum of the domestic "
                    "weights is "
                    "zero, unable to normalise\n", routine);
        goto RETURN;
    }

    for (i = 0; i < DomNbFactor; i++)
        DomA[i] /= sum;

    /* normalizing foreign weights*/

    sum = 0.;
    for (i = 0; i < ForNbFactor; i++)
        sum += ForA[i] * ForA[i];

    sum = sqrt(sum);
    if (sum < TINY)
    {
        DR_Error ("%s: Problem in Customised_Correl_Mtx: sum of the foreign "
                    "weights is "
                    "zero, unable to normalise\n", routine);
        goto RETURN;
    }

    for (i = 0; i < ForNbFactor; i++)
        ForA[i] /= sum;
    

    for (i = 0; i< DomNbFactor; i++)
    {
        for (j = 0; j < ForNbFactor; j++)
        {
            DummyRhoDomFor[j + i * ForNbFactor] = 
                                        DomA[i] * ForA[j] * RhoDswapFswap;
        }/* for j */

    } /* for i */

    
    if (choice == 1)
    {   
        /* correlations between domestic and foreign factors*/

        for (i = 0; i < DomNbFactor; i++)
        {
            for (j = 0; j < ForNbFactor; j++)
            {
                sum = 0.;
                for (k = 0; k < DomNbFactor; k++)
                {
                    cum = 0.;
                    for (l = 0; l < ForNbFactor; l++)
                    {
                        cum += ForTransform[j][l] * 
                               DummyRhoDomFor[l + k * ForNbFactor];
                    } /* for l */

                    sum += DomTransform[i][k] * cum;
                } /* for k */

                RhoDomFor[j + i * ForNbFactor] = sum;
            } /* for j */

        } /* for i */

        /* now test whether we recover input correlation */
        test_DVol = 0.;

        for (i = 0; i < DomNbFactor; i++)
        {
            test_DVol += DomAlpha[i] * DomAlpha[i]* DomB[i]* DomB[i];
        }

        if (DomNbFactor  > 1)
        {
            test_DVol += 2.0 * DomAlpha[0] * DomAlpha[1] * 
                         DomB[0] * DomB[1]* DomRho[0];
        }

        if (DomNbFactor > 2)
        {
            test_DVol += 2.0 * DomAlpha[0] * DomAlpha[2] * 
                         DomB[0] * DomB[2] * DomRho[1];
            test_DVol += 2.0 * DomAlpha[1] * DomAlpha[2] * 
                         DomB[1] * DomB[2] * DomRho[2];
        }

        if (test_DVol < TINY)
        {
            DR_Error ("Problem calculating test Vol for domestic rate "
                      "(%s)\n", routine);
            goto RETURN;
        }

        test_FVol = 0;
        for (i = 0; i < ForNbFactor; i++)
        {
            test_FVol += ForAlpha[i] * ForAlpha[i] * ForB[i] * ForB[i];
        }

        if (ForNbFactor  > 1)
        {
            test_FVol += 2.0 * ForAlpha[0] * ForAlpha[1] * 
                         ForB[0] * ForB[1] * ForRho[0];
        }

        if (ForNbFactor > 2)
        {
            test_FVol += 2.0 * ForAlpha[0] * ForAlpha[2] * 
                         ForB[0] * ForB[2] * ForRho[1];
            test_FVol += 2.0 * ForAlpha[1] * ForAlpha[2] * 
                         ForB[1] * ForB[2] * ForRho[2];
        }

        if (test_FVol < TINY)
        {
            DR_Error ("Problem calculating test Vol for foreign rate "
                      "(%s)\n", routine);
            goto RETURN;
        }

        
        /* Corr For versus Dom*/
        test_RhoDswapFswap = 0.;
        for (i = 0; i < DomNbFactor; i++)
        {
            for (j = 0; j < ForNbFactor; j++)
            {                
                test_RhoDswapFswap += DomAlpha[i] * DomB[i] *
                                      ForAlpha[j] * ForB[j] * 
                                      RhoDomFor[j + i * ForNbFactor];
            }/* for j*/

        }/* for i*/

        test_RhoDswapFswap /= sqrt(test_FVol);
        test_RhoDswapFswap /= sqrt(test_DVol);

        if ((fabs (test_RhoDswapFswap - RhoDswapFswap) > TINY))
        {
            DR_Error ("%s failed to reproduce"
                    "input correlation coefficient\n", routine);
            goto RETURN;
        }

    }/* end of choice=1*/

    else if (choice == 2 || choice == 3)
    {   
        /* the correlation in this case have already been obtained */

        for (i = 0; i < DomNbFactor; i++)
        {
            for (j = 0; j < ForNbFactor; j++)
            {
                RhoDomFor[j + i * ForNbFactor] = 
                                     DummyRhoDomFor[j + i * ForNbFactor];
            } /* for j */
        }/* for i */

    }/* end of choice=2 or choice=3*/

    status = SUCCESS;

RETURN:

    Free_DR_Array (DomA,DOUBLE,0,DomNbFactor - 1);
    Free_DR_Array (ForA,DOUBLE,0,ForNbFactor - 1);
    Free_DR_Array (DummyRhoDomFor,DOUBLE,0,DomNbFactor * ForNbFactor - 1);

    return (status);

}/* Hyb3_IR_IR_Corr */

