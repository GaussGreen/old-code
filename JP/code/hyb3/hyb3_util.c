/****************************************************************************/
/*      Utility routines.                                                   */
/****************************************************************************/
/*      UTIL.c                                                              */
/****************************************************************************/


/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "cupslib.h"


/*****  Hyb3_GetIndexStep  *******************************************************/
/**                                                                           
 *      Obtains the difference between values of an index at adjacent nodes.
 *      This value is then used in the smoothing algorithm (Smooth_Step).
 *      This function hides the dimension generality from  the  caller, 
 *      but it has no means to check that the adequate amount of space has 
 *      been allocated under the void * being passed.
 */
double   Hyb3_GetIndexStep(
		      TSLICE           Index,      /**< (I) Index pointer       */
                      int              Dim,        /**< (I) Index dimension     */
                      int              i,          /**< (I) Node indices        */
                      int              j,
                      int              k,
                      int              t,          /**< (I) Current time point  */
                      HYB3_TREE_DATA const* tree_data) /**< (I) Tree data structure */
{

    double  *IndexL;            /* Local pointer */

    double  IndexStep;          /* Output index step */
    double  IndexVal;           /* Index value at mid node */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */



    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    IndexStep = ERROR;  /* To avoid division by 0 */
                
    switch (Dim)
    {
        case 1:
        {
            IndexL = (double *)Index + Hyb3_Node_Offset(1, 0, 0, t, tree_data);

            if (i > Bottom1)
                IndexStep = MAX (IndexStep, fabs (IndexL[i-1] - IndexL[i]));
            if (i < Top1)                                                           
                IndexStep = MAX (IndexStep, fabs (IndexL[i+1] - IndexL[i]));

            break;
        }                    
        case 2:
        {
            IndexL = (double *)Index + Hyb3_Node_Offset(2, i, 0, t, tree_data);

            IndexVal = IndexL[j];

            if (j > Bottom2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j-1] - IndexVal));
            }
            if (j < Top2[i])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[j+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                IndexL = (double *)Index + Hyb3_Node_Offset(2, i-1, 0, t, tree_data);

                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }
            if (i < Top1)
            {         
                IndexL = (double *)Index + Hyb3_Node_Offset(2, i+1, 0, t, tree_data);

                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {                                                  
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }

            break;
        }                    
        case 3:
        {
            IndexL = (double *)Index + Hyb3_Node_Offset(3, i, j, t, tree_data);

            IndexVal = IndexL[k];

            if (k > Bottom3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k-1] - IndexVal));
            }
            if (k < Top3[i][j])
            {
                IndexStep = MAX (IndexStep, fabs (IndexL[k+1] - IndexVal));
            }

            if (i > Bottom1)
            {
                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexL = (double *)Index + Hyb3_Node_Offset(3, i-1, j, t, tree_data);

                    if ((k>=Bottom3[i-1][j])&&(k<=Top3[i-1][j]))
                    {
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (i < Top1)    
            {
                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {
                    IndexL = (double *)Index + Hyb3_Node_Offset(3, i+1, j, t, tree_data);

                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {                                                       
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (j > Bottom2[i])
            {
                IndexL = (double *)Index + Hyb3_Node_Offset(3, i, j-1, t, tree_data);

                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
            if (j < Top2[i])
            {
                IndexL = (double *)Index + Hyb3_Node_Offset(3, i, j+1, t, tree_data);

                if ((k>=Bottom3[i][j+1])&&(k<=Top3[i][j+1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
        }                    
        default:
        {
            break;
        }                    
    }  /* switch */


    return (IndexStep);

}  /* Hyb3_GetIndexStep */


/****************************************************************************/
/**Note: assumes that YieldMat and YieldFreq are such that they define
   an integer total number of coupons.
  
 ****************************************************************************/

int  Hyb3_CashAnnuity_t(
                   TSLICE     Annuity,      /**<(O)                     */
                   TSLICE     Yield,        /**<(I)                     */
                   int        YieldMat,     /**<(I) in months           */
                   char       YieldFreq,    /**<(I) A.S,Q,or M          */
                   int        YieldDim,     /**<(I) dim on tree         */
                   int         t,           /**<(I) Current time period */
                   HYB3_TREE_DATA   *tree_data)   
{
    TSLICE  AnnL = NULL;
    TSLICE  YieldL = NULL;
    int     i,j;
    int status = FAILURE;
    int n;
    int FreqInt;
    int offset;
    double  dummy;

    if (Annuity == NULL ||
        Yield   == NULL)
    {
        DR_Error("Invalid pointer input (NULL)\n");
        goto RETURN;
    }

    FreqInt = Conv_Freq(YieldFreq);

    n = (int) (YieldMat * FreqInt / 12);

    if (YieldDim == 1)
    {
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        YieldL = (double*) Yield + offset;
        AnnL =   (double*) Annuity + offset;
        for (i = tree_data->Bottom1[t] ; i <= tree_data->Top1[t]; i++)
        {    
            /* if Annuity is not defined, set it to 0.0 */

            if (fabs(1 + YieldL[i]/FreqInt) < TINY) AnnL[i] = 0.0;
            else /* Annuity is well defined  */
            {
                if(fabs(YieldL[i]) < TINY) AnnL[i] = n/FreqInt;
                else
                {   
                    /* n is integer:so pow well defined */
                    dummy = pow(1 + YieldL[i]/FreqInt , n) ;
                    AnnL[i] = (dummy - 1.0)/(YieldL[i] * dummy);
                }
            }
        }/* for i */
    }/* if dim = 1 */
    else if (YieldDim == 2)
    {
        for (i = tree_data->Bottom1[t] ; i <= tree_data->Top1[t]; i++)
        {
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            YieldL = (double*) Yield + offset;
            AnnL =   (double*) Annuity + offset;
            for (j = tree_data->Bottom2[t][i]; j <= tree_data->Top2[t][i]; j++)
            {
                if (fabs(1 + YieldL[j]/FreqInt) < TINY) AnnL[j] = 0.0;
                else 
                {
                    if(fabs(YieldL[j]) < TINY) AnnL[j] = n/FreqInt;
                    else
                    {
                        dummy = pow(1 + (YieldL[j]/FreqInt), n) ;
                        AnnL[j] = (dummy - 1.0)/(YieldL[j] * dummy);
                    }
                }

            }/* for j */
        }/* for i*/

    }/* of 2D */
    else
    {
        DR_Error("Cash Annuity function does not handle 3d slices \n ");
        goto RETURN;
    }


    status = SUCCESS;
RETURN:
    return(status);

}



/* ---------------------------------------------------------------------*/
/* PrintTree                                                            */
/*                                                                      */
/* "streaming" of the tree parameters into a file, the best possible    */
/* solution without having C++                                          */
/* */
int PrintTree( FILE *fp, HYB3_TREE_DATA  tree_data, int PrintTreeSizeInfo )
{

    int Status = FAILURE;

    int t, i, j; 


    /* validate the output file */
    if (fp == NULL) 
    {
        DR_Error("PrintTree: file-pointer is set to NULL!"); 
        goto RETURN;
    }

    /* print all the information into the file */
    fprintf( fp, "NbTP: %d\n", tree_data.NbTP); 
    fprintf( fp, "Ppy: %d\n", tree_data.Ppy); 

    for (i = 0; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d TPDate: %ld, Length: %12.10f LengthJ: %12.10f\n"
                , i, tree_data.TPDate[i], tree_data.Length[i], tree_data.LengthJ[i]); 
    }

    /* Critical dates */    
    for (i = 0; i < NBCRITDATE; ++i)
    {
        fprintf( fp, "i: %d CritType: %c, NbZeros: %d\n", i, tree_data.CritType[i], tree_data.NbZeros[i]);
    }
    /* not sure what to do with those */
    /* CRIT_DATE   *CritDate[NBCRITDATE];*/ /* Critical dates description         */  
    /* int         *TPType[NBCRITDATE];  */ /* Critical type of current node      */


    /* Express DEV dates */ 
    fprintf( fp, "NbEDevDates: %d\n", tree_data.NbEDevDates); 
    for (i = 0; i < tree_data.NbEDevDates; ++i)
    {
        fprintf( fp, "i: %d EDevDate: %ld\n", i, tree_data.EDevDate[i]); 
    }
    /* TSLICE      *EDevStPrice;       */       /* Corresp state prices           */

    
    /* Deterministic zero curve */
    fprintf(fp,"#ZeroCoupon\n");
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d ZeroCoupon: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",
                i, tree_data.ZeroCoupon[0][0][i], tree_data.ZeroCoupon[0][1][i], tree_data.ZeroCoupon[0][2][i], 
                tree_data.ZeroCoupon[1][0][i], tree_data.ZeroCoupon[1][1][i], tree_data.ZeroCoupon[1][2][i]); 
    }
    fprintf(fp,"#ZeroRate\n");
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d ZeroRate: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",
                i, tree_data.ZeroRate[0][0][i], tree_data.ZeroRate[0][1][i], tree_data.ZeroRate[0][2][i], 
                tree_data.ZeroRate[1][0][i], tree_data.ZeroRate[1][1][i], tree_data.ZeroRate[1][2][i]); 
    }
    fprintf(fp,"#FwdRate\n");
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d FwdRate: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n",
                i, tree_data.FwdRate[0][0][i], tree_data.FwdRate[0][1][i], tree_data.FwdRate[0][2][i], 
                tree_data.FwdRate[1][0][i], tree_data.FwdRate[1][1][i], tree_data.FwdRate[1][2][i]); 
    }

    /* Internal assigment of zero curve */
    fprintf(fp,"CvDiff: %d %d\n",tree_data.CvDiff[0],tree_data.CvDiff[1] );
    fprintf(fp,"CvIdx1: %d %d\n",tree_data.CvIdx1[0],tree_data.CvIdx1[1] );
    fprintf(fp,"CvIdx2: %d %d\n",tree_data.CvIdx2[0],tree_data.CvIdx2[1] );
    fprintf(fp,"CvDisc: %d %d\n",tree_data.CvDisc[0],tree_data.CvDisc[1] );


                                
    /* IR Volatility */
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {   
        fprintf( fp, "i: %d spotVol: %12.10f %12.10f Rho: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n"
            , i, tree_data.SpotVol[0][i], tree_data.SpotVol[1][i]
            , tree_data.Rho[0][i], tree_data.Rho[1][i], tree_data.Rho[2][i]
            , tree_data.Rho[3][i], tree_data.Rho[4][i], tree_data.Rho[5][i]);
    }
    for (i = 0; i < MAXINDEX; ++i)
    {
        fprintf( fp, "i: %d Index: %c %c\n", i, tree_data.Index[0][i], tree_data.Index[1][i]);
    }

    /* Equity */
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d EqStock: %12.10f EqMidNode: %12.10f NodeSettleTime: %12.10f "
                "NodeSettleDate: %ld SpotEqVol: %12.10f EqVol: %12.10f\n"
                , i, tree_data.FwdEq[i], tree_data.EqMidNode[i], tree_data.NodeSettleTime[i]
                , tree_data.NodeSettleDate[i], tree_data.SpotEqVol[i]
                , tree_data.EqVol[i]); 
    }


    /* FX */                    
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d FwdFx: %12.10f FxMidNode: %12.10f SpotFxVol: %12.10f "
                "FxVol: %12.10f A1/2/3: %12.10f %12.10f %12.10f\n"
                , i, tree_data.FwdFx[i], tree_data.FxMidNode[i], tree_data.SpotFxVol[i]
                , tree_data.FxVol[i], tree_data.A1C[i], tree_data.A2C[i], tree_data.A3C[i] );
    }

    for (i = 0 ; i < tree_data.SmileIndex[tree_data.NbTP]+1; ++i)
    {
        fprintf(fp, "i: %d A2: %12.10f A3: %12.10f A1: %12.10f\n", i, tree_data.A1[i], tree_data.A2[i], tree_data.A3[i]);
    }

    /* Forward FX mapping function */
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d tMin: %ld tMax: %ld SmileIndex: %d \n" 
                    , i, tree_data.tMin[i], tree_data.tMax[i], tree_data.SmileIndex[i] );
    }
    fprintf( fp, "SmileCached %d\n", tree_data.FXsmileCache.isCached);
    if( tree_data.FXsmileCache.isCached )
    {
        for (i = 0 ; i < tree_data.FXsmileCache.nbCachePts; ++i)
        {
            fprintf( fp, "i: %d nbPtSPL: %d\n",i,tree_data.FXsmileCache.nbPtSPL[i]);
            for (j = 0 ; j < tree_data.FXsmileCache.nbPtSPL[i]; ++j)
            {
                fprintf(fp, "i: %d j: %d X: %12.10f K: %12.10f SPL: %12.10f SPL_Inv: %12.10f"
                    " gd: %12.10f gd_SPL: %12.10f kdx: %12.10f kdx_SPL: %12.10f\n", 
                    i,j, tree_data.FXsmileCache.X[i][j], tree_data.FXsmileCache.K[i][j], 
                    tree_data.FXsmileCache.SPL[i][j],tree_data.FXsmileCache.SPL_Inv[i][j],
                    tree_data.FXsmileCache.gd[i][j],tree_data.FXsmileCache.gd_SPL[i][j],
                    tree_data.FXsmileCache.kdX[i][j],tree_data.FXsmileCache.kdX_SPL[i][j]);
            }
        }
    }

    fprintf( fp, "CalcCheckSlices: %d\n", tree_data.CalcCheckSlices); /* it's a helper state variable, though */

    /* NOT SURE IF WE NEED THIS ONE */
    /* $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/
    /*int     *ProxySkewIdx;   */   /* index for time dependent approxim'g skew func*/

    
    /* Model */ 
    fprintf( fp, "TreeType: %d\n", tree_data.TreeType) ;
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d  Aweight: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n"
                , i, tree_data.Aweight[0][i], tree_data.Aweight[1][i], tree_data.Aweight[2][i]
                , tree_data.Aweight[3][i], tree_data.Aweight[4][i], tree_data.Aweight[5][i] );
    }
    for (i = 0 ; i < tree_data.NbTP; ++i)
    {
        fprintf( fp, "i: %d DriftCUPS[0]: %12.10f DriftCUPS[1]: %12.10f IrZCenter: %12.10f  %12.10f\n"
            , i, tree_data.DriftCUPS[0][i], tree_data.DriftCUPS[1][i], tree_data.IrZCenter[0][i],tree_data.IrZCenter[1][i]);
    }
    
    /* Tree geometry */
    fprintf( fp , "NbSigmaMax: %d\n", tree_data.NbSigmaMax );
    fprintf( fp , "NbSigmaDBL: %12.10f\n", tree_data.NbSigmaDBL );
    fprintf( fp , "NbIRSigmaMax: %12.10f\n", tree_data.NbIRSigmaMax );
    fprintf( fp , "Width: %d %d %d\n", tree_data.Width[0], tree_data.Width[1], tree_data.Width[2]);
    fprintf( fp , "HalfWidth: %d %d %d\n", tree_data.HalfWidth[0], tree_data.HalfWidth[1], tree_data.HalfWidth[2]);

    fprintf( fp , "InnerTreeWidth: %d %d %d\n", 
                    tree_data.InnerTreeWidth[0], tree_data.InnerTreeWidth[1], tree_data.InnerTreeWidth[2]);
    fprintf( fp , "InnerTreeHalfWidth: %d %d %d\n", 
                    tree_data.InnerTreeHalfWidth[0], tree_data.InnerTreeHalfWidth[1], tree_data.InnerTreeHalfWidth[2]);

    if (PrintTreeSizeInfo)
    {

       for (t = 0 ; t < tree_data.NbTP; ++t)
        {
             for (i = tree_data.Bottom1[t]; i <= tree_data.Top1[t]; i++)
             {
                 for (j = tree_data.Bottom2[t][i]; j <= tree_data.Top2[t][i]; j++)
                 {   
                     fprintf(fp, "t: %d i: %d j: %d Bottom3: %d Top3: %d\n"
                                ,t,i,j, tree_data.Bottom3[t][i][j], tree_data.Top3[t][i][j]);
                 }  /* for j */
            }  /* for i */
        } /* t  */

        for (t = 0 ; t < tree_data.NbTP; ++t)
        {
             for (i = tree_data.InnerTreeBottom1[t]; i <= tree_data.InnerTreeTop1[t]; i++)
             {
                 for (j = tree_data.InnerTreeBottom2[t][i]; j <= tree_data.InnerTreeTop2[t][i]; j++)
                 {   
                     fprintf(fp, "t: %d i: %d j: %d InnerTreeBottom3: %d InnerTreeTop3: %d\n"
                                ,t,i,j, tree_data.InnerTreeBottom3[t][i][j], tree_data.InnerTreeTop3[t][i][j]);
                 }  /* for j */
            }  /* for i */
        } /* t  */
        /* These dimensions are obsolete and not currently used*/
        /*int     *OutTop1;       */    /* Outer limits of tree in its 1st dimension */
        /*int     *OutBottom1;    */    /* Outer limits of tree in its 1st dimension */
        /*int    **OutTop2;       */    /* Outer limits of tree in its 2nd dimension */
        /*int    **OutBottom2;    */    /* Outer limits of tree in its 2nd dimension */
        /*int   ***OutTop3;       */    /* Outer limits of tree in its 3rd dimension */
        /*int   ***OutBottom3;    */    /* Outer limits of tree in its 3rd dimension */
    }

    /* moment matching arrays */
    fprintf( fp, "FxMomentMatching: %d\n", tree_data.FxMomentMatching );


    Status = SUCCESS;

RETURN:

    return Status;

}


/* -------------------------------------------------------------------- */
/* PrintFxData                                                          */
/*                                                                      */
/* "streaming" of the fxdata parameters into a file, the best possible  */
/* solution without having C++                                          */
/* */
int PrintFXData( FILE *fp, FX_DATA  fx_data )
{

    int Status = FAILURE;

    long i;


    if (fp == NULL) 
    {
        DR_Error("PrintFXData: file-pointer is set to NULL!"); 
        goto RETURN;
    }


    /* print all the information into the file */
    fprintf(fp, "Spot: %6.5f\n", fx_data.Spot);
    fprintf(fp, "Today: %ld\n", fx_data.Today);
    fprintf(fp, "ValueDate: %ld\n", fx_data.ValueDate);

    fprintf(fp, "SpotDays: %d\n", fx_data.SpotDays);
    fprintf(fp, "NbVol: %d\n", fx_data.NbVol);    
    for (i = 0; i < fx_data.NbVol; ++i)
    {
        fprintf(fp, "i: %ld VolDate: %ld FXVol: %12.10f\n\n", i, fx_data.VolDate[i],  fx_data.FxVol[i] );
    }

    fprintf(fp, "NbInpSpotVol: %d\n", fx_data.NbInpSpotVol);
    for (i = 0; i < fx_data.NbInpSpotVol; ++i)
    {
        fprintf(fp, "i: %ld InpSpotVolDate: %ld  InpSpotVol: %12.10f\n"
                    , i, fx_data.InpSpotVolDate[i],fx_data.InpSpotVol[i] );
    }
                 

    /* Vols and correlations */
    for (i = 0; i < MAXNBDATE; ++i)
    {
        fprintf( fp, "i: %ld RhoIRIR: %12.10f RhoIRFX: %12.10f RhoIRFX: %12.10f\n"
            , i, fx_data.Rho[0][i], fx_data.Rho[1][i], fx_data.Rho[2][i] );
    }

    /* Calibration flags */
    fprintf( fp, "FxCutOffFlag: %d\n", fx_data.FxCutOffFlag);
    fprintf( fp, "FxCutOffLast: %d\n", fx_data.FxCutOffLast);
    fprintf( fp, "FxBootStrapMode: %d\n", fx_data.FxBootStrapMode);
    fprintf( fp, "FxCutOffLevel: %12.10f\n", fx_data.FxCutOffLevel);


    /* Input Smile Parameters */
    fprintf( fp, "NbSmilePt: %d\n", fx_data.NbSmilePt);
    for (i = 0; i < fx_data.NbSmilePt; ++i ) 
    {
        fprintf( fp, "i: %ld SmileDate: %ld  a1: %12.10f a2: %12.10f a3: %12.10f\n"
            , i, fx_data.SmileDate[i], fx_data.a1[i], fx_data.a2[i], fx_data.a3[i]); 
    }


    Status = SUCCESS;

RETURN:

    return Status;

}


/* ---------------------------------------------------------- */
/* PrintEqData                                                */
/*                                                            */
/* "streaming" of the eq data into a file, the best possible  */
/* solution without having C++                                */
/* */
/* ATTENTION: not tested yet (17 May 2005)*/
int PrintEqData( FILE *fp, EQ_DATA eq_data )
{

    int Status = FAILURE;

    long i;


    if (fp == NULL) 
    {
        DR_Error("PrintFXData: file-pointer is set to NULL!"); 
        goto RETURN;
    }


    /* print out the EQ-data */
    fprintf( fp , "Spot: %6.5f\n", eq_data.Spot );    
    fprintf( fp, "SettleType: %c\n", eq_data.SettleType);
    for (i = 0 ; i < MAXNBEQDIV; ++i )
    {
        fprintf( fp, "i: %ld FwdType: %d FwdDate: %ld\n", i, eq_data.FwdType[i], eq_data.FwdDate[i]); 
    }


    /* Fwds */
    fprintf( fp , "NbFwd: %d\n", eq_data.NbFwd);
    for (i = 0 ; i < eq_data.NbFwd; ++i )
    {
        fprintf( fp , "i: %ld FwdDate: %ld Fwd: %12.10f\n", i, eq_data.FwdDate[i], eq_data.Fwd[i]);
    }

    /* Vol */
    fprintf( fp , "NbVol: %d\n", eq_data.NbVol);
    for (i = 0 ; i < eq_data.NbVol; ++i )
    {
        fprintf( fp , "i: %ld VolDate: %ld Vol: %12.10f\n", i, eq_data.VolDate[i], eq_data.Vol[i]);
    }

    /* Borrow */
    fprintf( fp , "NbBorrow: %d\n", eq_data.NbBorrow);
    for (i = 0 ; i < eq_data.NbBorrow; ++i )
    {
        fprintf( fp , "i: %ld BorrowDate: %ld Borrow: %12.10f\n", i, eq_data.BorrowDate[i], eq_data.Borrow[i]);
    }

    /* Credit */
    fprintf( fp , "NbCredit: %d\n", eq_data.NbCredit);
    for (i = 0 ; i < eq_data.NbCredit; ++i )
    {
        fprintf( fp , "i: %ld CreditDate: %ld Credit: %12.10f\n", i, eq_data.CreditDate[i], eq_data.Credit[i]);
    }

    fprintf( fp , "NbSettle: %d\n", eq_data.NbSettle);
    for ( i = 0 ; i < eq_data.NbSettle; ++i )
    {
        fprintf( fp, "i: %ld LastTrading: %ld SettleDate: %ld\n", i, eq_data.LastTrading[i], eq_data.SettleDate[i] );
    }

    for( i = 0 ; i < MAXNBDATE; ++i )
    {
        fprintf( fp , "i: %ld  RhoEq_IRf: %12.10f  RhoEq_IRd: %12.10f  RhoEq_FX: %12.10f\n",
            i, eq_data.Rho[0][i], eq_data.Rho[1][i], eq_data.Rho[2][i] );
    }

    Status = SUCCESS;

RETURN:

    return Status;

}


/* ---------------------------------------------------------------- */
/* PrintTCurve                                                      */
/*                                                                  */
/* "streaming" of the t_curve data onto a file, the best possible   */
/* solution without having C++                                      */
/* */
/* ATTENTION: not tested yet: 17 May 2005 AlexK                     */

/* SK - passing T_CURVE by value is problematic
 */
int PrintTCurve( FILE *fp, T_CURVE  t_curve )
{

    int Status = FAILURE;

    if (fp == NULL) 
    {
        DR_Error("PrintFXData: file-pointer is set to NULL!"); 
        goto RETURN;
    }

    EslPrintZeroCurve(&t_curve, fp);
    Status = SUCCESS;

RETURN:

    return Status;

}


/* --------------------------------------------------------------- */
/* PrintMktVolData                                                 */
/*                                                                 */
/* "streaming" of the mktVol Data into a file, the best possible   */
/* solution without having C++                                     */
/* */
int PrintMktVolData( FILE *fp, MKTVOL_DATA  mktvol_data )
{

    int Status = FAILURE;

    long i;


    if (fp == NULL) 
    {
        DR_Error("PrintFXData: file-pointer is set to NULL!"); 
        goto RETURN;
    }

    /* Volatility curve */
    fprintf( fp, "BaseDate: %ld\n", mktvol_data.BaseDate);
    fprintf( fp, "NbVol: %d\n", mktvol_data.NbVol);

    for (i = 0 ; i < mktvol_data.NbVol; ++i)
    {
        fprintf( fp, "i: %ld VolDate: %ld Vol: %12.10f VolUsed: %d\n"
            , i, mktvol_data.VolDate[i], mktvol_data.Vol[i], mktvol_data.VolUsed[i]);
    }

    fprintf( fp, "CetNbIter: %d\n", mktvol_data.CetNbIter );
    fprintf( fp, "CetVegaError: %12.10f\n", mktvol_data.CetVegaError );
    for (i = 0 ; i < MAXNBDATE; ++i )
    {

        fprintf( fp , "i: %ld Aweight: %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n", i
                     , mktvol_data.Aweight[0][i], mktvol_data.Aweight[1][i], mktvol_data.Aweight[2][i]
                     , mktvol_data.Aweight[3][i], mktvol_data.Aweight[4][i], mktvol_data.Aweight[5][i] );
    }

    /* Forward swap information */
    fprintf( fp , "Freq: %c\n", mktvol_data.Freq );
    fprintf( fp , "DCC: %c\n", mktvol_data.DCC );

    for (i = 0 ; i < MAXNBDATE; ++i )
    {
        fprintf( fp , "i: %ld SwapSt: %ld SwapMat: %ld\n", i, mktvol_data.SwapSt[i], mktvol_data.SwapMat[i] );
    }

    /* Model parameters */
    fprintf( fp , "QLeft: %12.10f\n", mktvol_data.QLeft);
    fprintf( fp , "QRight: %12.10f\n", mktvol_data.QRight);
    fprintf( fp , "FwdShift: %12.10f\n", mktvol_data.FwdShift);
    fprintf( fp , "Alpha: %12.10f %12.10f %12.10f\n", mktvol_data.Alpha[0], mktvol_data.Alpha[1], mktvol_data.Alpha[2]);
    fprintf( fp , "Beta: %12.10f %12.10f %12.10f\n", mktvol_data.Beta[0], mktvol_data.Beta[1], mktvol_data.Beta[2]);
    fprintf( fp , "Rho: %12.10f %12.10f %12.10f\n", mktvol_data.Rho[0], mktvol_data.Rho[1],mktvol_data.Rho[2]);


    /* Calibration flags */
    fprintf( fp, "SkipFlag: %d\n", mktvol_data.SkipFlag);
    fprintf( fp, "CalibFlag: %d\n", mktvol_data.CalibFlag);


    Status = SUCCESS;

RETURN:

    return Status;

}

int PrintSlice_2d(HYB3_TREE_DATA *tree_data,
                  TSLICE          slice,
                  int             t,
                  FILE           *fp)
{
    int    Top1,      Bottom1;           /* Tree limits (IR foreign)       */
    int   *Top2,     *Bottom2;           /* Tree limits (IR domestic)      */

    int    offset;
    int    i,j;

    double *sliceL = NULL;

    Top1      = tree_data->OutTop1[t];
    Top2      = tree_data->OutTop2[t];
    Bottom1   = tree_data->OutBottom1[t];
    Bottom2   = tree_data->OutBottom2[t];

    for(i = Bottom1 ; i <= Top1; i++)
    {
        offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
        sliceL  = slice + offset;
            
        for(j = Bottom2[i]; j<= Top2[i]; j++)
        {
            fprintf(fp, "%12.10f\n", sliceL[j]);
        }/* for j */
    }/*for i */

    return SUCCESS;
}

int PrintSlice_3d(HYB3_TREE_DATA *tree_data,
               TSLICE          slice,
               int             t,
               FILE           *fp)
{
    int    Top1,      Bottom1;           /* Tree limits (IR foreign)       */
    int   *Top2,     *Bottom2;           /* Tree limits (IR domestic)      */
    int  **Top3,    **Bottom3;           /* Tree limits (third asset)      */

    int    offset;
    int    i,j,k;

    double *sliceL = NULL;

    Top1      = tree_data->OutTop1[t];
    Top2      = tree_data->OutTop2[t];
    Top3      = tree_data->OutTop3[t];
    Bottom1   = tree_data->OutBottom1[t];
    Bottom2   = tree_data->OutBottom2[t];
    Bottom3   = tree_data->OutBottom3[t];

    for(i = Bottom1 ; i <= Top1; i++)
    {
        for(j = Bottom2[i]; j<= Top2[i]; j++)
        {
            offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
            sliceL  = slice + offset;
            
            for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
            {
                fprintf(fp, "%12.10f\n", sliceL[k]);
            }
        }/* for j */
    }/*for i */

    return SUCCESS;
}

int PrintCriticalDates(FILE* stream, int nbCritDate, CRIT_DATE* critDate)
{
    int i;
    for (i = 0; i < nbCritDate; i++)
    {
        fprintf(stream, "%d %d %ld\n", i, critDate[i].Type, critDate[i].CritDate);
    }
    return SUCCESS;
}

/*****  Hyb3_DoubleArrayFloorIdx  ************************************************/
/*
*      Find the floor index "i" in a double array y[n] given search value x
*      such that y[i] <= x.
*      if x < y[0], i = -1
*/
int    Hyb3_DoubleArrayFloorIdx(double const* y, /* (I) Double ascending array */
                                int           n, /* (I) Size of array          */
                                double        x) /* (I) Search value           */
{
        register int jl, ju, jm;
        int    j;    /* return result */

        jl = (-1);
        ju = n;
        while (ju-jl > 1) {
                jm=(ju+jl) >> 1;
                if (x >= y[jm])
                        jl=jm;
                else
                        ju=jm;
        }
        j=jl;

        return(j);
}

/*****  Hyb3_DoubleQuadraticInterp  **********************************************/
/*
*      Special quadratic interpolation.  
*      Input X(i) must be in ascending order
*      and duplicated elements in both X(i) and Y(i) are allowed,
*/
int Hyb3_DoubleQuadraticInterp(
    double const* xa,     /* (I) array of X(i) (i=0,..,m-1) */
    double const* ya,     /* (I) array of Y(i) (i=0,..,m-1) */
    int           m,      /* (I) arrays size */
    double        x,      /* (I) point to intepolate */
    double*       y)      /* (O) interpolated value */
{
    static  char    routine[] = "Hyb3_DoubleQuadraticInterp";
    
    int    i;
    double d0, d1, d2;
    char   interpType;

    /* Input check
     * if X(i)=X(j), then Y(i)=Y(j)    
     * Duplication of X(i) is allowed 
     */
    for (i=1; i<m; i++)
    {
        if ( IS_EQUAL(xa[i], xa[i-1]) &&
            !IS_EQUAL(ya[i], ya[i-1]))
        {
            DR_Error("%s: X[%d] = X[%d], but Y[%d] != Y[%d]",
                     routine,
                     i-1, i, i-1, i);
            return FAILURE;
        }
    }

    
    /* Flat on both ends */
    if (x >= xa[m-1])
    {
        *y = ya[m-1];
        return SUCCESS;
    }
    
    if (x <= xa[0])
    {
        *y = ya[0];
        return SUCCESS;
    }

    /* Interpolate in between */

    /* 1. Two points only */
    if (m == 2)
    {
        linterp (x, 
                 y, 
                 xa[0], xa[1], 
                 ya[0], ya[1]);

        return SUCCESS;
    }
    
    /* Search the corresponding index */
    i = Hyb3_DoubleArrayFloorIdx(xa, m, x);

    /* 2. Between two flat points */
    if (IS_EQUAL(xa[i], xa[i+1])) 
    {
        *y = ya[i];
        return SUCCESS;
    }

    /* 3. Linear if both left and right buckets are flat,
     *    otherwise quadratic interpolation involving current period and
     *    either the one on the right or on the left
     */
    if (i == m-2)        /* Right-most bucket */
    {
        if (IS_EQUAL(xa[i-1], xa[i]))  /* left bucket is flat */
        {
            interpType = 'L';
        }
        else    /* quadratic involving i-1, i, i+1 */
        {
            interpType = 'Q';
            i = i-1;
        }
    }
    else if (i == 0)    /* Left-most bucket */
    {
        if (IS_EQUAL(xa[i+1], xa[i+2]))  /* right bucket is flat */
        {
            interpType = 'L';
        }
        else    /* quadratic involving i, i+1, i+2 */
        {
            interpType = 'Q';
        }
    }
    else                /* In the middle  */
    {
        if (!IS_EQUAL(xa[i+1], xa[i+2]))  /* right bucket is NOT flat */
        {
            interpType = 'Q';
        }
        else if (!IS_EQUAL(xa[i-1], xa[i]))  /* left bucket is NOT flat */
        {
            interpType = 'Q';
            i = i-1;
        }
        else        /* Both right and left buckets are flat */
        {
            interpType = 'L';
        }
    }


    /* Interpolate based on interp type */

    if (interpType == 'L')     /* Linear */
    {
        linterp (x, 
                 y, 
                 xa[i], xa[i+1], 
                 ya[i], ya[i+1]);

        return SUCCESS;
    }
    else                         /* Quadratic */
    {
        d0 = 1. / ((xa[i]-xa[i+1])
                  *(xa[i]-xa[i+2]));
        d1 = 1. / ((xa[i+1]-xa[i])
                  *(xa[i+1]-xa[i+2]));
        d2 = 1. / ((xa[i+2]-xa[i])
                  *(xa[i+2]-xa[i+1]));
 
        sqinterp(xa[i],
                 xa[i+1],
                 xa[i+2],
                 ya[i],  
                 ya[i+1],
                 ya[i+2],
                 d0,
                 d1,
                 d2,
                 x, 
                 y);

        return SUCCESS;
    }

}    /* Hyb3_DoubleQuadraticInterp */

/*************** Hyb3_Triangulation4D ***********************************************/
/**Produces the "usual" Lower triangular matrix for orthogonalising 4 corrolated 
   factors.  (Cholesky decomposition)                                            
   NOTE: If Nbfac < 4 then unused matrix elements are set to zero.             
                                                                               
 *****************************************************************************/

int Hyb3_Triangulation4D (
        double TriangMtx[4][4],/**<(O) Lower triangular matrix transformation*/
        int    Nbfac,          /**<(I) Number of factors                     */
        double Rho[4][4])           /**<(I) correlation coefficients              */
{

    int i,j,status = FAILURE;
   

    if ((Nbfac > 1) && (Rho == NULL))
    {
        DR_Error ("Invalid pointer input to Triangulation\n");
        return (status);
    }

    if (Nbfac != 1 && Nbfac != 2 && Nbfac != 3 && Nbfac != 4)
    {
        DR_Error ("Invalid number of factors in Triangulation: "
                    "should be 1,2,3 or 4\n");
        return (status);
    }

    /* initialise matrix*/
    for (i = 0 ; i < 4; i++)
    {
        for(j = 0; j < 4; j++)
        {
            TriangMtx[i][j] = 0.;

        }/* for j*/

    }/* for i */


    TriangMtx[0][0] = 1.;

    if  (Nbfac > 1)
    {
        TriangMtx[1][0] = Rho[1][0];

        if (1.- Rho[1][0] * Rho[1][0] < TINY)
        {
            DR_Error ("bad choice of correlation parameter Rho[1][0]\n");
            goto RETURN;
        }

        TriangMtx[1][1] = sqrt (1.0 - SQUARE(TriangMtx[1][0]));
    }

    if (Nbfac > 2)
    {
        TriangMtx[2][0] = Rho[2][0];
        TriangMtx[2][1] = (Rho[2][1] - TriangMtx[2][0] * TriangMtx[1][0]) 
                           / TriangMtx[1][1];

        if ( (1. - SQUARE( TriangMtx[2][0]) - SQUARE(TriangMtx[2][1])) < TINY)
        {
            DR_Error ("Bad choice of correlation parameters:\n"
                       "correlation matrix is not positive definite!!\n");
            goto RETURN;
        }

        TriangMtx[2][2] = sqrt (1. - SQUARE( TriangMtx[2][0]) 
                                   - SQUARE( TriangMtx[2][1]));
    }
    if (Nbfac > 3)
    {
        TriangMtx[3][0] = Rho[3][0] / TriangMtx[0][0];
        TriangMtx[3][1] = (Rho[3][1] - TriangMtx[3][0] * TriangMtx[1][0]) 
                           / TriangMtx[1][1];
        TriangMtx[3][2] = (Rho[3][2] 
                            - TriangMtx[3][0] * TriangMtx[2][0] 
                            - TriangMtx[3][1] * TriangMtx[2][1])
                           / TriangMtx[2][2];
        
        if ( (1 - SQUARE(TriangMtx[3][0]) 
                - SQUARE(TriangMtx[3][1]) 
                - SQUARE(TriangMtx[3][2])) < TINY)
        {
            DR_Error ("Bad choice of correlation parameters:\n"
                       "correlation matrix is not positive definite!!\n");
            goto RETURN;
        }
        TriangMtx[3][3] = sqrt(1 - SQUARE(TriangMtx[3][0]) 
                                 - SQUARE(TriangMtx[3][1]) 
                                 - SQUARE(TriangMtx[3][2]));
    }

    status = SUCCESS;

RETURN:
    return (status);

}/* Hyb3_Triangulation4D */


