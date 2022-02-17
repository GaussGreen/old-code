/****************************************************************************/
/*        Memory allocation.                                                */
/****************************************************************************/
/*        ALLOC.c                                                           */
/****************************************************************************/



/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "cupslib.h"

/*****  Hyb3_SwapRate_Init *******************************************************/
/**
*         Initialize SwapRate pointers to NULL.
*/
void     Hyb3_SwapRate_Init(SWAPRATE_DATA *sr)
{
    int i;

    for (i=0; i<MAXNBDATE; ++i)
    {
        sr->SwapRate[i] = NULL;
    }

    return;
}

/*****  Hyb3_ForDim  *************************************************************/
/**
*   For each tree mode returns the dimension of the slice corresponding to the 
*   foreign rate.
*   If there is no foreign slice then returns -1 and prints error message.
*/
int      Hyb3_ForDim(HYB3_TREE_DATA *tree_data)
{
    int fDim;
    if(tree_data->TreeType == TTYPE_2IR2F1D)
    {
        fDim = 2;        
    }
    else if((tree_data->TreeType == TTYPE_2IR)     ||
            (tree_data->TreeType == TTYPE_FX2IR)   ||
            (tree_data->TreeType == TTYPE_EQF2IR)  ||
            (tree_data->TreeType == TTYPE_EQC2IR)  ||
            (tree_data->TreeType == TTYPE_EQD2IR)  ||
            (tree_data->TreeType == TTYPE_EQDFX2IR)||
            (tree_data->TreeType == TTYPE_EQFFX2IR) )
    {
        fDim = 1;
    }
    else if(tree_data->TreeType == TTYPE_EQ1IR)
    {
        DR_Error("Hyb3_ForDim: tree_data->TreeType = %d not supported", tree_data->TreeType);
        fDim = -1;
    }
    else
    {
        DR_Error("Hyb3_ForDim: tree_data->TreeType = %d not supported", tree_data->TreeType);
        fDim = -1;
    }

    return fDim;
}

/*****  Hyb3_DomDim  *************************************************************/
/**
*   For each tree mode returns the dimension of the slice corresponding to the 
*   domestic rate.
*/
int      Hyb3_DomDim(HYB3_TREE_DATA *tree_data)
{
    int dDim;

    switch(tree_data->TreeType)
    {
    case TTYPE_2IR: 
        dDim = 2;
        break;
    case TTYPE_FX2IR: 
        dDim = 2;
        break;
    case TTYPE_EQD2IR:
        dDim = 2;
        break;
    case TTYPE_EQF2IR:
        dDim = 2;
        break;
    case TTYPE_EQC2IR:
        dDim = 2;
        break;
    case TTYPE_EQ1IR:
        dDim = 1;
        break;
    case TTYPE_2IR2F1D:
        dDim = 3;
        break;
    case TTYPE_EQDFX2IR:
        dDim = 2;
        break;
    case TTYPE_EQFFX2IR:
        dDim = 2;
        break;        
    default:
        DR_Error("Hyb3_DomDim: tree_data->TreeType = %d not supported", tree_data->TreeType);
        dDim = -1;
    }

    return dDim;
}

int      Hyb3_SwapRate_Alloc(SWAPRATE_DATA  *sr,
                             HYB3_TREE_DATA *tree_data)
{
    int status = FAILURE;

    int i;

    int dDim, fDim, dim;

    /* basic checks */
    if (sr->NbResetDate < 0)
    {
        DR_Error("NbResetDate in SwapRate structure must be non-negative.\n");
        goto RETURN;
    }

    if (sr->NbResetDate == 0)
    {
        return SUCCESS;
    }

    if (tree_data->TreeType == TTYPE_EQ1IR)
    {
        if (sr->Denom != 'D')
        {
            DR_Error("Swap rates must be domestic if tree mode is EQ1IR.\n");
            goto RETURN;
        }
    }
    
    /* set dim to the dimension of the slices to be allocated */
    dDim = Hyb3_DomDim(tree_data);
    fDim = Hyb3_ForDim(tree_data);
    dim  = (sr->Denom == 'D') ? dDim : fDim;

    for (i=0; i<sr->NbResetDate; ++i)
    {
        sr->SwapRate[i] = Hyb3_Alloc_Slice(tree_data, dim);

        if (sr->SwapRate[i] == NULL)
        {
            goto RETURN;
        }
    }
    
    status = SUCCESS;

  RETURN:
    
    if (status == FAILURE)
    {
        DR_Error("Hyb3_SwapRate_Alloc: "
                 "Unable to allocate memory for SwapRate structure.\n");
    }
    return(status);
}

/*****  Hyb3_SwapRate_Free *******************************************************/
/**
*         Free SwapRate structures.
*/
void     Hyb3_SwapRate_Free(SWAPRATE_DATA  *sr,
                            HYB3_TREE_DATA *tree_data)
{
    int i;

    int dDim, fDim, dim;

    /* set dim to the dimension of the slices to be allocated */
    dDim = Hyb3_DomDim(tree_data);
    fDim = Hyb3_ForDim(tree_data);
    dim  = (sr->Denom == 'D') ? dDim : fDim;

    for (i=0; i<sr->NbResetDate; ++i)
    {
        Hyb3_Free_Slice(sr->SwapRate[i], tree_data, dim);
        sr->SwapRate[i] = NULL;
    }
}

int Hyb3_TreeSim_Init(TREESIM_DATA *ts, int NbSwapRateData)
{
    int i;

    ts->NbPathDate = 0;

    ts->NbPathAll = 0;
    ts->NbPathSub = 0;

    ts->MaxNbState = 0;

    for (i=0; i<MAXNBDATE; ++i)
    {
        ts->NbState [i] = 0;

        ts->PathDate[i] = 0;

        ts->PathAll [i] = NULL;
        ts->PathSub [i] = NULL;
        ts->State   [i] = NULL;
    }

    ts->SwapRateData = NULL;

    ts->i = NULL;
    ts->j = NULL;
    ts->k = NULL;

    ts->seed = 0;

    ts->NbSwapRateData = NbSwapRateData;

    if (ts->NbSwapRateData < 0)
    {
        DR_Error("Nb of swap rate data structures must be >= 0.\n");
        return FAILURE;
    }

    if (ts->NbSwapRateData > 0)
    {
        ts->SwapRateData = (SWAPRATE_DATA *)
            DR_Array(TYPE_SWAPRATE_DATA, 0, ts->NbSwapRateData-1);

        if (ts->SwapRateData == NULL)
        {
            DR_Error("Unable to allocate memory for swap rate data.\n");
            return FAILURE;
        }
    }

    for (i=0; i<ts->NbSwapRateData-1; ++i)
    {
        Hyb3_SwapRate_Init(&(ts->SwapRateData[i]));
    }

    return SUCCESS;
}
    
int Hyb3_TreeSim_Alloc(TREESIM_DATA *ts, HYB3_TREE_DATA *tree_data)
{
    int i;

    if (ts->NbPathDate<0)
    {
        DR_Error("TreeSim number of sample dates must be >= 0.\n");
        return FAILURE;
    }

    if (ts->NbPathAll < ts->NbPathSub)
    {
        DR_Error("Size of full sample must be >= size of abridged sample.\n");
        return FAILURE;
    }

    if ((ts->NbPathSub < 4))
    {
        DR_Error("Size of abridged sample must be >= 4.\n");
        return FAILURE;
    }

    if (ts->NbPathDate>MAXNBDATE)
    {
        DR_Error("TreeSim number of path dates must be <= MAXNBDATE.\n");
        return FAILURE;
    }

    for (i=0; i<ts->NbSwapRateData; ++i)
    {
        if (Hyb3_SwapRate_Alloc(&(ts->SwapRateData[i]),
                                tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }

    ts->i = (int *)DR_Array(INT, 0, ts->NbPathAll-1);
    ts->j = (int *)DR_Array(INT, 0, ts->NbPathAll-1);
    ts->k = (int *)DR_Array(INT, 0, ts->NbPathAll-1);

    for (i=0; i<ts->NbPathDate; ++i)
    {
        ts->PathAll [i] = (double *)DR_Array(DOUBLE, 0, ts->NbPathAll-1);
        ts->PathSub [i] = (double *)DR_Array(DOUBLE, 0, ts->NbPathSub-1);
        ts->State   [i] = (double *)DR_Array(DOUBLE, 0, ts->NbPathSub-1);

        if ((ts->PathAll [i] == NULL) ||
            (ts->PathSub [i] == NULL) ||
            (ts->State   [i] == NULL))
        {
            DR_Error("Could not allocate memory for state variable simulation.\n");
            return FAILURE;
        }
    }

    return SUCCESS;
}

void Hyb3_TreeSim_FreeSR(TREESIM_DATA *ts, HYB3_TREE_DATA *tree_data)
{
    int i;

    for (i=0; i<ts->NbSwapRateData; ++i)
    {
        Hyb3_SwapRate_Free(&(ts->SwapRateData[i]),
                           tree_data);
    }

    if (ts->NbSwapRateData > 0)
    {
        Free_DR_Array(ts->SwapRateData,
                      TYPE_SWAPRATE_DATA, 0, ts->NbSwapRateData-1);
    }

    ts->SwapRateData = NULL;

    return;
}

void Hyb3_TreeSim_Free(TREESIM_DATA *ts, HYB3_TREE_DATA *tree_data)
{
    int i;

    if (ts->SwapRateData != NULL)
    {
        Hyb3_TreeSim_FreeSR(ts, tree_data);
    }

    Free_DR_Array(ts->i, INT, 0, ts->NbPathAll-1);
    Free_DR_Array(ts->j, INT, 0, ts->NbPathAll-1);
    Free_DR_Array(ts->k, INT, 0, ts->NbPathAll-1);

    for (i=0; i<ts->NbPathDate; ++i)
    {
        Free_DR_Array(ts->PathAll [i], DOUBLE, 0, ts->NbPathAll-1);
        Free_DR_Array(ts->PathSub [i], DOUBLE, 0, ts->NbPathSub-1);
        Free_DR_Array(ts->State   [i], DOUBLE, 0, ts->NbPathSub-1);
    }

    return;
}
    
/*****  Hyb3_Tree_Init  **********************************************************/
/**
*         Initialize tree pointers to NULL.
*/
void     Hyb3_Tree_Init(HYB3_TREE_DATA  *tree_data)   /**< Tree building data structure */
{
     int  i;
     int  j;

      /* The value -999 indicates that NbSigmaDBL has not been initialised */
      /* this value is tested for in Build tree                            */
      tree_data->NbSigmaDBL = -9999;

      tree_data->CalcCheckSlices = FALSE;

      tree_data->FxMomentMatching = TRUE;

      tree_data->debug = FALSE;

     /* This initialisation is important because some of the */
     /* freeing will rely on the number of timepoints.       */
     tree_data->NbTP = 0;

     for (i = 0; i < NBCRITDATE; i++)
     {
          tree_data->TPType[i] = NULL;
          tree_data->CritDate[i] = NULL;
     } 

     tree_data->Top1    = NULL;
     tree_data->Bottom1 = NULL;
     tree_data->Top2    = NULL;
     tree_data->Bottom2 = NULL;
     tree_data->Top3    = NULL;
     tree_data->Bottom3 = NULL;
     tree_data->Top4    = NULL;
     tree_data->Bottom4 = NULL;

     tree_data->InnerTreeTop1    = NULL;
     tree_data->InnerTreeBottom1 = NULL;
     tree_data->InnerTreeTop2    = NULL;
     tree_data->InnerTreeBottom2 = NULL;
     tree_data->InnerTreeTop3    = NULL;
     tree_data->InnerTreeBottom3 = NULL;
     tree_data->InnerTreeTop4    = NULL;
     tree_data->InnerTreeBottom4 = NULL;

     tree_data->OutTop1    = NULL;
     tree_data->OutBottom1 = NULL;
     tree_data->OutTop2    = NULL;
     tree_data->OutBottom2 = NULL;
     tree_data->OutTop3    = NULL;
     tree_data->OutBottom3 = NULL;
     tree_data->OutTop4    = NULL;
     tree_data->OutBottom4 = NULL;

     for (i = 0; i < 4; i++)
     {
        tree_data->Width[i] = 0;
        tree_data->HalfWidth[i] = 0;

        tree_data->xWidth[i] = 0;
        tree_data->xHalfWidth[i] = 0;
     }

     tree_data->xDate = 0;
     tree_data->xT    = 0;

     for (i = 0; i < 2; i++)
     {
         for (j=0; j<3; j++)
         {
             tree_data->ZeroCoupon[i][j] = NULL;
             tree_data->ZeroRate[i][j]   = NULL;
             tree_data->FwdRate[i][j]    = NULL;
         }
     } 

     
    for (i = 0; i<2; i++)
    {
        for (j = 0; j < 3; j++)
        {
            tree_data->IrAweight[i][j] = NULL;
        }
     }  

     tree_data->SpotVol[0]    = NULL;
     tree_data->SpotVol[1]    = NULL;
     tree_data->IrZCenter[0]  = NULL;
     tree_data->IrZCenter[1]  = NULL;
     tree_data->DriftCUPS[0]  = NULL;
     tree_data->DriftCUPS[1]  = NULL;
     tree_data->TPDate        = NULL;
     tree_data->Length        = NULL;
     tree_data->LengthJ       = NULL;

     tree_data->EDevStPrice = NULL;

     for (i = 0; i < 6; i++)
     {
          tree_data->Rho[i]     = NULL;
     }

     for (i = 0; i < 10; i++)
     {
          tree_data->Aweight[i] = NULL;
     }

     /* Equity */
     tree_data->FwdEq          = NULL;
     tree_data->EqMidNode      = NULL;
     tree_data->NodeSettleTime = NULL;
     tree_data->NodeSettleDate = NULL;
     tree_data->EqVol          = NULL;
     tree_data->SpotEqVol      = NULL;

     /* FX */
     tree_data->FwdFx      = NULL;
     tree_data->FxMidNode  = NULL;
     tree_data->FxVol      = NULL;
     tree_data->SpotFxVol  = NULL;

     /* Smile (Equity or FX) */
     tree_data->A1         = NULL;
     tree_data->A2         = NULL;
     tree_data->A3         = NULL;
     tree_data->A1C        = NULL;
     tree_data->A2C        = NULL;
     tree_data->A3C        = NULL;
     tree_data->SmileIndex = NULL;
     tree_data->tMin       = NULL;
     tree_data->tMax       = NULL;


     /*Set the initial values for the flags */
     tree_data->FXsmileCache.isCached  = FALSE;

     /* K function */
     tree_data->FXsmileCache.X          = NULL;
     tree_data->FXsmileCache.K          = NULL;
     tree_data->FXsmileCache.SPL        = NULL;
     tree_data->FXsmileCache.SPL_Inv    = NULL;
     tree_data->FXsmileCache.nbCachePts = 0;
     tree_data->FXsmileCache.nbPtSPL    = NULL;

     /* gDash and KdashTimesX function (generated together with K function) */
     tree_data->FXsmileCache.gd         = NULL;
     tree_data->FXsmileCache.gd_SPL     = NULL;
     tree_data->FXsmileCache.kdX        = NULL;
     tree_data->FXsmileCache.kdX_SPL    = NULL;


     /* Express DEV tool */
     tree_data->NbEDevDates = 0;
     tree_data->EDevDate    = NULL;
     tree_data->EDevStPrice = NULL;


     return;

}  /* Hyb3_Tree_Init */

/*****  Hyb3_Tree_Alloc  *********************************************************/
/**
*         Allocation of memory  for the arrays  constituting  the tree.
*         The limits of the tree are allocated directly in Hyb3_Tree_Limits.
*         We allocate two extra  nodes before the start and the end of 
*         the tree to deal with time difference (see dev.c for example).
*/
int     Hyb3_Tree_Alloc(HYB3_TREE_DATA   *tree_data)
{
    int    i;
    int    j;
    int    NbTP;  /* Total number of time points in the tree   */
    
    int    status = FAILURE;

    /* Total number of time points */
    NbTP = tree_data->NbTP;

    /* Critical dates */
    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->TPType[i]   = (int *)DR_Array(INT,      -1, NbTP+1);
        tree_data->CritDate[i] = (CRIT_DATE *)DR_Array(CRITDATE,-1, NbTP+1);

        if ( tree_data->TPType[i]   == NULL  ||
             tree_data->CritDate[i] == NULL  )
        {
            goto RETURN;
        }
    }

    /* All the ellipse limits */
    tree_data->Top1    =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->Bottom1 =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->Top2    =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->Bottom2 =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->Top3    = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->Bottom3 = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->Top4    = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);
    tree_data->Bottom4 = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);

    if ( (tree_data->Top1    == NULL)
      || (tree_data->Bottom1 == NULL)
      || (tree_data->Top2    == NULL)
      || (tree_data->Bottom2 == NULL)
      || (tree_data->Top3    == NULL)
      || (tree_data->Bottom3 == NULL)
      || (tree_data->Top4    == NULL)
      || (tree_data->Bottom4 == NULL))
    {
        goto RETURN;
    }

    /* All the limits for the outer ellipses */
    tree_data->OutTop1    =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->OutBottom1 =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->OutTop2    =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->OutBottom2 =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->OutTop3    = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->OutBottom3 = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->OutTop4    = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);
    tree_data->OutBottom4 = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);

    if ( (tree_data->OutTop1    == NULL)
      || (tree_data->OutBottom1 == NULL)
      || (tree_data->OutTop2    == NULL)
      || (tree_data->OutBottom2 == NULL)
      || (tree_data->OutTop3    == NULL)
      || (tree_data->OutBottom3 == NULL)
      || (tree_data->OutTop4    == NULL)
      || (tree_data->OutBottom4 == NULL))
    {
        goto RETURN;
    }

    /* All the Inner ellipse limits */
    tree_data->InnerTreeTop1    =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->InnerTreeBottom1 =   (int *)DR_Array(INT, -1, NbTP+1);
    tree_data->InnerTreeTop2    =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->InnerTreeBottom2 =  (int **)DR_Array(INT_PTR, -1, NbTP+1);
    tree_data->InnerTreeTop3    = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->InnerTreeBottom3 = (int ***)DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data->InnerTreeTop4    = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);
    tree_data->InnerTreeBottom4 = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);

    if ( (tree_data->InnerTreeTop1    == NULL)
      || (tree_data->InnerTreeBottom1 == NULL)
      || (tree_data->InnerTreeTop2    == NULL)
      || (tree_data->InnerTreeBottom2 == NULL)
      || (tree_data->InnerTreeTop3    == NULL)
      || (tree_data->InnerTreeBottom3 == NULL)
      || (tree_data->InnerTreeTop4    == NULL)
      || (tree_data->InnerTreeBottom4 == NULL))
    {
        goto RETURN;
    }

    /* Zero bond prices, zero rates and forwards */
    for (i=0; i < 2; i++)
    {
        for (j=0; j < 3; j++)
        {
            tree_data->ZeroCoupon[i][j] = (double *)DR_Array(DOUBLE,-1,NbTP+1);                  
            tree_data->ZeroRate[i][j]   = (double *)DR_Array(DOUBLE,-1,NbTP+1);                  
            tree_data->FwdRate[i][j]    = (double *)DR_Array(DOUBLE,-1,NbTP+1);                	

            if ( (tree_data->ZeroCoupon[i] == NULL)
             ||  (tree_data->ZeroRate[i]   == NULL)
             ||  (tree_data->FwdRate[i]    == NULL))
            {
                goto RETURN;
            }
        }  /* for j */
    } /* for i */

    tree_data->SpotVol[0] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->SpotVol[1] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    if ((tree_data->SpotVol[0] == NULL) ||
        (tree_data->SpotVol[1] == NULL))
    {
        goto RETURN;
    }
    for (i=0; i < 2; i++)
    {
        for (j = 0; j < 3; j++)
        {
            tree_data->IrAweight[i][j] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
            if ( tree_data->IrAweight[i][j] == NULL )
            {
                goto RETURN;
            }
        }
    }

    tree_data->IrZCenter[0] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->IrZCenter[1] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->DriftCUPS[0] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->DriftCUPS[1] = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    if ((tree_data->IrZCenter[0] == NULL) ||
        (tree_data->IrZCenter[1] == NULL) ||
        (tree_data->DriftCUPS[0] == NULL) ||
        (tree_data->DriftCUPS[1] == NULL))
    {
        goto RETURN;
    }

    tree_data->TPDate  = (long *)DR_Array (LONG, -1, NbTP+1);
    if (tree_data->TPDate == NULL)
    {
        goto RETURN;
    }

    tree_data->Length    = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->LengthJ   = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    if ((tree_data->Length  == NULL) ||
        (tree_data->LengthJ == NULL))
    {
        goto RETURN;
    }

    for (i = 0; i < 6; i++)
    {
        tree_data->Rho[i] = (double *)DR_Array(DOUBLE, -1, NbTP+1);

        if (tree_data->Rho[i] == NULL)
        {
            goto RETURN;
        }    

    }  /* for i */

    for (i = 0; i < 10; i++)
    {
        tree_data->Aweight[i] = (double *)DR_Array(DOUBLE, -1, NbTP+1);

        if (tree_data->Aweight[i] == NULL)
        {
            goto RETURN;
        }    

    }  /* for i */

    /* Equity */
    tree_data->FwdEq           = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->EqMidNode       = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->NodeSettleTime  = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->NodeSettleDate  = (long *)  DR_Array(LONG,   -1, NbTP+1);
    tree_data->EqVol           = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->SpotEqVol       = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    
    if (tree_data->FwdEq           == NULL  || 
        tree_data->EqMidNode       == NULL  || 
        tree_data->NodeSettleTime  == NULL  || 
        tree_data->NodeSettleDate  == NULL  || 
        tree_data->EqVol           == NULL  || 
        tree_data->SpotEqVol       == NULL)
    {
        DR_Error("Unable to allocate memory for eq data in tree struct.\n"); 
        goto RETURN;
    }

    /* FX */
    tree_data->FwdFx      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->FxMidNode  = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->FxVol      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->SpotFxVol  = (double *)DR_Array(DOUBLE, -1, NbTP+1);

    if (tree_data->FwdFx      == NULL   ||
        tree_data->FxMidNode  == NULL   ||
        tree_data->FxVol      == NULL   ||
        tree_data->SpotFxVol  == NULL)
    {
        DR_Error("Unable to allocate memory for FX data in tree struct.\n"); 
        goto RETURN;
    }

    /* Smile (Equity or FX) */
    tree_data->A1         = (double *)DR_Array(DOUBLE, -1, MAXNBDATE - 1);
    tree_data->A2         = (double *)DR_Array(DOUBLE, -1, MAXNBDATE - 1);
    tree_data->A3         = (double *)DR_Array(DOUBLE, -1, MAXNBDATE - 1);
    tree_data->A1C        = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->A2C        = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->A3C        = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    tree_data->tMin       = (long   *)DR_Array(LONG,   -1, NbTP+1);
    tree_data->tMax       = (long   *)DR_Array(LONG,   -1, NbTP+1);
    tree_data->SmileIndex = (int    *)DR_Array(DOUBLE, -1, NbTP+1);

    if (tree_data->A1         == NULL   ||
        tree_data->A2         == NULL   ||
        tree_data->A3         == NULL   || 
        tree_data->A1C        == NULL   ||
        tree_data->A2C        == NULL   ||
        tree_data->A3C        == NULL   || 
        tree_data->tMin       == NULL   ||
        tree_data->tMax       == NULL   ||
        tree_data->SmileIndex == NULL)
    {
        DR_Error("Unable to allocate memory for smile data in tree struct.\n"); 
        goto RETURN;
    }


    /* All is well up to here, we are done */
    status = SUCCESS;


  RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Hyb3_Tree_Alloc: could not allocate memory for the tree!\n");
    }

    return(status);

}  /* Hyb3_Tree_Alloc */


/*****  Hyb3_Tree_Free  ***********************************************************/
/**
*         Free the arrays constituting the tree.
*/
int     Hyb3_Tree_Free(HYB3_TREE_DATA    *tree_data)   /**< Tree building data structure */
{
    int   i;
    int   j;
    int   t;
    int   NbTP;  /* Total number of time points in the tree   */
    int   nCache; /* number of cached FX smile points */
                
    NbTP = tree_data->NbTP;

    /* First of all, free slices living in HYB3_TREE_DATA */
    /* NB: Assumes that CET always runs 1-D          */
    if (tree_data->EDevStPrice != NULL)
    {
        if(tree_data->TreeType == TTYPE_1IR2F) /*CET in 2 factor mode */
        {
        for (i = 0; i < tree_data->NbEDevDates; i++)
        {
                Hyb3_Free_Slice(tree_data->EDevStPrice[i], tree_data, 2);
            }
        }
        else
        {
            for (i = 0; i < tree_data->NbEDevDates; i++)
            {
            Hyb3_Free_Slice(tree_data->EDevStPrice[i], tree_data, 1);
        }
    }
    }
    Free_DR_Array(tree_data->EDevStPrice,DOUBLE_PTR,0,tree_data->NbEDevDates-1);
    Free_DR_Array(tree_data->EDevDate, LONG, 0, tree_data->NbEDevDates - 1);

    
    Free_DR_Array (tree_data->SpotVol[0], DOUBLE, -1, NbTP+1);
    
    Free_DR_Array (tree_data->SpotVol[1], DOUBLE, -1, NbTP+1);

    for (i = 0; i < 2; i++)
    {
        for (j =0; j<3;j++)
        {
            Free_DR_Array (tree_data->IrAweight[i][j], DOUBLE, -1, NbTP+1); 
        }
    } 
    
    Free_DR_Array (tree_data->IrZCenter[0], DOUBLE, -1, NbTP+1);
    
    Free_DR_Array (tree_data->IrZCenter[1], DOUBLE, -1, NbTP+1);
    
    Free_DR_Array (tree_data->DriftCUPS[0], DOUBLE, -1, NbTP+1);

    Free_DR_Array (tree_data->DriftCUPS[1], DOUBLE, -1, NbTP+1);

   
    Free_DR_Array (tree_data->Length, DOUBLE, -1, NbTP+1);
    
    Free_DR_Array (tree_data->LengthJ, DOUBLE, -1, NbTP+1);
   
    Free_DR_Array (tree_data->TPDate, LONG, -1, NbTP+1);
    
    for (i = 0; i < 6; i++)
    {
        Free_DR_Array (tree_data->Rho[i], DOUBLE, -1, NbTP+1);

    }  /* for i */

    for (i = 0; i < 10; i++)
    {
        Free_DR_Array (tree_data->Aweight[i], DOUBLE, -1, NbTP+1); 

    }  /* for i */

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 3; j++)
        {
            Free_DR_Array (tree_data->FwdRate[i][j], DOUBLE, -1, NbTP+1);
            Free_DR_Array (tree_data->ZeroRate[i][j], DOUBLE, -1, NbTP+1);
            Free_DR_Array(tree_data->ZeroCoupon[i][j], DOUBLE, -1, NbTP+1);
            
        } /* For j */

    }  /* for i */

    for (i = 0; i < NBCRITDATE; i++)
    {

        Free_DR_Array (tree_data->CritDate[i], CRITDATE, -1, NbTP+1);
        Free_DR_Array (tree_data->TPType[i], INT, -1, NbTP+1);
        

    }  /* for i */

    /*  Free the limits of the ellipse. */
    if (    (tree_data->TreeType != TTYPE_2IR) 
         && (tree_data->TreeType != TTYPE_EQ1IR)
         && (tree_data->TreeType != TTYPE_EQDFX2IR)
         && (tree_data->TreeType != TTYPE_EQFFX2IR) /* i.e. 3rd dimenstion ON, 4th OFF */
         && (tree_data->Top1 != NULL) && (tree_data->Bottom1 != NULL)
         && (tree_data->Top2 != NULL) && (tree_data->Bottom2 != NULL)
         && (tree_data->Top3 != NULL) && (tree_data->Bottom3 != NULL))
    {
        for (t = 0; t <= NbTP+1; t++)
        {    
            if((tree_data->Top2[t] != NULL) && (tree_data->Bottom2[t] != NULL)
            && (tree_data->Top3[t] != NULL) && (tree_data->Bottom3[t] != NULL))
            {
                for (i = tree_data->Bottom1[t]; i <= tree_data->Top1[t]; i++)
                {
                    if (tree_data->Top3[t][i] != NULL)
                    {
                        Free_DR_Array (tree_data->Top3[t][i], INT, 
                                       tree_data->Bottom2[t][i], 
                                       tree_data->Top2[t][i]);
                    }

                    if (tree_data->Bottom3[t][i] != NULL)
                    {
                        Free_DR_Array (tree_data->Bottom3[t][i], INT,
                                       tree_data->Bottom2[t][i], 
                                       tree_data->Top2[t][i]);
                    }

                }  
            }  

            if (tree_data->Top3[t] != NULL)
            {
                Free_DR_Array (tree_data->Top3[t], INT_PTR,
                               tree_data->Bottom1[t], tree_data->Top1[t]);
            }

            if (tree_data->Bottom3[t] != NULL)
            {
                 Free_DR_Array (tree_data->Bottom3[t], INT_PTR,
                                    tree_data->Bottom1[t], tree_data->Top1[t]);
            }

        }  /* for t */

    }  /* if */
    else if ((
                (tree_data->TreeType == TTYPE_EQDFX2IR) ||
                (tree_data->TreeType == TTYPE_EQFFX2IR)
             ) &&
             (tree_data->Top1    != NULL) &&
             (tree_data->Bottom1 != NULL) &&
             (tree_data->Top2    != NULL) &&
             (tree_data->Bottom2 != NULL) &&
             (tree_data->Top3    != NULL) &&
             (tree_data->Bottom3 != NULL) &&
             (tree_data->Top4    != NULL) &&
             (tree_data->Bottom4 != NULL))
    {
        for (t = 0; t <= NbTP+1; t++)
        {    
            if ((tree_data->Top2    [t] != NULL) &&
                (tree_data->Bottom2 [t] != NULL) &&
                (tree_data->Top3    [t] != NULL) &&
                (tree_data->Bottom3 [t] != NULL))
            {
                for (i = tree_data->Bottom1[t]; i <= tree_data->Top1[t]; i++)
                {
                    if ((tree_data->Top3    [t][i] != NULL) &&
                        (tree_data->Bottom3 [t][i] != NULL))
                    {
                        if ((tree_data->Top4 [t]    != NULL) &&
                            (tree_data->Top4 [t][i] != NULL))
                        {
                            for (j = tree_data->Bottom2[t][i]; j <= tree_data->Top2[t][i]; j++)
                            {
                                if (tree_data->Top4[t][i][j] != NULL)
                                {
                                    Free_DR_Array (tree_data->Top4    [t][i][j], INT, 
                                                   tree_data->Bottom3 [t][i][j], 
                                                   tree_data->Top3    [t][i][j]);
                                }

                                if (tree_data->Bottom4[t][i][j] != NULL)
                                {
                                    Free_DR_Array (tree_data->Bottom4 [t][i][j], INT, 
                                                   tree_data->Bottom3 [t][i][j], 
                                                   tree_data->Top3    [t][i][j]);
                                }
                            }

                            Free_DR_Array (tree_data->Top4    [t][i], INT_PTR, 
                                           tree_data->Bottom2 [t][i], 
                                           tree_data->Top2    [t][i]);

                            Free_DR_Array (tree_data->Bottom4 [t][i], INT_PTR, 
                                           tree_data->Bottom2 [t][i], 
                                           tree_data->Top2    [t][i]);
                        }

                        Free_DR_Array (tree_data->Top3    [t][i], INT, 
                                       tree_data->Bottom2 [t][i], 
                                       tree_data->Top2    [t][i]);

                        Free_DR_Array (tree_data->Bottom3 [t][i], INT, 
                                       tree_data->Bottom2 [t][i], 
                                       tree_data->Top2    [t][i]);
                    }
                }  

                Free_DR_Array (tree_data->Top3    [t], INT_PTR, 
                               tree_data->Bottom1 [t], 
                               tree_data->Top1    [t]);

                Free_DR_Array (tree_data->Bottom3 [t], INT_PTR, 
                               tree_data->Bottom1 [t], 
                               tree_data->Top1    [t]);

                Free_DR_Array (tree_data->Top4    [t], INT_D_PTR, 
                               tree_data->Bottom1 [t], 
                               tree_data->Top1    [t]);

                Free_DR_Array (tree_data->Bottom4 [t], INT_D_PTR, 
                               tree_data->Bottom1 [t], 
                               tree_data->Top1    [t]);
            }  
        }
    }

    if ((tree_data->Top1 != NULL) && (tree_data->Bottom1 != NULL) &&
        (tree_data->Top2 != NULL) && (tree_data->Bottom2 != NULL))
    {
        for (t = 0; t <= NbTP+1; t++)
        {    
            if (tree_data->Top2[t] != NULL)
            {
                Free_DR_Array(tree_data->Top2[t], INT,
                              tree_data->Bottom1[t], tree_data->Top1[t]);
            }

            if (tree_data->Bottom2[t] != NULL)
            {
                 Free_DR_Array(tree_data->Bottom2[t], INT,
                               tree_data->Bottom1[t], tree_data->Top1[t]);
            }

        }  /* for t */

    }  /* if */

    if (tree_data->Top1 != NULL)
    {
        Free_DR_Array(tree_data->Top1, INT, -1, NbTP+1);
    }

    if (tree_data->Bottom1 != NULL)
    {
        Free_DR_Array(tree_data->Bottom1, INT, -1, NbTP+1);
    }

    if (tree_data->Top2 != NULL)
    {
        Free_DR_Array(tree_data->Top2, INT_PTR, -1, NbTP+1);
    }

    if (tree_data->Bottom2 != NULL)
    {
        Free_DR_Array(tree_data->Bottom2, INT_PTR, -1, NbTP+1);
    }

    if (tree_data->Top3 != NULL)
    {
        Free_DR_Array(tree_data->Top3, INT_D_PTR, -1, NbTP+1);
    }

    if (tree_data->Bottom3 != NULL)
    {
        Free_DR_Array(tree_data->Bottom3, INT_D_PTR, -1, NbTP+1);
    }

    if (tree_data->Top4 != NULL)
    {
        Free_DR_Array(tree_data->Top4, INT_T_PTR, -1, NbTP+1);
    }

    if (tree_data->Bottom4 != NULL)
    {
        Free_DR_Array(tree_data->Bottom4, INT_T_PTR, -1, NbTP+1);
    }

    if (    (tree_data->TreeType != TTYPE_2IR) 
         && (tree_data->TreeType != TTYPE_EQ1IR)
         && (tree_data->TreeType != TTYPE_EQDFX2IR)
         && (tree_data->TreeType != TTYPE_EQFFX2IR) /* i.e. 3rd dimenstion ON, 4th OFF */
         && (tree_data->OutTop1 != NULL) && (tree_data->OutBottom1 != NULL)
         && (tree_data->OutTop2 != NULL) && (tree_data->OutBottom2 != NULL)
         && (tree_data->OutTop3 != NULL) && (tree_data->OutBottom3 != NULL))
    {
        for (t = 0; t <= NbTP+1; t++)
        {    
            if((tree_data->OutTop2[t] != NULL) && (tree_data->OutBottom2[t] != NULL)
            && (tree_data->OutTop3[t] != NULL) && (tree_data->OutBottom3[t] != NULL))
            {

                for (i=tree_data->OutBottom1[t]; i<=tree_data->OutTop1[t]; i++)
                {
                    if (tree_data->OutTop3[t][i] != NULL)
                    {
                        Free_DR_Array(tree_data->OutTop3[t][i], INT,
                                      tree_data->OutBottom2[t][i], 
                                      tree_data->OutTop2[t][i]);
                    }

                    if (tree_data->OutBottom3[t][i] != NULL)
                    {
                        Free_DR_Array(tree_data->OutBottom3[t][i], INT,
                                      tree_data->OutBottom2[t][i], 
                                      tree_data->OutTop2[t][i]);
                    }

                }  /* for i */
            }  /* if */

            if (tree_data->OutTop3[t] != NULL)
            {
                Free_DR_Array(tree_data->OutTop3[t], INT_PTR,
                              tree_data->OutBottom1[t], 
                              tree_data->OutTop1[t]);
            }

            if (tree_data->OutBottom3[t] != NULL)
            {
                Free_DR_Array(tree_data->OutBottom3[t], INT_PTR,
                              tree_data->OutBottom1[t], tree_data->OutTop1[t]);
            }

        }  /* for t */

    }  /* if */

    if ((tree_data->OutTop1 != NULL) && (tree_data->OutBottom1 != NULL) &&
        (tree_data->OutTop2 != NULL) && (tree_data->OutBottom2 != NULL))
    {
        for (t = 0; t <= NbTP+1; t++)
        {    
            if (tree_data->OutTop2[t] != NULL)
            {
                Free_DR_Array(tree_data->OutTop2[t], INT,
                              tree_data->OutBottom1[t], tree_data->OutTop1[t]);
            }

            if (tree_data->OutBottom2[t] != NULL)
            {
                Free_DR_Array(tree_data->OutBottom2[t], INT,
                              tree_data->OutBottom1[t], tree_data->OutTop1[t]);
            }

        }  /* for t */
    }  /* if */

    if (tree_data->OutTop1 != NULL)
        Free_DR_Array(tree_data->OutTop1, INT, -1, NbTP+1);

    if (tree_data->OutBottom1 != NULL)
        Free_DR_Array(tree_data->OutBottom1, INT, -1, NbTP+1);

    if (tree_data->OutTop2 != NULL)
        Free_DR_Array(tree_data->OutTop2, INT_PTR, -1, NbTP+1);

    if (tree_data->OutBottom2 != NULL)
        Free_DR_Array(tree_data->OutBottom2, INT_PTR, -1, NbTP+1);

    if (tree_data->OutTop3 != NULL)
        Free_DR_Array(tree_data->OutTop3, INT_D_PTR, -1, NbTP+1);

    if (tree_data->OutBottom3 != NULL)
        Free_DR_Array(tree_data->OutBottom3, INT_D_PTR, -1, NbTP+1);

    if (tree_data->OutTop4 != NULL)
        Free_DR_Array(tree_data->OutTop4, INT_T_PTR, -1, NbTP+1);

    if (tree_data->OutBottom4 != NULL)
        Free_DR_Array(tree_data->OutBottom4, INT_T_PTR, -1, NbTP+1);

    if (    (tree_data->TreeType != TTYPE_2IR) 
             && (tree_data->TreeType != TTYPE_EQ1IR)
             && (tree_data->TreeType != TTYPE_EQDFX2IR)
             && (tree_data->TreeType != TTYPE_EQFFX2IR) /* i.e. 3rd dimenstion ON, 4th OFF */
             && (tree_data->InnerTreeTop1 != NULL) && (tree_data->InnerTreeBottom1 != NULL)
             && (tree_data->InnerTreeTop2 != NULL) && (tree_data->InnerTreeBottom2 != NULL)
             && (tree_data->InnerTreeTop3 != NULL) && (tree_data->InnerTreeBottom3 != NULL))
        {
            for (t = 0; t <= NbTP+1; t++)
            {    
                if((tree_data->InnerTreeTop2[t] != NULL) && (tree_data->InnerTreeBottom2[t] != NULL)
                && (tree_data->InnerTreeTop3[t] != NULL) && (tree_data->InnerTreeBottom3[t] != NULL))
                {
    
                    for (i=tree_data->InnerTreeBottom1[t]; i<=tree_data->InnerTreeTop1[t]; i++)
                    {
                        if (tree_data->InnerTreeTop3[t][i] != NULL)
                        {
                            Free_DR_Array(tree_data->InnerTreeTop3[t][i], INT,
                                          tree_data->InnerTreeBottom2[t][i], 
                                          tree_data->InnerTreeTop2[t][i]);
                        }
    
                        if (tree_data->InnerTreeBottom3[t][i] != NULL)
                        {
                            Free_DR_Array(tree_data->InnerTreeBottom3[t][i], INT,
                                          tree_data->InnerTreeBottom2[t][i], 
                                          tree_data->InnerTreeTop2[t][i]);
                        }
    
                    }  /* for i */
                }  /* if */
    
                if (tree_data->InnerTreeTop3[t] != NULL)
                {
                    Free_DR_Array(tree_data->InnerTreeTop3[t], INT_PTR,
                                  tree_data->InnerTreeBottom1[t], 
                                  tree_data->InnerTreeTop1[t]);
                }
    
                if (tree_data->InnerTreeBottom3[t] != NULL)
                {
                    Free_DR_Array(tree_data->InnerTreeBottom3[t], INT_PTR,
                                  tree_data->InnerTreeBottom1[t], tree_data->InnerTreeTop1[t]);
                }
    
            }  /* for t */
    
        }  /* if */

        if ((tree_data->InnerTreeTop1 != NULL) && (tree_data->InnerTreeBottom1 != NULL) &&
            (tree_data->InnerTreeTop2 != NULL) && (tree_data->InnerTreeBottom2 != NULL))
        {
            for (t = 0; t <= NbTP+1; t++)
            {    
                if (tree_data->InnerTreeTop2[t] != NULL)
                {
                    Free_DR_Array(tree_data->InnerTreeTop2[t], INT,
                                  tree_data->InnerTreeBottom1[t], tree_data->InnerTreeTop1[t]);
                }
    
                if (tree_data->InnerTreeBottom2[t] != NULL)
                {
                    Free_DR_Array(tree_data->InnerTreeBottom2[t], INT,
                                  tree_data->InnerTreeBottom1[t], tree_data->InnerTreeTop1[t]);
                }
    
            }  /* for t */
        }  /* if */
    
    if (tree_data->InnerTreeTop1 != NULL)
        Free_DR_Array(tree_data->InnerTreeTop1, INT, -1, NbTP+1);
    
    if (tree_data->InnerTreeBottom1 != NULL)
        Free_DR_Array(tree_data->InnerTreeBottom1, INT, -1, NbTP+1);
    
    if (tree_data->InnerTreeTop2 != NULL)
        Free_DR_Array(tree_data->InnerTreeTop2, INT_PTR, -1, NbTP+1);
    
    if (tree_data->InnerTreeBottom2 != NULL)
        Free_DR_Array(tree_data->InnerTreeBottom2, INT_PTR, -1, NbTP+1);
    
    if (tree_data->InnerTreeTop3 != NULL)
        Free_DR_Array(tree_data->InnerTreeTop3, INT_D_PTR, -1, NbTP+1);
    
    if (tree_data->InnerTreeBottom3 != NULL)
        Free_DR_Array(tree_data->InnerTreeBottom3, INT_D_PTR, -1, NbTP+1);

    if (tree_data->InnerTreeTop4 != NULL)
        Free_DR_Array(tree_data->InnerTreeTop4, INT_T_PTR, -1, NbTP+1);
    
    if (tree_data->InnerTreeBottom4 != NULL)
        Free_DR_Array(tree_data->InnerTreeBottom4, INT_T_PTR, -1, NbTP+1);

    /* Memory related to the third asset */
    if (tree_data->FwdFx != NULL)
        Free_DR_Array(tree_data->FwdFx, DOUBLE, -1, NbTP+1);

    if (tree_data->FxMidNode != NULL)
        Free_DR_Array(tree_data->FxMidNode, DOUBLE, -1, NbTP+1);

    if (tree_data->SpotFxVol != NULL)
        Free_DR_Array(tree_data->SpotFxVol, DOUBLE, -1, NbTP+1);

    if (tree_data->FxVol != NULL)
        Free_DR_Array(tree_data->FxVol, DOUBLE, -1, NbTP+1);

    if (tree_data->SmileIndex != NULL)
        Free_DR_Array(tree_data->SmileIndex, DOUBLE, -1, NbTP+1);

    if (tree_data->A1 != NULL)
        Free_DR_Array(tree_data->A1, DOUBLE, -1, MAXNBDATE-1);

    if (tree_data->A2 != NULL)
        Free_DR_Array(tree_data->A2, DOUBLE, -1, MAXNBDATE-1);

    if (tree_data->A3 != NULL)
        Free_DR_Array(tree_data->A3, DOUBLE, -1, MAXNBDATE-1);

    if (tree_data->A1C != NULL)
        Free_DR_Array(tree_data->A1C, DOUBLE, -1, NbTP+1);

    if (tree_data->A2C != NULL)
        Free_DR_Array(tree_data->A2C, DOUBLE, -1, NbTP+1);

    if (tree_data->A3C != NULL)
        Free_DR_Array(tree_data->A3C, DOUBLE, -1, NbTP+1);

    if (tree_data->tMin != NULL)
        Free_DR_Array(tree_data->tMin, LONG, -1, NbTP+1);

    if (tree_data->tMax != NULL)
        Free_DR_Array(tree_data->tMax, LONG, -1, NbTP+1);

    /* FX smile caching memory: release as well */
    if (tree_data->FXsmileCache.isCached)
    {
        nCache = tree_data->FXsmileCache.nbCachePts;
        if (tree_data->FXsmileCache.X != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.X, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.K != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.K,DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.SPL != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.SPL, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.SPL_Inv != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.SPL_Inv, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.nbPtSPL != NULL)
            Free_DR_Array (tree_data->FXsmileCache.nbPtSPL, INT, 0, nCache);
        
        if (tree_data->FXsmileCache.gd != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.gd, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.gd_SPL != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.gd_SPL, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.kdX != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.kdX, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
        
        if (tree_data->FXsmileCache.kdX_SPL != NULL)
            Free_DR_Matrix (tree_data->FXsmileCache.kdX_SPL, DOUBLE, 0, nCache, 0 , NBSPLINE + 3);
    }

    /* other information */
    if (tree_data->FwdEq != NULL)
        Free_DR_Array(tree_data->FwdEq, DOUBLE, -1, NbTP+1);

    if (tree_data->EqMidNode != NULL)
        Free_DR_Array(tree_data->EqMidNode, DOUBLE, -1, NbTP+1);

    if (tree_data->NodeSettleTime != NULL)
        Free_DR_Array(tree_data->NodeSettleTime, DOUBLE, -1, NbTP+1);

    if (tree_data->SpotEqVol != NULL)
        Free_DR_Array(tree_data->SpotEqVol, DOUBLE, -1, NbTP+1);

    if (tree_data->EqVol != NULL)
        Free_DR_Array(tree_data->EqVol, DOUBLE, -1, NbTP+1);

    if (tree_data->NodeSettleDate != NULL)
        Free_DR_Array(tree_data->NodeSettleDate, LONG, -1, NbTP+1);

    Hyb3_Tree_Init(tree_data);
    return(SUCCESS);

}  /* Hyb3_Tree_Free */


/****** Hyb3_Dev_Init *************************************************/
/**
  
        Initialisation of memory for DEV structure
  
  
 ******************************************************************/


void    Hyb3_Dev_Init (HYB3_DEV_DATA  *dev_data)
{
    int i;

    for (i = 0; i < 3; i++)
    {
        dev_data->Discount_1D[i] = NULL;
        dev_data->Discount_2D[i] = NULL;
        dev_data->Discount_3D[i] = NULL;
    }


    dev_data->Aux1D = NULL;
    dev_data->Aux2D = NULL;
    dev_data->Aux3D = NULL;

    dev_data->DomZero = NULL;
    dev_data->FwdFX   = NULL;


    dev_data->Shift1 = NULL;
    dev_data->Shift2 = NULL;
    dev_data->Shift3 = NULL;
    dev_data->Shift4 = NULL;
    dev_data->Shift5 = NULL;

    dev_data->p = NULL;
    dev_data->s = NULL;
    dev_data->q = NULL;
    dev_data->t = NULL;
    dev_data->r = NULL;

    dev_data->quu = NULL;
    dev_data->qu0 = NULL;
    dev_data->qud = NULL;
    dev_data->q0u = NULL;
    dev_data->q00 = NULL;
    dev_data->q0d = NULL;
    dev_data->qdu = NULL;
    dev_data->qd0 = NULL;
    dev_data->qdd = NULL;

    dev_data->EqSpot = NULL;
    dev_data->NextEqSpot = NULL;
    dev_data->EqFwd = NULL;
    dev_data->FxSpot = NULL;
    dev_data->NextFxSpot = NULL;

    dev_data->gDash       = NULL;
    dev_data->kDashTimesX = NULL;
    dev_data->kVar        = NULL;

    dev_data->t2 = NULL;
    dev_data->t3 = NULL;

    return;
}/* Hyb3_Dev_Init */



/*****  Hyb3_Dev_Alloc  ***********************************************************/
/**
*         Allocation of memory for the dev structure.
*         We use rectangular matrices, not elliptic and we add one extra node
*         to be consistent with branching.
*/
int    Hyb3_Dev_Alloc(HYB3_DEV_DATA     *dev_data,     /**< (I) Hyb3_Dev data structure    */
                 HYB3_TREE_DATA    *tree_data)    /**< (I) Tree data structure   */
{
    
    
    
    int    Size1;    /* Number of positions needed to describe 1-D slice */
    int    Size2;    /* Number of positions needed to describe 2-D slice */
    int    Size3;    /* Number of positions needed to describe 3-D slice */


    int status = FAILURE;


    Size1 = tree_data->Width[0];
    if (tree_data->TreeType != TTYPE_1IR)
    {
        Size2 = tree_data->Width[0] * tree_data->Width[1];
    }
    else
    {
        Size2 = 0;
    }
   
    dev_data->Discount_1D[0] = (TSLICE)DR_Array(DOUBLE, 0, Size1);
    dev_data->Discount_1D[1] = (TSLICE)DR_Array(DOUBLE, 0, Size1);
    dev_data->Discount_1D[2] = (TSLICE)DR_Array(DOUBLE, 0, Size1);
    if (dev_data->Discount_1D[0] == NULL ||
        dev_data->Discount_1D[1] == NULL ||
        dev_data->Discount_1D[2] == NULL )
    {
        goto RETURN;
    }

    dev_data->Discount_2D[0] = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->Discount_2D[1] = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->Discount_2D[2] = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    if (dev_data->Discount_2D[0] == NULL ||
        dev_data->Discount_2D[1] == NULL ||
        dev_data->Discount_2D[2] == NULL )
    {
        goto RETURN;
    }

    dev_data->Aux1D   = (TSLICE)DR_Array(DOUBLE, 0, Size1);
    dev_data->Aux2D   = (TSLICE)DR_Array(DOUBLE, 0, Size2);

    if (dev_data->Aux1D == NULL ||
        dev_data->Aux2D == NULL )
    {
        goto RETURN;
    }

    dev_data->p = (TPROB_0 *)DR_Array(TYPE_TPROB_0, 0, Size1);
    if (dev_data->p == NULL)
    {
        goto RETURN;
    }

    dev_data->s = (TPROB_0 *)DR_Array(TYPE_TPROB_0, 0, Size1);
    if (dev_data->s == NULL)
    {
        goto RETURN;
    }

    dev_data->q = (TPROB_0 *)DR_Array(TYPE_TPROB_0, 0, Size2);
    if (dev_data->q == NULL)
    {
        goto RETURN;
    }
    
    dev_data->t = (TPROB_0 *)DR_Array(TYPE_TPROB_0, 0, Size2);
    if (dev_data->t == NULL)
    {
        goto RETURN;
    }

    dev_data->t2 = (TSLICE)DR_Array(DOUBLE, 0, tree_data->Width[1]);

    dev_data->Shift1 = (int *)DR_Array(INT, 0, Size1);
    dev_data->Shift4 = (int *)DR_Array(INT, 0, Size1);
    dev_data->Shift2 = (int *)DR_Array(INT, 0, Size2);
    dev_data->Shift5 = (int *)DR_Array(INT, 0, Size2);
    if (dev_data->Shift1 == NULL ||
        dev_data->Shift4 == NULL ||
        dev_data->Shift2 == NULL ||
        dev_data->Shift5 == NULL)
    {
        goto RETURN;
    }

    dev_data->quu = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->qu0 = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->qud = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->q0u = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->q00 = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->q0d = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->qdu = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->qd0 = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    dev_data->qdd = (TSLICE)DR_Array(DOUBLE, 0, Size2);
    if (dev_data->quu == NULL ||
        dev_data->qu0 == NULL ||
        dev_data->qud == NULL ||
        dev_data->q0u == NULL ||
        dev_data->q00 == NULL ||
        dev_data->q0d == NULL ||
        dev_data->qdu == NULL ||
        dev_data->qd0 == NULL ||
        dev_data->qdd == NULL )
    {
        goto RETURN;
    }

    if (tree_data->TreeType == TTYPE_EQ1IR)
    {
        dev_data->EqSpot     = (TSLICE)DR_Array(DOUBLE, 0, Size2);
        dev_data->NextEqSpot = (TSLICE)DR_Array(DOUBLE, 0, Size2);
        dev_data->EqFwd      = (TSLICE)DR_Array(DOUBLE, 0, Size2);

        dev_data->gDash       = (TSLICE)DR_Array(DOUBLE, 0, Size2);
        dev_data->kDashTimesX = (TSLICE)DR_Array(DOUBLE, 0, Size2);
        dev_data->kVar        = (TSLICE)DR_Array(DOUBLE, 0, Size2);

        if (dev_data->EqSpot     == NULL ||
            dev_data->NextEqSpot == NULL ||
            dev_data->EqFwd      == NULL ||
            dev_data->gDash      == NULL ||
            dev_data->kDashTimesX== NULL ||
            dev_data->kVar       == NULL )
        {
            goto RETURN;
        }

    }
    else if ((tree_data->TreeType != TTYPE_2IR) && /* i.e. must be 3-D */
             (tree_data->TreeType != TTYPE_1IR)&&
             (tree_data->TreeType != TTYPE_1IR2F))
    {
        Size3 = tree_data->Width[0] *tree_data->Width[1] *tree_data->Width[2];
         
        dev_data->Shift3 = (int *)DR_Array(INT, 0, Size3);
        dev_data->Aux3D = (TSLICE)DR_Array(DOUBLE, 0, Size3);
        dev_data->r   =   (TPROB_0 *)DR_Array(TYPE_TPROB_0, 0, Size3);
        dev_data->t3  =   (TSLICE)DR_Array(DOUBLE, 0, tree_data->Width[2]);

        dev_data->gDash       = (TSLICE)DR_Array(DOUBLE, 0, Size3);
        dev_data->kDashTimesX = (TSLICE)DR_Array(DOUBLE, 0, Size3);
        dev_data->kVar        = (TSLICE)DR_Array(DOUBLE, 0, Size3);

        if (dev_data->Shift3 == NULL ||
            dev_data->Aux3D  == NULL ||
            dev_data->r      == NULL ||
            dev_data->t3     == NULL ||
            dev_data->gDash      == NULL ||
            dev_data->kDashTimesX== NULL ||
            dev_data->kVar       == NULL )
        {
            goto RETURN;
        }

        if (tree_data->TreeType == TTYPE_FX2IR)
        {
            dev_data->FxSpot = (TSLICE)DR_Array(DOUBLE, 0, Size3);
            dev_data->NextFxSpot = (TSLICE)DR_Array(DOUBLE,0,Size3);
            dev_data->FwdFX      = (TSLICE)DR_Array(DOUBLE,0,Size3);
            dev_data->DomZero    = (TSLICE)DR_Array(DOUBLE,0,Size2);

            if (dev_data->FxSpot     == NULL || 
                dev_data->NextFxSpot == NULL ||
                dev_data->FwdFX      == NULL ||
                dev_data->DomZero    == NULL )
            {
                goto RETURN;
            }

        }
        else if (tree_data->TreeType == TTYPE_2IR2F1D)
        {
                dev_data->Discount_3D[0] = (TSLICE)DR_Array(DOUBLE, 0, Size3);
                dev_data->Discount_3D[1] = (TSLICE)DR_Array(DOUBLE, 0, Size3);
                dev_data->Discount_3D[2] = (TSLICE)DR_Array(DOUBLE, 0, Size3);

                if (dev_data->Discount_3D[0] == NULL ||
                    dev_data->Discount_3D[1] == NULL ||
                    dev_data->Discount_3D[2] == NULL )
                {
                    goto RETURN;
                }



        }
        else /* it must be one of the equity models */
        {
            dev_data->EqSpot     = (TSLICE)DR_Array(DOUBLE, 0, Size3);
            dev_data->NextEqSpot = (TSLICE)DR_Array(DOUBLE, 0, Size3);
            dev_data->EqFwd      = (TSLICE)DR_Array(DOUBLE, 0, Size3);
            if (dev_data->EqSpot     == NULL ||
                dev_data->NextEqSpot == NULL ||
                dev_data->EqFwd      == NULL )
            {
                goto RETURN;
            }

        }  /* if then else */

    }

    status = SUCCESS;

  RETURN:
    
    if (status == FAILURE)
    {
        DR_Error("Hyb3_Dev_Alloc: Unable to allocate memory for DEV structure.\n");
    }
    return(status);

}  /* Hyb3_Dev_Alloc */


/*****  Hyb3_Dev_Free  ************************************************************/
/**
 *         Free memory for the dev structure.
 */
int    Hyb3_Dev_Free(HYB3_DEV_DATA     *dev_data,           /**< (I) Hyb3_Dev data structure  */
                HYB3_TREE_DATA    *tree_data)          /**< (I) Tree data structure */
{
    int    Size1;    /* Number of positions needed to describe 1-D slice */
    int    Size2;    /* Number of positions needed to describe 2-D slice */
    int    Size3;    /* Number of positions needed to describe 3-D slice */

    Size1 = tree_data->Width[0];
    if (tree_data->TreeType != TTYPE_1IR)
    {
        Size2 = tree_data->Width[0] * tree_data->Width[1];
    }
    else
    {
        Size2 = 0;
    }

    Free_DR_Array(dev_data->Discount_1D[0], DOUBLE, 0, Size1);
    Free_DR_Array(dev_data->Discount_1D[1], DOUBLE, 0, Size1);
    Free_DR_Array(dev_data->Discount_1D[2], DOUBLE, 0, Size1);

    Free_DR_Array(dev_data->Discount_2D[0], DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->Discount_2D[1], DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->Discount_2D[2], DOUBLE, 0, Size2);

    Free_DR_Array(dev_data->Shift1, INT, 0, Size1);
    Free_DR_Array(dev_data->Shift4, INT, 0, Size1);
    Free_DR_Array(dev_data->Shift2, INT, 0, Size2);
    Free_DR_Array(dev_data->Shift5, INT, 0, Size2);

    Free_DR_Array(dev_data->Aux1D, DOUBLE, 0, Size1);
    Free_DR_Array(dev_data->Aux2D, DOUBLE, 0, Size2);

    Free_DR_Array(dev_data->p, TYPE_TPROB_0, 0, Size1);
    Free_DR_Array(dev_data->s, TYPE_TPROB_0, 0, Size1);
    Free_DR_Array(dev_data->q, TYPE_TPROB_0, 0, Size2);
    Free_DR_Array(dev_data->t, TYPE_TPROB_0, 0, Size2);
          
    Free_DR_Array(dev_data->t2, DOUBLE, 0, tree_data->Width[1]);

    Free_DR_Array(dev_data->quu, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->qu0, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->qud, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->q0u, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->q00, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->q0d, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->qdu, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->qd0, DOUBLE, 0, Size2);
    Free_DR_Array(dev_data->qdd, DOUBLE, 0, Size2);

    if (tree_data->TreeType == TTYPE_EQ1IR)
    {
        Free_DR_Array(dev_data->EqSpot,     DOUBLE, 0, Size2);
        Free_DR_Array(dev_data->NextEqSpot, DOUBLE, 0, Size2);
        Free_DR_Array(dev_data->EqFwd,      DOUBLE, 0, Size2);

        Free_DR_Array(dev_data->gDash,      DOUBLE, 0, Size2);
        Free_DR_Array(dev_data->kDashTimesX,DOUBLE, 0, Size2);
        Free_DR_Array(dev_data->kVar,       DOUBLE, 0, Size2);
    }
    else if ((tree_data->TreeType!=TTYPE_2IR) && /* i.e. running in 3-D mode */
             (tree_data->TreeType!=TTYPE_1IR) &&
             (tree_data->TreeType != TTYPE_1IR2F))
    {
        Size3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];

        Free_DR_Array(dev_data->Discount_3D[0], DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->Discount_3D[1], DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->Discount_3D[2], DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->Shift3, INT, 0, Size3);
        Free_DR_Array(dev_data->Aux3D, DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->r,     TYPE_TPROB_0, 0, Size3);
        Free_DR_Array(dev_data->t3,    DOUBLE, 0, tree_data->Width[2]);

        Free_DR_Array(dev_data->gDash,      DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->kDashTimesX,DOUBLE, 0, Size3);
        Free_DR_Array(dev_data->kVar,       DOUBLE, 0, Size3);

        if (tree_data->TreeType == TTYPE_FX2IR)
        {
            Free_DR_Array(dev_data->FxSpot, DOUBLE, 0, Size3);
            Free_DR_Array(dev_data->NextFxSpot, DOUBLE,0,Size3);
            Free_DR_Array(dev_data->FwdFX,DOUBLE,0,Size3);
            Free_DR_Array(dev_data->DomZero,DOUBLE,0,Size2);
        }
        else
        {
            Free_DR_Array(dev_data->EqSpot,     DOUBLE, 0, Size3);
            Free_DR_Array(dev_data->NextEqSpot, DOUBLE, 0, Size3);
            Free_DR_Array(dev_data->EqFwd,      DOUBLE, 0, Size3);
        }
    }  /* if then else */

    Hyb3_Dev_Init(dev_data);
    return(SUCCESS);
}  /* Hyb3_Dev_Free */


/*****  Hyb3_Alloc_Slice  ********************************************************/
/**
*         Allocate memory for one time slice.
*/
TSLICE    Hyb3_Alloc_Slice (HYB3_TREE_DATA    *tree_data, /**< (I) Tree data structure  */
                       int           dimension) /**< (I) Dim of desired slice */

{
    TSLICE    Slice;

    int     Size1;    /* Number of positions needed to describe 1-D slice */
    int     Size2;    /* Number of positions needed to describe 2-D slice */
    int     Size3;    /* Number of positions needed to describe 3-D slice */
    
    Slice = NULL;

    switch (dimension)
    {
        case 1:
          Size1 = tree_data->Width[0];          
          Slice = (TSLICE) DR_Array(DOUBLE, 0, Size1);
          break;

        case 2:
          Size2 = tree_data->Width[0] * tree_data->Width[1];
          Slice = (TSLICE) DR_Array(DOUBLE, 0, Size2);
          break;

        case 3:
          Size3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
          Slice = (TSLICE) DR_Array(DOUBLE, 0, Size3);
          break;

        default:

          DR_Error("Invalid slice dimension! (Hyb3_Alloc_Slice)\n");
     }

     return (Slice);
}  /* Hyb3_Alloc_Slice */


/*****  Hyb3_Free_Slice  *********************************************************/
/**
*         Free memory for time slice.
*/
void     Hyb3_Free_Slice(TSLICE       Slice,        /**< (I) Time slice            */
                    HYB3_TREE_DATA   *tree_data,    /**< (I) Tree data structure   */
                    int          dimension)
{
    int     Size1;    /* Number of positions needed to describe 1-D slice */
    int     Size2;    /* Number of positions needed to describe 2-D slice */
    int     Size3;    /* Number of positions needed to describe 3-D slice */

    if (Slice == NULL)
    {
        return;
    }

    switch(dimension)
    {
        case 1:
            Size1 = tree_data->Width[0];
            Free_DR_Array(Slice, DOUBLE, 0, Size1);
            break;
        case 2:
            Size2 = tree_data->Width[0] * tree_data->Width[1];
            Free_DR_Array(Slice, DOUBLE, 0, Size2);
            break;
        case 3:
            Size3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
            Free_DR_Array(Slice, DOUBLE, 0, Size3);
            break;
        default:

            DR_Error("Invalid slice dimension! (Hyb3_Free_Slice)\n");

    } /* switch */

     return;
}  /* Hyb3_Free_Slice */

/*****  Hyb3_Tree_Dim  *********************************************************/
/**
*        Returns the tree dimension for a given tree-mode.
*/
int Hyb3_Tree_Dim(HYB3_TREE_DATA* tree_data)
{

    if (tree_data->TreeType == TTYPE_1IR)     return 1;
    if (tree_data->TreeType == TTYPE_2IR)     return 2;
    if (tree_data->TreeType == TTYPE_EQ1IR)   return 2;
    if (tree_data->TreeType == TTYPE_1IR2F)   return 2;
    if (tree_data->TreeType == TTYPE_FX2IR)   return 3;
    if (tree_data->TreeType == TTYPE_EQD2IR)  return 3;
    if (tree_data->TreeType == TTYPE_EQF2IR)  return 3;
    if (tree_data->TreeType == TTYPE_EQC2IR)  return 3;
    if (tree_data->TreeType == TTYPE_2IR2F1D) return 3;
    if (tree_data->TreeType == TTYPE_EQDFX2IR) return 4;
    if (tree_data->TreeType == TTYPE_EQFFX2IR) return 4;
    
    return -1; /* not a valid tree type */
}  /* Hyb3_Tree_Dim */

