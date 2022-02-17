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

/*****  Fix3_SwapRate_Init *******************************************************/
/**
*         Initialize SwapRate pointers to NULL.
*/
void     Fix3_SwapRate_Init(SWAPRATE_DATA *sr)
{
    int i;

    for (i=0; i<MAXNBDATE; ++i)
    {
        sr->SwapRate[i] = NULL;
    }

    return;
}

/*****  Fix3_SwapRate_Alloc *******************************************************/
/**
*         Allocate more SwapRate structures.
*/
int      Fix3_SwapRate_Alloc(SWAPRATE_DATA  *sr,
                             FIX3_TREE_DATA *tree_data)
{
    int i;

    int Size, Size1, Size2;

    int status = FAILURE;

    Size1 = tree_data->Width[0];
    Size2 = tree_data->Width[1] * Size1;

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
        if (sr->Denom == 'D')
        {
            Size = Size1;
        }
        else
        {
            DR_Error("Swap rates must be domestic if tree mode is EQ1IR.\n");
            goto RETURN;
        }
    }
    else
    {
        switch (sr->Denom)
        {
            case 'F':
                Size = Size1; break;

            case 'D':
                Size = Size2; break;

            default:
                DR_Error("SwapRate->Denom not properly initialized.\n");
                goto RETURN;
        }
    }

    for (i=0; i<sr->NbResetDate; ++i)
    {
        sr->SwapRate[i] = (TSLICE)DR_Array(DOUBLE, 0, Size);

        if (sr->SwapRate[i] == NULL)
        {
            goto RETURN;
        }
    }
    
    status = SUCCESS;

  RETURN:
    
    if (status == FAILURE)
    {
        DR_Error("Fix3_SwapRate_Alloc: "
                 "Unable to allocate memory for SwapRate structure.\n");
    }
    return(status);
}
    
/*****  Fix3_SwapRate_Free *******************************************************/
/**
*         Free SwapRate structures.
*/
void     Fix3_SwapRate_Free(SWAPRATE_DATA  *sr,
                            FIX3_TREE_DATA *tree_data)
{
    int i;

    int Size, Size1, Size2;

    Size1 = tree_data->Width[0];
    Size2 = tree_data->Width[1] * Size1;

    if (sr->NbResetDate > 0)
    {
        if (tree_data->TreeType == TTYPE_EQ1IR)
        {
            if (sr->Denom == 'D')
            {
                Size = Size1;
            }
            else
            {
                return;
            }
        }
        else
        {
            switch (sr->Denom)
            {
                case 'F':
                    Size = Size1; break;

                case 'D':
                    Size = Size2; break;

                default:
                    return;
            }
        }

        for (i=0; i<sr->NbResetDate; ++i)
        {
            Free_DR_Array(sr->SwapRate[i], DOUBLE, 0, Size);
        }
    }

    Fix3_SwapRate_Init(sr);
    return;
}

int Fix3_TreeSim_Init(TREESIM_DATA *ts, int NbSwapRateData)
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
        Fix3_SwapRate_Init(&(ts->SwapRateData[i]));
    }

    return SUCCESS;
}
    
int Fix3_TreeSim_Alloc(TREESIM_DATA *ts, FIX3_TREE_DATA *tree_data)
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
        if (Fix3_SwapRate_Alloc(&(ts->SwapRateData[i]),
                                tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }

    ts->i = (int *)DR_Array(INT, 0, ts->NbPathAll-1);
    ts->j = (int *)DR_Array(INT, 0, ts->NbPathAll-1);
    ts->k = (int *)DR_Array(INT, 0, ts->NbPathAll-1);

    /* ts->State[0] is used to store the initial state */ 
    ts->State [0] = (double *)DR_Array(DOUBLE, 0, ts->NbPathSub-1);

    if (ts->State[0] == NULL)
    {
        DR_Error("Could not allocate memory for state variable simulation.\n");
        return FAILURE;
    }

    for (i=0; i<ts->NbPathDate; ++i)
    {
        ts->PathAll [i] = (double *)DR_Array(DOUBLE, 0, ts->NbPathAll-1);
        ts->PathSub [i] = (double *)DR_Array(DOUBLE, 0, ts->NbPathSub-1);
        ts->State [i+1] = (double *)DR_Array(DOUBLE, 0, ts->NbPathSub-1);

        if ((ts->PathAll [i] == NULL) ||
            (ts->PathSub [i] == NULL) ||
            (ts->State [i+1] == NULL))
        {
            DR_Error("Could not allocate memory for state variable simulation.\n");
            return FAILURE;
        }
    }

    return SUCCESS;
}

void Fix3_TreeSim_FreeSR(TREESIM_DATA *ts, FIX3_TREE_DATA *tree_data)
{
    int i;

    for (i=0; i<ts->NbSwapRateData; ++i)
    {
        Fix3_SwapRate_Free(&(ts->SwapRateData[i]),
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

void Fix3_TreeSim_Free(TREESIM_DATA *ts, FIX3_TREE_DATA *tree_data)
{
    int i;

    if (ts->SwapRateData != NULL)
    {
        Fix3_TreeSim_FreeSR(ts, tree_data);
    }

    Free_DR_Array(ts->i, INT, 0, ts->NbPathAll-1);
    Free_DR_Array(ts->j, INT, 0, ts->NbPathAll-1);
    Free_DR_Array(ts->k, INT, 0, ts->NbPathAll-1);

    Free_DR_Array(ts->State[0], DOUBLE, 0, ts->NbPathSub-1);

    for (i=0; i<ts->NbPathDate; ++i)
    {
        Free_DR_Array(ts->PathAll [i], DOUBLE, 0, ts->NbPathAll-1);
        Free_DR_Array(ts->PathSub [i], DOUBLE, 0, ts->NbPathSub-1);
        Free_DR_Array(ts->State [i+1], DOUBLE, 0, ts->NbPathSub-1);
    }

    return;
}
