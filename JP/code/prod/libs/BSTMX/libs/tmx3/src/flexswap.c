/****************************************************************************/
/*      Flexi swap.                                                         */
/****************************************************************************/
/*      FLEXSWAP.c                                                          */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"



/*****  Flexswap_t  **********************************************************/
/*
*       Flexi swap price.
*/
int     Flexswap_t (double      **Swap,      /* (I/O) Flexi swaps            */
                    double      **Coupon,    /* (I) Recv and paid coupons    */
                    long        *ResetFlag,  /* (I) Recv and paid reset flag */
                    char        *Arrears,    /* (I) Recv and paid arrears    */        
                    long        ExerFlag,    /* (I) Exercise flag            */
                    double      Strike,      /* (I) Strike for penalty       */
                    int         NbOuts,      /* (I) Maximum number of swaps  */
                    double      *OutsLevel,  /* (I) Array of outstandings    */
                    double      LastLowOuts, /* (I) Last lowest outs         */
                    double      LastHighOuts,/* (I) Last highest outs        */
                    double      NewLowOuts,  /* (I) New lowest outs          */
                    double      NewHighOuts, /* (I) New highest outs         */
                    long        OptExerFlag, /* (I) Optimal exer search flag */
                    double      *OptimalOuts,/* (I) Optimal outs on vald exer*/
                    int         t,           /* (I) Current time point       */
                    int         T,           /* (I) Last time point          */
                    int         DCurve,      /* (I) Discount curve           */
                    DEV_DATA    *dev_data,   /* (I) Dev data structure       */
                    TREE_DATA   *tree_data)  /* (I) Tree data structure      */
{

    double  *CombCoupon=NULL;        /* Combined coupon variable          */
    double  **SwapL=NULL;            /* Local slice pointer               */

    double  currValue;               /* Processed swap value              */
    double  maxValue;                /* Maximal swap value for given outs */

    int     optimal;                 /* Optimal outs level index          */
    int     idx;                     /* Outstanding index                 */
    int     New, Last;               /* Last and new indices              */
    int     LastLow=0, LastHigh=0;   /* Last low and high positions       */
    int     NewLow=0, NewHigh=0;     /* New low and high positions        */
    int     Weight[2];               /* Each coupon weight                */

    int     MinEff;                  /* Effective min amortisation        */

    int     Top1, Bottom1;           /* Tree limits (1rst dim)            */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)             */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)             */

    int     i, j, k;                 /* Node indices                      */
    int     offset;                  /* Node offset                       */
    int     l;                       /* Chooser cap index                 */
    int     status = FAILURE;        /* Error status                      */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    SwapL = (double **) DR_Array (DOUBLE_PTR, 0, NbOuts-1);
    if (SwapL == NULL)
    {
        DR_Error ("Flexswap_t: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }

    /* Allocate combined coupon variable if reset */
    if (ResetFlag[0] || ResetFlag[1])
    {
        CombCoupon = Alloc_Slice (tree_data);
        if (CombCoupon == NULL)
        {
            DR_Error ("Flexswap_t: could not allocate memory!");
            goto FREE_MEM_AND_RETURN;        
        }
    }


    /* Find last and new, low and high, position indices */
    for (idx = 0; idx < NbOuts; idx ++)
    {
        if (fabs(OutsLevel[idx]-LastLowOuts)  < ERROR)  LastLow  = idx;
        if (fabs(OutsLevel[idx]-LastHighOuts) < ERROR)  LastHigh  = idx;
        if (fabs(OutsLevel[idx]-NewLowOuts)   < ERROR)  NewLow  = idx;
        if (fabs(OutsLevel[idx]-NewHighOuts)  < ERROR)  NewHigh  = idx;
    }


    /* Dev swaps for relevant outstanding levels */
    for (idx = LastLow; idx <= LastHigh; idx++)
    {
        if (Dev (Swap[idx],
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }  
    } /* for idx */


    /* Add coupon to swap if reset-in-advance. Relevant outstandings are in [LastLow,LastHigh] */
    for (l = 0; l < 2; l++) Weight[l] = (ResetFlag[l] && (Arrears[l] == 'N')) ? 1 : 0;

    if (Weight[0] || Weight[1])
    {
        for (idx = LastLow; idx <= LastHigh; idx++)
        {
            if (LCombTwoSlices(CombCoupon,
                               Coupon[0],
                               +1. * Weight[0] * OutsLevel[idx],
                               Coupon[1],
                               -1. * Weight[1] * OutsLevel[idx],
                               t,
                               tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
            
            if (AddTwoSlices(Swap[idx],
                             Swap[idx],
                             CombCoupon,
                             t,
                             tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
        } /* for ic */
    } /* if reset */


    /* Determine optimal amortization level */
    if (OptExerFlag)
    {
        /* Loop over outstanding between last low and last high to determine the index   */
        /* of the most expensive swap. The given outstanding is the new high outstanding */
        New      = NewHigh; 
        Last     = LastLow; 
        optimal  = Last; 
        MinEff = MIN(LastHigh, New);                
        
        if (tree_data->NbFactor == 1)
        {
            for (idx = 0; idx < NbOuts; idx ++)
            {
                SwapL[idx] = Swap[idx] + Node_Offset(1, 0, 0, t, tree_data);
            }
            
            
            maxValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
            for (Last = LastLow + 1; Last <= LastHigh; Last++)
            {
                currValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                if (currValue > maxValue)
                {
                    maxValue = currValue;
                    optimal  = Last;
                }
            }  /* for Last */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (idx = 0; idx < NbOuts; idx ++)
            {
                SwapL[idx] = Swap[idx] + Node_Offset(2, 0, 0, t, tree_data);
            }
                    
            maxValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);            
            for (Last = LastLow + 1; Last <= LastHigh; Last++)
            {
                currValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                if (currValue > maxValue)
                {
                    maxValue = currValue;
                    optimal  = Last;
                }
            }  /* for Last */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (idx = 0; idx < NbOuts; idx ++)
            {
                SwapL[idx] = Swap[idx] + Node_Offset(3, 0, 0, t, tree_data);
            }
            
            maxValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);            
            for (Last = LastLow + 1; Last <= LastHigh; Last++)
            {
                currValue = SwapL[Last][0] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                if (currValue > maxValue)
                {
                    maxValue = currValue;
                    optimal  = Last;
                }
             }  /* for Last */
        }  /* if then else */

        /* Record optimal outstanding */
        *OptimalOuts = OutsLevel[optimal];
        
    }  /* if */


    /* Flexi swap action: amortization of notional on exercise date */
    if (ExerFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (idx = 0; idx < NbOuts; idx ++)
            {
                SwapL[idx] = Swap[idx] + offset;
            }
            
            /* Need reverse order on new outstandings here: start with highest, end with lowest */

             for (i = Bottom1; i <= Top1; i ++)
             {
                for (New = NewHigh; New >= NewLow; New--)
                {
                    Last = LastLow; 
                    MinEff = MIN(LastHigh, New);
                    maxValue = SwapL[Last][i] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);

                    for (Last = LastLow + 1; Last <= MinEff; Last++)
                    {
                        currValue = SwapL[Last][i] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                        maxValue  = MAX(currValue, maxValue); 
                    }
                    SwapL[New][i] = maxValue;
                 
                }  /* for New */        
            }  /* for i */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (idx = 0; idx < NbOuts; idx ++)
                {
                    SwapL[idx] = Swap[idx] + offset;
                }
            
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (New = NewHigh; New >= NewLow; New--)
                    {
                        Last = LastLow; 
                        MinEff = MIN(LastHigh, New);
                        maxValue = SwapL[Last][j] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                        
                        for (Last = LastLow + 1; Last <= MinEff; Last++)
                        {
                            currValue = SwapL[Last][j] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                            maxValue  = MAX(currValue, maxValue); 
                        }
                        SwapL[New][j] = maxValue;

                    }  /* for New */        
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (idx = 0; idx < NbOuts; idx ++)
                    {
                        SwapL[idx] = Swap[idx] + offset;
                    }
            
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (New = NewHigh; New >= NewLow; New--)
                        {
                            Last = LastLow;
                            MinEff = MIN(LastHigh, New); 
                            maxValue = SwapL[Last][k] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                            
                            for (Last = LastLow + 1; Last <= MinEff; Last++)
                            {
                                currValue = SwapL[Last][k] - Strike * (OutsLevel[MinEff] - OutsLevel[Last]);
                                maxValue  = MAX(currValue, maxValue); 
                            }
                            SwapL[New][k] = maxValue;
                            
                        }  /* for New */        
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if */


    /* Add coupon to swap if reset-in-arrears. Relevant outstandings are in [NewLow,NewHigh] */
    for (l = 0; l < 2; l++) Weight[l] = (ResetFlag[l] && (Arrears[l] == 'Y')) ? 1 : 0;

    if (Weight[0] || Weight[1])
    {
        for (idx = NewLow; idx <= NewHigh; idx++)
        {
            if (LCombTwoSlices(CombCoupon,
                               Coupon[0],
                               +1. * Weight[0] * OutsLevel[idx],
                               Coupon[1],
                               -1. * Weight[1] * OutsLevel[idx],
                               t,
                               tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
            
            if (AddTwoSlices(Swap[idx],
                             Swap[idx],
                             CombCoupon,
                             t,
                             tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }  
        } /* for ic */
    } /* if reset */


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (SwapL, DOUBLE_PTR, 0, NbOuts-1);

    if (ResetFlag[0] || ResetFlag[1])
    {
        Free_Slice (CombCoupon, tree_data);
    }

    return (status);

}  /* Flexswap_t */
