/****************************************************************************/
/*      Utility routines.                                                   */
/****************************************************************************/
/*      HYB4_UTIL.c                                                         */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>  
#include <math.h>
#include <string.h>
#include "hyb4_lib.h"

/*****  Hyb4_GetIndexStep  *******************************************************/

double   Hyb4_GetIndexStep(
		              TSLICE           Index,      /**< (I) Index pointer       */
                      int              Dim,        /**< (I) Index dimension     */
                      int              i,          /**< (I) Node indices        */
                      int              j,
                      int              k,
                      int              L,
                      int              t,          /**< (I) Current time point  */
                      HYB4_TREE_DATA  *tree_data)  /**< (I) Tree data structure */
{

    double  *IndexL;            /* Local pointer */

    double  IndexStep;          /* Output index step */
    double  IndexVal;           /* Index value at mid node */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */
    int     ***Top4, ***Bottom4;  /* Tree limits (3rd dim)  */

    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    IndexStep = ERROR;  /* To avoid division by 0 */
                
    switch (Dim)
    {
        case 1:
        {
            IndexL = (double *)Index + tree_data->NodeOffset0[t];

            if (i > Bottom1)
                IndexStep = MAX (IndexStep, fabs (IndexL[i-1] - IndexL[i]));
            if (i < Top1)                                                           
                IndexStep = MAX (IndexStep, fabs (IndexL[i+1] - IndexL[i]));

            break;
        }                    
        case 2:
        {
            IndexL = (double *)Index + tree_data->NodeOffset1[t][i];

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
                IndexL = (double *)Index + tree_data->NodeOffset1[t][i-1];

                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }
            if (i < Top1)
            {         
                IndexL = (double *)Index + tree_data->NodeOffset1[t][i+1];

                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {                                                  
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }

            break;
        }                    
        case 3:
        {
            IndexL = (double *)Index + tree_data->NodeOffset2[t][i][j];

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
                    IndexL = (double *)Index + tree_data->NodeOffset2[t][i-1][j];

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
                    IndexL = (double *)Index + tree_data->NodeOffset2[t][i+1][j];

                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {                                                       
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (j > Bottom2[i])
            {
                IndexL = (double *)Index + tree_data->NodeOffset2[t][i][j-1];

                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
            if (j < Top2[i])
            {
                IndexL = (double *)Index + tree_data->NodeOffset2[t][i][j+1];

                if ((k>=Bottom3[i][j+1])&&(k<=Top3[i][j+1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }

            break;
        }                    
        case 4:
        {
            IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j][k];

            IndexVal = IndexL[L];

            if (i-1 >= Bottom1)
            {
                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    if ((k>=Bottom3[i-1][j])&&(k<=Top3[i-1][j]))
                    {
                        if ((L>=Bottom4[i-1][j][k])&&(L<=Top4[i-1][j][k]))
                        {
                            IndexL = (double *)Index + tree_data->NodeOffset3[t][i-1][j][k];

                            IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                        }
                    }
                }
            }
            if (i+1 <= Top1)
            {
                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {
                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {
                        if ((L>=Bottom4[i+1][j][k])&&(L<=Top4[i+1][j][k]))
                        {
                            IndexL = (double *)Index + tree_data->NodeOffset3[t][i+1][j][k];

                            IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                        }
                    }
                }
            }

            if (j-1 >= Bottom2[i])
            {
                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    if ((L>=Bottom4[i][j-1][k])&&(L<=Top4[i][j-1][k]))
                    {
                        IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j-1][k];

                        IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                    }
                }
            }
            if (j+1 <= Top2[i])
            {
                if ((k>=Bottom3[i][j+1])&&(k<=Top3[i][j+1]))
                {
                    if ((L>=Bottom4[i][j+1][k])&&(L<=Top4[i][j+1][k]))
                    {
                        IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j+1][k];

                        IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                    }
                }
            }

            if (k-1 >= Bottom3[i][j])
            {
                if ((L>=Bottom4[i][j][k-1])&&(L<=Top4[i][j][k-1]))
                {
                    IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j][k-1];

                    IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                }
            }
            if (k+1 <= Top3[i][j])
            {
                if ((L>=Bottom4[i][j][k+1])&&(L<=Top4[i][j][k+1]))
                {
                    IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j][k+1];

                    IndexStep = MAX (IndexStep, fabs (IndexL[L] - IndexVal));
                }
            }

            if (L-1 >= Bottom4[i][j][k])
            {
                IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j][k];

                IndexStep = MAX (IndexStep, fabs (IndexL[L-1] - IndexVal));
            }
            if (L+1 <= Top4[i][j][k])
            {
                IndexL = (double *)Index + tree_data->NodeOffset3[t][i][j][k];

                IndexStep = MAX (IndexStep, fabs (IndexL[L+1] - IndexVal));
            }

            break;
        }                    
        default:
        {
            break;
        }                    
    }

    return (IndexStep);

}  /* Hyb4_GetIndexStep */
