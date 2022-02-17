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
#include "fix123head.h"


/*****  Fix3_GetIndexStep  *******************************************************/
/*                                                                           
 *      Obtains the difference between values of an index at adjacent nodes.
 *      This value is then used in the smoothing algorithm (Smooth_Step).
 *      This function hides the dimension generality from  the  caller, 
 *      but it has no means to check that the adequate amount of space has 
 *      been allocated under the void * being passed.
 */
double   Fix3_GetIndexStep ( double const*    Index,      /* (I) Index pointer       */
                        int              Dim,        /* (I) Index dimension     */
                        int              i,          /* (I) Node indices        */
                        int              j,
                        int              k,
                        int              t,          /* (I) Current time point  */
                        FIX3_TREE_DATA const* tree_data)	/* (I) Tree data structure */
{

    double const*   IndexL;            /* Local pointer */

    double          IndexStep;          /* Output index step */
    double          IndexVal;           /* Index value at mid node */

    int             Top1, Bottom1;      /* Tree limits (1rst dim) */
    int             *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int             **Top3, **Bottom3;  /* Tree limits (3rd dim)  */



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
            IndexL = Index + Fix3_Node_Offset(1, 0, 0, t, tree_data);

            if (i > Bottom1)
                IndexStep = MAX (IndexStep, fabs (IndexL[i-1] - IndexL[i]));
            if (i < Top1)                                                           
                IndexStep = MAX (IndexStep, fabs (IndexL[i+1] - IndexL[i]));

            break;
        }                    
        case 2:
        {
            IndexL = Index + Fix3_Node_Offset(2, i, 0, t, tree_data);

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
                IndexL = Index + Fix3_Node_Offset(2, i-1, 0, t, tree_data);

                if ((j>=Bottom2[i-1])&&(j<=Top2[i-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }
            if (i < Top1)
            {         
                IndexL = Index + Fix3_Node_Offset(2, i+1, 0, t, tree_data);

                if ((j>=Bottom2[i+1])&&(j<=Top2[i+1]))
                {                                                  
                    IndexStep = MAX (IndexStep, fabs (IndexL[j] - IndexVal));
                }
            }

            break;
        }                    
        case 3:
        {
            IndexL = Index + Fix3_Node_Offset(3, i, j, t, tree_data);

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
                    IndexL = Index + Fix3_Node_Offset(3, i-1, j, t, tree_data);

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
                    IndexL = Index + Fix3_Node_Offset(3, i+1, j, t, tree_data);

                    if ((k>=Bottom3[i+1][j])&&(k<=Top3[i+1][j]))
                    {                                                       
                        IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                    }
                }
            }
            if (j > Bottom2[i])
            {
                IndexL = Index + Fix3_Node_Offset(3, i, j-1, t, tree_data);

                if ((k>=Bottom3[i][j-1])&&(k<=Top3[i][j-1]))
                {
                    IndexStep = MAX (IndexStep, fabs (IndexL[k] - IndexVal));
                }
            }
            if (j < Top2[i])
            {
                IndexL = Index + Fix3_Node_Offset(3, i, j+1, t, tree_data);

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

}  /* Fix3_GetIndexStep */


