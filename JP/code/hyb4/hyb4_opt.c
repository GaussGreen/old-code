/****************************************************************************/
/*      Calculation of option price in the lattice.                         */
/****************************************************************************/
/*      HYB4_OPT.c                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

int Hyb4_Option_t(TSLICE     Opt,             /**< (I/O) Option prices         */
                  TSLICE      Under,           /**< (I) Underlying prices       */
                  double      Notional,        /**< (I) Notional                */
                  double      Strike,          /**< (I) Strike                  */
                  long        ExerFlag,        /**< (I) Exercise flag           */
                  int         CoP,             /**< (I) =1 for call, -1 for put */
                  int         t,               /**< (I) Current time point      */
                  int         T,               /**< (I) Total number of point   */
                  int         DCurve,          /**< (I) Discount curve          */
                  int         DMode,           /**< (I) Dim of DEV and slices   */
                  HYB4_DEV_DATA  *dev_data,       /**< (I) Hyb4_Dev data structure */
                  HYB4_TREE_DATA *tree_data)      /**< (I) Tree data structure     */
{

    double  *OptL; 
    double  *UnderL;

    int	                                
            Top1,    Bottom1,       /* Tree limits (1st dim) */
           *Top2,   *Bottom2,       /* Tree limits (2nd dim) */
          **Top3,  **Bottom3,       /* Tree limits (3rd dim) */
         ***Top4, ***Bottom4,       /* Tree limits (4th dim) */

            i, j, k, L,             /* Node indices          */
            offset,
            status = FAILURE;       /* Error status	         */


    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    if (Hyb4_Dev(Opt,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */


    if (ExerFlag)
    {
 
        if (DMode == DISC_1D_NOCUPS)
        {
            offset = tree_data->NodeOffset0[t];

            OptL   = (double *)Opt   + offset;
            UnderL = (double *)Under + offset;

            /* 
            *   If this is the last exercise date OptL has been initialized to
            *   zero so the exercise decision is Max (intrinsic value, 0).
            *   If this is not the last exercise date we have an American 
            *   exercise decision Max (intrinsic value, live option).
            */

            for (i = Bottom1; i <= Top1; i ++)                                      
            {
                OptL[i] = MAX (Notional * CoP * (UnderL[i] - Strike), OptL[i]);

            }  /* for i */	
        }
        else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
        {


            for (i = Bottom1; i <= Top1; i ++)
            {

                offset = tree_data->NodeOffset1[t][i];
                OptL   = (double *)Opt   + offset;
                UnderL = (double *)Under + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    OptL[j] = MAX (Notional * CoP * 
                         (UnderL[j] - Strike), OptL[j]);

                }  /* for j */	
            }
        }
        else if (DMode == DISC_3D_CUPS)
        {
 
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = tree_data->NodeOffset2[t][i][j];
                    OptL   = (double *)Opt   + offset;
                    UnderL = (double *)Under + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        OptL[k] = MAX (Notional * CoP * 
                             (UnderL[k] - Strike), OptL[k]);

                    }  /* for k */	
                }
            }
        }  /* if then else */
        else if (DMode == DISC_4D_CUPS)
        {
 
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        OptL   = (double *)Opt   + offset;
                        UnderL = (double *)Under + offset;

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {
                            OptL[L] = MAX (Notional * CoP * 
                                 (UnderL[L] - Strike), OptL[L]);
                        }

                    }  /* for k */	
                }
            }
        }  /* if then else */

    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb4_Option_t */

