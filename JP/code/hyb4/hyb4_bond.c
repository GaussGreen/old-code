/*****************************************************************************/
/*        Calculation of bond price in the lattice.                          */
/*****************************************************************************/
/*        HYB4_BOND.C                                                        */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

/*****  Hyb4_Bond_Price  **********************************************************/
/*
*         Calculate the bond price in the lattice using the new dev routine.
*/
int  Hyb4_Bond_Price(
         TSLICE      BondPtr,       /* (I/O) Array of bond prices at t       */
         double      Principal,     /* (I) Principal of the bond             */
         long        PrincipalFlag, /* (I) TRUE if principal paid at t       */
         double      Coupon,        /* (I) Coupon of the bond                */
         long        CouponFlag,    /* (I) TRUE if a coupon is paid at t     */
         int         t,             /* (I) Current time period               */
         int         T,             /* (I) Total number of periods in tree   */
         int         DCurve,        /* (I) Discount curve index (0,1,2)      */
         int         DMode,         /* (I) Defines the dim of the DEV        */
         HYB4_DEV_DATA  *dev_data,      /* (I) Data required for DEV'ing         */
         HYB4_TREE_DATA *tree_data)     /* (I) Structure of tree data            */
{


    double   *BondL;      /* Local for convenience in addressing slice */

    int
        Top1,    Bottom1,       /* Limits of the tree (1rst dimension) */
       *Top2,   *Bottom2,
      **Top3,  **Bottom3,
        i,  j,  k,              /* Node indices                        */

    status = FAILURE;       /* Error status = FAILURE initially    */

    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];

    /* Discounted expected value function */           
    if (Hyb4_Dev(BondPtr,   
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) != SUCCESS)
    {
        goto RETURN;
    }
         
    if (DMode == DISC_1D_NOCUPS)
    {
        BondL = (double *)(BondPtr) + tree_data->NodeOffset0[t];

        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL[i] += Principal;

            }  /* for i */
        }  /* if */
    

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL[i] += Coupon;

            }  /* for i */
        }  /* if */
    }

    else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
    {
        
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL = (double *)(BondPtr) + tree_data->NodeOffset1[t][i];

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL[j] += Principal;

                }  /* for j */
            }
        }  /* if */
    

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {

                BondL = (double *)(BondPtr) + tree_data->NodeOffset1[t][i];

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL[j] += Coupon;

                }  /* for j */
            }
        }  /* if */

    }
    else if (DMode == DISC_3D_CUPS)
    {
        

        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL = (double *)(BondPtr) + tree_data->NodeOffset2[t][i][j];

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        BondL[k] += Principal;

                    }  /* for k */
                }
            }
        }  /* if */
    

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL = (double *)(BondPtr) + tree_data->NodeOffset2[t][i][j];

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        BondL[k] += Coupon;

                    }  /* for k */
                }
            }
        }  /* if */
    }  /* if then else */

    status = SUCCESS;

  RETURN:

    return (status);

}  /* Hyb4_Bond_Price */
