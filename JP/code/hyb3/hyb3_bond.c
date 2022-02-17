/*****************************************************************************/
/*        Calculation of bond price in the lattice.                          */
/*****************************************************************************/
/*        BOND.C                                                             */
/*****************************************************************************/


/*
$Header$
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

     



/*****  Hyb3_Bond_Price  **********************************************************/
/*
*         Calculate the bond price in the lattice using the dev routine.
*/
int  Hyb3_Bond_Price(
         TSLICE      BondPtr,       /* (I/O) Array of bond prices at t       */
         double      Principal,     /* (I) Principal of the bond             */
         long        PrincipalFlag, /* (I) TRUE if principal paid at t       */
         double      Coupon,        /* (I) Coupon of the bond                */
         long        CouponFlag,    /* (I) TRUE if a coupon is paid at t     */
         int         t,             /* (I) Current time period               */
         int         T,             /* (I) Total number of periods in tree   */
         int         DCurve,        /* (I) Discount curve index (0,1,2)      */
         int         DMode,         /* (I) Defines the dim of the DEV        */
         HYB3_DEV_DATA   *dev_data,      /* (I) Data required for DEV'ing         */
         HYB3_TREE_DATA  *tree_data)     /* (I) Structure of tree data            */
{


    double   *BondL;      /* Local for convenience in addressing slice */

    int
        Top1,    Bottom1,       /* Limits of the tree (1rst dimension) */
       *Top2,   *Bottom2,
      **Top3,  **Bottom3,
        i,  j,  k,              /* Node indices                        */

        status = FAILURE;       /* Error status = FAILURE initially    */

          
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];


        /* Discounted expected value function */           
        if (Hyb3_Dev(BondPtr,   
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
        BondL = (double *)BondPtr + Hyb3_Node_Offset(1, 0, 0, t, tree_data);

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

    else if (DMode == DISC_2D_CUPS   || 
             DMode == DISC_2D_NOCUPS ||
             DMode == DISC_2D_1IR2F_NOCUPS)
    {
        
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL = (double *)BondPtr + Hyb3_Node_Offset(2,i,0,t,tree_data);

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

                BondL = (double *)BondPtr + Hyb3_Node_Offset(2,i,0,t,tree_data);

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL[j] += Coupon;

                }  /* for j */
            }
        }  /* if */

    }
    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        

        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL = (double *)BondPtr + Hyb3_Node_Offset(3, i, j, t, tree_data);

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
                    BondL = (double *)BondPtr + Hyb3_Node_Offset(3, i, j, t, tree_data);

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


}  /* Hyb3_Bond_Price */



/*****  Hyb3_FwdBond_t  ***********************************************************/
/*
*         Calculate the fwd bond price in the lattice using the dev routine.
*   
*         This routine does not deal with principal exchanges.
*/
int  Hyb3_FwdBond_t(
         TSLICE      FwdBond,
         TSLICE      Bond,          /* (I/O) Array of bond prices at t      */
         TSLICE      ZeroToPmt,     /* (I) Zero bond to cp pmt date         */
         double      Coupon,        /* (I) Coupon of the bond               */
         long        CouponFlag,    /* (I) TRUE if a coupon is paid at t    */
         double      CouponRate,    /* (I) Fixed rate paid                  */
         long        CouponAccSt,   /* (I) Accrual start for current coupon */
         char        CouponDCConv,  /* (I) Day count convention used        */
         double      Outstanding,   /* (I) Outstanding for cp calculation   */
         long        ExerFlag,      /* (I) Exercise flag                    */
         char        OptStubConv,   /* (I) Option stub rule                 */
         long        CurrentDate,   /* (I) Current date                     */
         int         t,             /* (I) Current time period              */
         int         T,             /* (I) Total number of periods in tree  */
         int         DCurve,        /* (I) Discount curve index (0,1,2)     */
         int         DMode,         /* (I) Defines the dim of the DEV       */
         HYB3_DEV_DATA   *dev_data,      /* (I) Data required for DEV'ing        */
         HYB3_TREE_DATA  *tree_data)     /* (I) Structure of tree data           */

{


    double   *FwdBondL;
    double   *BondL;        /* Local for convenience in addressing */
    double    ZeroAtNode;

    double    Accrued = 0.0;

    int
        Top1,    Bottom1,   /* Limits of the tree (1rst dimension) */
       *Top2,   *Bottom2,
      **Top3,  **Bottom3,
        i,  j,  k,          /* Node indices                        */
        offset,

    status = FAILURE;       /* Error status = FAILURE initially    */

      
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /*  If this is not an exercise date discount and return. */
    if (!ExerFlag)                                              
    {

        /* Discounted expected value function */           
        if (Hyb3_Dev(FwdBond,   
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) != SUCCESS)
        {
            goto RETURN;
        }
         
        status = SUCCESS;

        return (status);
        
    } 


    /*   Start from "basic" bond as passed in  */
    if (DMode == DISC_1D_NOCUPS)
    {
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        FwdBondL = (double *)FwdBond + offset;
        BondL    = (double *)Bond + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            FwdBondL[i] = BondL[i];
        }  

    }
    else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);  
            FwdBondL  = (double *)FwdBond + offset; 
            BondL     = (double *)Bond + offset;    

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                FwdBondL[j] = BondL[j];
            }  
        }
    }
    else if (DMode == DISC_3D_CUPS)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                                                            
                FwdBondL = (double *)FwdBond + offset;          
                BondL    = (double *)Bond    + offset;          
          
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    FwdBondL[k] = BondL[k];
                }
            }
        }
    }



    /* Now do the stubs (i.e. pmt of accruals) */

    /*  Calculate accrued interests if we are within accrual period. */
    if (CouponAccSt <= CurrentDate)
    {
        if (DrDayCountFraction (CouponAccSt,
                                CurrentDate,                                 
                                CouponDCConv,
                                &Accrued) == FAILURE)
        {
            DR_Error("Unable to calculate accrued dcf for fwd bond (Hyb3_FwdBond_t)");
            goto RETURN;
        }        
    }

    if (DMode == DISC_1D_NOCUPS)
    {
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        FwdBondL = (double *)FwdBond + offset;

        if (CouponFlag) /* if just paid cp, then take it away */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                FwdBondL[i] -= Coupon;
            }
        } 
        else if (OptStubConv == 'S' || OptStubConv == 'B')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                ZeroAtNode  = Hyb3_GetValueAtNode(1,ZeroToPmt,i,0,0,t,tree_data);
                FwdBondL[i] -= CouponRate * ZeroAtNode * Accrued * Outstanding ;    
            }           
        }
    }
    else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
    {
        
        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);  
                FwdBondL  = (double *)FwdBond + offset; 
                
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    FwdBondL[j] -= Coupon;

                } 
            }
        }  
        else if (OptStubConv == 'S' || OptStubConv == 'B')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                FwdBondL = (double *)FwdBond + Hyb3_Node_Offset(2,i,0,t,tree_data);
                
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZeroAtNode  = Hyb3_GetValueAtNode(2,ZeroToPmt,i,j,0,t,tree_data);
                    FwdBondL[j] -= CouponRate*ZeroAtNode*Accrued*Outstanding;
                } 
            }
        }
    }
    else if (DMode == DISC_3D_CUPS)
    {

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    FwdBondL =(double *)FwdBond+Hyb3_Node_Offset(3,i,j,t,tree_data);

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        FwdBondL[k] -= Coupon;
                    }  
                }
            }
        }  
        else if (OptStubConv == 'S' || OptStubConv == 'B')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    FwdBondL = (double *)FwdBond + offset;                 
          
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        ZeroAtNode   = Hyb3_GetValueAtNode(3,ZeroToPmt,i,j,k,t,tree_data);
                        FwdBondL[k] -= CouponRate * ZeroAtNode * Accrued * Outstanding ;
                    }
                }
            }
        }
    }  /* if then else */




    status = SUCCESS;

  RETURN:

    return (status);


}  /* Hyb3_FwdBond_t */

