/****************************************************************************/
/*      Standard bond.                                                      */
/****************************************************************************/
/*      BOND.c                                                              */
/****************************************************************************/

/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/bond.c,v 1.5 2001/04/09 09:46:57 lcooper Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"




/*****  Bond_t  *************************************************************/
/*
*       Bond price. 
*/
int     Bond_t (double     *Bond,           /* (I/O) Bond prices             */
                double      Principal,      /* (I) Principal payment         */
                long        PrincipalFlag,  /* (I) TRUE if principal is paid */
                double      Coupon,         /* (I) Coupon payment            */
                long        CouponFlag,     /* (I) TRUE if a coupon is paid  */
                int         t,              /* (I) Current time point        */
                int         T,              /* (I) Last time point           */
                int         DCurve,         /* (I) Discount curve            */
                DEV_DATA    *dev_data,      /* (I) Dev data structure        */
                TREE_DATA   *tree_data)     /* (I) Tree data structure       */
{

    double  *BondL;                         /* Local slice pointer    */

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */
    int     i, j, k;                        /* Node indices           */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   Bond,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;                    
    }

    
    if (tree_data->NbFactor == 1)
    {
        BondL = Bond + Node_Offset(1, 0, 0, t, tree_data);
    
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL[i] += Principal;
            }
        }  /* if */

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL[i] += Coupon;
            }
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL = Bond + Node_Offset(2, i, 0, t, tree_data);

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL[j] += Principal;
                }
            }  /* for i */
        }  /* if */

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                BondL = Bond + Node_Offset(2, i, 0, t, tree_data);

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL[j] += Coupon;
                }
            }  /* for i */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL = Bond + Node_Offset(3, i, j, t, tree_data);

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        BondL[k] += Principal;
                    }
                }  /* for j */
        }  /* if */    

        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    BondL = Bond + Node_Offset(3, i, j, t, tree_data);

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        BondL[k] += Coupon;
                    }
                }  /* for j */
        }  /* if */
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Bond_t */



/*****  FwdBond_t  **********************************************************/
/*
*       Returns the price of a fwd bond contract against a given strike.
*       Input: Underlying bond price (including initial principle payment if
*              it occured on/after current date)
*
*/
int     FwdBond_t (double      *FwdBond,    /* (O) Forward bond              */
                   double      *Bond,       /* (I) Underlying bond           */
                   double      *ZeroToPmt,  /* (I) Zero to next payment      */
                   long        CouponFlag,  /* (I) Coupon flag               */
                   double      CpnOS,       /* (I) OS Prn for the cpn rate   */
                   double      CpnRate,     /* (I) Coupon rate               */
                   char        DayCount,    /* (I) Day count convention      */
                   long        CpAccStDate, /* (I) Accrual start of last cpn */
                   double      LastCpn,     /* (I) Last coupon payment       */
                   double      TotalPrn,    /* (I) Tot Prn paid/rec to date  */
                   long        ExerFlag,    /* (I) Exercise flag             */
                   double      Strike,      /* (I) Current strike            */
                   char        OptStub,     /* (I) Option stub rule          */
                   long        CurrentDate, /* (I) Current date              */
                   int         t,           /* (I) Current time period       */
                   int         T,           /* (I) Total number of period    */
                   int         DCurve,      /* (I) Discount curve            */
                   DEV_DATA    *dev_data,   /* (I) Dev data structure        */
                   TREE_DATA   *tree_data)  /* (I) Tree data structure       */
{

    double  *FwdBondL;               /* Local slice pointer     */
    double  *BondL;
    double  *ZeroToPmtL;

    double  InitialPmt;              /* Amt to pay at exercise  */
    double  Accrued =0.0;            /* Accrued amt at exercise */

    int     Top1, Bottom1;           /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)   */

    int     i, j, k;                 /* Node indices            */
    int     offset;                  /* Node offset             */
    int     status = FAILURE;        /* Error status            */
    int     StubAdjRequired;         /* TRUE if stub adj */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /*
    *   If this is not an exercise date discount and return.
    */
    if (!ExerFlag)                                              
    {
        if (Dev (   FwdBond,
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
        }
    
        return (SUCCESS);
        
    }  /* if */


    /* 
    *   Calculate accrued interest if necessary
    */              

    if (CurrentDate > CpAccStDate)
    {
        if (DrDayCountFraction (CpAccStDate,
                                CurrentDate,                                 
                                DayCount,
                                &Accrued) == FAILURE)
        {
            DR_Error("FwdBond_t: unable to calculate dcf for accrued!");
            goto RETURN;
        }
        
        Accrued *= CpnOS * CpnRate;
        StubAdjRequired = TRUE;
    }
    else
    {
        StubAdjRequired = FALSE;
    }
    
                                                
    /* 
    *   Calculate the initial payment at exercise
    */

    InitialPmt = TotalPrn + Strike;


    /* 
    *   Calculate the fwd bond contract price at exercise
    */

    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        FwdBondL   = FwdBond   + offset;
        BondL      = Bond      + offset;
        ZeroToPmtL = ZeroToPmt + offset;
    
        for (i = Bottom1; i <= Top1; i ++)
        {
            FwdBondL[i] = BondL[i] - InitialPmt;
        }

        if (StubAdjRequired)
        {

            /* Subtract coupon payment just been made */
            if (CouponFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    FwdBondL[i] -= LastCpn;
                }
            }
            /* Bond stub convention: subtract accrued interests */
            else if (OptStub == 'B')
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    FwdBondL[i] -= Accrued;
                }
            }
            /* FwdBond convention: subtract last coupon and replace by stub */
            else if (OptStub == 'S')
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    FwdBondL[i] -= Accrued * ZeroToPmtL[i];
                }
            }  /* if then else */

        } /* StubAdjRequired */
    }
    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Node_Offset(2, i, 0, t, tree_data);

            FwdBondL = FwdBond + offset;
            BondL    = Bond    + offset;
    
            for (j = Bottom2[i]; j <= Top2[i]; j++)             
            {
                FwdBondL[j] = BondL[j] - InitialPmt;
            }
        }  /* for i */


        if (StubAdjRequired)
        {
        
            if (CouponFlag)
            {                   
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);
        
                    FwdBondL = FwdBond + offset;
        
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        FwdBondL[j] -= LastCpn;
                    }
                }  /* for i */
            }                                                           
            else if (OptStub == 'B')
            {                                                                   
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);
        
                    FwdBondL = FwdBond + offset;
        
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        FwdBondL[j] -= Accrued;
                    }
                }  /* for i */
            }
            else if (OptStub == 'S')
            {                                                                       
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);
        
                    FwdBondL      = FwdBond      + offset;
                    ZeroToPmtL = ZeroToPmt + offset;
        
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        FwdBondL[j] -= Accrued * ZeroToPmtL[j];
                    }
                }  /* for i */
            }  /* if then else */                        

        } /* StubAdjRequired */
    }
    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Node_Offset(3, i, j, t, tree_data);

                FwdBondL = FwdBond + offset;
                BondL    = Bond    + offset;
    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                {
                    FwdBondL[k] = BondL[k] - InitialPmt;
                }
            }  /* for j */


        if (StubAdjRequired)
        {
        
            if (CouponFlag)
            {                   
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
        
                        FwdBondL = FwdBond + offset;
        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            FwdBondL[k] -= LastCpn;
                        }
                    }  /* for j */
            }                                                           
            else if (OptStub == 'B')
            {                                                                   
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
        
                        FwdBondL = FwdBond + offset;
        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            FwdBondL[k] -= Accrued;
                        }
                    }  /* for j */
            }
            else if (OptStub == 'S')
            {                                                                       
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)             
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
        
                        FwdBondL      = FwdBond      + offset;
                        ZeroToPmtL = ZeroToPmt + offset;
        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            FwdBondL[k] -= Accrued * ZeroToPmtL[k];
                        }
                    }  /* for j */
            }  /* if then else */

        } /* StubAdjRequired */

    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* FwdBond_t */
