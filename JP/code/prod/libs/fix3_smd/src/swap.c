 /****************************************************************************/
/*      Calculation of forward swap in the lattice.                         */
/****************************************************************************/
/*      SWAP.c                                                              */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"


                        
/*****  Fix3_ParSwap_t  **********************************************************/
/*
*       Forward par swap price. Floating leg is a fraction of par.
*/
int     Fix3_ParSwap_t (double      *Swap,       /* (O) Forward swap              */
                   double      *Bond,       /* (I) Underlying bond           */
                   double      *ZeroToPmt,  /* (I) Zero to next payment      */
                   double      Outstanding, /* (I) Outstanding notional      */
                   long        CouponFlag,  /* (I) Coupon flag               */
                   double      CpnRate,     /* (I) Coupon rate               */
                   char        DayCount,    /* (I) Day count convention      */
                   long        CpAccStDate, /* (I) Accrual start of last cpn */
                   double      LastCpn,     /* (I) Last coupon payment       */
                   long        ExerFlag,    /* (I) Exercise flag             */
                   double      Strike,      /* (I) Current strike            */
                   char        OptStub,     /* (I) Option stub rule          */
                   long        CurrentDate, /* (I) Current date              */
                   int         t,           /* (I) Current time period       */
                   int         T,           /* (I) Total number of period    */
                   int         DCurve,      /* (I) Discount curve            */
                   FIX3_DEV_DATA    *dev_data,   /* (I) Fix3_Dev data structure        */
                   FIX3_TREE_DATA   *tree_data)  /* (I) Tree data structure       */
{

    double  *SwapL;                  /* Local slice pointer    */
    double  *BondL;
    double  *ZeroToPmtL;

    double  FloatLeg;                /* Value of the floating leg        */
    double  Accrued;                 /* Accrued interest at current date */

    int     Top1, Bottom1;           /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)  */

    int     i, j, k;                 /* Node indices           */
    int     offset;                  /* Node offset            */
    int     status = FAILURE;        /* Error status           */


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
        if (Fix3_Dev (   Swap,
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
    *   Calculate accrued interest if we are within accrual period.
    */              
    if (CpAccStDate <= CurrentDate)
    {
        if (DrDayCountFraction (CpAccStDate,
                                CurrentDate,                                 
                                DayCount,
                                &Accrued) == FAILURE)
        {
            DR_Error("Fix3_ParSwap_t: unable to calculate dcc for accrued!");
            goto RETURN;
        }
        
        Accrued *= Outstanding * CpnRate;        
    }
    else
    {
        DR_Error("Fix3_ParSwap_t: version does not allow for fwd starting swap!");
        goto RETURN;

    }  /* else if */
    
                                                
    /* 
    *   Value of floating leg: par + possible strike.
    */
    FloatLeg = Outstanding + Strike;


    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL      = Swap      + offset;
        BondL      = Bond      + offset;
        ZeroToPmtL = ZeroToPmt + offset;
    
        for (i = Bottom1; i <= Top1; i ++)
        {
            SwapL[i] = BondL[i] - FloatLeg;
        }

        /* Subtract coupon payment just been made */
        if (CouponFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= LastCpn;
            }
        }
        /* Bond stub convention: subtract accrued interests */
        else if (OptStub == 'B')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= Accrued;
            }
        }
        /* Swap convention: subtract last coupon and replace by stub */
        else if (OptStub == 'S')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= Accrued * ZeroToPmtL[i];
            }
        }  /* if then else */                        
    }
    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            SwapL = Swap + offset;
            BondL = Bond + offset;
    
            for (j = Bottom2[i]; j <= Top2[i]; j++)             
            {
                SwapL[j] = BondL[j] - FloatLeg;
            }
        }  /* for i */
        
        if (CouponFlag)
        {                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    SwapL[j] -= LastCpn;
                }
            }  /* for i */
        }                                                           
        else if (OptStub == 'B')
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    SwapL[j] -= Accrued;
                }
            }  /* for i */
        }
        else if (OptStub == 'S')
        {                                                                       
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL      = Swap      + offset;
                ZeroToPmtL = ZeroToPmt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    SwapL[j] -= Accrued * ZeroToPmtL[j];
                }
            }  /* for i */
        }  /* if then else */                        
    }
    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                SwapL = Swap + offset;
                BondL = Bond + offset;
    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                {
                    SwapL[k] = BondL[k] - FloatLeg;
                }
            }  /* for j */
        
        if (CouponFlag)
        {                   
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        SwapL[k] -= LastCpn;
                    }
                }  /* for j */
        }                                                           
        else if (OptStub == 'B')
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        SwapL[k] -= Accrued;
                    }
                }  /* for j */
        }
        else if (OptStub == 'S')
        {                                                                       
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)             
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL      = Swap      + offset;
                    ZeroToPmtL = ZeroToPmt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        SwapL[k] -= Accrued * ZeroToPmtL[k];
                    }
                }  /* for j */
        }  /* if then else */                        
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_ParSwap_t */



/*****  Fix3_AmortSwap_t  ********************************************************/
/*
*       Forward swap price. Floating leg is arbitrary.
*/
int     Fix3_AmortSwap_t (
             double      *Swap,         /* (O) Forward swap                  */
             double      *Bond,         /* (I) Underlying bond               */
             double      *Floater,      /* (I) Underlying floater            */
             double      *FltIndex,     /* (I) Last/current index rate       */
             double      *FixZeroToPmt, /* (I) Zero w/mat on next fix coupon */
             double      *FltZeroToPmt, /* (I) Zero w/mat on next flt coupon */
             long        FixCpnFlag,    /* (I) Fix coupon flag               */
             double      FixCpnRate,    /* (I) Fix coupon rate               */
             double      FixCpnOuts,    /* (I) Fix coupon outstanding        */
             double      FixCpnDcf,     /* (I) Fix coupon dcf                */
             long        FixCpnAccSt,   /* (I) Fix coupon accrued start      */
             char        FixDCConv,     /* (I) Fix coupon day count conv     */
             long        FltResetFlag,  /* (I) Flt coupon flag               */
             double      FltCpnOuts,    /* (I) Flt coupon outstanding        */
             double      FltCpnDcf,     /* (I) Flt coupon dcf                */
             long        FltCpnAccSt,   /* (I) Flt coupon accrued start      */
             char        FltDCConv,     /* (I) Flt coupon day count conv     */
             double      FloatSpd,      /* (I) Floating leg spread           */
             char        CoS,           /* (I) Flt coupon compound flag      */
             long        ExerFlag,      /* (I) Exercise flag                 */
             double      Strike,        /* (I) Current strike                */
             char        OptStub,       /* (I) Option stub rule              */
             char        ArreasReset,   /* (I) 'Y' if set-in-arreas          */
             long        CurrentDate,   /* (I) Current date                  */
             int         t,             /* (I) Current time period           */
             int         T,             /* (I) Total number of period        */
             int         DCurve,        /* (I) Discount curve                */
             FIX3_DEV_DATA    *dev_data,     /* (I) Fix3_Dev data structure            */
             FIX3_TREE_DATA   *tree_data)    /* (I) Tree data structure           */
{

    double  *SwapL;                  /* Local slice pointer    */
    double  *BondL;
    double  *FloaterL;
    double  *FltIndexL;
    double  *FixZeroToPmtL;
    double  *FltZeroToPmtL;

    double  FixAccrued = 0.0;        /* Fix coupon accrued at current date */
    double  FltAccrued = 0.0;        /* Flt coupon accrued at current date */

    int     Top1, Bottom1;           /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)  */

    int     i, j, k;                 /* Node indices           */
    int     offset;                  /* Node offset            */
    int     status = FAILURE;        /* Error status           */


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
        if (Fix3_Dev (Swap,
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
     * Calculate dirty (without accrueds) price of the swap
     */
    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL    = Swap    + offset;
        BondL    = Bond    + offset;
        FloaterL = Floater + offset;
    
        for (i = Bottom1; i <= Top1; i ++)
        {
            SwapL[i] = BondL[i] - FloaterL[i] - Strike;
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            SwapL    = Swap    + offset;
            BondL    = Bond    + offset;
            FloaterL = Floater + offset;
    
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                SwapL[j] = BondL[j] - FloaterL[j] - Strike;
            }
        }  /* for j */
    }
    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                SwapL    = Swap    + offset;
                BondL    = Bond    + offset;
                FloaterL = Floater + offset;
    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    SwapL[k] = BondL[k] - FloaterL[k] - Strike;
                }
            }  /* for j */
    }  /* if then else */

    /* 
     *  Calculate fix accrued interests if we are within accrual period.
     */             
    if (FixCpnAccSt <= CurrentDate)    /* no accureds for fwd starting swap */
    {   
        if (DrDayCountFraction (FixCpnAccSt,
                                CurrentDate,                                 
                                FixDCConv,
                                &FixAccrued) == FAILURE)
        {
            DR_Error("Fix3_AmortSwap_t: unable to calculate dcf for fix accrued!");
            goto RETURN;
        }        
    }  /* if */
  
    if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL         = Swap         + offset;
        FixZeroToPmtL = FixZeroToPmt + offset;
    
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
            }
        }                                                           
        else if (OptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixAccrued * FixCpnOuts;
            }
        }
        else if (OptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixAccrued * FixCpnOuts 
                          * FixZeroToPmtL[i];
            }
        }  /* if then else */                        
    }
    else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
    {
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
                }
            }  /* for i */
        }                                                           
        else if (OptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts;
                }
            }  /* for i */
        }
        else if (OptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL         = Swap         + offset;
                FixZeroToPmtL = FixZeroToPmt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts * FixZeroToPmtL[j];
                }
            }  /* for i */
        }  /* if then else */                        
    }
    else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
    {
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
                    }
                }  /* for j */
        }                                                           
        else if (OptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixAccrued * FixCpnOuts;
                    }
                }  /* for j */
        }
        else if (OptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL         = Swap         + offset;
                    FixZeroToPmtL = FixZeroToPmt + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixAccrued * FixCpnOuts 
                                        * FixZeroToPmtL[k];
                    }
                }  /* for j */
        }  /* if then else */                        
    }  /* if then else on nb factors*/

    /*
     *  Calculate flt accrued interests if we are within accrual period.
     */
    if (FltCpnAccSt <= CurrentDate)
    {
        if (DrDayCountFraction (FltCpnAccSt,
                                CurrentDate,                                 
                                FltDCConv,
                                &FltAccrued) == FAILURE)
        {
            DR_Error("Fix3_AmortSwap_t: unable to calculate dcf for flt accrued!");
            goto RETURN;
        }        
    }  /* if */
 

    if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 1))
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL         = Swap         + offset;
        FltIndexL     = FltIndex     + offset;
        FltZeroToPmtL = FltZeroToPmt + offset;
    
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += (FltIndexL[i] + FloatSpd) * FltCpnOuts * FltCpnDcf;
                    }
                }                          /* Swap,Bond stub : add PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += (FltIndexL[i] + FloatSpd) * FltCpnOuts * FltAccrued;
                    }
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += (FltIndexL[i] + FloatSpd) * FltCpnOuts 
                            * (FltAccrued - FltCpnDcf * FltZeroToPmtL[i]);
                    }
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] -= (FltIndexL[i] + FloatSpd) * FltCpnOuts 
                            * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[i];
                    }
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] -= (FltIndexL[i] + FloatSpd) * FltCpnOuts 
                            * FltCpnDcf * FltZeroToPmtL[i];
                    }
                }  /* if then else */ 
            }  /* if then else */
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += (pow(1. + (FltIndexL[i] + FloatSpd), FltCpnDcf) -1.)
                                  * FltCpnOuts;
                    }
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 2))
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += (FltIndexL[j] + FloatSpd) * FltCpnOuts * FltCpnDcf;
                        }
                    }  /* for i */
                }                          /* Swap, Bond stub : add PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += (FltIndexL[j] + FloatSpd) * FltCpnOuts * FltAccrued;
                        }
                    }  /* for i */
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += (FltIndexL[j] + FloatSpd) * FltCpnOuts 
                                 * (FltAccrued - FltCpnDcf * FltZeroToPmtL[j]);
                        }
                    }  /* for i */
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] -= (FltIndexL[j] + FloatSpd) * FltCpnOuts 
                                * (FltCpnDcf - FltAccrued) *  FltZeroToPmtL[j];
                        }
                    }  /* for i */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] -= (FltIndexL[j] + FloatSpd) * FltCpnOuts 
                                            * FltCpnDcf * FltZeroToPmtL[j];
                        }
                    }  /* for i */
                }  /* if then else */ 
            } 
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += (pow(1. + (FltIndexL[j] + FloatSpd), FltCpnDcf) -1.)
                                         * FltCpnOuts;
                        }
                    }  /* for i */
                }       
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 3))
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += (FltIndexL[k] + FloatSpd) * FltCpnOuts * FltCpnDcf;
                            }
                        }  /* for j */
                }                          /* Swap, Bond stub : add PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += (FltIndexL[k] + FloatSpd) * FltCpnOuts * FltAccrued;
                            }
                        }  /* for j */
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += (FltIndexL[k] + FloatSpd) * FltCpnOuts 
                                    * (FltAccrued - FltCpnDcf * FltZeroToPmtL[k]);
                            }
                        }  /* for j */
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] -= (FltIndexL[k] + FloatSpd) * FltCpnOuts 
                                    * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[k];
                            }
                        }  /* for j */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] -= (FltIndexL[k] + FloatSpd) * FltCpnOuts 
                                                * FltCpnDcf * FltZeroToPmtL[k];
                            }
                        }  /* for j */
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += (pow(1. + (FltIndexL[k] + FloatSpd), FltCpnDcf) -1.)
                                                * FltCpnOuts;
                            }
                        }  /* for j */
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_AmortSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound */
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_AmortSwap_t */



/*****  Fix3_Swaplet_t  ***********************************************************/
/*
*       Single swaplet.
*/
int   Fix3_Swaplet_t (double      *Swaplet,    /* (I/O) Swaplet                   */
                 double      *Index,      /* (I) Floating index              */
                 double      *Zero,       /* (I) Zero to next payment date   */
                 long        SwapletFlag, /* (I) Caplet reset flag           */
                 double      FixAmount,   /* (I) Fix payment                 */
                 double      DayCntFtn,   /* (I) Day count fraction          */
                 char        CoS,         /* (I) 'C'ompound or 'S'imple rate */
                 double      FloatSpd,    /* (I) Floating leg spread         */
                 char        Arrears,     /* (I) 'Y' if reset in arrears     */
                 double      Notional,    /* (I) Caplet notional             */
                 int         t,           /* (I) Current time period         */
                 int         T,           /* (I) Total number of period      */
                 int         DCurve,      /* (I) Discount curve              */
                 FIX3_DEV_DATA    *dev_data,   /* (I) Fix3_Dev data structure          */
                 FIX3_TREE_DATA   *tree_data)  /* (I) Tree data structure         */
{

    double  *SwapletL;               /* Local slice pointer    */
    double  *IndexL;
    double  *ZeroL;        

    int     Top1, Bottom1;           /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)  */

    int     i, j, k;                 /* Node indices           */
    int     offset;                  /* Node offset            */
    int     status = FAILURE;        /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Fix3_Dev (   Swaplet,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
                    
    }  /* if */

    
    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapletL = Swaplet + offset;
        IndexL   = Index   + offset;
        ZeroL    = Zero    + offset;

        /* If there is a fixing we add a swaplet */
        if (SwapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapletL[i] = FixAmount - Notional * DayCntFtn * (IndexL[i] + FloatSpd);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapletL[i] = ZeroL[i] * (FixAmount - Notional * DayCntFtn * (IndexL[i] + FloatSpd));
                    }
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapletL[i] = FixAmount - Notional * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) -1.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapletL[i] = ZeroL[i] * (FixAmount - Notional * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) - 1.));
                    }
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (SwapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapletL = Swaplet + offset;
                        IndexL   = Index   + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapletL[j] = FixAmount - Notional * DayCntFtn * (IndexL[j] + FloatSpd);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapletL = Swaplet + offset;
                        IndexL   = Index   + offset;
                        ZeroL    = Zero    + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapletL[j] = ZeroL[j] * (FixAmount - Notional * DayCntFtn * (IndexL[j] + FloatSpd));
                        }
                    }  /* for i */
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapletL = Swaplet + offset;
                        IndexL   = Index   + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapletL[j] = FixAmount - Notional * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) - 1.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapletL = Swaplet + offset;
                        IndexL   = Index   + offset;
                        ZeroL    = Zero    + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapletL[j] = ZeroL[j] * (FixAmount - Notional * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) - 1.));
                        }
                    }  /* for i */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (SwapletFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapletL = Swaplet + offset;
                            IndexL   = Index   + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                SwapletL[k] = FixAmount - Notional * DayCntFtn * (IndexL[k] + FloatSpd);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapletL = Swaplet + offset;
                            IndexL   = Index   + offset;
                            ZeroL    = Zero    + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                SwapletL[k] = ZeroL[k] * (FixAmount - Notional * DayCntFtn * (IndexL[k] + FloatSpd));
                            }
                        }  /* for j */
                }  /* if then else */
            }
            else
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapletL = Swaplet + offset;
                            IndexL   = Index   + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                SwapletL[k] = FixAmount - Notional * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) - 1.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapletL = Swaplet + offset;
                            IndexL   = Index   + offset;
                            ZeroL    = Zero    + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                SwapletL[k] = ZeroL[k] * (FixAmount - Notional * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) - 1.));
                            }
                        }  /* for j */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Swaplet_t */




/*****  Fix3_FwdSwap_t  ***********************************************************/
/*
*       Forward swap price. Floating leg is arbitrary. Same as Fix3_AmortSwap_t
*       except it allows different stub conventions for the fix and the flt
*       leg. Flt stub can be 'P'ar as well.
*/
int     Fix3_FwdSwap_t (
             double      *Swap,         /* (O) Forward swap                  */
             double      *Bond,         /* (I) Underlying bond               */
             double      *Floater,      /* (I) Underlying floater            */
             double      *FltIndex,     /* (I) Last/current index rate       */
             double      *FixZeroToPmt, /* (I) Zero w/mat on next fix coupon */
             double      *FltZeroToPmt, /* (I) Zero w/mat on next flt coupon */
             double      *IFltZeroToPmt,/* (I) Zero for par stub calculation */
             long        FixCpnFlag,    /* (I) Fix coupon flag               */
             double      FixCpnRate,    /* (I) Fix coupon rate               */
             double      FixCpnOuts,    /* (I) Fix coupon outstanding        */
             double      FixCpnDcf,     /* (I) Fix coupon dcf                */
             long        FixCpnAccSt,   /* (I) Fix coupon accrued start      */
             char        FixDCConv,     /* (I) Fix coupon day count conv     */
             long        FltResetFlag,  /* (I) Flt coupon flag               */
             double      FltCpnOuts,    /* (I) Flt coupon outstanding        */
             double      FltCpnDcf,     /* (I) Flt coupon dcf                */
             long        FltCpnAccSt,   /* (I) Flt coupon accrued start      */
             char        FltDCConv,     /* (I) Flt coupon day count conv     */
             char        CoS,           /* (I) Flt coupon compound flag      */
             long        ExerFlag,      /* (I) Exercise flag                 */
             double      Strike,        /* (I) Current strike                */
             char        FixOptStub,    /* (I) Fix option stub rule          */
             char        FltOptStub,    /* (I) Flt option stub rule          */
             char        ArreasReset,   /* (I) 'Y' if set-in-arreas          */
             long        CurrentDate,   /* (I) Current date                  */
             int         t,             /* (I) Current time period           */
             int         T,             /* (I) Total number of period        */
             int         DCurve,        /* (I) Discount curve                */
             FIX3_DEV_DATA    *dev_data,     /* (I) Fix3_Dev data structure            */
             FIX3_TREE_DATA   *tree_data)    /* (I) Tree data structure           */
{

    double  *SwapL;                  /* Local slice pointer    */
    double  *BondL;
    double  *FloaterL;
    double  *FltIndexL;
    double  *FixZeroToPmtL;
    double  *FltZeroToPmtL;
    double  *IFltZeroToPmtL;

    double  FixAccrued = 0.0;        /* Fix coupon accrued at current date */
    double  FltAccrued = 0.0;        /* Flt coupon accrued at current date */

    int     Top1, Bottom1;           /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;         /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;       /* Tree limits (3rd dim)  */

    int     i, j, k;                 /* Node indices           */
    int     offset;                  /* Node offset            */
    int     status = FAILURE;        /* Error status           */


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
        if (Fix3_Dev (Swap,
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

    /* basic checks */
    if ((Swap == NULL) || (Bond == NULL) || (Floater == NULL)) goto RETURN;

    if ((FixCpnAccSt < CurrentDate) && (!FixCpnFlag) && /* has fix stub and */
        (FixOptStub == 'S'))                            /* is swap stub     */
    {
        if (FixZeroToPmt == NULL) goto RETURN;
    }

    if (FltCpnAccSt < CurrentDate)
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)
                {
                    if (FltIndex == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Flt Index is not supplied.");
                        goto RETURN;
                    }
                }
                else if (FltOptStub == 'P')
                {
                    if (FltIndex      == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Flt Index is not supplied.");
                        goto RETURN;
                    }
                    if (IFltZeroToPmt == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Idx Zero is not supplied.");
                        goto RETURN;
                    }
                    if (FltZeroToPmt  == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Pmt Zero is not supplied.");
                        goto RETURN;
                    }
                }
                else if (FltOptStub == 'B')
                {;
                }                             
                else if (FltOptStub == 'S')
                {                          
                    if (FltIndex      == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Flt Index is not supplied.");
                        goto RETURN;
                    }
                }
            }
            else 
            {
                if (FltResetFlag)
                {;
                }
                else if (FltOptStub == 'P')
                {
                    if (IFltZeroToPmt == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Idx Zero is not supplied.");
                        goto RETURN;
                    }
                    if (FltZeroToPmt  == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Pmt Zero is not supplied.");
                        goto RETURN;
                    }
                }
                else
                {
                    if (FltIndex      == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Flt Index is not supplied.");
                        goto RETURN;
                    }
                    if (FltZeroToPmt  == NULL)
                    {
                        DR_Error("Fix3_FwdSwap_t: Pmt Zero is not supplied.");
                        goto RETURN;
                    }
                }
            }
        }
        else
        {
            if ((ArreasReset == 'Y') && (FltResetFlag))
            {
                if (FltIndex == NULL)
                {
                    DR_Error("Fix3_FwdSwap_t: Flt Index is not supplied.");
                    goto RETURN;
                }
            }
        } /* if simple or compound flag */
    }


    /*
     * Calculate dirty (without accrued) price of the swap
     */
    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL    = Swap    + offset;
        BondL    = Bond    + offset;
        FloaterL = Floater + offset;
    
        for (i = Bottom1; i <= Top1; i ++)
        {
            SwapL[i] = BondL[i] - FloaterL[i] - Strike;
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            SwapL    = Swap    + offset;
            BondL    = Bond    + offset;
            FloaterL = Floater + offset;
    
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                SwapL[j] = BondL[j] - FloaterL[j] - Strike;
            }
        }  /* for j */
    }
    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                SwapL    = Swap    + offset;
                BondL    = Bond    + offset;
                FloaterL = Floater + offset;
    
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    SwapL[k] = BondL[k] - FloaterL[k] - Strike;
                }
            }  /* for j */
    }  /* if then else */

    /* 
     *  Calculate fix accrued interests if we are within accrual period.
     */             
    if (FixCpnAccSt < CurrentDate)    /* no accrued for fwd starting swap */
    {   
        if (DrDayCountFraction (FixCpnAccSt,
                                CurrentDate,                                 
                                FixDCConv,
                                &FixAccrued) == FAILURE)
        {
            DR_Error("Fix3_FwdSwap_t: unable to calculate dcf for fix accrued!");
            goto RETURN;
        }        
    }  /* if */
  
    if ((FixCpnAccSt < CurrentDate) && (tree_data->NbFactor == 1))
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL         = Swap         + offset;
        FixZeroToPmtL = FixZeroToPmt + offset;
    
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
            }
        }                                                           
        else if (FixOptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixAccrued * FixCpnOuts;
            }
        }
        else if (FixOptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
            {
                SwapL[i] -= FixCpnRate * FixAccrued * FixCpnOuts 
                          * FixZeroToPmtL[i];
            }
        }  /* if then else */                        
    }
    else if ((FixCpnAccSt < CurrentDate) && (tree_data->NbFactor == 2))
    {
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
                }
            }  /* for i */
        }                                                           
        else if (FixOptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL = Swap + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts;
                }
            }  /* for i */
        }
        else if (FixOptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                SwapL         = Swap         + offset;
                FixZeroToPmtL = FixZeroToPmt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts * FixZeroToPmtL[j];
                }
            }  /* for i */
        }  /* if then else */                        
    }
    else if ((FixCpnAccSt < CurrentDate) && (tree_data->NbFactor == 3))
    {
        if (FixCpnFlag)           /* Subtract fix cpn just made  */
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixCpnOuts * FixCpnDcf;
                    }
                }  /* for j */
        }                                                           
        else if (FixOptStub == 'B')  /* Bond stub : subtract accrued */
        {                                                                   
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL = Swap + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixAccrued * FixCpnOuts;
                    }
                }  /* for j */
        }
        else if (FixOptStub == 'S')  /* Swap stub: subtract PVed accrued  */
        {                                                            
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    SwapL         = Swap         + offset;
                    FixZeroToPmtL = FixZeroToPmt + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        SwapL[k] -= FixCpnRate * FixAccrued * FixCpnOuts 
                                        * FixZeroToPmtL[k];
                    }
                }  /* for j */
        }  /* if then else */                        
    }  /* if then else on nb factors*/

    /*
     *  Calculate flt accrued interests if we are within accrual period.
     */
    if (FltCpnAccSt < CurrentDate)
    {
        if (DrDayCountFraction (FltCpnAccSt,
                                CurrentDate,                                 
                                FltDCConv,
                                &FltAccrued) == FAILURE)
        {
            DR_Error("Fix3_FwdSwap_t: unable to calculate dcf for flt accrued!");
            goto RETURN;
        }        
    }  /* if */
 

    if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 1))
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        SwapL          = Swap          + offset;
        FltIndexL      = FltIndex      + offset;
        FltZeroToPmtL  = FltZeroToPmt  + offset;
        IFltZeroToPmtL = IFltZeroToPmt + offset;
    
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += FltIndexL[i] * FltCpnOuts * FltCpnDcf;
                    }
                }
                else if (FltOptStub == 'P')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += FltCpnOuts * 
                                    (FltIndexL[i] * FltCpnDcf -
                                     (1/IFltZeroToPmtL[i]-1) * 
                                     FltZeroToPmtL[i]);
                    }                    
                }
                else if (FltOptStub == 'B')
                {
                    DR_Error("Fix3_FwdSwap_t: Bond stub not supported "
                             "for arrears reset!");
                    goto RETURN;
                }                             
                else if (FltOptStub == 'S')  /* Swap stub : add PVed accrued */
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += FltIndexL[i] * FltCpnOuts * FltAccrued;
                    }
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }
                else if (FltOptStub == 'P')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] -= FltCpnOuts * 
                                    (1/IFltZeroToPmtL[i]-1) * 
                                    FltZeroToPmtL[i];
                    }
                }
                else if (FltOptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += FltIndexL[i] * FltCpnOuts 
                            * (FltAccrued - FltCpnDcf * FltZeroToPmtL[i]);
                    }
                }
                else if (FltOptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] -= FltIndexL[i] * FltCpnOuts 
                            * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[i];
                    }
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] -= FltIndexL[i] * FltCpnOuts 
                            * FltCpnDcf * FltZeroToPmtL[i];
                    }
                }  /* if then else */ 
            }  /* if then else */
        }
        else /* pmt convention is compounding */
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        SwapL[i] += (pow(1. + FltIndexL[i], FltCpnDcf) -1.)
                                  * FltCpnOuts;
                    }
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 2))
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += FltIndexL[j] * FltCpnOuts * FltCpnDcf;
                        }
                    }  /* for i */
                }
                else if (FltOptStub == 'P')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL          = Swap          + offset;
                        FltIndexL      = FltIndex      + offset;
                        FltZeroToPmtL  = FltZeroToPmt  + offset;
                        IFltZeroToPmtL = IFltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += FltCpnOuts * 
                                        (FltIndexL[j] * FltCpnDcf -
                                         (1/IFltZeroToPmtL[j]-1) *
                                         FltZeroToPmtL[j]);
                        }
                    }  /* for i */
                }
                else if (FltOptStub == 'B')
                {
                    DR_Error("Fix3_FwdSwap_t: Bond stub not supported "
                             "for arrears reset!");
                    goto RETURN;
                }                             
                else if (FltOptStub == 'S') /* Swap stub : add PVed accrued */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += FltIndexL[j] * FltCpnOuts * FltAccrued;
                        }
                    }  /* for i */
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }
                else if (FltOptStub == 'P') 
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL          = Swap          + offset;
                        FltIndexL      = FltIndex      + offset;
                        FltZeroToPmtL  = FltZeroToPmt  + offset;
                        IFltZeroToPmtL = IFltZeroToPmt + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] -= FltCpnOuts * 
                                        (1/IFltZeroToPmtL[j]-1) * 
                                        FltZeroToPmtL[j];
                        }
                    }  /* for i */
                }  /* if then else */ 
                else if (FltOptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += FltIndexL[j] * FltCpnOuts 
                                 * (FltAccrued - FltCpnDcf * FltZeroToPmtL[j]);
                        }
                    }  /* for i */
                }
                else if (FltOptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] -= FltIndexL[j] * FltCpnOuts 
                                * (FltCpnDcf - FltAccrued) *  FltZeroToPmtL[j];
                        }
                    }  /* for i */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL         = Swap         + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] -= FltIndexL[j] * FltCpnOuts 
                                            * FltCpnDcf * FltZeroToPmtL[j];
                        }
                    }  /* for i */
                }  /* if then else */ 
            } 
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        SwapL     = Swap     + offset;
                        FltIndexL = FltIndex + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            SwapL[j] += (pow(1. + FltIndexL[j], FltCpnDcf) -1.)
                                         * FltCpnOuts;
                        }
                    }  /* for i */
                }       
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 3))
    {
        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += FltIndexL[k] * FltCpnOuts * FltCpnDcf;
                            }
                        }  /* for j */
                }                          
                else if (FltOptStub == 'P')  
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL          = Swap          + offset;
                            FltIndexL      = FltIndex      + offset;
                            FltZeroToPmtL  = FltZeroToPmt  + offset;
                            IFltZeroToPmtL = IFltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += FltCpnOuts * 
                                            (FltIndexL[k] * FltCpnDcf -
                                             (1/IFltZeroToPmtL[k]-1) *
                                             FltZeroToPmtL[k]);
                            }
                        }  /* for j */
                }
                else if (FltOptStub == 'B')
                {
                    DR_Error("Fix3_FwdSwap_t: Bond stub not supported "
                             "for arrears reset!");
                    goto RETURN;
                }                             
                else if (FltOptStub == 'S')  /* Swap stub : add PVed accrued */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += FltIndexL[k] * FltCpnOuts * FltAccrued;
                            }
                        }  /* for j */
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }
                else if (FltOptStub == 'P')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL          = Swap          + offset;
                            FltIndexL      = FltIndex      + offset;
                            FltZeroToPmtL  = FltZeroToPmt  + offset;
                            IFltZeroToPmtL = IFltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] -= FltCpnOuts * 
                                            (1/IFltZeroToPmtL[k]-1) * 
                                            FltZeroToPmtL[k];
                            }
                        }  /* for j */
                }
                else if (FltOptStub == 'B')  /* Bond stub : add accrued     */
                {                         /*             pay PVed coupon */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += FltIndexL[k] * FltCpnOuts 
                                    * (FltAccrued - FltCpnDcf * FltZeroToPmtL[k]);
                            }
                        }  /* for j */
                }
                else if (FltOptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] -= FltIndexL[k] * FltCpnOuts 
                                    * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[k];
                            }
                        }  /* for j */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL         = Swap         + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] -= FltIndexL[k] * FltCpnOuts 
                                                * FltCpnDcf * FltZeroToPmtL[k];
                            }
                        }  /* for j */
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Add flt cpn just made */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            SwapL     = Swap     + offset;
                            FltIndexL = FltIndex + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                SwapL[k] += (pow(1. + FltIndexL[k], FltCpnDcf) -1.)
                                                * FltCpnOuts;
                            }
                        }  /* for j */
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }        
            }
            else
            {
                if (FltResetFlag)          /* Floater is clean priced on reset date */ 
                {
                    ;
                }
                else                       /* Stub : Not supported */
                {
                    DR_Error("Fix3_FwdSwap_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound */
    }  /* if then else */


    status = SUCCESS;

RETURN:

    return (status);

}  /* Fix3_FwdSwap_t */



/*****  Fix3_FwdSwapFlows_t  ***********************************************/
/*
 *      Calculates a general forward swap value. Used in the version of
 *      fwd swap with arbitrary resets and payments.
 *      NOTE: This is only a payoff function
 */

int     Fix3_FwdSwapFlows_t(
                double     *Swap,             /* (O) Forward swap            */
                /* zero banks */
                CLAIM_BANK *IdxZBK,           /* (I) ZeroBk for pmt idx      */
                CLAIM_BANK *DiscZBK,          /* (I) ZeroBk for discounting  */
                /* fix leg stub details */
                double     *FullBond,         /* (I) Underlying full bond    */
                double      FixStubNotional,  /* (I) prn for fix stub prd    */
                double      FixStubRate,      /* (I) cpm for fix stub prd    */
                double      FixStubDcf,       /* (I) dcf of fix stub prd     */
                long        FixStubAccStart,  /* (I) fix stub prd acc start  */
                long        FixStubAccEnd,    /* (I) fix stub prd acc end    */
                long        FixStubPmtDate,   /* (I) fix stub prd pmt date   */
                char        FixDCC,           /* (I) fix stub DCC            */
                char        FixStubConv,      /* (I) fix stub convention     */
                /* flt leg stub details */
                double     *FullFloater,      /* (I) Underlying full floater */
                double      FltStubNotional,  /* (I) prn for flt stub prd    */
                double     *FltIndex,         /* (I) flt idex for stub prd   */
                double      FltStubSpread,    /* (I) spread for flt stub prd */
                double     *FltCoupon,        /* (I) pv of flt stub prd cpn  */
                double      FltStubDcf,       /* (I) flt stub prd dcf        */
                long        FltStubAccStart,  /* (I) flt stub prd acc start  */
                long        FltStubAccEnd,    /* (I) flt stub prd acc end    */
                long        FltStubResetEff,  /* (I) flt stub prd reset date */
                long        FltStubPmtDate,   /* (I) flt stub prd pmt date   */
                char        FltDCC,           /* (I) flt stub DCC            */
                char        CompFloat,        /* (I) Simple or Compounding   */
                char        FltStubConv,      /* (I) flt stub convention     */
                /* exercise details */
                double      Strike,           /* (I) Strike amount           */
                long        CurrentDate,      /* (I) date of current timept  */
                int         t,                /* (I) current timept          */
                FIX3_TREE_DATA  *tree_data)        /* (I) tree data               */
{
    /* --------------------------------------------- */

#undef  FULL_SWAP
#define FULL_SWAP(x) (SwapL[x] = FullBondL[x] - FullFloaterL[x] - Strike)

    /* --------------------------------------------- */

#undef  FIXSTUB
#define FIXSTUB(x)                                                          \
    {                                                                       \
        switch (FixStubConv)                                                \
        {                                                                   \
        case 'B':                                                           \
            SwapL[x] -= FixStubNotional * FixStubRate * FixAccruedDCF;      \
            break;                                                          \
        case 'S':                                                           \
            SwapL[x] += FixStubNotional * FixStubRate *                     \
                        (FixRemainDCF - FixStubDcf) * FixZeroL[x];          \
            break;                                                          \
        case 'N':                                                           \
            break;                                                          \
        default:                                                            \
            goto RETURN;                                                    \
        }                                                                   \
    }

    /* --------------------------------------------- */

#undef  FLTSTUB
#define FLTSTUB(x)                                                          \
    {                                                                       \
        switch (FltStubConv)                                                \
        {                                                                   \
        case 'B':                                                           \
            SwapL[x] += FltStubNotional * (FltIndexL[x] + FltStubSpread) *  \
                        FltAccruedDCF;                                      \
            break;                                                          \
        case 'S':                                                           \
            SwapL[x] -= FltCouponL[x] * (FltRemainDCF/FltStubDcf - 1.0);    \
            break;                                                          \
        case 'N':                                                           \
            break;                                                          \
        case 'P':                                                           \
            SwapL[x] -= FltStubNotional * (1.0/PZeroL[x]-1.0) *             \
                        FltZeroL[x] - FltCouponL[x];                        \
            break;                                                          \
        default:                                                            \
            goto RETURN;                                                    \
        }                                                                   \
    }

    /* --------------------------------------------- */

    int      status = FAILURE;
    int      hasFixStub;
    int      hasFltStub;
    double   FixAccruedDCF = 0.0; /* initialized to avoid compiler warning */
    double   FixRemainDCF  = 0.0;
    double   FltAccruedDCF = 0.0;
    double   FltRemainDCF  = 0.0;
    int      offset;
    int      i, j, k;

    double  *FixZero = NULL;
    double  *FltZero = NULL;
    double  *PZero   = NULL;

    double  *SwapL;
    double  *FullBondL;
    double  *FullFloaterL;
    double  *FixZeroL;
    double  *FltZeroL;
    double  *FltIndexL;
    double  *FltCouponL;
    double  *PZeroL;

    int      Top1,   Bottom1;       /* Tree limits (1st dim)  */
    int     *Top2,  *Bottom2;       /* Tree limits (2nd dim)  */
    int    **Top3, **Bottom3;       /* Tree limits (3rd dim)  */

    /* basic check */
    if ((Swap      == NULL) || (FullFloater == NULL) || 
        (FullBond  == NULL) || (tree_data   == NULL)) goto RETURN;

    /* set tree parameters */
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /* set stub flags */
    hasFixStub = (FixStubAccStart < CurrentDate);
    hasFltStub = (FltStubAccStart < CurrentDate);

    /* perform stub checks */
    if (hasFltStub)
    {
        if ((CompFloat == 'C') && (!(FltStubConv == 'N')))
        {
            DR_Error("Fix3_FwdSwapFlows_t: "
                     "cannot handle float stub with compounding payment");
            goto RETURN;
        }

        if ((FltStubConv == 'B') && (FltStubResetEff > CurrentDate))
        {
            DR_Error("Fix3_FwdSwapFlows_t: "
                     "bond stub not allowed when reset occurs "
                     "after exercise date");
            goto RETURN;
        }

        if ((FltStubConv == 'B') && (FltIndex == NULL)) goto RETURN;

        if ((FltStubConv == 'S') || (FltStubConv == 'P'))
        {
            if (FltCoupon == NULL) goto RETURN;
        }
    }


    /* find the stub quantities */
    if (hasFixStub)
    {
        if (FixStubConv == 'S')
        {
            FixZero = Fix3_ZbkReadZero(DiscZBK,
                                  FixStubPmtDate,
                                  TRUE,
                                  CurrentDate,
                                  t,
                                  tree_data);
            if (FixZero == NULL) goto RETURN;
        }

        if (DrDayCountFraction(FixStubAccStart,
                               CurrentDate,
                               FixDCC,
                               &FixAccruedDCF) == FAILURE) goto RETURN;

        if (DrDayCountFraction(CurrentDate,
                               FixStubAccEnd,
                               FixDCC,
                               &FixRemainDCF) == FAILURE) goto RETURN;
    }

    if (hasFltStub)
    {
        if (FltStubConv == 'P')
        {
            FltZero = Fix3_ZbkReadZero(DiscZBK,
                                  FltStubPmtDate,
                                  TRUE,
                                  CurrentDate,
                                  t,
                                  tree_data);

            PZero   = Fix3_ZbkReadZero(IdxZBK,
                                  FltStubAccEnd,
                                  TRUE,
                                  CurrentDate,
                                  t,
                                  tree_data);

            if ((FltZero == NULL) || (PZero == NULL)) goto RETURN;
        }

        if (DrDayCountFraction(FltStubAccStart,
                               CurrentDate,
                               FltDCC,
                               &FltAccruedDCF) == FAILURE) goto RETURN;

        if (DrDayCountFraction(CurrentDate,
                               FltStubAccEnd,
                               FltDCC,
                               &FltRemainDCF) == FAILURE) goto RETURN;

    }

    /*************************** 1 Factor *****************************/

    if (tree_data->NbFactor == 1)
    {
        offset       = Fix3_Node_Offset(1, 0, 0, t, tree_data);
        SwapL        = Swap         + offset;
        FullBondL    = FullBond     + offset;
        FullFloaterL = FullFloater  + offset;
        FixZeroL     = FixZero      + offset;
        FltZeroL     = FltZero      + offset;
        FltIndexL    = FltIndex     + offset;
        FltCouponL   = FltCoupon    + offset;
        PZeroL       = PZero        + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
                FULL_SWAP(i);
                if (hasFixStub) FIXSTUB(i);
                if (hasFltStub) FLTSTUB(i);
        }  /* for i */
    }

    /*************************** 2 Factor *****************************/

    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset       = Fix3_Node_Offset(2, i, 0, t, tree_data);
            SwapL        = Swap         + offset;
            FullBondL    = FullBond     + offset;
            FullFloaterL = FullFloater  + offset;
            FixZeroL     = FixZero      + offset;
            FltZeroL     = FltZero      + offset;
            FltIndexL    = FltIndex     + offset;
            FltCouponL   = FltCoupon    + offset;
            PZeroL       = PZero        + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                FULL_SWAP(j);
                if (hasFixStub) FIXSTUB(j);
                if (hasFltStub) FLTSTUB(j);
            }  /* for j */
        }  /* for i */
    }

    /*************************** 3 Factor *****************************/

    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset       = Fix3_Node_Offset(3, i, j, t, tree_data);
                SwapL        = Swap         + offset;
                FullBondL    = FullBond     + offset;
                FullFloaterL = FullFloater  + offset;
                FixZeroL     = FixZero      + offset;
                FltZeroL     = FltZero      + offset;
                FltIndexL    = FltIndex     + offset;
                FltCouponL   = FltCoupon    + offset;
                PZeroL       = PZero        + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    FULL_SWAP(k);
                    if (hasFixStub) FIXSTUB(k);
                    if (hasFltStub) FLTSTUB(k);
                }  /* for k */
            }  /* for j */
        }  /* for i */
    }  /* if then else */

    status = SUCCESS;

RETURN:

#undef FULL_SWAP
#undef FIXSTUB
#undef FLTSTUB

    return (status);

} /* Fix3_FwdSwapFlows_t */

