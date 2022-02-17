/****************************************************************************/
/*      Calculation of floater price in the lattice.                        */
/****************************************************************************/
/*      FLOATER.c                                                           */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"




/*****  Floater_t  **********************************************************/
/*
*       Floater price with arbitrary index.
*/
int     Floater_t ( double      *Floater,      /* (I/O) Floater              */
                    double      *Index,        /* (I) Floating index         */
                    double      *Zero,         /* (I) Zero to next pay date  */
                    long        FixingFlag,    /* (I) TRUE if index refixes  */
                    double      DayCntFtn,     /* (I) Day count fraction     */
                    double      FloatSpd,      /* (I) Spread                 */
                    char        CoS,           /* (I) 'C'ompound or 'S'imple */
                    char        Arrears,       /* (I) 'Y' for arrears reset  */
                    double      Outstanding,   /* (I) Outstanding notional   */
                    double      Principal,     /* (I) Principal payment      */
                    long        PrincipalFlag, /* (I) TRUE if principal paid */
                    int         t,             /* (I) Current time point     */
                    int         T,             /* (I) Last time point        */
                    int         DCurve,        /* (I) Discount curve         */
                    DEV_DATA    *dev_data,     /* (I) Dev data structure     */
                    TREE_DATA   *tree_data)    /* (I) Tree data structure    */
{

    double  *FloaterL;                  /* Local slice pointers */
    double  *IndexL;
    double  *ZeroL;
        
    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   Floater,
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
        offset = Node_Offset(1, 0, 0, t, tree_data);

        FloaterL = Floater + offset;
        IndexL   = Index   + offset;
        ZeroL    = Zero    + offset;
    
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                FloaterL[i] += Principal;
            }
        }  /* if */


        if (FixingFlag)                                                         /* Add fixing */
        {
            if (CoS == 'S')                                                     /* Simple rate */
            {
                if (Arrears == 'Y')
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FloaterL[i] += Outstanding * DayCntFtn * (IndexL[i] + FloatSpd);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FloaterL[i] += Outstanding * ZeroL[i] * DayCntFtn * (IndexL[i] + FloatSpd);
                    }
                }  /* if then else */
            }
            else                                                                /* Compounded floating rate */
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FloaterL[i] += Outstanding * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) - 1.);
                    }
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FloaterL[i] += Outstanding * ZeroL[i] * (pow (1. + (IndexL[i] + FloatSpd), DayCntFtn) - 1.);
                    }
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                FloaterL = Floater + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    FloaterL[j] += Principal;
                }
            }  /* for i */
        }  /* if */


        if (FixingFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FloaterL = Floater + offset;
                        IndexL   = Index   + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FloaterL[j] += Outstanding * DayCntFtn * (IndexL[j] + FloatSpd);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FloaterL = Floater + offset;
                        IndexL   = Index   + offset;
                        ZeroL    = Zero    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FloaterL[j] += Outstanding * ZeroL[j] * DayCntFtn * (IndexL[j] + FloatSpd);
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
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FloaterL = Floater + offset;
                        IndexL   = Index   + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FloaterL[j] += Outstanding * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) - 1.);
                        }
                    }  /* for i */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FloaterL = Floater + offset;
                        IndexL   = Index   + offset;
                        ZeroL    = Zero    + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FloaterL[j] += Outstanding * ZeroL[j] * (pow (1. + (IndexL[j] + FloatSpd), DayCntFtn) - 1.);
                        }
                    }  /* for i */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (PrincipalFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    FloaterL = Floater + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        FloaterL[k] += Principal;
                    }
                }  /* for j */
        }  /* if */


        if (FixingFlag)
        {
            if (CoS == 'S')
            {
                if (Arrears == 'Y')                                         
                {                                                       
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FloaterL = Floater + offset;
                            IndexL   = Index   + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                FloaterL[k] += Outstanding * DayCntFtn * (IndexL[k] + FloatSpd);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FloaterL = Floater + offset;
                            IndexL   = Index   + offset;
                            ZeroL    = Zero    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                FloaterL[k] += Outstanding * ZeroL[k] * DayCntFtn * (IndexL[k] + FloatSpd);
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
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FloaterL = Floater + offset;
                            IndexL   = Index   + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                FloaterL[k] += Outstanding * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) - 1.);
                            }
                        }  /* for j */
                } 
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FloaterL = Floater + offset;
                            IndexL   = Index   + offset;
                            ZeroL    = Zero    + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                FloaterL[k] += Outstanding * ZeroL[k] * (pow (1. + (IndexL[k] + FloatSpd), DayCntFtn) - 1.);
                            }
                        }  /* for j */
                }  /* if then else */
            }  /* if then else */
        }  /* if */
    }  /* if then else */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Floater_t */



/*****  FwdFloater_t  ********************************************************/
/*
*       Forward floater price. Floating leg is arbitrary.
*/
int     FwdFloater_t (  
            double      *FwdFloater,   /* (O) Forward floater               */
            double      *Floater,      /* (I) Underlying floater            */
            double      *FltIndex,     /* (I) Last/current index rate       */
            double      *FltZeroToPmt, /* (I) Zero w/mat on next flt coupon */
            long        FltResetFlag,  /* (I) Flt coupon flag               */
            double      Outstanding,   /* (I) Flt coupon outstanding        */
            double      FltCpnDcf,     /* (I) Flt coupon dcf                */
            long        FltCpnAccSt,   /* (I) Flt coupon accrued start      */
            char        FltDCConv,     /* (I) Flt coupon day count conv     */
            double      FloatSpd,      /* (I) Floating leg spread           */
            char        CoS,           /* (I) Flt coupon compound flag      */
            long        ExerFlag,      /* (I) Exercise flag                 */
            char        OptStub,       /* (I) Option stub rule              */
            char        ArreasReset,   /* (I) 'Y' if set-in-arreas          */
            long        CurrentDate,   /* (I) Current date                  */
            int         t,             /* (I) Current time point            */
            int         T,             /* (I) Last time point               */
            int         DCurve,        /* (I) Discount curve                */
            DEV_DATA    *dev_data,     /* (I) Dev data structure            */
            TREE_DATA   *tree_data)    /* (I) Tree data structure           */
{

    double  *FwdFloaterL;           /* Local slice pointers               */
    double  *FloaterL;
    double  *FltIndexL;
    double  *FltZeroToPmtL;
        
    double  FltAccrued =0.0;         /* Flt coupon accrued at current date */

    int     Top1, Bottom1;          /* Tree limits (1rst dim)             */
    int     *Top2, *Bottom2;        /* Tree limits (2nd dim)              */
    int     **Top3, **Bottom3;      /* Tree limits (3rd dim)              */

    int     i, j, k;                /* Node indices                       */
    int     offset;                 /* Node offset                        */
    int     status = FAILURE;       /* Error status                       */


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
        if (Dev (FwdFloater,/* Disc expd value of fwd flaoter starting on */
                 t,         /* last exercise date                         */
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
     * Calculate fwd floater
     */
    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        FwdFloaterL = FwdFloater + offset;
        FloaterL    = Floater    + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            FwdFloaterL[i] = FloaterL[i];
        }  
    }
    else if (tree_data->NbFactor == 2)
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Node_Offset(2, i, 0, t, tree_data);

            FwdFloaterL = FwdFloater + offset;
            FloaterL    = Floater    + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                FwdFloaterL[j] = FloaterL[j];
            }  
        }  /* for i */
    }
    else if (tree_data->NbFactor == 3)
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Node_Offset(3, i, j, t, tree_data);

                FwdFloaterL = FwdFloater + offset;
                FloaterL    = Floater    + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    FwdFloaterL[k] = FloaterL[k];
                } 
            }  /* for j */
    }

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
            DR_Error("FwdFloater_t: Unable to calculate dcf for flt accrued!");
            goto RETURN;
        }        
    }  /* if */
 

    if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        FwdFloaterL   = FwdFloater   + offset;
        FloaterL      = Floater      + offset;
        FltIndexL     = FltIndex     + offset;
        FltZeroToPmtL = FltZeroToPmt + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            FwdFloaterL[i] -= Outstanding;
        }    

        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Subtract flt cpn just made             */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] -= (FltIndexL[i] + FloatSpd) * FltCpnDcf
                                        * Outstanding;
                    }          
                }                          /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] -= (FltIndexL[i] + FloatSpd * FltZeroToPmtL[i]) 
                                        * FltAccrued * Outstanding ;
                    }  
                }  /* if then else */      /* Nothing to do for No-Accrued stub */
            }
            else 
            {
                if (FltResetFlag)          /* Nothing to do, next cpn just made */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub : subtract accrued       */
                {                         /*             add paid PVed coupon   */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] += (FltIndexL[i] + FloatSpd) * Outstanding 
                                        * (FltCpnDcf * FltZeroToPmtL[i] - FltAccrued);
                    }  
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] += (FltIndexL[i] + FloatSpd) * Outstanding 
                                        * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[i];
                    }  
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] += (FltIndexL[i] + FloatSpd) * Outstanding 
                                        * FltCpnDcf * FltZeroToPmtL[i];
                    }  
                }  /* if then else */ 
            }
        }
        else
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)         /* Subtract flt cpn just made          */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        FwdFloaterL[i] -= (pow(1. + (FltIndexL[i] + FloatSpd), FltCpnDcf) -1.)
                                        * Outstanding;
                    }  
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
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
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 2))
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Node_Offset(2, i, 0, t, tree_data);

            FwdFloaterL = FwdFloater + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                FwdFloaterL[j] -= Outstanding;
            }
        }  /* for i */

        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Subtract flt cpn just made             */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL = FwdFloater + offset;
                        FltIndexL   = FltIndex   + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] -= (FltIndexL[j] + FloatSpd) * FltCpnDcf
                                            * Outstanding;
                        }          
                    }  /* for i */
                }                          /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL   = FwdFloater   + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] -= (FltIndexL[j] + FloatSpd * FltZeroToPmtL[j]) 
                                            * FltAccrued * Outstanding ;
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
                else if (OptStub == 'B')  /* Bond stub : subtract accrued       */
                {                         /*             add paid PVed coupon   */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL   = FwdFloater   + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] += (FltIndexL[j] + FloatSpd) * Outstanding 
                                            * (FltCpnDcf * FltZeroToPmtL[j] - FltAccrued);
                        }  
                    }  /* for i */
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL   = FwdFloater   + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] += (FltIndexL[j] + FloatSpd) * Outstanding 
                                            * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[j];
                        }  
                    }  /* for i */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL   = FwdFloater   + offset;
                        FltIndexL     = FltIndex     + offset;
                        FltZeroToPmtL = FltZeroToPmt + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] += (FltIndexL[j] + FloatSpd) * Outstanding 
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
                if (FltResetFlag)         /* Subtract flt cpn just made          */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FwdFloaterL = FwdFloater + offset;
                        FltIndexL   = FltIndex   + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            FwdFloaterL[j] -= (pow(1. + (FltIndexL[j] + FloatSpd), FltCpnDcf) -1.)
                                            * Outstanding;
                        }  
                    }  /* for i */
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
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
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }
    else if ((FltCpnAccSt < CurrentDate) && (tree_data->NbFactor == 3))
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Node_Offset(3, i, j, t, tree_data);

                FwdFloaterL = FwdFloater + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    FwdFloaterL[k] -= Outstanding;
                }
            }  /* for j */

        if (CoS == 'S')
        {
            if (ArreasReset == 'Y')
            {
                if (FltResetFlag)          /* Subtract flt cpn just made             */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL = FwdFloater + offset;
                            FltIndexL   = FltIndex   + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] -= (FltIndexL[k] + FloatSpd) * FltCpnDcf
                                                * Outstanding;
                            }          
                        }  /* for j */
                }                          /* Swap,Bond stub : subtract PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL   = FwdFloater   + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] -= (FltIndexL[k] + FloatSpd * FltZeroToPmtL[k]) 
                                                * FltAccrued * Outstanding ;
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
                else if (OptStub == 'B')  /* Bond stub : subtract accrued       */
                {                         /*             add paid PVed coupon   */
                    
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL   = FwdFloater   + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] += (FltIndexL[k] + FloatSpd) * Outstanding 
                                                * (FltCpnDcf * FltZeroToPmtL[k] - FltAccrued);
                            }  
                        }  /* for j */
                }
                else if (OptStub == 'S')  /* Swap stub : add PVed (cpn-accrued) */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL   = FwdFloater   + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] += (FltIndexL[k] + FloatSpd) * Outstanding 
                                                * (FltCpnDcf - FltAccrued) * FltZeroToPmtL[k];
                            }  
                        }  /* for j */
                }
                else                      /* No-Accrued stub: add est.& PVed cpn */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL   = FwdFloater   + offset;
                            FltIndexL     = FltIndex     + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] += (FltIndexL[k] + FloatSpd) * Outstanding 
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
                if (FltResetFlag)         /* Subtract flt cpn just made          */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FwdFloaterL = FwdFloater + offset;
                            FltIndexL   = FltIndex   + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                FwdFloaterL[k] -= (pow(1. + (FltIndexL[k] + FloatSpd), FltCpnDcf) -1.)
                                                * Outstanding;
                            }  
                        }  /* for j */
                } 
                else                       /* Stub : Not supported */
                {
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
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
                    DR_Error("FwdFloater_t: accrued not supported for comp flt payment!");
                    goto RETURN;
                }
            }        
        } /* if simple or compound flag */
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* FwdFloater_t */


/*****  CmpFloater_t  *******************************************************/
/*
*       Compounded floater price with arbitrary index.
*/
int     CmpFloater_t (double      *CmpFloater,   /* (I/O) Floater           */
                      double      *Index,        /* (I) Floating index      */
                      double      *CmpInt,       /* (I) Cmp interest        */
                      long        AmortAbs,      /* (I) 'Y' if absolute amo */
                      long        FixingFlag,    /* (I) TRUE if refix       */
                      double      DayCntFtn,     /* (I) Day count fraction  */
                      double      FloatSpd,      /* (I) Spread              */
                      char        CoS,           /* (I) 'C'mp or 'S'imple   */
                      char        Arrears,       /* (I) 'Y' for arrears     */
                      double      Principal,     /* (I) Principal payment   */
                      long        PrincipalFlag, /* (I) TRUE if princ paid  */
                      int         t,             /* (I) Current time point  */
                      int         T,             /* (I) Last time point     */
                      int         DCurve,        /* (I) Discount curve      */
                      DEV_DATA    *dev_data,     /* (I) Dev data structure  */
                      TREE_DATA   *tree_data)    /* (I) Tree data structure */
{

    double  *CmpFloaterL;               /* Local slice pointers */
    double  *IndexL;
    double  *CmpIntL;
        
    double  PrincCI;

    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];



    if (Dev (   CmpFloater,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;                    
    }

    PrincCI = (AmortAbs == 'Y') ? Principal : 0.;

    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);
     
        CmpFloaterL = CmpFloater + offset;
        IndexL      = Index      + offset;
        CmpIntL     = CmpInt     + offset;
        
        if (Arrears == 'Y')
        {                                                       
            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CmpFloaterL[i] += Principal - PrincCI * CmpIntL[i];
                }
            }  /* if */
     
            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CmpIntL[i] *= 1. + DayCntFtn * (IndexL[i] + FloatSpd);
                    }
                }
                else                                                     /* Compounded floating rate */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CmpIntL[i] *= pow (1. + (IndexL[i] + FloatSpd), DayCntFtn);
                    }
                } 
            }  /* if Fixing ... */
        }
        else
        {
            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CmpIntL[i] *= 1. + DayCntFtn * (IndexL[i] + FloatSpd);
                    }
                }  /* if then else */
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CmpIntL[i] *= pow (1. + (IndexL[i] + FloatSpd), DayCntFtn);
                    }
                }  /* if then else */
            }
     
            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CmpFloaterL[i] += Principal - PrincCI * CmpIntL[i];
                }
            }  /* if */
        } /* if Arrears ... */
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Arrears == 'Y')
        {                                                       
            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);
                    
                    CmpFloaterL = CmpFloater + offset;
                    CmpIntL     = CmpInt     + offset;
    
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CmpFloaterL[j] += Principal - PrincCI * CmpIntL[j];
                    }
                }  /* for i */
            }  /* if */

            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CmpIntL = CmpInt + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CmpIntL[j] *= 1. + DayCntFtn * (IndexL[j] + FloatSpd);
                        }
                    }  /* for i */
                }
                else                                                     /* Compounded floating rate */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CmpIntL = CmpInt + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CmpIntL[j] *= pow (1. + (IndexL[j] + FloatSpd), DayCntFtn);
                        }
                    }  /* for i */
                } 
            }  /* if Fixing ... */
        }
        else
        {
            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CmpIntL = CmpInt + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CmpIntL[j] *= 1. + DayCntFtn * (IndexL[j] + FloatSpd);
                        }
                    }
                }  /* if then else */
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CmpIntL = CmpInt + offset;
                        IndexL  = Index  + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CmpIntL[j] *= pow (1. + (IndexL[j] + FloatSpd), DayCntFtn);
                        }
                    }
                }  /* if then else */
            }

            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);
                    
                    CmpFloaterL = CmpFloater + offset;
                    CmpIntL     = CmpInt     + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CmpFloaterL[j] += Principal - PrincCI * CmpIntL[j];
                    }
                }  /* for i */
            }  /* if */
        } /* if Arrears ... */
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Arrears == 'Y')
        {                                                       
            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
                        
                        CmpFloaterL = CmpFloater + offset;
                        CmpIntL     = CmpInt     + offset;
                        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CmpFloaterL[k] += Principal - PrincCI * CmpIntL[k];
                        }
                    }  /* for j */
            }  /* if */
            
            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CmpIntL = CmpInt + offset;
                            IndexL  = Index  + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CmpIntL[k] *= 1. + DayCntFtn * (IndexL[k] + FloatSpd);
                            }
                        }  /* for j */
                }
                else                                                     /* Compounded floating rate */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CmpIntL = CmpInt + offset;
                            IndexL  = Index  + offset;
                            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CmpIntL[k] *= pow (1. + (IndexL[k] + FloatSpd), DayCntFtn);
                            }
                        }  /* for j */
                } 
            }  /* if Fixing ... */
        }
        else
        {
            if (FixingFlag)                                              /* Account for fixing */
            {
                if (CoS == 'S')                                          /* Simple rate        */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CmpIntL = CmpInt + offset;
                            IndexL  = Index  + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CmpIntL[k] *= 1. + DayCntFtn * (IndexL[k] + FloatSpd);
                            }
                        }  /* for j */
                }
                else                                                     /* Compounded floating rate */
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)             
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CmpIntL = CmpInt + offset;
                            IndexL  = Index  + offset;
                            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CmpIntL[k] *= pow (1. + (IndexL[k] + FloatSpd), DayCntFtn);
                            }
                        }  /* for j */
                } 
            }  /* if Fixing ... */

            if (PrincipalFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
                        
                        CmpFloaterL = CmpFloater + offset;
                        CmpIntL     = CmpInt     + offset;
                        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CmpFloaterL[k] += Principal - PrincCI * CmpIntL[k];
                        }
                    }  /* for j */
            }  /* if */
        } /* if Arrears ... */
    }  /* if NbFactor ... */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* CmpFloater_t */
