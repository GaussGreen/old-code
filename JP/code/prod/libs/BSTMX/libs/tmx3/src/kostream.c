/****************************************************************************/
/*      Calculation of a knock-out underlying price in the lattice.         */
/****************************************************************************/
/*      KOSTREAM.c                                                          */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"

/*****  KoCap_t  ********************************************************/
/*
 *  Knock out cap function. 
 */
int     KoCap_t (double      *KoCap,       /* (O) Knock out cap              */
                 double      *KoPway,      /* (O) Knock out partway          */
                 double      *Caplet,      /* (I) Current caplet             */
                 double      *FltZeroToPmt,/* (I) Zero w/mat on next flt cpn */
                 double      *KoIndex,     /* (I) Knock out index            */
                 long        KoFlag,       /* (I) Knock out flag             */
                 double      LoBarrier,    /* (I) Lower barrier              */
                 double      HiBarrier,    /* (I) Higher barrier             */
                 double      Rebate,       /* (I) Rebate                     */
                 char        IoO,          /* (I) Knock-out 'I'n or 'O'ut    */
                 char        Smoothing,    /* (I) Smoothing ('Y' or 'N')     */
                 char        OptStub,      /* (I) Option stub rule           */
                 long        FltResetFlag, /* (I) Flt reset flag             */
                 long        FltCpnFlag,   /* (I) Flt payment flag           */
                 double      FltCpnDcf,    /* (I) Flt coupon dcf             */
                 long        FltCpnAccSt,  /* (I) Flt coupon accrued start   */
                 char        FltDCConv,    /* (I) Flt coupon day count conv  */
                 char        ArrearsReset, /* (I) 'Y' if set-in-arreas       */
                 long        CurrentDate,  /* (I) Current date               */
                 int         t,            /* (I) Current time point         */
                 int         T,            /* (I) Last time point            */
                 int         DCurve,       /* (I) Discount curve             */
                 DEV_DATA    *dev_data,    /* (I) Dev data structure         */
                 TREE_DATA   *tree_data)   /* (I) Tree data structure        */
{

    double  *KoCapRebateSl=NULL;   /* Rebate slice for KoCap     */
    double  *KoPwayRebateSl=NULL;  /* Rebate slice for KoPartway */

    double  *KoCapL;               /* Local slice pointers */
    double  *KoPwayL;
    double  *CapletL;              
    double  *FltZeroToPmtL;
    double  *KoCapRebateSlL;
    double  *KoPwayRebateSlL;
            
    double  *KoCapRebatePtr=NULL;  /* Ptr to rebate var for KoCap        */
    double  *KoPwayRebatePtr=NULL; /* Ptr to rebate var for KoPartway    */

    double  FltAccrued=0.;         /* Flt coupon accrued at current date */
    double  FltAccRatio=0.;        /* CurrAccrued/FullPeriodAccrued      */
    double  KoCapRebate;           /* Rebate number for KoCap            */
    double  KoPwayRebate;          /* Rebate number for KoPartway        */

    int     KoCapRebateDim=0;      /* Dimension of rebate for KoCap  */
    int     KoPwayRebateDim=0;     /* Dimension of rebate for KoPway */

    int     Top1, Bottom1;         /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;       /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;     /* Tree limits (3rd dim)  */

    int     i, j, k;               /* Node indices           */
    int     offset;                /* Node offset            */
    int     status = FAILURE;      /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* Allways discount */

    if (Dev (KoCap,    /* Disc expd value of a cap starting */
             t,        /* on the  last ko date              */
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;        
    }
    
    if (ArrearsReset == 'N')
    {
        if (Dev (KoPway,
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;            
        }
    }
    
    /*
     * Create rebates if Ko date 
     */
    if (KoFlag)
    {
        KoCapRebateSl  = Alloc_Slice (tree_data);
        KoPwayRebateSl = Alloc_Slice (tree_data);

        if ((KoCapRebateSl == NULL) || (KoPwayRebateSl == NULL))
        {
            DR_Error("KoCap_t: could not allocate memory for rebate slices!");
            goto FREE_MEM_AND_RETURN;
        }

        /* Calculate flt accrued interests if we are within accrual period. */

        if (FltCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FltCpnAccSt,
                                    CurrentDate,                                 
                                    FltDCConv,
                                    &FltAccrued) == FAILURE)
            {
                DR_Error("KoCap_t: unable to calculate dcf for flt accrued");
                goto FREE_MEM_AND_RETURN;
            }        
            FltAccRatio = FltAccrued / FltCpnDcf;
            
        }  /* if */
 
        /* Construct rebate variable for KoCap and KoPartway. There is       */
        /* always a rebate for fwd starting cap, but not for KoPway variable */

        KoCapRebate    = Rebate;
        KoCapRebateDim = 0;
        KoCapRebatePtr = &KoCapRebate;
        
        if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            CapletL         = Caplet         + offset;
            FltZeroToPmtL   = FltZeroToPmt   + offset;
            KoCapRebateSlL  = KoCapRebateSl  + offset;
            KoPwayRebateSlL = KoPwayRebateSl + offset;

            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)        /* Standard rebate for KoCap: done! */
                {
                    ;
                }                       /* Swap,Bond stub : add PVed accrued */
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoCapRebateSlL[i] = Rebate + FltAccRatio * CapletL[i];
                        
                    }  /* for i */
                    KoCapRebateDim = 1;
                    KoCapRebatePtr = KoCapRebateSl;
                }
                else                       /* NoStub: st KoCap rebate: done! */
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;                        /* Standard rebate for KoCap: done! */
                
                KoPwayRebate    = 0.;
                KoPwayRebateDim = 0;
                KoPwayRebatePtr = &KoPwayRebate;
                
                if (FltResetFlag)            /* KoPart rebate is zero: done! */
                {
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub: KoPart rebate is dcf */
                {                       
                    KoPwayRebate    = FltAccRatio;
                    KoPwayRebateDim = 0;
                    KoPwayRebatePtr = &KoPwayRebate;
                }
                else if (OptStub == 'S')    /* Swap stub: rebate is dcf*PVed */
                {                        
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoPwayRebateSlL[i] = FltAccRatio * FltZeroToPmtL[i];
                        
                    }  /* for i */
                    KoPwayRebateDim = 1;
                    KoPwayRebatePtr = KoPwayRebateSl;
                }
                else                  /* NoStub: KoPart rebate is zero:done! */
                {
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {
                    ;
                }
                else if((OptStub == 'B') || (OptStub == 'S'))
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        CapletL        = Caplet        + offset;
                        KoCapRebateSlL = KoCapRebateSl + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoCapRebateSlL[j] = Rebate+FltAccRatio*CapletL[j];
                        }
                    }
                    KoCapRebateDim = 1;
                    KoCapRebatePtr = KoCapRebateSl;
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;
                
                KoPwayRebate    = 0.;
                KoPwayRebateDim = 0;
                KoPwayRebatePtr = &KoPwayRebate;
                
                if (FltResetFlag)
                {
                    ;
                }                                                           
                else if (OptStub == 'B')
                {                       
                    KoPwayRebate    = FltAccRatio;
                    KoPwayRebateDim = 0;
                    KoPwayRebatePtr = &KoPwayRebate;
                }
                else if (OptStub == 'S')
                {                        
                    for (i = Bottom1; i <= Top1; i ++)          
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FltZeroToPmtL   = FltZeroToPmt   + offset;
                        KoPwayRebateSlL = KoPwayRebateSl + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoPwayRebateSlL[j] = FltAccRatio*FltZeroToPmtL[j];
                        }
                    }
                    KoPwayRebateDim = 1;
                    KoPwayRebatePtr = KoPwayRebateSl;
                }
                else
                {
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {
                    ;
                }
                else if((OptStub == 'B') || (OptStub == 'S'))
                {                          
                    for (i = Bottom1; i <= Top1; i ++)          
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CapletL        = Caplet        + offset;
                            KoCapRebateSlL = KoCapRebateSl + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                               KoCapRebateSlL[k]=Rebate+FltAccRatio*CapletL[k];
                            }
                        }
                    KoCapRebateDim = 1;
                    KoCapRebatePtr = KoCapRebateSl;
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;
                
                KoPwayRebate    = 0.;
                KoPwayRebateDim = 0;
                KoPwayRebatePtr = &KoPwayRebate;
                
                if (FltResetFlag)
                {
                    ;
                }                                                           
                else if (OptStub == 'B')
                {                       
                    KoPwayRebate    = FltAccRatio;
                    KoPwayRebateDim = 0;
                    KoPwayRebatePtr = &KoPwayRebate;
                }
                else if (OptStub == 'S')
                {                        
                    for (i = Bottom1; i <= Top1; i ++)          
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FltZeroToPmtL   = FltZeroToPmt   + offset;
                            KoPwayRebateSlL = KoPwayRebateSl + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                               KoPwayRebateSlL[k]=FltAccRatio*FltZeroToPmtL[k];
                            }
                        }
                    KoPwayRebateDim = 1;
                    KoPwayRebatePtr = KoPwayRebateSl;
                }
                else
                {
                    ;
                }  /* if then else */ 
            } /* if Arrears */
            
            
        }  /* if then else */

        /* Knockout KoCap and KoPway (if reset in advance ) */
    
        if (KoOptionVarRebate_t(KoCap,
                                KoIndex,
                                KoCapRebatePtr,
                                1,          /* This is Ko date */
                                LoBarrier,
                                HiBarrier,
                                IoO,                 
                                Smoothing,
                                KoCapRebateDim,
                                t,
                                tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
        
        if ((FltCpnAccSt <= CurrentDate) && (ArrearsReset == 'N'))
        {
            if (KoOptionVarRebate_t(KoPway,
                                    KoIndex,
                                    KoPwayRebatePtr,
                                    1,
                                    LoBarrier,
                                    HiBarrier,
                                    IoO,                 
                                    Smoothing,
                                    KoPwayRebateDim,
                                    t,
                                    tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
    } /* if KoFlag */

    /* Add caplet (arrears) or caplet*KoPway (advance) on reset day */
    
    if (FltResetFlag && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoCapL        = KoCap        + offset;
        KoPwayL       = KoPway       + offset;
        CapletL       = Caplet       + offset;
        FltZeroToPmtL = FltZeroToPmt + offset;

        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoCapL[i] += CapletL[i];
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoCapL[i] += KoPwayL[i] * CapletL[i] / FltZeroToPmtL[i];
            }
        }                                               
    }
    else if (FltResetFlag && (tree_data->NbFactor == 2))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoCapL  = KoCap  + offset;
                CapletL = Caplet + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoCapL[j] += CapletL[j];
                }
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoCapL        = KoCap        + offset;
                KoPwayL       = KoPway       + offset;
                CapletL       = Caplet       + offset;
                FltZeroToPmtL = FltZeroToPmt + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoCapL[j] += KoPwayL[j] * CapletL[j] / FltZeroToPmtL[j];
                }
            }
        } /* if Arrears */                                               
    }  
    else if (FltResetFlag && (tree_data->NbFactor == 3))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)          
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoCapL  = KoCap  + offset;
                    CapletL = Caplet + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        KoCapL[k] += CapletL[k];
                    }
                }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)          
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoCapL        = KoCap        + offset;
                    KoPwayL       = KoPway       + offset;
                    CapletL       = Caplet       + offset;
                    FltZeroToPmtL = FltZeroToPmt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        KoCapL[k] += KoPwayL[k]*CapletL[k] / FltZeroToPmtL[k];
                    }
                }
        } /* if Arrears */                                               
    }  /* if then else */

    /* Initialize the KoPartway if on reset date */

    if (FltCpnFlag && (ArrearsReset == 'N'))
    {
        if (Set_Slice(KoPway,
                      1.,
                      t,
                      tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    if (KoFlag)
    {
        Free_Slice (KoCapRebateSl, tree_data);
        Free_Slice (KoPwayRebateSl, tree_data);
    }

    return (status);

}  /* KoCap_t */


    
/*****  KoSwap_t  *******************************************************/
/*
 *  Knock out swap function. 
 */
int  KoSwap_t (double      *KoSwap,       /* (O) Knock out swap              */
               double      *KoPway,       /* (O) Knock out partway           */
               double      *FltSwaplet,   /* (I) Flt part of current swaplet */
               double      *FltZeroToPmt, /* (I) Zero w/ mat on next flt cpn */
               double      *FixZeroToPmt, /* (I) Zero w/ mat on next fix cpn */
               double      *KoIndex,      /* (I) Knock out index             */
               long        KoFlag,        /* (I) Knock out flag              */
               double      LoBarrier,     /* (I) Lower barrier               */
               double      HiBarrier,     /* (I) Higher barrier              */
               double      Rebate,        /* (I) Rebate                      */
               char        IoO,           /* (I) Knock-out 'I'n or 'O'ut     */
               char        Smoothing,     /* (I) Smoothing ('Y' or 'N')      */
               char        OptStub,       /* (I) Option stub rule            */
               long        FltResetFlag,  /* (I) Flt reset flag              */
               long        FltCpnFlag,    /* (I) Flt payment flag            */
               double      FltCpnDcf,     /* (I) Flt coupon dcf              */
               long        FltCpnAccSt,   /* (I) Flt coupon accrued start    */
               char        FltDCConv,     /* (I) Flt coupon day count conv   */
               long        FixCpnFlag,    /* (I) Fix coupon flag             */
               double      FixCpnAmount,  /* (I) Fix coupon amount           */
               double      FixCpnDcf,     /* (I) Fix coupon dcf              */
               long        FixCpnAccSt,   /* (I) Fix coupon accrued start    */
               char        FixDCConv,     /* (I) Fix coupon day count conv   */
               char        ArrearsReset,  /* (I) 'Y' if set-in-arreas        */
               long        CurrentDate,   /* (I) Current date                */
               int         t,             /* (I) Current time point          */
               int         T,             /* (I) Last time point             */
               int         DCurve,        /* (I) Discount curve              */
               DEV_DATA    *dev_data,     /* (I) Dev data structure          */
               TREE_DATA   *tree_data)    /* (I) Tree data structure         */
{

    double  *KoSwapRebate=NULL;     /* Rebate slice for KoCap     */
    double  *KoPwayRebate=NULL;     /* Rebate slice for KoPartway */

    double  *KoSwapL;               /* Local slice pointers */
    double  *KoPwayL;
    double  *FltSwapletL;
    double  *FixZeroToPmtL;
    double  *FltZeroToPmtL;
    double  *KoSwapRebateL;
    double  *KoPwayRebateL;
            
    double  FltAccrued=0.;         /* Flt coupon accrued at current date */
    double  FixAccrued=0.;         /* Fix coupon accrued at current date */
    double  FltAccRatio=0.;        /* CurrAccrued/FullPeriodAccrued      */
    double  FixAccRatio=0.;        /* CurrAccr/FullPeriodAccr for fix    */

    int     Top1, Bottom1;         /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;       /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;     /* Tree limits (3rd dim)  */

    int     i, j, k;               /* Node indices           */
    int     offset;                /* Node offset            */
    int     status = FAILURE;      /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* Allways discount */

    if (Dev (KoSwap,   /* Disc expd value of a swap starting */
             t,        /* on the  last ko date               */
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }
    
    if (ArrearsReset == 'N')
    {
        if (Dev (KoPway,
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }
    
    /*
     * Create rebates if Ko date 
     */
    if (KoFlag)
    {
        KoSwapRebate = Alloc_Slice (tree_data);
        KoPwayRebate = Alloc_Slice (tree_data);

        if ((KoSwapRebate == NULL) || (KoPwayRebate == NULL))
        {
            DR_Error("KoSwap_t: could not allocate memory for rebate slices!");
            goto FREE_MEM_AND_RETURN;
        }

        /* Calculate flt and fix accrued if we are within accrual period. */

        if (FltCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FltCpnAccSt,
                                    CurrentDate,                                 
                                    FltDCConv,
                                    &FltAccrued) == FAILURE)
            {
                DR_Error("KoSwap_t: unable to calculate dcf for flt accrued!");
                goto FREE_MEM_AND_RETURN;
            }        
            FltAccRatio = FltAccrued / FltCpnDcf;            
        }  
        else
        {
            FltAccRatio = 0.;
        } 

        if (FixCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FixCpnAccSt,
                                    CurrentDate,                                 
                                    FixDCConv,
                                    &FixAccrued) == FAILURE)
            {
                DR_Error("KoSwap_t: unable to calculate dcf for fix accrued!");
                goto FREE_MEM_AND_RETURN;
            }        
            FixAccRatio = FixAccrued / FixCpnDcf;
        } 
        else
        {
            FixAccRatio = 0.;
        } 


        /* Construct rebate variable for KoSwap and KoPartway. There is    */
        /* always rebate for fwd starting cap, but not for KoPway variable */

        /* Standard rebate */
        if (KoFlag && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            KoSwapRebateL = KoSwapRebate + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapRebateL[i] = Rebate;
            }
        }
        else if (KoFlag && (tree_data->NbFactor == 2))
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapRebateL = KoSwapRebate + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapRebateL[j] = Rebate;
                }
            }
        }
        else if (KoFlag && (tree_data->NbFactor == 3))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapRebateL[k] = Rebate;
                    }
                }
        }  /* if then else */


        /* Fix leg rebate */
        if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            KoSwapRebateL = KoSwapRebate + offset;
            FixZeroToPmtL = FixZeroToPmt + offset;

            if (FixCpnFlag)             /* Standard rebate for KoSwap: done! */
            {                
                ;
            }                                                           
            else if (OptStub == 'B')  /* Bond stub: KoSwap rebate is accrued */
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    KoSwapRebateL[i] += FixLegRebate;
                }
            }
            else if (OptStub == 'S')    /* Swap stub: rebate is PVed accrued */
            {             
                for (i = Bottom1; i <= Top1; i ++)
                {
                   KoSwapRebateL[i]+=FixAccRatio*FixCpnAmount*FixZeroToPmtL[i];
                }
            }
            else                 /* NoStub: KoSwap rebate is standard: done! */
            {                  
                ;
            }  /* if then else */ 
        } 
        else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
        {
            if (FixCpnFlag)
            {                
                ;
            }                                                           
            else if (OptStub == 'B')
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        KoSwapRebateL[j] += FixLegRebate;
                    }
                }
            }
            else if (OptStub == 'S')
            {             
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;
                    FixZeroToPmtL = FixZeroToPmt + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        KoSwapRebateL[j] += FixAccRatio * FixCpnAmount 
                                                            * FixZeroToPmtL[j];
                    }
                }
            }
            else
            {                  
                ;
            }  /* if then else */ 
        }
        else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
        {
            if (FixCpnFlag)
            {                
                ;
            }                                                           
            else if (OptStub == 'B')
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoSwapRebateL = KoSwapRebate + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                        {
                            KoSwapRebateL[k] += FixLegRebate;
                        }
                    }
            }
            else if (OptStub == 'S')
            {             
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoSwapRebateL = KoSwapRebate + offset;
                        FixZeroToPmtL = FixZeroToPmt + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                        {
                            KoSwapRebateL[k] += FixAccRatio * FixCpnAmount 
                                                            * FixZeroToPmtL[k];
                        }
                    }
            }
            else
            {                  
                ;
            }  /* if then else */ 
        }


        /* Floating leg rebate */

        if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            FltSwapletL   = FltSwaplet   + offset;
            FltZeroToPmtL = FltZeroToPmt + offset;
            KoSwapRebateL = KoSwapRebate + offset;
            KoPwayRebateL = KoPwayRebate + offset;

            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)     /* Std rebate for KoSwap: done in fix! */
                {        
                    ;
                }                          
                else if(OptStub == 'B')         /* Bond stub : NOT SUPPORTED */
                {
                    DR_Error("KoSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(OptStub == 'S')      /* Swap stub : add PVed accrued */
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoSwapRebateL[i] += FltAccRatio * FltSwapletL[i];
                    }
                }
                else                      /* NoStub: st KoSwap rebate: done! */
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;                       /* Rebate for KoSwap: done !         */
                                        /* Standard rebate for KoPway: zero! */

                if (FltResetFlag)       /* KoPart rebate is zero: done !     */
                {                      
                    ;
                }                                                           
                else if (OptStub == 'B')  /* Bond stub: KoPart rebate is dcf */
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoPwayRebateL[i] = FltAccRatio;
                    }
                }
                else if (OptStub == 'S')    /* Swap stub: rebate is dcf*PVed */
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoPwayRebateL[i] = FltAccRatio * FltZeroToPmtL[i];
                    }
                }
                else                 /* NoStub: KoPart rebate is zero: done! */
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {        
                    ;
                }                          
                else if(OptStub == 'B')
                {
                    DR_Error("KoSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(OptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FltSwapletL   = FltSwaplet   + offset;
                        KoSwapRebateL = KoSwapRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoSwapRebateL[j] += FltAccRatio * FltSwapletL[j];
                        }
                    }
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;


                if (FltResetFlag)
                {                      
                    ;
                }                                                           
                else if (OptStub == 'B')
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        KoPwayRebateL = KoPwayRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoPwayRebateL[j] = FltAccRatio;
                        }
                    }
                }
                else if (OptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FltZeroToPmtL = FltZeroToPmt + offset;
                        KoPwayRebateL = KoPwayRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoPwayRebateL[j] = FltAccRatio * FltZeroToPmtL[j];
                        }
                    }
                }
                else
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {        
                    ;
                }                          
                else if(OptStub == 'B')
                {
                    DR_Error("KoSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(OptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FltSwapletL   = FltSwaplet   + offset;
                            KoSwapRebateL = KoSwapRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoSwapRebateL[k] += FltAccRatio*FltSwapletL[k];
                            }
                        }
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;


                if (FltResetFlag)
                {                      
                    ;
                }                                                           
                else if (OptStub == 'B')
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            KoPwayRebateL = KoPwayRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoPwayRebateL[k] = FltAccRatio;
                            }
                        }
                }
                else if (OptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FltZeroToPmtL = FltZeroToPmt + offset;
                            KoPwayRebateL = KoPwayRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoPwayRebateL[k]=FltAccRatio*FltZeroToPmtL[k];
                            }
                        }
                }
                else
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }  /* if then else */


        /* Knockout KoSwap and KoPway (if reset in advance ) */
    
        if (KoOptionVarRebate_t(KoSwap,
                                KoIndex,
                                KoSwapRebate,
                                1,          /* This is Ko date */
                                LoBarrier,
                                HiBarrier,
                                IoO,                 
                                Smoothing,
                                1,         /* There is always a slice rebate */
                                t,
                                tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
        
        if ((FltCpnAccSt <= CurrentDate) && (ArrearsReset == 'N'))
        {
            if (KoOptionVarRebate_t(KoPway,
                                    KoIndex,
                                    KoPwayRebate,
                                    1,
                                    LoBarrier,
                                    HiBarrier,
                                    IoO,                 
                                    Smoothing,
                                    1,     /* There is always a slice rebate */
                                    t,
                                    tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
    } /* if KoFlag */


    /* Add fix swaplet */

    if (FixCpnFlag && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoSwapL = KoSwap + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            KoSwapL[i] += FixCpnAmount;
        }
    }
    else if (FixCpnFlag && (tree_data->NbFactor == 2))
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Node_Offset(2, i, 0, t, tree_data);

            KoSwapL = KoSwap + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                KoSwapL[j] += FixCpnAmount;
            }
        }
    }  
    else if (FixCpnFlag && (tree_data->NbFactor == 3))
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Node_Offset(3, i, j, t, tree_data);

                KoSwapL = KoSwap + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                {
                    KoSwapL[k] += FixCpnAmount;
                }
            }
    }  /* if then else */


    /* Add flt swaplet (arrears) or flt swaplet*KoPway (adv) on reset day */

    if (FltResetFlag && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoSwapL       = KoSwap       + offset;
        KoPwayL       = KoPway       + offset;
        FltSwapletL   = FltSwaplet   + offset;
        FltZeroToPmtL = FltZeroToPmt + offset;

        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapL[i] += FltSwapletL[i];
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapL[i] += KoPwayL[i] * FltSwapletL[i] / FltZeroToPmtL[i];
            }
        }                                               
    }
    else if (FltResetFlag && (tree_data->NbFactor == 2))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapL     = KoSwap     + offset;
                FltSwapletL = FltSwaplet + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapL[j] += FltSwapletL[j];
                }
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapL       = KoSwap       + offset;
                KoPwayL       = KoPway       + offset;
                FltSwapletL   = FltSwaplet   + offset;
                FltZeroToPmtL = FltZeroToPmt + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapL[j] += KoPwayL[j] * FltSwapletL[j]/FltZeroToPmtL[j];
                }
            }
        }                                               
    }  
    else if (FltResetFlag && (tree_data->NbFactor == 3))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapL       = KoSwap     + offset;
                    FltSwapletL   = FltSwaplet + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapL[k] += FltSwapletL[k];
                    }
                }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapL       = KoSwap       + offset;
                    KoPwayL       = KoPway       + offset;
                    FltSwapletL   = FltSwaplet   + offset;
                    FltZeroToPmtL = FltZeroToPmt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapL[k]+=KoPwayL[k]*FltSwapletL[k]/FltZeroToPmtL[k];
                    }
                }
        }                                               
    }  /* if then else */

    /* Initialize the KoPartway if on reset date */

    if (FltCpnFlag && (ArrearsReset == 'N'))
    {
        if (Set_Slice(KoPway,
                      1.,
                      t,
                      tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    if (KoFlag)
    {
        Free_Slice (KoSwapRebate, tree_data);
        Free_Slice (KoPwayRebate, tree_data);
    }

    return (status);

}  /* KoSwap_t */



/*****  KoOptionVarRebate_t   ***********************************************/
/*                                                                           
*       Knock-out an option price without discounting it. It uses Arnon's
*       smoothing algorithm. This version allows for node-dependent rebate.
*/
int     KoOptionVarRebate_t (
             double      *KoOpt,      /* (I/O) Knock-out option prices       */
             double      *Index,      /* (I) Index prices                    */
             double      *Rebate,     /* (I) Rebate                          */
             long        KoFlag,      /* (I) Knock-out flag                  */
             double      LowBarrier,  /* (I) Lower barrier                   */
             double      HighBarrier, /* (I) Higher barrier                  */
             char        IoO,         /* (I) Knock-out 'I'nside or 'O'utside */
             char        Smoothing,   /* (I) Smoothing ('Y' or 'N')          */
             int         RebateDim,   /* (I) Dimension of rebate slice       */
             int         t,           /* (I) Current time point              */
             TREE_DATA   *tree_data)  /* (I) Structure of tree data          */
{

    double  *RebateN=NULL;        /* Numerical value of the rabate at node   */

    double  *KoOptL;              /* Local slice pointers */
    double  *IndexL;
    double  *RebateL=NULL;
        
    double  x;                    /* Intermediate smoothed value             */
    double  IndexStep;            /* Maximum difference between index values */

    int     Top1, Bottom1;        /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;      /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;    /* Tree limits (3rd dim)  */

    int     i, j, k;              /* Node indices           */
    int     offset;               /* Node offset            */
    int     status = FAILURE;     /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (!KoFlag)    /* Nothing to do */
    {
        return (SUCCESS);
    }


    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoOptL = KoOpt + offset;
        IndexL = Index + offset;

        if (RebateDim == 0) 
            RebateN = Rebate;
        else 
            RebateL = Rebate + offset;


        if (Smoothing == 'N')
        {
            if (IoO == 'I')     /* Knock-out inside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if (   (IndexL[i] > LowBarrier * (1. - BARRIER_TOL))
                        && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i] = (RebateDim == 0) ? *RebateN :  RebateL[i];
                    }
                }  /* for i */  
            }           
            else                /* Knock-out outside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if (   (IndexL[i] < LowBarrier * (1. + BARRIER_TOL))
                        || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KoOptL[i] = (RebateDim == 0) ? *RebateN :  RebateL[i];
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                IndexStep = GetIndexStep (Index,
                                          1,                                                  
                                          i, 0, 0,
                                          t,
                                          tree_data);
                    
                if (IoO == 'I')
                {                                               
                    if (Smooth_Step (&x,
                                     (RebateDim == 0) ? *RebateN :  RebateL[i],
                                     KoOptL[i],
                                     IndexL[i],
                                     LowBarrier,
                                     IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                    
                    if (Smooth_Step (&(KoOptL[i]),
                                     KoOptL[i],
                                     x,
                                     IndexL[i],
                                     HighBarrier,
                                     IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                }   
                else
                {   
                    if (Smooth_Step (&x,
                                     KoOptL[i],
                                     (RebateDim == 0) ? *RebateN :  RebateL[i],
                                     IndexL[i],
                                     LowBarrier,
                                     IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                    
                    if (Smooth_Step (&(KoOptL[i]),
                                     (RebateDim == 0) ? *RebateN :  RebateL[i],
                                     x,
                                     IndexL[i],
                                     HighBarrier,
                                     IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                }  /* if then else */   
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    if (RebateDim == 0) 
                        RebateN = Rebate;
                    else 
                        RebateL = Rebate + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if (   (IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                            && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KoOptL[j]=(RebateDim == 0) ? *RebateN : RebateL[j];
                        }
                    }  /* for j */
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    if (RebateDim == 0) 
                        RebateN = Rebate;
                    else 
                        RebateL = Rebate + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if (   (IndexL[j] < LowBarrier * (1. + BARRIER_TOL))
                            || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KoOptL[j]=(RebateDim == 0) ? *RebateN : RebateL[j];
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoOptL = KoOpt + offset;
                IndexL = Index + offset;

                if (RebateDim == 0) 
                    RebateN = Rebate;
                else 
                    RebateL = Rebate + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexStep = GetIndexStep (Index,
                                              2,                                                  
                                              i, j, 0,
                                              t,
                                              tree_data);
                    
                    if (IoO == 'I')
                    {                                               
                        if (Smooth_Step(&x,
                                        (RebateDim==0) ? *RebateN : RebateL[j],
                                        KoOptL[j],
                                        IndexL[j],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        if (Smooth_Step (&(KoOptL[j]),
                                         KoOptL[j],
                                         x,
                                         IndexL[j],
                                         HighBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }   
                    else
                    {   
                        if (Smooth_Step(&x,
                                        KoOptL[j],
                                        (RebateDim==0) ? *RebateN : RebateL[j],
                                        IndexL[j],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Smooth_Step(&(KoOptL[j]),
                                        (RebateDim==0) ? *RebateN : RebateL[j],
                                        x,
                                        IndexL[j],
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }  /* if then else */   
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        if (RebateDim == 0) 
                            RebateN = Rebate;
                        else 
                            RebateL = Rebate + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL))
                             && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                              KoOptL[k]=(RebateDim==0) ? *RebateN : RebateL[k];
                            }
                        }  /* for k */  
                    }  /* for j */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        if (RebateDim == 0) 
                            RebateN = Rebate;
                        else 
                            RebateL = Rebate + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if (  (IndexL[k] < LowBarrier * (1. + BARRIER_TOL))
                               || (IndexL[k] > HighBarrier * (1.-BARRIER_TOL)))
                            {        
                                KoOptL[k] = (RebateDim==0) ? *RebateN : RebateL[k];
                            }
                        }  /* for k */  
                    }  /* for j */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    if (RebateDim == 0) 
                        RebateN = Rebate;
                    else 
                        RebateL = Rebate + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        IndexStep = GetIndexStep (Index,
                                                  3,                                                  
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                        
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (&x,
                                             (RebateDim == 0) ? *RebateN : RebateL[k],
                                             KoOptL[k],
                                             IndexL[k],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                    
                            if (Smooth_Step (&(KoOptL[k]),
                                             KoOptL[k],
                                             x,
                                             IndexL[k],
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }   
                        else
                        {   
                            if (Smooth_Step (&x,
                                             KoOptL[k],
                                             (RebateDim == 0) ? *RebateN :  RebateL[k],
                                             IndexL[k],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step (&(KoOptL[k]),
                                             (RebateDim == 0) ? *RebateN :  RebateL[k],
                                             x,
                                             IndexL[k],
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }  /* if then else */   
                    }  /* for k */  
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* KoOptionVarRebate_t */


/*****  KoFwdSwap_t  *******************************************************/
/*
 *  Knock out fwd swap function. This is the same as the KoSwap_t except
 *  it supports separate fixed and floating stub conventions
 */
int  KoFwdSwap_t(double    *KoSwap,       /* (O) Knock out swap              */
                 double    *KoPway,       /* (O) Knock out partway           */
                 double    *FltSwaplet,   /* (I) Flt part of current swaplet */
                 double    *FltZeroToPmt, /* (I) Zero w/ mat on next flt cpn */
                 double    *FixZeroToPmt, /* (I) Zero w/ mat on next fix cpn */
                 double    *IFltZeroToPmt,/* (I) as above but on Idx curve   */
                 double    *KoIndex,      /* (I) Knock out index             */
                 long      KoFlag,        /* (I) Knock out flag              */
                 double    LoBarrier,     /* (I) Lower barrier               */
                 double    HiBarrier,     /* (I) Higher barrier              */
                 double    Rebate,        /* (I) Rebate                      */
                 char      IoO,           /* (I) Knock-out 'I'n or 'O'ut     */
                 char      Smoothing,     /* (I) Smoothing ('Y' or 'N')      */
                 char      FixOptStub,    /* (I) Fix stub rule               */
                 char      FltOptStub,    /* (I) Flt stub rule               */
                 long      FltResetFlag,  /* (I) Flt reset flag              */
                 long      FltCpnFlag,    /* (I) Flt payment flag            */
                 double    FltCpnOuts,    /* (I) Flt coupon o/s notional     */
                 double    FltCpnDcf,     /* (I) Flt coupon dcf              */
                 long      FltCpnAccSt,   /* (I) Flt coupon accrued start    */
                 char      FltDCConv,     /* (I) Flt coupon day count conv   */
                 long      FixCpnFlag,    /* (I) Fix coupon flag             */
                 double    FixCpnAmount,  /* (I) Fix coupon amount           */
                 double    FixCpnDcf,     /* (I) Fix coupon dcf              */
                 long      FixCpnAccSt,   /* (I) Fix coupon accrued start    */
                 char      FixDCConv,     /* (I) Fix coupon day count conv   */
                 char      ArrearsReset,  /* (I) 'Y' if set-in-arreas        */
                 long      CurrentDate,   /* (I) Current date                */
                 int       t,             /* (I) Current time point          */
                 int       T,             /* (I) Last time point             */
                 int       DCurve,        /* (I) Discount curve              */
                 DEV_DATA  *dev_data,     /* (I) Dev data structure          */
                 TREE_DATA *tree_data)    /* (I) Tree data structure         */
{

    double  *KoSwapRebate=NULL;     /* Rebate slice for KoSwap    */
    double  *KoPwayRebate=NULL;     /* Rebate slice for KoPartway */

    double  *KoSwapL;               /* Local slice pointers */
    double  *KoPwayL;
    double  *FltSwapletL;
    double  *FixZeroToPmtL;
    double  *FltZeroToPmtL;
    double  *KoSwapRebateL;
    double  *KoPwayRebateL;
    double  *IFltZeroToPmtL;
            
    double  FltAccrued=0.;         /* Flt coupon accrued at current date */
    double  FixAccrued=0.;         /* Fix coupon accrued at current date */
    double  FltAccRatio=0.;        /* CurrAccrued/FullPeriodAccrued      */
    double  FixAccRatio=0.;        /* CurrAccr/FullPeriodAccr for fix    */

    int     Top1, Bottom1;         /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;       /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;     /* Tree limits (3rd dim)  */

    int     i, j, k;               /* Node indices           */
    int     offset;                /* Node offset            */
    int     status = FAILURE;      /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* Allways discount */

    if (Dev (KoSwap,   /* Disc expd value of a swap starting */
             t,        /* on the  last ko date               */
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }
    
    if (ArrearsReset == 'N')
    {
        if (Dev (KoPway,
                 t,
                 T,
                 DCurve,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }
    
    /*
     * Create rebates if Ko date 
     */
    if (KoFlag)
    {
        KoSwapRebate = Alloc_Slice (tree_data);
        KoPwayRebate = Alloc_Slice (tree_data);

        if ((KoSwapRebate == NULL) || (KoPwayRebate == NULL))
        {
            DR_Error("KoFwdSwap_t: could not allocate memory "
                     "for rebate slices!");
            goto FREE_MEM_AND_RETURN;
        }

        /* Calculate flt and fix accrued if we are within accrual period. */

        if (FltCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FltCpnAccSt,
                                    CurrentDate,
                                    FltDCConv,
                                    &FltAccrued) == FAILURE)
            {
                DR_Error("KoFwdSwap_t: unable to calculate dcf "
                         "for flt accrued!");
                goto FREE_MEM_AND_RETURN;
            }        
            FltAccRatio = FltAccrued / FltCpnDcf;            
        }  
        else
        {
            FltAccRatio = 0.;
        } 

        if (FixCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FixCpnAccSt,
                                    CurrentDate,
                                    FixDCConv,
                                    &FixAccrued) == FAILURE)
            {
                DR_Error("KoFwdSwap_t: unable to calculate dcf "
                         "for fix accrued!");
                goto FREE_MEM_AND_RETURN;
            }        
            FixAccRatio = FixAccrued / FixCpnDcf;
        } 
        else
        {
            FixAccRatio = 0.;
        } 


        /* Construct rebate variable for KoSwap and KoPartway. There is    */
        /* always rebate for fwd starting swap, but not for KoPway variable */

        /* Standard rebate */
        if (KoFlag && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            KoSwapRebateL = KoSwapRebate + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapRebateL[i] = Rebate;
            }
        }
        else if (KoFlag && (tree_data->NbFactor == 2))
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapRebateL = KoSwapRebate + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapRebateL[j] = Rebate;
                }
            }
        }
        else if (KoFlag && (tree_data->NbFactor == 3))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapRebateL[k] = Rebate;
                    }
                }
        }  /* if then else */


        /* Fix leg rebate */
        if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            KoSwapRebateL = KoSwapRebate + offset;
            FixZeroToPmtL = FixZeroToPmt + offset;

            if (FixCpnFlag)             /* Standard rebate for KoSwap: done! */
            {                
                ;
            }                                                           
            else if (FixOptStub == 'B') /* Bond stub: */
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    KoSwapRebateL[i] += FixLegRebate;
                }
            }
            else if (FixOptStub == 'S')    /* Swap stub: */
            {             
                for (i = Bottom1; i <= Top1; i ++)
                {
                   KoSwapRebateL[i]+=FixAccRatio*FixCpnAmount*FixZeroToPmtL[i];
                }
            }
            else                 /* NoStub: KoSwap rebate is standard: done! */
            {                  
                ;
            }  /* if then else */ 
        } 
        else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
        {
            if (FixCpnFlag)
            {                
                ;
            }                                                           
            else if (FixOptStub == 'B')
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        KoSwapRebateL[j] += FixLegRebate;
                    }
                }
            }
            else if (FixOptStub == 'S')
            {             
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    KoSwapRebateL = KoSwapRebate + offset;
                    FixZeroToPmtL = FixZeroToPmt + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        KoSwapRebateL[j] += FixAccRatio * FixCpnAmount 
                                                            * FixZeroToPmtL[j];
                    }
                }
            }
            else
            {                  
                ;
            }  /* if then else */ 
        }
        else if ((FixCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
        {
            if (FixCpnFlag)
            {                
                ;
            }                                                           
            else if (FixOptStub == 'B')
            {   
                double FixLegRebate;

                FixLegRebate = FixAccRatio * FixCpnAmount;

                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoSwapRebateL = KoSwapRebate + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                        {
                            KoSwapRebateL[k] += FixLegRebate;
                        }
                    }
            }
            else if (FixOptStub == 'S')
            {             
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

                        KoSwapRebateL = KoSwapRebate + offset;
                        FixZeroToPmtL = FixZeroToPmt + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                        {
                            KoSwapRebateL[k] += FixAccRatio * FixCpnAmount 
                                                            * FixZeroToPmtL[k];
                        }
                    }
            }
            else
            {                  
                ;
            }  /* if then else */ 
        }


        /* Floating leg rebate */

        if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 1))
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            FltSwapletL   = FltSwaplet   + offset;
            FltZeroToPmtL = FltZeroToPmt + offset;
            KoSwapRebateL = KoSwapRebate + offset;
            KoPwayRebateL = KoPwayRebate + offset;
            IFltZeroToPmtL = IFltZeroToPmt + offset;

            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)     /* Std rebate for KoSwap: done in fix! */
                {        
                    ;
                }                          
                else if(FltOptStub == 'B')   /* Bond stub : NOT SUPPORTED */
                {
                    DR_Error("KoFwdSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(FltOptStub == 'S')   /* Swap stub : add PVed accrued */
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoSwapRebateL[i] += FltAccRatio * FltSwapletL[i];
                    }
                }
                else if(FltOptStub == 'P') /* Par stub */
                {
                    /* if strictly within accrual period */
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {
                        for (i = Bottom1; i <= Top1; i ++)
                        {
                            KoSwapRebateL[i] += FltSwapletL[i] + 
                                                FltCpnOuts * 
                                                (1/IFltZeroToPmtL[i]-1) * 
                                                FltZeroToPmtL[i];
                        }
                    }
                }
                else                      /* NoStub: st KoSwap rebate: done! */
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;                       /* Rebate for KoSwap: done !         */
                                        /* Standard rebate for KoPway: zero! */

                if (FltResetFlag)       /* KoPart rebate is zero: done !     */
                {                      
                    ;
                }                                                           
                else if (FltOptStub == 'B')  /* Bond stub: */
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoPwayRebateL[i] = FltAccRatio;
                    }
                }
                else if (FltOptStub == 'S')    /* Swap stub: */
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        KoPwayRebateL[i] = FltAccRatio * FltZeroToPmtL[i];
                    }
                }
                else if (FltOptStub == 'P')    /* Par stub: */
                {
                    /* if strictly within accrual period */
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {                                          
                        for (i = Bottom1; i <= Top1; i ++)
                        {
                            KoPwayRebateL[i] = FltZeroToPmtL[i];
                            KoSwapRebateL[i] += FltCpnOuts *
                                                (1/IFltZeroToPmtL[i]-1) * 
                                                FltZeroToPmtL[i];
                        }
                    }
                }
                else                 /* NoStub: KoPart rebate is zero: done! */
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 2))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {        
                    ;
                }                          
                else if(FltOptStub == 'B')
                {
                    DR_Error("KoFwdSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(FltOptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FltSwapletL   = FltSwaplet   + offset;
                        KoSwapRebateL = KoSwapRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoSwapRebateL[j] += FltAccRatio * FltSwapletL[j];
                        }
                    }
                }
                else if(FltOptStub == 'P')
                {
                    /* if strictly within accrual period */
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {
                        for (i = Bottom1; i <= Top1; i ++)
                        {
                            offset = Node_Offset(2, i, 0, t, tree_data);
                    
                            FltSwapletL   = FltSwaplet   + offset;
                            KoSwapRebateL = KoSwapRebate + offset;
                            FltZeroToPmtL = FltZeroToPmt + offset;
                            IFltZeroToPmtL = IFltZeroToPmt + offset;
                    
                            for (j = Bottom2[i]; j <= Top2[i]; j++)
                            {
                                KoSwapRebateL[j] += FltSwapletL[j] + 
                                                    FltCpnOuts * 
                                                    (1/IFltZeroToPmtL[j]-1) * 
                                                    FltZeroToPmtL[j];
                            }
                        }
                    }
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;


                if (FltResetFlag)
                {                      
                    ;
                }                                                           
                else if (FltOptStub == 'B')
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        KoPwayRebateL = KoPwayRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoPwayRebateL[j] = FltAccRatio;
                        }
                    }
                }
                else if (FltOptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Node_Offset(2, i, 0, t, tree_data);

                        FltZeroToPmtL = FltZeroToPmt + offset;
                        KoPwayRebateL = KoPwayRebate + offset;

                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            KoPwayRebateL[j] = FltAccRatio * FltZeroToPmtL[j];
                        }
                    }
                }
                else if (FltOptStub == 'P')    /* Par stub: */
                {                          
                    /* if strictly within accrual period */
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {
                        for (i = Bottom1; i <= Top1; i ++)
                        {
                            offset = Node_Offset(2, i, 0, t, tree_data);
                    
                            FltZeroToPmtL = FltZeroToPmt + offset;
                            IFltZeroToPmtL = IFltZeroToPmt + offset;
                            KoPwayRebateL = KoPwayRebate + offset;
                            KoSwapRebateL = KoSwapRebate + offset;
                    
                            for (j = Bottom2[i]; j <= Top2[i]; j++)
                            {
                                KoPwayRebateL[j] = FltZeroToPmtL[j];
                                KoSwapRebateL[j] += FltCpnOuts * 
                                                    (1/IFltZeroToPmtL[j]-1) * 
                                                    FltZeroToPmtL[j];
                            }
                        }
                    }
                }
                else
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }
        else if ((FltCpnAccSt <= CurrentDate) && (tree_data->NbFactor == 3))
        {
            if (ArrearsReset == 'Y')
            {
                if (FltResetFlag)
                {        
                    ;
                }                          
                else if(FltOptStub == 'B')
                {
                    DR_Error("KoFwdSwap_t: bond stub is not supported for "
                                "reset-in-arrears swap!");
                    goto FREE_MEM_AND_RETURN;
                }
                else if(FltOptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FltSwapletL   = FltSwaplet   + offset;
                            KoSwapRebateL = KoSwapRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoSwapRebateL[k] += FltAccRatio*FltSwapletL[k];
                            }
                        }
                }
                else if(FltOptStub == 'P')
                {
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {
                        for (i = Bottom1; i <= Top1; i ++)
                            for (j = Bottom2[i]; j <= Top2[i]; j++)
                            {
                                offset = Node_Offset(3, i, j, t, tree_data);
                    
                                FltSwapletL   = FltSwaplet   + offset;
                                KoSwapRebateL = KoSwapRebate + offset;
                                FltZeroToPmtL = FltZeroToPmt + offset;
                                IFltZeroToPmtL = IFltZeroToPmt + offset;
                    
                                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                                {
                                    KoSwapRebateL[k]+=FltSwapletL[k] + 
                                                      FltCpnOuts * 
                                                      (1/IFltZeroToPmtL[k]-1) * 
                                                      FltZeroToPmtL[k];
                                } /* for k */
                            } /* for j */
                    }
                }
                else
                {
                    ;
                } /* if then else */
            }
            else 
            {
                ;


                if (FltResetFlag)
                {                      
                    ;
                }                                                           
                else if (FltOptStub == 'B')
                {                         
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            KoPwayRebateL = KoPwayRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoPwayRebateL[k] = FltAccRatio;
                            }
                        }
                }
                else if (FltOptStub == 'S')
                {                          
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Node_Offset(3, i, j, t, tree_data);

                            FltZeroToPmtL = FltZeroToPmt + offset;
                            KoPwayRebateL = KoPwayRebate + offset;

                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                            {
                                KoPwayRebateL[k]=FltAccRatio*FltZeroToPmtL[k];
                            }
                        }
                }
                else if (FltOptStub == 'P')    /* Par stub: */
                {
                    /* if strictly within accrual period */
                    if ((FltCpnAccSt < CurrentDate) && (!FltCpnFlag))
                    {
                        for (i = Bottom1; i <= Top1; i ++)
                            for (j = Bottom2[i]; j <= Top2[i]; j++)
                            {
                                offset = Node_Offset(3, i, j, t, tree_data);
                    
                                FltZeroToPmtL = FltZeroToPmt + offset;
                                IFltZeroToPmtL = IFltZeroToPmt + offset;
                                KoPwayRebateL = KoPwayRebate + offset;
                                KoSwapRebateL = KoSwapRebate + offset;
                    
                                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                                {
                                    KoPwayRebateL[k] = FltZeroToPmtL[k];
                                    KoSwapRebateL[k]+=FltCpnOuts * 
                                                      (1/IFltZeroToPmtL[k]-1) * 
                                                      FltZeroToPmtL[k];
                                }
                            }
                    }
                }
                else
                {                      
                    ;
                }  /* if then else */ 
            } /* if Arrears */
        }  /* if then else */


        /* Knockout KoSwap and KoPway (if reset in advance ) */
    
        if (KoOptionVarRebate_t(KoSwap,
                                KoIndex,
                                KoSwapRebate,
                                1,          /* This is Ko date */
                                LoBarrier,
                                HiBarrier,
                                IoO,                 
                                Smoothing,
                                1,         /* There is always a slice rebate */
                                t,
                                tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
        
        if ((FltCpnAccSt <= CurrentDate) && (ArrearsReset == 'N'))
        {
            if (KoOptionVarRebate_t(KoPway,
                                    KoIndex,
                                    KoPwayRebate,
                                    1,
                                    LoBarrier,
                                    HiBarrier,
                                    IoO,                 
                                    Smoothing,
                                    1,     /* There is always a slice rebate */
                                    t,
                                    tree_data) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
    } /* if KoFlag */


    /* Add fix swaplet */

    if (FixCpnFlag && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoSwapL = KoSwap + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            KoSwapL[i] += FixCpnAmount;
        }
    }
    else if (FixCpnFlag && (tree_data->NbFactor == 2))
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Node_Offset(2, i, 0, t, tree_data);

            KoSwapL = KoSwap + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                KoSwapL[j] += FixCpnAmount;
            }
        }
    }  
    else if (FixCpnFlag && (tree_data->NbFactor == 3))
    {
        for (i = Bottom1; i <= Top1; i ++)
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Node_Offset(3, i, j, t, tree_data);

                KoSwapL = KoSwap + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                {
                    KoSwapL[k] += FixCpnAmount;
                }
            }
    }  /* if then else */


    /* Add flt swaplet (arrears) or flt swaplet*KoPway (adv) on reset day */

    if (FltResetFlag && (tree_data->NbFactor == 1))
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        KoSwapL       = KoSwap       + offset;
        KoPwayL       = KoPway       + offset;
        FltSwapletL   = FltSwaplet   + offset;
        FltZeroToPmtL = FltZeroToPmt + offset;

        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapL[i] += FltSwapletL[i];
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                KoSwapL[i] += KoPwayL[i] * FltSwapletL[i] / FltZeroToPmtL[i];
            }
        }                                               
    }
    else if (FltResetFlag && (tree_data->NbFactor == 2))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapL     = KoSwap     + offset;
                FltSwapletL = FltSwaplet + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapL[j] += FltSwapletL[j];
                }
            }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                KoSwapL       = KoSwap       + offset;
                KoPwayL       = KoPway       + offset;
                FltSwapletL   = FltSwaplet   + offset;
                FltZeroToPmtL = FltZeroToPmt + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    KoSwapL[j] += KoPwayL[j] * FltSwapletL[j]/FltZeroToPmtL[j];
                }
            }
        }                                               
    }  
    else if (FltResetFlag && (tree_data->NbFactor == 3))
    {
        if (ArrearsReset == 'Y')
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapL       = KoSwap     + offset;
                    FltSwapletL   = FltSwaplet + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapL[k] += FltSwapletL[k];
                    }
                }
        }   
        else
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    KoSwapL       = KoSwap       + offset;
                    KoPwayL       = KoPway       + offset;
                    FltSwapletL   = FltSwaplet   + offset;
                    FltZeroToPmtL = FltZeroToPmt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)   
                    {
                        KoSwapL[k]+=KoPwayL[k]*FltSwapletL[k]/FltZeroToPmtL[k];
                    }
                }
        }                                               
    }  /* if then else */

    /* Initialize the KoPartway if on reset date */

    if (FltCpnFlag && (ArrearsReset == 'N'))
    {
        if (Set_Slice(KoPway,
                      1.,
                      t,
                      tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    if (KoFlag)
    {
        Free_Slice (KoSwapRebate, tree_data);
        Free_Slice (KoPwayRebate, tree_data);
    }

    return (status);

}  /* KoFwdSwap_t */



