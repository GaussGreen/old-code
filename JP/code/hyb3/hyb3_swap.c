/****************************************************************************/
/*      Calculation of forward swap in the lattice.                         */
/****************************************************************************/
/*      SWAP.C                                                              */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


                        
 

/*****  Hyb3_FwdSwap_t  ***********************************************************/
/*
 *      Calculates the  forward swap price on the lattice minus an offset
 *      level.  In other  words, this  function calculates  the intrinsic
 *      value of an option on the forward swap struck at the offset level.
 *
 *      The  swap  is described by a  fixed and a floating  legs and  the 
 *      forward value is computed BY ADDING the two legs.The direction of
 *      the flows (if relevant) must be set while constructing the legs.
 *
 *      In addition, the function takes a flag to indicate whether the flt
 *      leg is ON. If FALSE, then the value of the floater is still priced
 *      in, but no processing of accruals is done.
 *
 *
 */
int    Hyb3_FwdSwap_t(TSLICE      Swap,           /* (O) Forward swap             */
                 TSLICE      Bond,           /* (I) Underlying bond          */
                 TSLICE      Floater,        /* (I) Underlying floater       */
                 int         FltLegOn,       /* (I) True if flt being priced */
                 TSLICE      FltIndex,       /* (I) Last/current index rate  */
                 TSLICE      FixZeroToPmt,   /* (I) Zero w/mat on next fix cp*/
                 TSLICE      FltZeroToPmt,   /* (I) Zero w/mat on next flt cp*/
                 long        FixCpnFlag,     /* (I) Fix coupon flag          */
                 double      FixCpnRate,     /* (I) Fix coupon rate          */
                 double      FixCpnOuts,     /* (I) Fix coupon outstanding   */
                 double      FixCpnDcf,      /* (I) Fix coupon dcf           */
                 long        FixCpnAccSt,    /* (I) Fix coupon accrued start */
                 char        FixDCConv,      /* (I) Fix coupon day count conv*/
                 long        FltResetFlag,   /* (I) Flt coupon flag          */
                 double      FltCpnOuts,     /* (I) Flt coupon outstanding   */
                 double      FltCpnDcf,      /* (I) Flt coupon dcf           */
                 double      FltCpnSpd,      /* (I) Spread to be added to flt*/
                 long        FltCpnAccSt,    /* (I) Flt coupon accrued start */
                 char        FltDCConv,      /* (I) Flt coupon day count conv*/
                 char        SimpOrComp,     /* (I) Flt rt Simple/Cpounding  */
                 long        ExerFlag,       /* (I) Exercise flag            */
                 double      Strike,         /* (I) Current strike           */
                 char        OptStub,        /* (I) Option stub rule         */
                 char        ArrearsReset,   /* (I) 'Y' if set-in-arreas     */
                 long        CurrentDate,    /* (I) Current date             */
                 int         t,              /* (I) Current time period      */
                 int         T,              /* (I) Total number of period   */
                 int         DCurve,         /* (I) Discount curve           */
                 int         DMode,          /* (I) DEV mode                 */
                 HYB3_DEV_DATA   *dev_data,       /* (I) Hyb3_Dev data structure       */
                 HYB3_TREE_DATA  *tree_data)      /* (I) Tree data structure      */
{



    
    double      *SwapL = NULL;   
    double      *BondL = NULL;   
    double      *FloaterL = NULL;

    double      *FltIndexL = NULL;    
    double      *FixZeroToPmtL = NULL;
    double      *FltZeroToPmtL = NULL;
    
    
    double       FixAccrued;     /* Fix coupon accrued at current date */
    double       FltAccrued;     /* Flt coupon accrued at current date */
    
    int       Top1,   Bottom1;   /* Tree limits (1rst dim)             */
    int      *Top2,  *Bottom2;   /* Tree limits (2nd dim)              */
    int     **Top3, **Bottom3;   /* Tree limits (3rd dim)              */

    int     i, j, k;             /* Node indices                       */
    int     offset;
    int     status = FAILURE;    /* Error status                       */



    /* Local assignment of tree limits */   
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* If this is not an exercise date discount and return */
    if (!ExerFlag)                                              
    {
        if (Hyb3_Dev (Swap,     /* Disc expd value of a swap starting on the  */
                 t,    /* last exercise date                         */
                 T,
                 DCurve,
                 DMode,
                 dev_data,
                 tree_data) == FAILURE)
        {
            goto RETURN;
            
        }  
    
        status = SUCCESS;
        goto RETURN;
        
    }  


    /* Must check if DMode is compatible with underlying tree */
    if (DMode == DISC_3D_CUPS)
    {
        if ((tree_data->TreeType != TTYPE_FX2IR)   &&
            (tree_data->TreeType != TTYPE_EQF2IR)  &&
            (tree_data->TreeType != TTYPE_EQD2IR))
        {
            DR_Error("Slice dimension chosen would need 3-D tree type! "
                    "(Hyb3_FwdSwap_t)\n");
            goto RETURN;
        }
    }
    if (DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if ( tree_data->TreeType != TTYPE_2IR2F1D )
        {
            DR_Error("Slice dimension chosen would need 3-D tree type! "
                    "(Hyb3_FwdSwap_t)\n");
            goto RETURN;
        }
    }


    /* Calculate dirty price of the swap, i.e. no accruals */
    if (DMode == DISC_1D_NOCUPS)
    {

        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        SwapL    = (double *) Swap + offset;
        BondL    = (double *) Bond + offset;
        FloaterL = (double *) Floater + offset;
        
        for (i = Bottom1; i <= Top1; i ++)
        {
            SwapL[i] = BondL[i] + FloaterL[i] - Strike;                                     
        }  
          
    }
    else if (DMode == DISC_2D_CUPS        ||
             DMode == DISC_2D_NOCUPS      ||
             DMode == DISC_2D_1IR2F_NOCUPS  )
    {
         
        for (i = Bottom1; i <= Top1; i ++)
        {
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            SwapL    = (double *) Swap + offset;   
            BondL    = (double *) Bond + offset;   
            FloaterL = (double *) Floater + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                SwapL[j] = BondL[j] + FloaterL[j] - Strike;  
            }
        } 
                  
    }
    else if (DMode == DISC_3D_CUPS || DISC_3D_2IR2F1D_CUPS)
    {

        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                SwapL    = (double *) Swap + offset;    
                BondL    = (double *) Bond + offset;    
                FloaterL = (double *) Floater + offset; 

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {
                    SwapL[k] = BondL[k] + FloaterL[k] 
                                   - Strike;
                } 
            }
        }     
    }


    /*  Calculate FIXED accrual if within accrual period */            	
    if (FixCpnAccSt <= CurrentDate)
    {   
        if (DrDayCountFraction (FixCpnAccSt,
                                CurrentDate,                                 
                                FixDCConv,
                                &FixAccrued) == FAILURE)
        {
            DR_Error("Unable to calculate dcf for fix accrued. (Hyb3_FwdSwap_t)");
            goto RETURN;
        }        

  
        if (DMode == DISC_1D_NOCUPS)
        {

            offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
            SwapL   = (double *) Swap + offset;
            FixZeroToPmtL = (double *) FixZeroToPmt + offset;
 

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
        else if (DMode == DISC_2D_CUPS       || 
                 DMode == DISC_2D_NOCUPS     ||
                DMode == DISC_2D_1IR2F_NOCUPS  )
        {

            if (FixCpnFlag)           /* Subtract fix cpn just made  */
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);          
                    SwapL   = (double *) Swap + offset;               
                    
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        SwapL[j] -= FixCpnRate * FixCpnOuts * FixCpnDcf;                 
                    }  
                }
            }                                                           
            else if (OptStub == 'B')  /* Bond stub : subtract accrued */
            {                                                                   
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);          
                    SwapL   = (double *) Swap + offset;               
                    
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts;                  
                    } 
                }
            }
            else if (OptStub == 'S')  /* Swap stub: subtract PVed accrued  */
            {                                                            
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2,i,0,t,tree_data);          
                    SwapL   = (double *) Swap + offset;               
                    FixZeroToPmtL = (double *) FixZeroToPmt + offset; 
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        SwapL[j] -= FixCpnRate * FixAccrued * FixCpnOuts 
                                     * FixZeroToPmtL[j];                  
                    } 
                }
            }  /* if then else */                        
        }
        else if(DMode == DISC_3D_CUPS || DISC_3D_2IR2F1D_CUPS)
        {

            if (FixCpnFlag)           /* Subtract fix cpn just made  */
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);            
                        SwapL   = (double *) Swap + offset;                 
                        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            SwapL[k] -= FixCpnRate*FixCpnOuts*FixCpnDcf;
                        }
                    }
                }  
            }                                                           
            else if (OptStub == 'B')  /* Bond stub : subtract accrued */
            {                                                                   
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                        SwapL   = (double *) Swap + offset;              
                        
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            SwapL[k] -= FixCpnRate*FixAccrued*FixCpnOuts;
                        } 
                    }
                } 
            }
            else if (OptStub == 'S')  /* Swap stub: subtract PVed accrued  */
            {                                                            
                for (i = Bottom1; i <= Top1; i ++)
                {
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                        SwapL   = (double *) Swap + offset;              
                        FixZeroToPmtL = (double *) FixZeroToPmt + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                        {
                            SwapL[k] -= FixCpnRate * FixAccrued * 
                                            FixCpnOuts *FixZeroToPmtL[k];
                        }  
                    }
               }
            }  /* if then else */   
                                 
        }  /* if then else on nb factors*/

    } /* if within accrual period */



    /* Calculate FLT accrual if within accrual period and floating leg ON */
    if (FltLegOn)
    {
        if (FltCpnAccSt <= CurrentDate)
        {
            if (DrDayCountFraction (FltCpnAccSt,
                                    CurrentDate,                                 
                                    FltDCConv,
                                    &FltAccrued) == FAILURE)
            {
                DR_Error("Unable to calculate dcf for flt accrued.(Hyb3_FwdSwap_t)");
                goto RETURN;
            }        
 

            if (DMode == DISC_1D_NOCUPS)
            {
                
                offset = Hyb3_Node_Offset(1,0,0,t,tree_data);           
                SwapL         = (double *) Swap + offset;
                FltIndexL     = (double *) FltIndex + offset;                
                FltZeroToPmtL = (double *) FltZeroToPmt + offset;  

                if (SimpOrComp == 'S')
                {
                    if (ArrearsReset == 'Y')
                    {
                        if (FltResetFlag) /* Subtract flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] -= (FltIndexL[i]+FltCpnSpd)
                                            *FltCpnOuts*FltCpnDcf;
                            }  
                        }/* Swap,Bond stub : subtract PVed accrued */
                        else if((OptStub == 'B') || (OptStub == 'S'))
                        {                          
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] -= (FltIndexL[i]+FltCpnSpd*FltZeroToPmtL[i])
                                            *FltCpnOuts*FltAccrued;
                            } 
                        }  /* if then else */ 
                    
                       /* Nothing to do for No-Accrued stub */
                    }
                    else 
                    {
                        if (FltResetFlag) /* Nothing to do, next cp just paid*/
                        {
                        ;
                        }                                                           
                        else if (OptStub == 'B')  /* Bond stub : sub accrued */
                        {                         /*         pay PVed coupon */
                
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] -= (FltIndexL[i]+FltCpnSpd)*FltCpnOuts 
                                    *(FltAccrued - FltCpnDcf*FltZeroToPmtL[i]);
                            } 
                        }
                        else if (OptStub == 'S')  /* Swap stub : add PVed  */
                        {                        
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] += (FltIndexL[i]+FltCpnSpd)*FltCpnOuts 
                                    *(FltCpnDcf - FltAccrued)*FltZeroToPmtL[i];
                            } 
                        }
                        else  /* No-Accrued stub: add est.& PVed cpn */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] += (FltIndexL[i]+FltCpnSpd)*FltCpnOuts 
                                          * FltCpnDcf * FltZeroToPmtL[i];
                            }  
                        }  /* if then else */ 
                    }
                }
                else /* SimpOrComp */
                {
                    if (FltResetFlag)
                    {
                        if (ArrearsReset == 'Y') /* Subtract flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                SwapL[i] -= (pow(1.0 + FltIndexL[i] + FltCpnSpd,
                                                 FltCpnDcf)-1.0)* FltCpnOuts;
                            }  
                        }
                        else /* Nothing to do, next cp just paid*/
                        {                          
                            ;
                        }   
                    
                    }
                    else 
                    {
                        DR_Error("Accrual calculations are not supported \n"
                                "when flt pmt is on a compounding rate."
                                "(Hyb3_FwdSwap_t)\n");
                        goto RETURN;
                    }
                } /* If SimpOrComp */
            }
            else if (DMode == DISC_2D_CUPS        ||
                     DMode == DISC_2D_NOCUPS      ||
                     DMode == DISC_2D_1IR2F_NOCUPS )
            {
                


                if (SimpOrComp == 'S')
                {
                    if (ArrearsReset == 'Y')
                    {
                        if (FltResetFlag)          /* Add flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);                    
                                SwapL         = (double *) Swap + offset;                   
                                FltIndexL     = (double *) FltIndex + offset;               
                                FltZeroToPmtL = (double *) FltZeroToPmt + offset;           

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] -= (FltIndexL[j]+FltCpnSpd) 
                                         * FltCpnOuts * FltCpnDcf;
                                }  
                            }
                        }   /* Swap, Bond stub : add PVed accrued */
                        else if((OptStub == 'B') || (OptStub == 'S'))
                        {                        
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);         
                                SwapL         = (double *) Swap + offset;        
                                FltIndexL     = (double *) FltIndex + offset;    
                                FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] -= (FltIndexL[j]+FltCpnSpd*FltZeroToPmtL[j])
                                             * FltCpnOuts * FltAccrued;
                                } 
                            }
                        }  /* if then else */      
                        
                        /* Nothing to do for No-Accrued stub */
                    }
                    else 
                    {
                        if (FltResetFlag) /* Nothing to do, next cp just paid*/
                        {
                            ;
                        }                                                           
                        else if (OptStub == 'B')  /* Bond stub : add accrued */
                        {                         /*         pay PVed coupon */
                
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);         
                                SwapL         = (double *) Swap + offset;        
                                FltIndexL     = (double *) FltIndex + offset;    
                                FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] -= (FltIndexL[j]+FltCpnSpd)
                                     * FltCpnOuts * (FltAccrued  - 
                                             FltCpnDcf * FltZeroToPmtL[j]);
                                } 
                            }
                        }
                        else if (OptStub == 'S')  /* Swap stub : add PVed */
                        {                        
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);         
                                SwapL         = (double *) Swap + offset;        
                                FltIndexL     = (double *) FltIndex + offset;    
                                FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] += (FltIndexL[j]+FltCpnSpd)
                                         * (FltCpnDcf - FltAccrued) *FltCpnOuts
                                         *  FltZeroToPmtL[j];
                                } 
                            }
                        }
                        else  /* No-Accrued stub: add est.& PVed cpn */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                               offset = Hyb3_Node_Offset(2,i,0,t,tree_data);         
                               SwapL         = (double *) Swap + offset;        
                               FltIndexL     = (double *) FltIndex + offset;    
                               FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] += (FltIndexL[j]+FltCpnSpd)
                                             * FltCpnOuts
                                             * FltCpnDcf * FltZeroToPmtL[j];
                                } 
                            }
                        }  /* if then else */ 
                    }
                }
                else /* SimpOrComp */
                {
                    if (FltResetFlag)
                    {
                        if (ArrearsReset == 'Y')          /* Add flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);         
                                SwapL         = (double *) Swap + offset;        
                                FltIndexL     = (double *) FltIndex + offset;    
                                FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    SwapL[j] -= (pow(1.+FltIndexL[j]
                                                    +FltCpnSpd,FltCpnDcf)-1.)
                                                  * FltCpnOuts;
                                } 
                            } 
                        }
                        else /* Nothing to do, next cp just paid*/
                        {                        
                            ;
                        }  /* if then else */      
                    }
                    else 
                    {
                        DR_Error("Accrual calculations are not supported \n"
                                "when flt pmt is on a compounding rate."
                                "(Hyb3_FwdSwap_t)\n");
                        goto RETURN;
                    }
                }


            }
            else if (DMode == DISC_3D_CUPS || DISC_3D_2IR2F1D_CUPS)
            {
                
        
                if (SimpOrComp)
                {
                    if (ArrearsReset == 'Y')
                    {
                        if (FltResetFlag)          /* Add flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                                    SwapL         = (double *) Swap + offset;        
                                    FltIndexL     = (double *) FltIndex + offset;    
                                    FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] -= (FltIndexL[k] +
                                            FltCpnSpd) * FltCpnOuts * FltCpnDcf;
                                    } 
                                }
                            } 
                        }  /* Swap, Bond stub : add PVed accrued */
                        else if((OptStub == 'B') || (OptStub == 'S'))
                        {                        
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                                    SwapL         = (double *) Swap + offset;        
                                    FltIndexL     = (double *) FltIndex + offset;    
                                    FltZeroToPmtL = (double *) FltZeroToPmt + offset;


                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] -= (FltIndexL[k] 
                                            + FltCpnSpd*FltZeroToPmtL[k])
                                            * FltCpnOuts * FltAccrued;
                                    }
                                }
                            } 
                        }  /* if then else */      
                        
                        /* Nothing to do for No-Accrued stub */
                    }
                    else 
                    {
                        if (FltResetFlag) /* Nothing to do, next cp just made*/
                        {
                        ;
                        }                                                           
                        else if (OptStub == 'B')  /* Bond stub : add accrued */
                        {                         /*         pay PVed coupon */
                
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                                    SwapL         = (double *) Swap + offset;        
                                    FltIndexL     = (double *) FltIndex + offset;    
                                    FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] -= (FltIndexL[k]
                                                           +FltCpnSpd)
                                                * FltCpnOuts * (FltAccrued 
                                                - FltCpnDcf 
                                                * FltZeroToPmtL[k]);
                                    }  
                                }
                            }
                        }
                        else if (OptStub == 'S')  /* Swap stub : add PVed cp */
                        {                        
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                                    SwapL         = (double *) Swap + offset;        
                                    FltIndexL     = (double *) FltIndex + offset;    
                                    FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] += (FltIndexL[k]
                                                           + FltCpnSpd)
                                                * FltCpnOuts 
                                                * (FltCpnDcf - FltAccrued) 
                                                *  FltZeroToPmtL[k];
                                    } 
                                }
                            }
                        }
                        else  /* No-Accrued stub: add est.& PVed cpn */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] += (FltIndexL[k] + FltCpnSpd)
                                                * FltCpnOuts * FltCpnDcf 
                                                * FltZeroToPmtL[k];
                                    } 
                                }
                            }
                        }  /* if then else */ 

                    } /* Arrears or advance */
                }
                else /* SimpOrComp */
                {
                    if (FltResetFlag)
                    {
                        if (ArrearsReset == 'Y')          /* Add flt cpn just made */
                        {
                            for (i = Bottom1; i <= Top1; i ++)
                            {
                                for (j = Bottom2[i]; j <= Top2[i]; j++)
                                {
                                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);         
                                    SwapL         = (double *) Swap + offset;        
                                    FltIndexL     = (double *) FltIndex + offset;    
                                    FltZeroToPmtL = (double *) FltZeroToPmt + offset;

                                    for (k=Bottom3[i][j]; k <= Top3[i][j]; k++)
                                    {
                                        SwapL[k] -= (pow(1.+ FltIndexL[k]
                                                               +FltCpnSpd ,
                                             FltCpnDcf) - 1.0)* FltCpnOuts;
                                    }  
                                }
                            }
                        } 
                        else /* Nothing to do, next cp just made*/
                        {                        
                            ;
                        }  
                    }
                    else 
                    {
                        DR_Error("Accrual calculations are not supported \n"
                                "when flt pmt is on a compounding rate."
                                "(Hyb3_FwdSwap_t)\n");
                        goto RETURN;

                    } /* Arrears or advance */

                } /* SimpOrComp */
        
            }  /* if then else DMode type */

        }

    }




    status = SUCCESS;

  RETURN:

    return (status);

}  /* Hyb3_FwdSwap_t */
