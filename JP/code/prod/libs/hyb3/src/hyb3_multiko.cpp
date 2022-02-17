/****************************************************************************/
/*      Calculation of a knock-out option price in the lattice.             */
/****************************************************************************/
/*     MULTIKO.c                                                            */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"




/*****  Hyb3_MultiKo_t  ***********************************************************/
/*                                                                           
 *    Calculate the knock-out option price and smoothed knock-out option
 *    price in the lattice using multi-smoothing. It does not include 
 *    discounting.
 *
 *    The price slice and the smoothed price slice must be the
 *    same dimension, but the rebate and the index can be of 
 *    variable dimension.
 */


int    Hyb3_MultiKo_t(
           TSLICE       SmoothKo,    /* (I/O) Smoothed ko option prices      */
           TSLICE       KoOpt,       /* (I/O) Unsmoothed ko option prices    */
           int          IndexDim,    /* (I) Dim of var in IndexPtr (1 or 2)  */
           TSLICE       IndexPtr,    /* (I) Index prices at time pint        */
           long         KoFlag,      /* (I) TRUE if time point is a KO date  */
           double       LowBarrier,  /* (I) Lower barrier levels             */
           double       HighBarrier, /* (I) Higher barrier levels            */
           int          RebateDim,   /* (I) Dim of var pointed by RebatePtr  */
           TSLICE       RebatePtr,   /* (I) Rebate variable pointer          */
           char         IoO,         /* (I) KO 'I'nside or 'O'utside barrier */
           char         Smoothing,   /* (I) Smoothing ('Y' or 'N')           */
           int          t,           /* (I) Current time period              */
           int          DMode,       /* (I) Also used for dim of slices      */
           HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure              */
{



    double  *KoOptL    = NULL;    /* Local for slice addressing              */
    double  *SmoothKoL = NULL;
    double   x;                   /* Smoothed value                          */
    double   IndexStep;           /* Maximum difference between index values */
    double   RebateAtNode;        /* Rebate value at the node (i,j)          */
    double   IndexAtNode;         /* Barrier index value at the node (i,j)   */

    int     Top1,  Bottom1;       /* Limits of the tree (1rst dimension)     */
    int    *Top2, *Bottom2;       /* Limits of the tree (2nd dimension)      */
    int   **Top3,**Bottom3;       /* Limits of the tree (3rd dimension)      */
    int     i, j, k;              /* Node indices                            */
    int     offset;
    int     status = FAILURE;     /* Error status = FAILURE initially        */

        
    Top1    = tree_data->Top1[t];                   
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];                   
    Bottom1 = tree_data->Bottom1[t];                    
    Bottom2 = tree_data->Bottom2[t];    
    Bottom3 = tree_data->Bottom3[t];                      


    if (!KoFlag)  /* Nothing to do */
    {
        return (SUCCESS);
    }  


    /* Quick check */
    if (LowBarrier > HighBarrier * (1. - BARRIER_TOL) )
    {
        DR_Error("Low barrier must be less than high barrier! "
                "(Hyb3_MultiKo_t)\n");
        goto RETURN;
    }


    if (DMode == DISC_1D_NOCUPS)
    {        
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        KoOptL = (double *)KoOpt + offset;
        SmoothKoL = (double *)SmoothKo + offset;


        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                              RebatePtr,
                                              i,0,0,t,tree_data);

                IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,t,tree_data);

                if (IoO == 'I')
                {
                   /* Knock-out KoOpt AND SmoothKo with no smoothing */        
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i]    = RebateAtNode;
                        SmoothKoL[i] = RebateAtNode;
                    }
                }
                else /* IoO is 'O' */
                {
                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KoOptL[i]    = RebateAtNode;
                        SmoothKoL[i] = RebateAtNode;
                    }
                } /* if then else */
            } /* for i */
        }
        else /* smoothing is ON */
        {
            for (i = Bottom1; i <= Top1; i ++)
            {

                RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                              RebatePtr,
                                              i,0,0,t,tree_data);

                IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,t,tree_data);

                IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                          IndexDim,
                                          i, 0, 0,
                                          t,
                                          tree_data);

                if (IoO == 'I')
                {                    
                    /* Knock-out SmoothKo using new smoothing */ 
                    if (Hyb3_Multi_Smooth(
                                &(x),
                                RebateAtNode, /* up value            */
                                RebateAtNode, /* unsmooth up value   */
                                SmoothKoL[i], /* down value          */
                                KoOptL[i],    /* unsmooth down value */
                                IndexAtNode,
                                LowBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                    if (Hyb3_Multi_Smooth(
                                &(SmoothKoL[i]), 
                                SmoothKoL[i], /* up value            */
                                KoOptL[i],    /* unsmooth up value   */
                                x,            /* down value          */
                                x,            /* unsmooth down value */
                                IndexAtNode,
                                HighBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;                                   
                    }

                   /* Knock-out KoOpt with no smoothing */        
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i] = RebateAtNode;               
                    }
                }
                else /* IoO is 'O' */
                {
                    if (Hyb3_Multi_Smooth(
                                &(x),
                                SmoothKoL[i], /* up value            */     
                                KoOptL[i],    /* unsmooth up value   */
                                RebateAtNode, /* down value          */  
                                RebateAtNode, /* unsmooth down value */ 
                                IndexAtNode,
                                LowBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }

                    if (Hyb3_Multi_Smooth(
                                &(SmoothKoL[i]),
                                RebateAtNode, /* up value            */ 
                                RebateAtNode, /* unsmooth up value   */
                                x,            /* down value          */
                                KoOptL[i],    /* unsmooth down value */
                                IndexAtNode,
                                HighBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }

                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KoOptL[i] = RebateAtNode;
                    }                
                }  /* if then else */
            }  /* for i */
        } /* smoothing ON/OFF */
    }
    if (DMode == DISC_2D_CUPS   || 
        DMode == DISC_2D_NOCUPS ||
        DMode == DISC_2D_1IR2F_NOCUPS)
    {
        if (Smoothing == 'N')
        {
             for (i = Bottom1; i <= Top1; i ++)
             {    
                 offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                 KoOptL = (double *)KoOpt + offset;
                 SmoothKoL = (double *)SmoothKo + offset;

                 for (j = Bottom2[i]; j <= Top2[i]; j++)
                 {	
                    RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                                  RebatePtr,
                                                  i,j,0,t,tree_data);

                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,t,tree_data);

                    if (IoO == 'I')
                    {
                        if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                            (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KoOptL[j]    = RebateAtNode;
                            SmoothKoL[j] = RebateAtNode;
                        }  
                    }
                    else /* IoO IS 'O' */
                    {
                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KoOptL[j]    = RebateAtNode;
                            SmoothKoL[j] = RebateAtNode;
                        }                            
                    }  /* if then else */
                }  /* for j */
             }  /* for i */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                KoOptL = (double *)KoOpt + offset;
                SmoothKoL = (double *)SmoothKo + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                                  RebatePtr,
                                                  i,j,0,t,tree_data);

                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,t,tree_data);

                    IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                              IndexDim,
                                              i, j, 0,
                                              t,
                                              tree_data);

                    if (IoO == 'I')
                    {
                        if (Hyb3_Multi_Smooth(
                                    &(x),
                                    RebateAtNode, /* up value            */
                                    RebateAtNode, /* unsmooth up value   */
                                    SmoothKoL[j], /* down value          */
                                    KoOptL[j],    /* unsmooth down value */
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        if (Hyb3_Multi_Smooth(
                                    &(SmoothKoL[j]), /* up value            */
                                    SmoothKoL[j],    /* unsmooth up value   */
                                    KoOptL[j],       /* down value          */
                                    x,               /* unsmooth down value */
                                    x,
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;                                   
                        }

                        /* Knock-out KoOpt with no smoothing */        
                        if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                            (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KoOptL[j] = RebateAtNode;               
                        }                        
                    }	
                    else /* IoO is 'O' */
                    {	
                        if (Hyb3_Multi_Smooth(
                                    &(x),
                                    SmoothKoL[j], /* up value            */
                                    KoOptL[j],    /* unsmooth up value   */
                                    RebateAtNode, /* down value          */
                                    RebateAtNode, /* unsmooth down value */
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Hyb3_Multi_Smooth(
                                    &(SmoothKoL[j]),
                                    RebateAtNode, /* up value            */
                                    RebateAtNode, /* unsmooth up value   */
                                    x,            /* down value          */
                                    KoOptL[j],    /* unsmooth down value */
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KoOptL[j] = RebateAtNode;
                        }
                    } /* if then else */	
                } /* for j */	
            } /* for i */	
        } /* smoothing ON/OFF */	
    } /* if */
    if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    KoOptL = (double *)KoOpt + offset;
                    SmoothKoL = (double *)SmoothKo + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                                      RebatePtr,
                                                      i,j,k,t,tree_data);

                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      IndexPtr,
                                                      i,j,k,t,tree_data);

                        if (IoO == 'I')
                        {
                            /* Knock-out KoOpt AND SmoothKo with no smoothing */
                            if ((IndexAtNode > LowBarrier *(1.-BARRIER_TOL)) &&
                                (IndexAtNode < HighBarrier*(1.+BARRIER_TOL)))
                            {        
                                KoOptL[k]    = RebateAtNode;               
                                SmoothKoL[k] = RebateAtNode;
                            }  
                        }
                        else /* IoO IS 'O' */
                        {
                            if ((IndexAtNode < LowBarrier *(1.+BARRIER_TOL)) ||
                                (IndexAtNode > HighBarrier*(1.-BARRIER_TOL)))
                            {        
                                KoOptL[k]    = RebateAtNode;
                                SmoothKoL[k] = RebateAtNode; 
                            }   
                        } /* if then else */
                    } /* for k */
                } /* for j */    
            } /* for i */
        }           
        else /* smoothing is ON */
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    KoOptL = (double *)KoOpt + offset;
                    SmoothKoL = (double *)SmoothKo + offset; 
                    
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {

                        RebateAtNode = Hyb3_GetValueAtNode(RebateDim, 
                                                  RebatePtr,
                                                  i,j,k,t,tree_data);

                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,k,t,tree_data);

                        IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                                  IndexDim,
                                                  i, j, k,
                                                  t,
                                                  tree_data);

                        if (IoO == 'I')
                        {
                            /* Knock-out SmoothKo using new smoothing */ 
                            if (Hyb3_Multi_Smooth(
                                        &(x),
                                        RebateAtNode, /* up value            */
                                        RebateAtNode, /* unsmooth up value   */
                                        SmoothKoL[k], /* down value          */
                                        KoOptL[k],    /* unsmooth down value */
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            if (Hyb3_Multi_Smooth(
                                        &(SmoothKoL[k]), 
                                        SmoothKoL[k], /* up value            */
                                        KoOptL[k],    /* unsmooth up value   */
                                        x,            /* down value          */
                                        x,            /* unsmooth down value */
                                        IndexAtNode,
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;                                   
                            }
                        
                            /* Knock-out KoOpt with no smoothing */        
                            if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) &&
                                (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                KoOptL[k] = RebateAtNode;               
                            }                        
                        }
                        else /* IoO is 'O' */
                        {
                            if (Hyb3_Multi_Smooth(
                                        &(x),
                                        SmoothKoL[k], /* up value            */
                                        KoOptL[k],    /* unsmooth up value   */
                                        RebateAtNode, /* down value          */
                                        RebateAtNode, /* unsmooth down value */
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            if (Hyb3_Multi_Smooth(
                                        &(SmoothKoL[k]),
                                        RebateAtNode, /* up value            */
                                        RebateAtNode, /* unsmooth up value   */
                                        x,            /* down value          */
                                        KoOptL[k],    /* unsmooth down value */
                                        IndexAtNode,
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        
                            if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                                (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                KoOptL[k] = RebateAtNode;
                            }
                        } /* if then else */
                    } /* for k */	
                } /* for j */	
            } /* for i */	
        } /* smoothing ON/OFF */	
    }  /* if */    
    

    status = SUCCESS;

    RETURN:


    return (status);

}  /* Hyb3_MultiKo_t */







/*****  Hyb3_MultiKi_t  ***********************************************************/
/*                                                                           
 *  Calculate the knock-in  option price and smoothed knock-in option price
 *  in the lattice using multi-smoothing.
 *
 *  This is used for path-dependent knock-in structures when the in/out parity 
 *  does not hold. It includes discounting.
 */


int    Hyb3_MultiKi_t(
           TSLICE       SmoothKi,    /* (I/O) Smoothed ki option prices      */
           TSLICE       KiOpt,       /* (I/O) Unsmoothed ki option prices    */
           TSLICE       Opt,         /* (I) Standard option                  */
           int          IndexDim,    /* (I) Dim of var pointed by IndexPtr   */
           TSLICE       IndexPtr,    /* (I) Index prices at current period   */
           long         KiFlag,      /* (I) TRUE if time point is a KI date  */
           double       LowBarrier,  /* (I) Lower barrier levels             */
           double       HighBarrier, /* (I) Higher barrier levels            */
           char         IoO,         /* (I) KO 'I'nside or 'O'utside barrier */
           char         Smoothing,   /* (I) Smoothing ('Y' or 'N')           */
           int          t,           /* (I) Current time period              */
           int          T,           /* (I) Total number of period           */
           int          DCurve,      /* (I) Discount curve                   */
           int          DMode,       /* (I) Also used for dim of slices      */
           HYB3_DEV_DATA    *dev_data,    /* (I) Hyb3_Dev data structure               */
           HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure              */
{


    double  *SmoothKiL = NULL;       /* Local slice pointers                 */
    double  *KiOptL = NULL;
    double  *OptL   = NULL;

    double  x;                    /* Smoothed value                          */
    double  IndexStep;            /* Maximum difference between index values */
    double  IndexAtNode;          /* Barrier index value at the node (i,j)   */

    int     Top1,  Bottom1;       /* Limits of the tree (1rst dimension)     */
    int    *Top2, *Bottom2;       /* Limits of the tree (2nd dimension)      */
    int   **Top3,**Bottom3;       /* Limits of the tree (3rd dimension)      */
    int     i, j, k;              /* Node indices                            */
    int     offset;
    int     status = FAILURE;     /* Error status = FAILURE initially        */

        
    Top1    = tree_data->Top1[t];                   
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];                   
    Bottom1 = tree_data->Bottom1[t];                    
    Bottom2 = tree_data->Bottom2[t];    
    Bottom3 = tree_data->Bottom3[t];                      


    if (Hyb3_Dev (SmoothKi,
             t,
             T,
             DCurve,
             DMode,
             dev_data,
             tree_data) == FAILURE)
    {
        goto RETURN;
        
    }  


    if (Hyb3_Dev (KiOpt,
             t,
             T,
             DCurve,
             DMode,
             dev_data,
             tree_data) == FAILURE)
    {
        goto RETURN;
        
    }  


    if (!KiFlag)  /* Nothing to do */
    {
        return (SUCCESS);
    }  


    /* Quick check */
    if (LowBarrier > HighBarrier * (1. - BARRIER_TOL) )
    {
        DR_Error("Low barrier must be less than high barrier! "
                "(Hyb3_MultiKi_t)\n");
        goto RETURN;
    }


    if (DMode == DISC_1D_NOCUPS)
    {

        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        SmoothKiL = (double *)SmoothKi + offset;
        KiOptL = (double *)KiOpt + offset;
        OptL = (double *)Opt   + offset;
       

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {                   
                IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,t,tree_data);

                if (IoO == 'I')
                {
                    /* Knock-in KiOpt AND SmoothKi with no smoothing */
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KiOptL[i]    = OptL[i];
                        SmoothKiL[i] = OptL[i];
                    }  
                }
                else /* IoO IS 'O' */
                {
                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KiOptL[i]    = OptL[i];
                        SmoothKiL[i] = OptL[i];
                    }        
                } /* if then else */
            } /* for i */
        }	            
        else /* smoothing is ON */
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,t,tree_data);

                IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                          IndexDim,
                                          i, 0, 0,
                                          t,
                                          tree_data);

                if (IoO == 'I') 
                {
                    /* Knock-in SmoothKi using new smoothing */
                    if (Hyb3_Multi_Smooth(
                                &(x),
                                OptL[i],      /* up value            */
                                OptL[i],      /* unsmooth up value   */
                                SmoothKiL[i], /* down value          */
                                KiOptL[i],    /* unsmooth down value */
                                IndexAtNode,
                                LowBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                    if (Hyb3_Multi_Smooth(
                                &(SmoothKiL[i]),
                                SmoothKiL[i], /* up value            */
                                KiOptL[i],    /* unsmooth up value   */
                                x,            /* down value          */
                                x,            /* unsmooth down value */
                                IndexAtNode,
                                HighBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }

                    /* Knock-in KiOpt with no smoothing */
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];
                    }
                }	
                else /* IoO is 'O' */
                {	
                    if (Hyb3_Multi_Smooth(
                                &(x),
                                SmoothKiL[i], /* up value            */
                                KiOptL[i],    /* unsmooth up value   */
                                OptL[i],      /* down value          */
                                OptL[i],      /* unsmooth down value */
                                IndexAtNode,   
                                LowBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }

                    if (Hyb3_Multi_Smooth(
                                &(SmoothKiL[i]),
                                OptL[i],      /* up value            */
                                OptL[i],      /* unsmooth up value   */
                                x,            /* down value          */
                                KiOptL[i],    /* unsmooth down value */
                                IndexAtNode,
                                HighBarrier,
                                IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }

                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];
                    }

                }  /* if then else */	
            }  /* for i */	
        }  /* smoothing ON/OFF */	
    }
    if (DMode == DISC_2D_CUPS   || 
        DMode == DISC_2D_NOCUPS ||
        DMode == DISC_2D_1IR2F_NOCUPS)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {     
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                
                SmoothKiL = (double *)SmoothKi + offset;                
                KiOptL = (double *)KiOpt + offset;
                OptL = (double *)Opt   + offset;
                  
              
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {	
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,t,tree_data);

                    if (IoO == 'I')
                    {
                        /* Knock-in KiOpt AND SmoothKi with no smoothing */
                        if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) &&
                            (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KiOptL[j]    = OptL[j];
                            SmoothKiL[j] = OptL[j];
                        }  
                    }
                    else /* IoO IS 'O' */
                    {
                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KiOptL[j]    = OptL[j];
                            SmoothKiL[j] = OptL[j];
                        }                            
                    }  /* if then else */
                }  /* for j */
            }  /* for i */
        }           
        else /* smoothing is ON */
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
             
                SmoothKiL = (double *)SmoothKi + offset;
                KiOptL = (double *)KiOpt + offset;
                OptL = (double *)Opt   + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,t,tree_data);

                    IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                              IndexDim,
                                              i, j, 0,
                                              t,
                                              tree_data);
                                             
                    if (IoO == 'I')
                    {
                        /* Knock-in SmoothKi using new smoothing */
                        if (Hyb3_Multi_Smooth(
                                    &(x),
                                    OptL[j],      /* up value            */
                                    OptL[j],      /* unsmooth up value   */
                                    SmoothKiL[j], /* down value          */
                                    KiOptL[j],    /* unsmooth down value */
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        if (Hyb3_Multi_Smooth(
                                    &(SmoothKiL[j]),
                                    SmoothKiL[j], /* up value            */
                                    KiOptL[j],    /* unsmooth up value   */
                                    x,            /* down value          */
                                    x,            /* unsmooth down value */
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        /* Knock-in KiOpt with no smoothing */
                        if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) &&
                            (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];
                        }
                    }
                    else /* IoO is 'O' */
                    {
                        if (Hyb3_Multi_Smooth(
                                    &(x),
                                    SmoothKiL[j], /* up value            */
                                    KiOptL[j],    /* unsmooth up value   */
                                    OptL[j],      /* down value          */
                                    OptL[j],      /* unsmooth down value */
                                    IndexAtNode,   
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        if (Hyb3_Multi_Smooth(
                                    &(SmoothKiL[j]),
                                    OptL[j],      /* up value            */
                                    OptL[j],      /* unsmooth up value   */
                                    x,            /* down value          */
                                    KiOptL[j],    /* unsmooth down value */
                                    IndexAtNode,   
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }

                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];
                        }

                    }  /* if then else */	
                }  /* for j */	
            }  /* for i */	
        }  /* if then else */	
    }
    if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    
                    SmoothKiL = (double *)SmoothKi + offset;
                    KiOptL = (double *)KiOpt + offset;
                    OptL = (double *)Opt   + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      IndexPtr,
                                                      i,j,k,t,tree_data);

                        if (IoO == 'I')
                        {
                            /* Knock-in KiOpt AND SmoothKi with no smoothing */
                            if((IndexAtNode > LowBarrier *(1.- BARRIER_TOL)) &&
                               (IndexAtNode < HighBarrier*(1.+BARRIER_TOL)))
                            {        
                                KiOptL[k]    = OptL[k];               
                                SmoothKiL[k] = OptL[k];
                            }  
                        }
                        else /* IoO IS 'O' */
                        {
                            if ((IndexAtNode < LowBarrier*(1.+BARRIER_TOL)) ||
                                (IndexAtNode > HighBarrier*(1.-BARRIER_TOL)))
                            {        
                                KiOptL[k]    = OptL[k];
                                SmoothKiL[k] = OptL[k];
                            }
                            
                        }  /* if then else */
                    }  /* for k */
                }  /* for j */    
            }  /* for i */
        }           
        else /* smoothing is ON */
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    SmoothKiL = (double *)SmoothKi + offset;
                    KiOptL = (double *)KiOpt + offset;
                    OptL = (double *)Opt   + offset;
                 
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      IndexPtr,
                                                      i,j,k,t,tree_data);

                        IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                                  IndexDim,
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                                             
                        if (IoO == 'I')
                        {                                                                                           
                            /* Knock-in SmoothKi using new smoothing */
                            if (Hyb3_Multi_Smooth(
                                        &(x),
                                        OptL[k],      /* up value            */
                                        OptL[k],      /* unsmooth up value   */
                                        SmoothKiL[k], /* down value          */
                                        KiOptL[k],    /* unsmooth down value */
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            if (Hyb3_Multi_Smooth(
                                        &(SmoothKiL[k]),
                                        SmoothKiL[k], /* up value            */
                                        KiOptL[k],    /* unsmooth up value   */
                                        x,            /* down value          */
                                        x,            /* unsmooth down value */
                                        IndexAtNode,
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            /* Knock-in KiOpt with no smoothing */
                            if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                                (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                KiOptL[k] = OptL[k];
                            }
                        }	
                        else /* IoO is 'O' */
                        {	
                            if (Hyb3_Multi_Smooth(
                                        &(x),
                                        SmoothKiL[k], /* up value            */
                                        KiOptL[k],    /* unsmooth up value   */
                                        OptL[k],      /* down value          */
                                        OptL[k],      /* unsmooth down value */
                                        IndexAtNode,   
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            if (Hyb3_Multi_Smooth(
                                        &(SmoothKiL[k]),
                                        OptL[k],      /* up value            */
                                        OptL[k],      /* unsmooth up value   */
                                        x,            /* down value          */
                                        KiOptL[k],    /* unsmooth down value */
                                        IndexAtNode,   
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }

                            if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                                (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                KiOptL[k] = OptL[k];
                            }

                        } /* if then else */
                    } /* for k */	
                } /* for j */	
            } /* for i */	
         } /* if then else */	
    } /* if */
        

    status = SUCCESS;

    RETURN:

    return (status);


}  /* Hyb3_MultiKi_t */











/*****  Hyb3_Multi_Smooth  ********************************************************/
/*
*       Smoothing for use with multi-knockout/in (cf. corresponding memo).
*/
int     Hyb3_Multi_Smooth(double  *SmoothValue,   /* (O) Output smoothed value   */
                     double  UpValue,        /* (I) Up value                */
                     double  USUpValue,      /* (I) Unsmoothed up value     */
                     double  DownValue,      /* (I) Down value              */
                     double  USDownValue,    /* (I) Unsmoothed down value   */
                     double  Index,          /* (I) Index value             */
                     double  Barrier,        /* (I) Barrier level           */
                     double  MaxDiff)        /* (I) Maximum index spread    */
{
    double
            x, y;
    int
            status = FAILURE;   /* Error status = FAILURE initially */
                     
                     
    x = (Index - Barrier) / MaxDiff;
        
    if (x < -1.)        /* Step function = 0 below step */
    { 
        *SmoothValue = DownValue;
    }
    else if (x > 1.)    /* Step function = 1 above step */
    { 
        *SmoothValue = UpValue;
    }
    else	
    {	
        y = 0.5 + x / 16. * (15. + x * x * (-10. + 3. * x * x));
        
        *SmoothValue = (USUpValue - USDownValue) * y + USDownValue;
                
    }  /* if then else */	


    status = SUCCESS;

    return (status);

}  /* Hyb3_Multi_Smooth */



