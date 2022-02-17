/****************************************************************************/
/*      Calculation of a knock-out option price in the lattice.             */
/****************************************************************************/
/*      KOOPT.c                                                             */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

/*****  Hyb4_KoOption_t  *********************************************************/
/*                                                                           
 *    Calculate the knock-out option price in the lattice using
 *    Arnon's smoothing. It does not include discounting.
 *
 *    The price slice must be 2-D, but the rebate and the index
 *    can be of variable dimension.
 */
int    Hyb4_KoOption_t(
           TSLICE       KoOpt,       /* (I/O) KO option prices at t point    */
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
           HYB4_TREE_DATA *tree_data)   /* (I) Tree data structure              */
{
    double  *KoOptL = NULL;       /* Local for slice addressing              */
    double   x;                   /* Smoothed value                          */
    double   IndexStep;           /* Maximum difference between index values */
    double   RebateAtNode = 0;    /* Rebate value at the node (i,j)          */
    double   IndexAtNode;         /* Barrier index value at the node (i,j)   */
    double  *RebateAtNodeSlice = NULL, 
            *IndexAtNodeSlice = NULL;

    int     Top1,   Bottom1;       /* Limits of the tree (1rst dimension)     */
    int    *Top2,  *Bottom2;       /* Limits of the tree (2nd dimension)      */
    int   **Top3, **Bottom3;       /* Limits of the tree (3rd dimension)      */
    int  ***Top4,***Bottom4;       /* Limits of the tree (3rd dimension)      */

    int     i, j, k, L;            /* Node indices                            */
    int     offset;
    int     status = FAILURE;     /* Error status = FAILURE initially        */

    Top1    = tree_data->iMax[t];
    Top2    = tree_data->jMax[t];
    Top3    = tree_data->kMax[t];
    Top4    = tree_data->LMax[t];

    Bottom1 = tree_data->iMin[t];
    Bottom2 = tree_data->jMin[t];
    Bottom3 = tree_data->kMin[t];
    Bottom4 = tree_data->LMin[t];

    if (!KoFlag)  /* Nothing to do */
    {
        return (SUCCESS);
    }  



    if (DMode == DISC_1D_NOCUPS)
    {


        offset = tree_data->NodeOffset0[t];
        KoOptL = (double *)KoOpt + offset;

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                RebateAtNode = Hyb4_GetValueAtNode(RebateDim, 
                                              RebatePtr,
                                              i,0,0,0,t,tree_data);

                IndexAtNode  = Hyb4_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,0,t,tree_data);

                if (IoO == 'I')
                {
                    if ((IndexAtNode > LowBarrier  * (1. - HYB4_BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + HYB4_BARRIER_TOL)))
                    {        
                        KoOptL[i] = RebateAtNode;               
                    }  
                }
                else /* IoO IS 'O' */
                {
                    if ((IndexAtNode < LowBarrier  * (1. + HYB4_BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - HYB4_BARRIER_TOL)))
                    {        
                        KoOptL[i] = RebateAtNode;
                    }        
                }  /* if then else */
            } /* For i */
        }                
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {


                RebateAtNode = Hyb4_GetValueAtNode(RebateDim, 
                                              RebatePtr,
                                              i,0,0,0,t,tree_data);

                IndexAtNode  = Hyb4_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,0,t,tree_data);

                IndexStep = Hyb4_GetIndexStep (IndexPtr,
                                          IndexDim,
                                          i, 0, 0, 0,
                                          t,
                                          tree_data);
                                
                /* Knock-out inside the two barriers */
                if (IoO == 'I') 
                {                                  
                    /* 'up' value = rebate        */
                    /* 'down' value = live option */
                    if (Smooth_Step(&(x),
                                    RebateAtNode,  
                                    KoOptL[i],     
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }                                                 
                                
                    /* 'up' value = live option      */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step(&(KoOptL[i]),
                                    KoOptL[i],   
                                    x,   
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;                                   
                    }  
                }    
                else                /* Knock-out outside the two barriers */
                {    
                    /* 'up' value = live option */
                    /* 'down' value = rebate    */
                    if (Smooth_Step(&(x),
                                    KoOptL[i],      
                                    RebateAtNode,   
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }                                             
                                
                    /* 'up' value = rebate           */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step(&(KoOptL[i]),
                                    RebateAtNode,  
                                    x,   
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                }  /* if then else */    
            }  /* for i */    
        }  /* if then else on Smoothing ON */    
    }
    if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
    {
        

        if (Smoothing == 'N')
        {
             for (i = Bottom1; i <= Top1; i ++)
             {    
             
                 offset = tree_data->NodeOffset1[t][i];
                 KoOptL = (double *)KoOpt + offset;
                                   
                 for (j = Bottom2[i]; j <= Top2[i]; j++)
                 {    
                    RebateAtNode = Hyb4_GetValueAtNode(RebateDim, 
                                                  RebatePtr,
                                                  i,j,0,0,t,tree_data);

                    IndexAtNode  = Hyb4_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,0,t,tree_data);

                    if (IoO == 'I')
                    {
                        if ((IndexAtNode > LowBarrier  * (1. - HYB4_BARRIER_TOL)) && 
                            (IndexAtNode < HighBarrier * (1. + HYB4_BARRIER_TOL)))
                        {        
                                KoOptL[j] = RebateAtNode;               
                        }  
                    }
                    else /* IoO IS 'O' */
                    {
                        if ((IndexAtNode < LowBarrier  * (1. + HYB4_BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - HYB4_BARRIER_TOL)))
                        {        
                            KoOptL[j] = RebateAtNode;
                        }                            
                    }  /* if then else */
                }  /* for j */
             }  /* for i */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {

                offset = tree_data->NodeOffset1[t][i];
                KoOptL = (double *)KoOpt + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    RebateAtNode = Hyb4_GetValueAtNode(RebateDim, 
                                                  RebatePtr,
                                                  i,j,0,0,t,tree_data);

                    IndexAtNode  = Hyb4_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,0,t,tree_data);

                    IndexStep = Hyb4_GetIndexStep (IndexPtr,
                                              IndexDim,
                                              i, j, 0, 0,
                                              t,
                                              tree_data);
                                                                                               
                    if (IoO == 'I')
                    {                                                                                           
                        if (Smooth_Step(&(x),
                                        RebateAtNode,
                                        KoOptL[j],
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }  /* if */                                                
                                
                        if (Smooth_Step(&(KoOptL[j]),
                                        KoOptL[j],
                                        x,
                                        IndexAtNode,
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;                                   
                        }
                    }    
                    else
                    {    
                        if (Smooth_Step(&(x),
                                        KoOptL[j],
                                        RebateAtNode,
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Smooth_Step(&(KoOptL[j]),
                                        RebateAtNode,
                                        x,
                                        IndexAtNode,
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
    if (DMode == DISC_3D_CUPS)
    {
        if (Smoothing == 'N')
        {
            double lowBar, hiBar;

            if (RebateDim == 0)
                RebateAtNode = (double)RebatePtr[0];
            
            if (IoO == 'I')
            {
                lowBar = LowBarrier*(1.-HYB4_BARRIER_TOL);
                hiBar = HighBarrier*(1.+HYB4_BARRIER_TOL);
            }
            else
            {
                lowBar = LowBarrier *(1.+HYB4_BARRIER_TOL);
                hiBar = HighBarrier*(1.-HYB4_BARRIER_TOL);
            }

            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = tree_data->NodeOffset2[t][i][j];
                    KoOptL = (double *)KoOpt + offset;

                    IndexAtNodeSlice  = (double *)IndexPtr + offset;

                    /* it is very verbose code to repeat the inner loop, but we 
                       want to maximise performance from this loop and thus remove
                       all the code possible from it */
                    if (RebateDim != 0)
                    {
                        RebateAtNodeSlice = (double *)RebatePtr + offset;

                        if (IoO == 'I')
                        {
                            for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                            {    
                                IndexAtNode  = IndexAtNodeSlice[k];

                                if ((IndexAtNode > lowBar) && 
                                    (IndexAtNode < hiBar))
                                {        
                                    KoOptL[k] = RebateAtNodeSlice[k];
                                }
                            }
                        }
                        else 
                        {
                            for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                            {    
                                IndexAtNode  = IndexAtNodeSlice[k];

                                if ((IndexAtNode < lowBar) ||
                                    (IndexAtNode > hiBar))
                                {        
                                    KoOptL[k] = RebateAtNodeSlice[k];
                                }
                            }
                        }
                    }
                    else  /* rebate dimension = 0 -ie. single value rebate */
                    {
                        if (IoO == 'I'  )
                        {
                            for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                            {    
                                IndexAtNode  = IndexAtNodeSlice[k];

                                if ((IndexAtNode > lowBar) && 
                                    (IndexAtNode < hiBar))
                                {        
                                    KoOptL[k] = RebateAtNode;
                                }
                            }
                        }
                        else 
                        {
                            for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                            {    
                                IndexAtNode  = IndexAtNodeSlice[k];

                                if ((IndexAtNode < lowBar) ||
                                    (IndexAtNode > hiBar))
                                {        
                                    KoOptL[k] = RebateAtNode;
                                }
                            }
                        }
                    } /* rebate dim */
                }/* For j */    
            } /* For i */
        }           
        else /* Smoothing required */
        {           
            if (RebateDim == 0)
                RebateAtNode = (double) RebatePtr[0];

            for (i = Bottom1; i <= Top1; i ++)                  
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = tree_data->NodeOffset2[t][i][j];
                    KoOptL = (double *)KoOpt + offset;

                    if (RebateDim != 0)
                        RebateAtNodeSlice = (double *)RebatePtr + offset;

                    IndexAtNodeSlice  = (double *)IndexPtr + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {

                        if (RebateDim != 0)
                            RebateAtNode = RebateAtNodeSlice[k];

                        IndexAtNode  = IndexAtNodeSlice[k];

                        IndexStep = Hyb4_GetIndexStep (IndexPtr,
                                                  IndexDim,
                                                  i, j, k, 0,
                                                  t,
                                                  tree_data);
                                             
                        if (IoO == 'I')
                        {                                                                                           
                            if (Smooth_Step(&(x),
                                            RebateAtNode,
                                            KoOptL[k],
                                            IndexAtNode,
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step(&(KoOptL[k]),
                                            KoOptL[k],
                                            x,
                                            IndexAtNode,
                                            HighBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;                                   
                            }
                        }    
                        else
                        {    
                            if (Smooth_Step(&(x),
                                            KoOptL[k],
                                            RebateAtNode,
                                            IndexAtNode,
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step(&(KoOptL[k]),
                                            RebateAtNode,
                                            x,
                                            IndexAtNode,
                                            HighBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }  /* if then else */
                    } /* For k */    
                }  /* for j */    
            }  /* for i */    
         }  /* if then else */    
    }  /* if */    
    if (DMode == DISC_4D_CUPS)
    {
        if (Smoothing == 'N')
        {
            double lowBar, hiBar;

            if (RebateDim == 0)
                RebateAtNode = (double)RebatePtr[0];
            
            if (IoO == 'I')
            {
                lowBar = LowBarrier*(1.-HYB4_BARRIER_TOL);
                hiBar = HighBarrier*(1.+HYB4_BARRIER_TOL);
            }
            else
            {
                lowBar = LowBarrier *(1.+HYB4_BARRIER_TOL);
                hiBar = HighBarrier*(1.-HYB4_BARRIER_TOL);
            }

            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        KoOptL = (double *)KoOpt + offset;

                        IndexAtNodeSlice  = (double *)IndexPtr + offset;

                        /* it is very verbose code to repeat the inner loop, but we 
                           want to maximise performance from this loop and thus remove
                           all the code possible from it */
                        if (RebateDim != 0)
                        {
                            RebateAtNodeSlice = (double *)RebatePtr + offset;

                            if (IoO == 'I')
                            {
                                for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                                {    
                                    IndexAtNode  = IndexAtNodeSlice[L];

                                    if ((IndexAtNode > lowBar) && 
                                        (IndexAtNode < hiBar))
                                    {        
                                        KoOptL[L] = RebateAtNodeSlice[L];
                                    }
                                }
                            }
                            else 
                            {
                                for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                                {    
                                    IndexAtNode  = IndexAtNodeSlice[L];

                                    if ((IndexAtNode < lowBar) ||
                                        (IndexAtNode > hiBar))
                                    {        
                                        KoOptL[L] = RebateAtNodeSlice[L];
                                    }
                                }
                            }
                        }
                        else  /* rebate dimension = 0 -ie. single value rebate */
                        {
                            if (IoO == 'I')
                            {
                                for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                                {    
                                    IndexAtNode  = IndexAtNodeSlice[L];

                                    if ((IndexAtNode > lowBar) && 
                                        (IndexAtNode < hiBar))
                                    {        
                                        KoOptL[L] = RebateAtNode;
                                    }
                                }
                            }
                            else 
                            {
                                for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                                {    
                                    IndexAtNode  = IndexAtNodeSlice[L];

                                    if ((IndexAtNode < lowBar) ||
                                        (IndexAtNode > hiBar))
                                    {        
                                        KoOptL[L] = RebateAtNode;
                                    }
                                }
                            }
                        } /* rebate dim */
                    } /* For k */
                }/* For j */    
            } /* For i */
        }           
        else /* Smoothing required */
        {           
            if (RebateDim == 0)
                RebateAtNode = (double) RebatePtr[0];

            for (i = Bottom1; i <= Top1; i ++)                  
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        offset = tree_data->NodeOffset3[t][i][j][k];
                        KoOptL = (double *)KoOpt + offset;

                        if (RebateDim != 0)
                            RebateAtNodeSlice = (double *)RebatePtr + offset;

                        IndexAtNodeSlice  = (double *)IndexPtr + offset;

                        for (L=Bottom4[i][j][k]; L<=Top4[i][j][k]; L++)
                        {

                            if (RebateDim != 0)
                                RebateAtNode = RebateAtNodeSlice[L];

                            IndexAtNode  = IndexAtNodeSlice[L];

                            IndexStep = Hyb4_GetIndexStep (IndexPtr,
                                                      IndexDim,
                                                      i, j, k, L,
                                                      t,
                                                      tree_data);
                                             
                            if (IoO == 'I')
                            {                                                                                           
                                if (Smooth_Step(&(x),
                                                RebateAtNode,
                                                KoOptL[L],
                                                IndexAtNode,
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step(&(KoOptL[L]),
                                                KoOptL[L],
                                                x,
                                                IndexAtNode,
                                                HighBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;                                   
                                }
                            }    
                            else
                            {    
                                if (Smooth_Step(&(x),
                                                KoOptL[L],
                                                RebateAtNode,
                                                IndexAtNode,
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step(&(KoOptL[L]),
                                                RebateAtNode,
                                                x,
                                                IndexAtNode,
                                                HighBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }  /* if then else */
                        } /* For L */
                    } /* For k */    
                }  /* for j */    
            }  /* for i */    
         }  /* if then else */    
    }  /* if */    
    

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb4_KoOption_t */

