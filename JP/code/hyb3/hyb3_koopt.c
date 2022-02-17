/****************************************************************************/
/*      Calculation of a knock-out option price in the lattice.             */
/****************************************************************************/
/*      KOOPT.c                                                             */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"







/*****  Hyb3_KoOption_t  ***************************************************/
/* 
 *   wrapper around new function with additional parameter                                                                                
 */
int    Hyb3_KoOption_t(
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
           HYB3_TREE_DATA   *tree_data) /* (I) Tree data structure         */
{

    return Hyb3_KoOption_Am_t( KoOpt, IndexDim, IndexPtr, KoFlag, LowBarrier,
           HighBarrier, RebateDim, RebatePtr, IoO, Smoothing, t, DMode, tree_data, FALSE );
}

/*****  Hyb3_KoOption_Am_t  ********************************************************/
/*                                                                           
 *    Calculate the knock-out option price in the lattice using
 *    Arnon's smoothing. It does not include discounting.
 *
 *    The price slice must be 2-D, but the rebate and the index
 *    can be of variable dimension.
 *
 *    Smoothing Flag changes stepsize in smoothing so that American/Continuous
 *    knock out is repriced more in line with analytical prices (+ better Greeks)
 * 
 *    Only switched on in the 3rd dimension for the time being 
 */
int    Hyb3_KoOption_Am_t(
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
           HYB3_TREE_DATA   *tree_data,   /* (I) Tree data structure         */
           int          useAmSmoothing)  /* (I) American Smoothing Flag      */
{





    double  *KoOptL = NULL;       /* Local for slice addressing              */
    double   x;                   /* Smoothed value                          */
    double   IndexStep;           /* Maximum difference between index values */
    double   RebateAtNode = 0;    /* Rebate value at the node (i,j)          */
    double   IndexAtNode;         /* Barrier index value at the node (i,j)   */
    double  *RebateAtNodeSlice = NULL, 
            *IndexAtNodeSlice = NULL;

    int     Top1,  Bottom1;       /* Limits of the tree (1rst dimension)     */
    int    *Top2, *Bottom2;       /* Limits of the tree (2nd dimension)      */
    int   **Top3,**Bottom3;       /* Limits of the tree (3rd dimension)      */
    int     i, j, k;              /* Node indices                            */
    int     offset;
    int     status = FAILURE;     /* Error status = FAILURE initially        */

    /* new smoothing factor for american/continuous exercise */
    double  americanSmoothFactor = (useAmSmoothing ? 1.5 : 1.0);


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



    if (DMode == DISC_1D_NOCUPS)
    {


        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        KoOptL = (double *)KoOpt + offset;

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
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i] = RebateAtNode;               
                    }  
                }
                else /* IoO IS 'O' */
                {
                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
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
                                KoOptL[j] = RebateAtNode;               
                        }  
                    }
                    else /* IoO IS 'O' */
                    {
                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
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

                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                KoOptL = (double *)KoOpt + offset;

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
    if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        if (Smoothing == 'N')
        {
            double lowBar, hiBar;

            if (RebateDim == 0)
                RebateAtNode = (double)RebatePtr[0];
            
            if (IoO == 'I')
            {
                lowBar = LowBarrier*(1.-BARRIER_TOL);
                hiBar = HighBarrier*(1.+BARRIER_TOL);
            }
            else
            {
                lowBar = LowBarrier *(1.+BARRIER_TOL);
                hiBar = HighBarrier*(1.-BARRIER_TOL);
            }

            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
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
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    KoOptL = (double *)KoOpt + offset;

                    if (RebateDim != 0)
                        RebateAtNodeSlice = (double *)RebatePtr + offset;

                    IndexAtNodeSlice  = (double *)IndexPtr + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {

                        if (RebateDim != 0)
                            RebateAtNode = RebateAtNodeSlice[k];

                        IndexAtNode  = IndexAtNodeSlice[k];

                        IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                                  IndexDim,
                                                  i, j, k,
                                                  t,
                                                  tree_data);

                        IndexStep *= americanSmoothFactor;
                                             
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

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_KoOption_t */


/*****  Hyb3_KoOption_t  ***************************************************/
/* 
 *   wrapper around new function with additional parameter                                                                                
 */
int    Hyb3_KiOption_t(
           TSLICE       KiOpt,       /* (I/O) KO option prices at time point */
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
           HYB3_DEV_DATA    *dev_data,  /* (I) Hyb3_Dev data structure       */
           HYB3_TREE_DATA   *tree_data) /* (I) Tree data structure           */
{

    return Hyb3_KiOption_Am_t( 
            KiOpt, Opt, IndexDim, IndexPtr, KiFlag, LowBarrier, HighBarrier, 
            IoO, Smoothing,t,  T, DCurve, DMode, dev_data, tree_data, FALSE) ;


}

/*****  Hyb3_KiOption_t  **********************************************************/
/*                                                                           
 *  Calculate the knock-in  option price in the lattice usingArnon's smoothing.
 *  This is used for path-dependent knock-in structures when the in/out parity 
 *  does not hold. It includes discounting.
 *
 *  Smoothing Flag changes stepsize in smoothing so that American/Continuous
 *  knock out is repriced more in line with analytical prices (+ better Greeks)
 *
 *  Feature only switched on in case of 3rd dimension for time being
 */
int    Hyb3_KiOption_Am_t(
           TSLICE       KiOpt,       /* (I/O) KO option prices at time point */
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
           HYB3_DEV_DATA    *dev_data,  /* (I) Hyb3_Dev data structure       */
           HYB3_TREE_DATA   *tree_data, /* (I) Tree data structure           */
           int         useAmSmoothing)  /* (I) American Smoothing Flag      */
{



    double  *KiOptL;              /* Local slice pointers                    */
    double  *OptL;

    double  x;                    /* Smoothed value                          */
    double  IndexStep;            /* Maximum difference between index values */
    double  IndexAtNode;          /* Barrier index value at the node (i,j)   */
    double  *IndexAtNodeSlice;

    int     Top1,  Bottom1;       /* Limits of the tree (1rst dimension)     */
    int    *Top2, *Bottom2;       /* Limits of the tree (2nd dimension)      */
    int   **Top3,**Bottom3;       /* Limits of the tree (3rd dimension)      */
    int     i, j, k;              /* Node indices                            */
    int     offset;
    int     status = FAILURE;     /* Error status = FAILURE initially        */

    /* new smoothing factor for american/continuous exercise */
    double  americanSmoothFactor = (useAmSmoothing ? 1.5 : 1.0);

    
    Top1    = tree_data->Top1[t];                   
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];                   
    Bottom1 = tree_data->Bottom1[t];                    
    Bottom2 = tree_data->Bottom2[t];    
    Bottom3 = tree_data->Bottom3[t];                      


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


    if (DMode == DISC_1D_NOCUPS)
    {

        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        KiOptL = (double *)KiOpt + offset;
        OptL   = (double *)Opt   + offset;
       

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {                   
                IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                              IndexPtr,
                                              i,0,0,t,tree_data);

                if (IoO == 'I')
                {
                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];               
                    }  
                }
                else /* IoO IS 'O' */
                {
                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];
                    }        
                }
            } /* For i */
        }                
        else
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
                                
                /* Knock-out inside the two barriers */
                if (IoO == 'I') 
                {                              
                    /* 'up' value = std option    */
                    /* 'down' value = live option */
                    if (Smooth_Step(&(x),
                                    OptL[i],   
                                    KiOptL[i], 
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }                                                 
                                
                    /* 'up' value = live option      */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step(&(KiOptL[i]),
                                    KiOptL[i],   
                                    x,      
                                    IndexAtNode,
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;                                   
                    }  
                }    
                else            /* Knock-out outside the two barriers */
                {    
                    /* 'up' value = live option  */
                    /* 'down' value = std option */
                    if (Smooth_Step(&(x),
                                    KiOptL[i],   
                                    OptL[i],     
                                    IndexAtNode,
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }                                             
                                
                    /* 'up' value = std option       */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step(&(KiOptL[i]),
                                    OptL[i],   
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
    if (DMode == DISC_2D_CUPS   || 
        DMode == DISC_2D_NOCUPS ||
        DMode == DISC_2D_1IR2F_NOCUPS)
    {
        

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {     
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                KiOptL = (double *)KiOpt + offset;
                OptL   = (double *)Opt   + offset;
                  
              
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {    
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  IndexPtr,
                                                  i,j,0,t,tree_data);

                    if (IoO == 'I')
                    {
                        if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                           (IndexAtNode < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];               
                        }  
                    }
                    else
                    {
                        if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) ||
                            (IndexAtNode > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];
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
                KiOptL = (double *)KiOpt + offset;
                OptL   = (double *)Opt   + offset;

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
                        if (Smooth_Step(&(x),
                                        OptL[j],
                                        KiOptL[j],
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Smooth_Step(&(KiOptL[j]),
                                        KiOptL[j],
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
                                        KiOptL[j],
                                        OptL[j],
                                        IndexAtNode,
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Smooth_Step(&(KiOptL[j]),
                                        OptL[j],
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
    if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {

        if (Smoothing == 'N')
        {
            double lowBar, hiBar;

            if (IoO == 'I')
            {
                lowBar = LowBarrier*(1.- BARRIER_TOL);
                hiBar = HighBarrier*(1.+BARRIER_TOL);
            }
            else
            {
                lowBar = LowBarrier*(1.+BARRIER_TOL);
                hiBar = HighBarrier*(1.-BARRIER_TOL);
            }


            for (i = Bottom1; i <= Top1; i ++)
            {                   
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    KiOptL = (double *)KiOpt + offset;
                    OptL   = (double *)Opt   + offset;
                    IndexAtNodeSlice  = (double *)IndexPtr + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {    
                        IndexAtNode  = IndexAtNodeSlice[k];

                        if (IoO == 'I')
                        {
                            if((IndexAtNode > lowBar) && 
                               (IndexAtNode < hiBar))
                            {        
                                KiOptL[k] = OptL[k];               
                            }  
                        }
                        else
                        {
                            if ((IndexAtNode < lowBar) ||
                                (IndexAtNode > hiBar))
                            {        
                                KiOptL[k] = OptL[k];
                            }
                            
                        }
                    } /* For k */
                }/* For j */    
            } /* For i */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)                  
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    KiOptL = (double *)KiOpt + offset;
                    OptL   = (double *)Opt   + offset;
                    IndexAtNodeSlice  = (double *)IndexPtr + offset;
                 
                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        IndexAtNode  = IndexAtNodeSlice[k];

                        IndexStep = Hyb3_GetIndexStep (IndexPtr,
                                                  IndexDim,
                                                  i, j, k,
                                                  t,
                                                  tree_data);

                                    
                        IndexStep *= americanSmoothFactor;

                        if (IoO == 'I')
                        {                                                                                           
                            if (Smooth_Step(&(x),
                                            OptL[k],
                                            KiOptL[k],
                                            IndexAtNode,
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step(&(KiOptL[k]),
                                            KiOptL[k],
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
                                            KiOptL[k],
                                            OptL[k],
                                            IndexAtNode,
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step(&(KiOptL[k]),
                                            OptL[k],
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
        

    status = SUCCESS;

    RETURN:

    return (status);


}  /* Hyb3_KiOption_t */



/*****  Hyb3_Trigger_t  *******************************************************/
/*                                                            
 *       Option price with exercise triggered by an index. It uses Arnon's
 *       smoothing algorithm.
 *
 *       WARNING: This function uses the CoP (call or put) flag to calc-
 *       ulate  the price of 'long the call' or 'short the put'. This is
 *       intended to extend  the functionality  of Hyb3_Trigger_t() for usage
 *       with the 'top' options (i.e. callable or puttable).  The caller
 *       of this function  must therefore be  aware that, when pricing a
 *       simple put, the value will always be <0.
 *
 *
 */
int     Hyb3_Trigger_t  
           (TSLICE      TriggerOpt,  /* (I/O) Trigger option prices     */
            TSLICE      Under,       /* (I) Underlying swap prices      */
            int         IndexDim,     /* (I) Dimension of index slice    */
            TSLICE      Index,       /* (I) Trigger values              */
            double      Notional,     /* (I) Notional of the option      */
            double      Strike,       /* (I) Strike of the option        */
            long        TriggerFlag,  /* (I) Trigger flag                */
            int         CoP,          /* (I) =1 for a call, -1 for a put */
            double      LowBarrier,   /* (I) Lower barrier               */
            double      HighBarrier,  /* (I) Higher barrier              */
            char        IoO,          /* (I) Trigger 'I'n or 'O'ut side  */
            char        Smoothing,    /* (I) Smoothing ('Y' or 'N')      */
            int         t,            /* (I) Current time period         */
            int         T,            /* (I) Total number of period      */
            int         DCurve,       /* (I) Discount curve              */
            int         DMode,        /* (I) Also used for dim of slices */
            HYB3_DEV_DATA    *dev_data,    /* (I) Hyb3_Dev data structure          */
            HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure         */
{


    double  *TriggerL;                /* Local slice pointers            */
    double  *UnderL;

    double
            CoPd,
            x,                        /* Intermediate smoothed value     */
            Intrinsic,                /* Intrinsic value of option       */
            IndexAtNode,              /* Current value of index          */
            IndexStep;                /* Max diff between index values   */

    double  *IndexAtNodeSlice;

    int          
            Top1, Bottom1,            /* Tree limits (1rst dim)          */
            *Top2, *Bottom2,          /* Tree limits (2nd dim)           */
            **Top3, **Bottom3,        /* Tree limits (3rd dim)           */
            i, j, k,                  /* Node indices                    */
            offset,
            status = FAILURE;         /* Error status                    */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    
    if (Hyb3_Dev (TriggerOpt,
             t,
             T,
             DCurve,
             DMode,
             dev_data,
             tree_data) == FAILURE)
    {
        goto RETURN;
        
    }  
    
    
    if (!TriggerFlag)                  /* Nothing to do */
    {
        status = SUCCESS;

        return (status);

    }  /* if */


    /* Check that index dimension is in line with DEV mode */
    if (IndexDim!=1 && IndexDim!=2 && IndexDim!=3)
    {
        DR_Error("Invalid index dimension (must be 1, 2 or 3).(Hyb3_Trigger_t)\n");
        goto RETURN;
    }

    if (DMode == DISC_1D_NOCUPS)
    {
        if (IndexDim != 1)
        {
           DR_Error("DEV Mode (1D) incompatible with index dim!(Hyb3_Trigger_t)\n");
           goto RETURN;
        }
    }
    else if (DMode == DISC_2D_CUPS   || 
             DMode == DISC_2D_NOCUPS ||
             DMode == DISC_2D_1IR2F_NOCUPS)
    {
        if ((IndexDim != 1) && (IndexDim != 2))
        {
           DR_Error("DEV Mode (2D) incompatible with index dim!(Hyb3_Trigger_t)\n");
           goto RETURN;
        }
    }
    

    if (CoP == 1)
    {
        CoPd = 1.000;
    }
    else
    {
        CoPd = -1.000;
    }

        

    if (DMode == DISC_1D_NOCUPS)
    {
        
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);
        TriggerL = (double *)TriggerOpt + offset;
        UnderL   = (double *)Under   + offset;

        if (Smoothing == 'N')
        {
            if (IoO == 'I')                                                
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  Index,
                                                  i,0,0,t,tree_data);

                    if ((IndexAtNode > LowBarrier  * (1. - BARRIER_TOL)) && 
                        (IndexAtNode < HighBarrier * (1. + BARRIER_TOL))  )
                    {        
                        TriggerL[i] = CoPd *  MAX (CoPd*TriggerL[i],
                                         Notional*CoPd*(UnderL[i] - Strike));              
                    }  /* if */                 
                }  /* for i */    
            }           
            else                                                           
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  Index,
                                                  i,0,0,t,tree_data);

                    if ((IndexAtNode < LowBarrier  * (1. + BARRIER_TOL)) || 
                        (IndexAtNode > HighBarrier * (1. - BARRIER_TOL))  )
                    {        
                        TriggerL[i] = CoPd* MAX (CoPd*TriggerL[i], 
                                       Notional * CoPd * (UnderL[i] - Strike));                      
                    }  /* if */                 
                }  /* for i */    
            }  /* if then else */
        }           
        else /* Smoothing ON */
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {   
                Intrinsic = Notional * CoPd * (UnderL[i] - Strike);
                TriggerL[i] *= CoPd;
                
                if(TriggerL[i] >= Intrinsic) /* No value to am exer */
                {
                    TriggerL[i] *= CoPd;
                }
                else
                {
                    /* There is value to american exercise              */
                    
                    /* Calculate max difference between values of index */
                    /* at adjacent nodes (check for boundaries).        */
                    /* Avoid division by 0                              */
                    IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                  Index,
                                                  i,0,0,t,tree_data);
                
                    IndexStep = Hyb3_GetIndexStep (Index,
                                              IndexDim,
                                              i, 0, 0,
                                              t,
                                              tree_data);

                    if (IoO == 'I')     
                    {                                               
                        if (Smooth_Step (&x,
                                         Intrinsic,     /* up=intrinsic   */
                                         TriggerL[i],/* down=option    */
                                         IndexAtNode,
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                            
                        }  /* if */                             
                        
                        if (Smooth_Step (&(TriggerL[i]),
                                         TriggerL[i], /* up=option     */
                                         x,              /* down=smoothed */
                                         IndexAtNode,
                                         HighBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                            
                        }  /* if */
                    }    
                    else                                              
                    {    
                        if (Smooth_Step (&x,
                                         TriggerL[i],  /* up=option     */
                                         Intrinsic,       /* down=intrinsic*/
                                         IndexAtNode,
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                            
                        }  /* if */                                       
                        
                        if (Smooth_Step (&(TriggerL[i]),
                                         Intrinsic,       /* up=intrinsic */
                                         x,               /* down=smoothed*/
                                         IndexAtNode,
                                         HighBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                            
                        }  /* if */
                    }  /* if then else */    
                
                    TriggerL[i] *= CoPd;
                
                } /* if there is value to amer exer */   
            }  /* for i over nodes */    
        }  /* if then else for Smoothing */    
    }

    else if (DMode == DISC_2D_CUPS   || 
             DMode == DISC_2D_NOCUPS ||
             DMode == DISC_2D_1IR2F_NOCUPS)
    {    

        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++) 
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    TriggerL = (double *)TriggerOpt + offset;
                    UnderL   = (double *)Under   + offset;
                  
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      Index,
                                                      i,j,0,t,tree_data);

                        if ((IndexAtNode > LowBarrier  * (1.- BARRIER_TOL)) && 
                            (IndexAtNode < HighBarrier * (1.+ BARRIER_TOL)))
                        {        
                            TriggerL[j] = CoPd*MAX(CoPd*TriggerL[j],
                                         Notional* CoPd*(UnderL[j]-Strike));
                            
                        }  /* if */                 
                    }  /* for j */
                }    
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++) 
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                    TriggerL = (double *)TriggerOpt + offset;
                    UnderL   = (double *)Under   + offset;
         
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      Index,
                                                      i,j,0,t,tree_data);

                        if ((IndexAtNode < LowBarrier  * (1.+BARRIER_TOL)) || 
                            (IndexAtNode > HighBarrier * (1.-BARRIER_TOL)))
                        {        
                            TriggerL[j] = CoPd*MAX(CoPd*TriggerL[j],
                                    Notional * CoPd * (UnderL[j] - Strike));
                            
                        }  /* if */                 
                    }  /* for j */    
                }
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {

                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
                TriggerL = (double *)TriggerOpt + offset;
                UnderL   = (double *)Under   + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Intrinsic = Notional * CoPd * (UnderL[j] - Strike);
                    TriggerL[j] *= CoPd;
                    
                    if(TriggerL[j] >= Intrinsic) /* No value */
                    {
                        TriggerL[j] *= CoPd;
                    }
                    else
                    {
                        /* There is value to american exericse        */
                        IndexAtNode  = Hyb3_GetValueAtNode(IndexDim,
                                                      Index,
                                                      i,j,0,t,tree_data);

                        IndexStep = Hyb3_GetIndexStep (Index,
                                                  IndexDim,
                                                  i, j, 0,
                                                  t,
                                                  tree_data);
                        
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (&x,
                                             Intrinsic,
                                             TriggerL[j],
                                             IndexAtNode,
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                                
                            }  /* if */                        
                            
                            if (Smooth_Step (&(TriggerL[j]),
                                             TriggerL[j],
                                             x,
                                             IndexAtNode,
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                                
                            }  /* if */
                        }    
                        else
                        {    
                            if (Smooth_Step (&x,
                                             TriggerL[j],
                                             Intrinsic,
                                             IndexAtNode,
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                                
                            }  /* if */                                                
                            
                            if (Smooth_Step (&(TriggerL[j]),
                                             Intrinsic,
                                             x,
                                             IndexAtNode,
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                                
                            }  /* if */
                        }  /* if then else */

                        TriggerL[j] *= CoPd;

                    } /* if amer exer has value */
                }  /* for j */    
            }  /* for i */    
        }  /* if then else smoothing ON */
    }  /* if then else */
    else if (DMode == DISC_3D_CUPS || DMode == DISC_3D_2IR2F1D_CUPS)
    {
        

        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++)     
                {     
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                        TriggerL = (double *)TriggerOpt + offset;
                        UnderL   = (double *)Under   + offset;
                        IndexAtNodeSlice  = (double *)Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            IndexAtNode  = IndexAtNodeSlice[k];

                            if ((IndexAtNode > LowBarrier *(1.-BARRIER_TOL)) && 
                                (IndexAtNode < HighBarrier*(1.+ BARRIER_TOL)))
                            {        
                                TriggerL[k] = CoPd * 
                                  MAX (CoPd*TriggerL[k],
                                  Notional * CoPd * (UnderL[k] - Strike));
                                       
                            }  /* if */                 
                        }  /* for k */
                    }
                }    
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                        TriggerL = (double *)TriggerOpt + offset;
                        UnderL   = (double *)Under   + offset;
                        IndexAtNodeSlice  = (double *)Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            IndexAtNode  = IndexAtNodeSlice[k];

                            if ((IndexAtNode < LowBarrier *(1.+BARRIER_TOL)) || 
                                (IndexAtNode > HighBarrier*(1.-BARRIER_TOL)))
                            {        
                                TriggerL[k] = CoPd *
                                    MAX (CoPd*TriggerL[k], 
                                    Notional * CoPd * (UnderL[k]-Strike));
                                       
                            }  /* if */                 
                        }  /* for k */
                    }
                }    
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
                    TriggerL = (double *)TriggerOpt + offset;
                    UnderL   = (double *)Under   + offset;
                    IndexAtNodeSlice  = (double *)Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)          
                    {
                        Intrinsic = Notional * CoPd * (UnderL[k] - Strike);
                        TriggerL[k] *= CoPd;
                
                        if(TriggerL[k] >= Intrinsic) /* No value */
                        {
                            TriggerL[k] *= CoPd;
                        }
                        else
                        {
                            /* There is value to american exericse        */
                        
                            /* Calculate max diff between values of index */
                            /* at adjacent nodes (check for boundaries).  */
                            /* Avoid division by 0                        */
                    
                            IndexAtNode  = IndexAtNodeSlice[k];

                            IndexStep = Hyb3_GetIndexStep (Index,
                                                      IndexDim,
                                                      i, j, k,
                                                      t,
                                                      tree_data);
                            
                            if (IoO == 'I')
                            {                                               
                                if (Smooth_Step (&x,
                                                 Intrinsic,
                                                 TriggerL[k],
                                                 IndexAtNode,
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                    
                                }  /* if */                        
                                
                                if (Smooth_Step (&(TriggerL[k]),
                                                 TriggerL[k],
                                                 x,
                                                 IndexAtNode,
                                                 HighBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                    
                                }  /* if */
                            }    
                            else
                            {    
                                if (Smooth_Step (&x,
                                                 TriggerL[k],
                                                 Intrinsic,
                                                 IndexAtNode,
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                    
                                }  /* if */                                                
                                
                                if (Smooth_Step (&(TriggerL[k]),
                                                 Intrinsic,
                                                 x,
                                                 IndexAtNode,
                                                 HighBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                    
                                }  /* if */
                            }  /* if then else */
                        
                            TriggerL[k] *= CoPd;
                        
                        } /* if american exercise has value */
                    }  /* for k */    
                }  /* for j */    
            }  /* for i */    
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Trigger_t */




/*******************************************************************************
** These should be in a separate file. Since new .cpp modules are combustible,
** in the interest of fire safety they are kept here.
*******************************************************************************/
int Hyb3_KOIndicator(
         double*                Indicator,    /**< (O) Slice to be modified   */
         double const*          Index,        /**< (I) observation index      */
         double                 low,          /**< (I) Low barrier            */
         double                 high,         /**< (I) High barrier           */
         char                   io,           /**< (I) 'I'nside or 'O'utside  */
         int                    t,            /**< (I) Current time point     */
         int                    tree_dim,     /**< (I) dimension of the slice */
         HYB3_TREE_DATA const*  tree_data)    /**< (I) Tree data              */
{


    double       *IndicatorL;
    double const *IndexL;

    int     Top1, Bottom1;                     /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                   /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;                 /* Tree limits (3rd dim)   */

    int     i, j, k;                           /* Node indices            */
    int     offset;                            /* Node offset             */

    double  inside = ((io == 'I') ? 1.0 : -1.0);


    Top1    = tree_data->Top1[t];   
    Bottom1 = tree_data->Bottom1[t];
    Top2    = tree_data->Top2[t];    
    Bottom2 = tree_data->Bottom2[t];
    Top3    = tree_data->Top3[t];   
    Bottom3 = tree_data->Bottom3[t];
    
    switch (tree_dim)
    {
        case 1:
        {
            offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

            IndicatorL = Indicator + offset;
            IndexL     = Index + offset;

            for (i=Bottom1; i<=Top1; i++)
            {
                IndicatorL[i] = ( ( IndexL[i] < 0.5*(low+high) ) 
                            ? inside*(IndexL[i] - low) : inside*(high-IndexL[i]) ) ;
            }

            break;
        }        
        case 2:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                IndicatorL = Indicator + offset;
                IndexL     = Index + offset;

                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    IndicatorL[j] = ( ( IndexL[j] < 0.5*(low+high) ) 
                            ? inside*(IndexL[j] - low) : inside*(high-IndexL[j]) ) ;
                }
            }

            break;
        }        
        case 3:
        {
            for (i=Bottom1; i<=Top1; i++)
            {
                for (j=Bottom2[i]; j<=Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                    IndicatorL = Indicator + offset;
                    IndexL     = Index + offset;

                    for (k=Bottom3[i][j]; k<=Top3[i][j]; k++)
                    {
                        IndicatorL[k] = ( ( IndexL[k] < 0.5*(low+high) ) 
                            ? inside*(IndexL[k] - low) : inside*(high-IndexL[k]) ) ;
                    }
                }
            }

            break;
        }        
        default:
        {
            DR_Error("Hyb3_KoIndicator: invalid slice dimension "
                        "(must be 1, 2 or 3)!");
            return(FAILURE);
        }
    } /* End of switch() */
    
    return (SUCCESS);

}
