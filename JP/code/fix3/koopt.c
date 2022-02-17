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
#include "fix123head.h"

#define  EPS_TOL       1e-6           /* Tol for smoothing epsilon      */


/*****  Fix3_KoOption_t  *********************************************************/
/*                                                                           
*       Knock-out an option price without discounting it. It uses Arnon's
*       smoothing algorithm.
*/
int     Fix3_KoOption_t (double      *KoOpt,      /* (I/O) Knock-out option       */
                    double      *Index,      /* (I) Index prices             */
                    long        KoFlag,      /* (I) Knock-out flag           */
                    double      LowBarrier,  /* (I) Lower barrier            */
                    double      HighBarrier, /* (I) Higher barrier           */
                    double      Rebate,      /* (I) Rebate                   */
                    char        IoO,         /* (I) Ko 'I'nside or 'O'utside */
                    char        Smoothing,   /* (I) Smoothing ('Y' or 'N')   */
                    int         t,           /* (I) Current time point       */
                    FIX3_TREE_DATA   *tree_data)  /* (I) Structure of tree data   */
{

    double  *KoOptL;            /* Local slice pointers */
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  IndexStep;          /* Diff between index values   */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (!KoFlag) /* Nothing to do */
    {
        return (SUCCESS);
    }

    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        KoOptL = KoOpt + offset;
        IndexL = Index + offset;

        if (Smoothing == 'N')
        {
            if (IoO == 'I') /* Knock-out inside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] > LowBarrier * (1. - BARRIER_TOL)) 
                     && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i] = Rebate;
                    }
                }  /* for i */  
            }           
            else            /* Knock-out outside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] < LowBarrier * (1. + BARRIER_TOL)) 
                     || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KoOptL[i] = Rebate;
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                IndexStep = Fix3_GetIndexStep (Index,
                                          1,                                                  
                                          i, 0, 0,
                                          t,
                                          tree_data);
                                                  
                if (IoO == 'I')
                {                                               
                    /* 'up' value = rebate        */
                    /* 'down' value = live option */
                    if (Smooth_Step (	&x,
                                        Rebate,                                 
                                        KoOptL[i],                              
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = live option      */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (	&(KoOptL[i]),
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
                    /* 'up' value = live option */
                    /* 'down' value = rebate    */
                    if (Smooth_Step (	&x,
                                        KoOptL[i],                              
                                        Rebate,                                 
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = rebate           */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (	&(KoOptL[i]),
                                        Rebate,                                 
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
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                         && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KoOptL[j] = Rebate;
                        }
                    }  /* for j */  
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] < LowBarrier * (1. + BARRIER_TOL)) 
                         || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KoOptL[j] = Rebate;
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                KoOptL = KoOpt + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexStep = Fix3_GetIndexStep (Index,
                                              2,                                                  
                                              i, j, 0,
                                              t,
                                              tree_data);
                                                  
                    if (IoO == 'I')
                    {                                               
                        if (Smooth_Step (	&x,
                                            Rebate,
                                            KoOptL[j],
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        if (Smooth_Step (	&(KoOptL[j]),
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
                        if (Smooth_Step (	&x,
                                            KoOptL[j],
                                            Rebate,
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                
                        if (Smooth_Step (	&(KoOptL[j]),
                                            Rebate,
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
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL)) 
                             && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                KoOptL[k] = Rebate;
                            }
                        }  /* for k */  
                    }  /* for j */  
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] < LowBarrier * (1. + BARRIER_TOL)) 
                             || (IndexL[k] > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                KoOptL[k] = Rebate;
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  3,                                                  
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                                                  
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (	&x,
                                                Rebate,
                                                KoOptL[k],
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                    
                            if (Smooth_Step (	&(KoOptL[k]),
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
                            if (Smooth_Step (	&x,
                                                KoOptL[k],
                                                Rebate,
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                
                            if (Smooth_Step (	&(KoOptL[k]),
                                                Rebate,
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

}  /* Fix3_KoOption_t */




/*****  Fix3_KoOptionEps_t  *********************************************************/
/*                                                                           
*  Knock-out an option price without discounting it. Use Epsilon rule 
* for smoothing: same convetion as in bivariate for pos, neg epsilons
*/
int     Fix3_KoOptionEps_t (
                    double      *KoOpt,      /* (I/O) Knock-out option       */
                    double      *Index,      /* (I) Index prices             */
                    long        KoFlag,      /* (I) Knock-out flag           */
                    double      LowBarrier,  /* (I) Lower barrier            */
                    double      HighBarrier, /* (I) Higher barrier           */
                    double      EpsLo,       /* (I) Epsilon Lo               */
                    double      EpsHi,       /* (I) Epsilon Hi               */
                    double      Rebate,      /* (I) Rebate                   */
                    char        IoO,         /* (I) Ko 'I'nside or 'O'utside */
                    char        Smoothing,   /* (I) Smoothing ('Y' or 'N')   */
                    int         t,           /* (I) Current time point       */
                    FIX3_TREE_DATA   *tree_data) /* (I) Structure of tree data   */
{

    double  *KoOptL;            /* Local slice pointers */
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  IndexStepLo;
    double  IndexStepHi;
    double  IndexStep;          /* Diff between index values   */
    int     IsZeroEpsLo = FALSE;
    int     IsZeroEpsHi = FALSE;
    

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    IndexStepLo = EpsLo / 2.0;
    IndexStepHi = EpsHi / 2.0;  

    /* if Epsilon smaller than tolerance then use old smoothing alg */
    if (fabs(EpsLo) < EPS_TOL)
    {
        IsZeroEpsLo = TRUE;
    }
    if (fabs(EpsHi) < EPS_TOL)
    {
        IsZeroEpsHi = TRUE;
    }


    if (!KoFlag) /* Nothing to do */
    {
        return (SUCCESS);
    }

    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0,  0, t, tree_data);

        KoOptL = KoOpt + offset;
        IndexL = Index + offset;


        if (Smoothing == 'N' ) 
        {
            if (IoO == 'I') /* Knock-out inside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] > LowBarrier * (1. - BARRIER_TOL)) 
                     && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KoOptL[i] = Rebate;
                    }
                }  /* for i */  
            }           
            else            /* Knock-out outside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] < LowBarrier * (1. + BARRIER_TOL)) 
                     || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KoOptL[i] = Rebate;
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {    
                IndexStep = Fix3_GetIndexStep (Index,
                                          1,                                                  
                                          i, 0, 0,
                                          t,
                                          tree_data);
                                                  
                if (IoO == 'I')
                {                                               
                    /* 'up' value = rebate        */
                    /* 'down' value = live option */
                    /* is zero epsilon use old smoothing algorithm */
                    if (IsZeroEpsLo)
                    {
                        if (Smooth_Step (&x,
                                        Rebate,                                 
                                        KoOptL[i],                              
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }
                    else
                    {
                        if (Smooth_Step   ( &x,
                                            Rebate,                                 
                                            KoOptL[i],                              
                                            IndexL[i],
                                            LowBarrier - IndexStepLo,
                                            fabs(IndexStepLo)) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }
                    
                    
                    /* 'up' value = live option      */
                    /* 'down' value = smoothed value */
                    if (IsZeroEpsHi)
                    {
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
                        if (Smooth_Step( &(KoOptL[i]),
                                        KoOptL[i],                              
                                        x,                                      
                                        IndexL[i],
                                        HighBarrier + IndexStepHi,
                                        fabs(IndexStepHi)) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }
                    
                }
                else           
                {
                    /* 'up' value = live option */
                    /* 'down' value = rebate    */  
                    if (IsZeroEpsLo)
                    {
                         if (Smooth_Step (&x,
                                        KoOptL[i],                              
                                        Rebate,                                 
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                        {
                            goto RETURN;
                         }            
                    }
                    else
                    {
                        if (Smooth_Step   ( &x,
                                        KoOptL[i],                              
                                        Rebate,                                 
                                        IndexL[i],
                                        LowBarrier - IndexStepLo,
                                        fabs(IndexStepLo)) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }
                    /* 'up' value = rebate           */
                    /* 'down' value = smoothed value */
                    if (IsZeroEpsHi)
                    {
                        if (Smooth_Step (&(KoOptL[i]),
                                        Rebate,                                 
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
                        if (Smooth_Step   ( &(KoOptL[i]),
                                            Rebate,                                 
                                            x,                                      
                                            IndexL[i],
                                            HighBarrier + IndexStepHi,
                                            fabs(IndexStepHi)) == FAILURE)
                        {
                            goto RETURN;
                        }
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
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                         && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KoOptL[j] = Rebate;
                        }
                    }  /* for j */  
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] < LowBarrier * (1. + BARRIER_TOL)) 
                         || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KoOptL[j] = Rebate;
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                KoOptL = KoOpt + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    
                    IndexStep = Fix3_GetIndexStep (Index,
                                              2,                                                  
                                              i, j, 0,
                                              t,
                                              tree_data);
                    if (IoO == 'I')
                    {     
                        if (IsZeroEpsLo)
                        {
                            if (Smooth_Step (&x,
                                            Rebate,
                                            KoOptL[j],
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        else
                        {
                            if (Smooth_Step    (&x,
                                            Rebate,
                                            KoOptL[j],
                                            IndexL[j],
                                            LowBarrier - IndexStepLo,
                                            fabs(IndexStepLo)) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                          
                        if (IsZeroEpsHi)
                        {
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
                            if (Smooth_Step    (&(KoOptL[j]),
                                                KoOptL[j],
                                                x,
                                                IndexL[j],
                                                HighBarrier + IndexStepHi,
                                                fabs(IndexStepHi)) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                    }
                    else
                    {
                        if (IsZeroEpsLo)
                        {
                            if (Smooth_Step (&x,
                                            KoOptL[j],
                                            Rebate,
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        else
                        {
                            if (Smooth_Step    (&x,
                                                KoOptL[j],
                                                Rebate,
                                                IndexL[j],
                                                LowBarrier - IndexStepLo,
                                                fabs(IndexStepLo)) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }
                        if (IsZeroEpsHi)
                        {
                            if (Smooth_Step (&(KoOptL[j]),
                                            Rebate,
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
                            if (Smooth_Step    (&(KoOptL[j]),
                                                Rebate,
                                                x,
                                                IndexL[j],
                                                HighBarrier + IndexStepHi,
                                                fabs(IndexStepHi)) == FAILURE)
                            {
                                goto RETURN;
                            }
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
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL)) 
                             && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                KoOptL[k] = Rebate;
                            }
                        }  /* for k */  
                    }  /* for j */  
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KoOptL = KoOpt + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] < LowBarrier * (1. + BARRIER_TOL)) 
                             || (IndexL[k] > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                KoOptL[k] = Rebate;
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    KoOptL = KoOpt + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {    
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  3,                                                  
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                        if (IoO == 'I')
                        {    
                            if (IsZeroEpsLo)
                            {
                                 if (Smooth_Step (&x,
                                                Rebate,
                                                KoOptL[k],
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }    
                            }
                            else
                            {
                                if (Smooth_Step (&x,
                                                Rebate,
                                                KoOptL[k],
                                                IndexL[k],
                                                LowBarrier - IndexStepLo,
                                                fabs(IndexStepLo)) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }
                            if(IsZeroEpsHi)
                            {
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
                                if (Smooth_Step (&(KoOptL[k]),
                                                    KoOptL[k],
                                                    x,
                                                    IndexL[k],
                                                    HighBarrier + IndexStepHi,
                                                    fabs(IndexStepHi)) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }
                            
                        }
                        else
                        {
                            if(IsZeroEpsLo)
                            {
                                 if (Smooth_Step (&x,
                                                KoOptL[k],
                                                Rebate,
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }
                            else
                            {
                                if (Smooth_Step (&x,
                                                KoOptL[k],
                                                Rebate,
                                                IndexL[k],
                                                LowBarrier - IndexStepLo,
                                                fabs(IndexStepLo)) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }

                            if(IsZeroEpsHi)
                            {
                                if (Smooth_Step (&(KoOptL[k]),
                                                Rebate,
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
                                if (Smooth_Step  (&(KoOptL[k]),
                                                Rebate,
                                                x,
                                                IndexL[k],
                                                HighBarrier + IndexStepHi,
                                                fabs(IndexStepHi)) == FAILURE)
                                {
                                    goto RETURN;
                                }
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

}  /* Fix3_KoOptionEpsilon_t */


/*****  Fix3_KiOption_t  *********************************************************/
/*                                                                           
*       Calculate the knock-in option price in the lattice using
*       Arnon's smoothing.
*       This is used for path-dependent knock-in structures when the in / out
*       parity does not hold. It includes discounting.
*/
int     Fix3_KiOption_t (double      *KiOpt,      /* (I/O) Knock-in option prices */
                    double      *Opt,        /* (I) Standard option prices   */
                    double      *Index,      /* (I) Index prices             */
                    long        KiFlag,      /* (I) Knock-in flag            */
                    double      LowBarrier,  /* (I) Lower barrier            */
                    double      HighBarrier, /* (I) Higher barrier           */
                    char        IoO,         /* (I) Ki 'I'nside or 'O'utside */
                    char        Smoothing,   /* (I) Smoothing ('Y' or 'N')   */
                    int         t,           /* (I) Current time point       */
                    int         T,           /* (I) Last time point          */
                    int         DCurve,      /* (I) Discount curve           */
                    FIX3_DEV_DATA    *dev_data,   /* (I) Fix3_Dev data structure       */
                    FIX3_TREE_DATA   *tree_data)  /* (I) Tree data structure      */
{

    double  *KiOptL;            /* Local slice pointers */
    double  *OptL;
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  IndexStep;          /* Diff between index values   */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Fix3_Dev (   KiOpt,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
    }


    if (!KiFlag)    /* Nothing to do */
    {
        return (SUCCESS);
    }

    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        KiOptL = KiOpt + offset;
        OptL   = Opt   + offset;
        IndexL = Index + offset;

        if (Smoothing == 'N')
        {
            if (IoO == 'I') /* Knock-in inside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] > LowBarrier * (1. - BARRIER_TOL)) 
                     && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];
                    }
                }  /* for i */  
            }           
            else            /* Knock-in outside the two levels */
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] < LowBarrier * (1. + BARRIER_TOL)) 
                     || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        KiOptL[i] = OptL[i];
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                IndexStep = Fix3_GetIndexStep (Index,
                                          1,                                                  
                                          i, 0, 0,
                                          t,
                                          tree_data);
                                                  
                if (IoO == 'I')
                {              
                    /* 'up' value = standard option */
                    /* 'down' value = live option   */
                    if (Smooth_Step (	&x,
                                        OptL[i],                                
                                        KiOptL[i],                              
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = live option      */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (	&(KiOptL[i]),
                                        KiOptL[i],                              
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
                    /* 'up' value = live option */
                    /* 'down' value = standard option */
                    if (Smooth_Step (	&x,
                                        KiOptL[i],                              
                                        OptL[i],                                
                                        IndexL[i],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = standard option  */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (	&(KiOptL[i]),
                                        OptL[i],                                
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
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KiOptL = KiOpt + offset;
                    OptL   = Opt   + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                         && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];
                        }
                    }  /* for j */
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    KiOptL = KiOpt + offset;
                    OptL   = Opt   + offset;
                    IndexL = Index + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] < LowBarrier * (1. + BARRIER_TOL)) 
                         || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            KiOptL[j] = OptL[j];
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                KiOptL = KiOpt + offset;
                OptL   = Opt   + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexStep = Fix3_GetIndexStep (Index,
                                              2,                                                  
                                              i, j, 0,
                                              t,
                                              tree_data);
                                                  
                    if (IoO == 'I')
                    {                                               
                        if (Smooth_Step (	&x,
                                            OptL[j],
                                            KiOptL[j],
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        if (Smooth_Step (	&(KiOptL[j]),
                                            KiOptL[j],
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
                        if (Smooth_Step (	&x,
                                            KiOptL[j],
                                            OptL[j],
                                            IndexL[j],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        if (Smooth_Step (	&(KiOptL[j]),
                                            OptL[j],
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
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KiOptL = KiOpt + offset;
                        OptL   = Opt   + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL)) 
                             && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                KiOptL[k] = OptL[k];
                            }
                        }  /* for k */  
                    }  /* for j */  
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        KiOptL = KiOpt + offset;
                        OptL   = Opt   + offset;
                        IndexL = Index + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] < LowBarrier * (1. + BARRIER_TOL)) 
                             || (IndexL[k] > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                KiOptL[k] = OptL[k];
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    KiOptL = KiOpt + offset;
                    OptL   = Opt   + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  3,                                                  
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                                                  
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (	&x,
                                                OptL[k],
                                                KiOptL[k],
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                        
                            if (Smooth_Step (	&(KiOptL[k]),
                                                KiOptL[k],
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
                            if (Smooth_Step (	&x,
                                                KiOptL[k],
                                                OptL[k],
                                                IndexL[k],
                                                LowBarrier,
                                                IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                                    
                            if (Smooth_Step (	&(KiOptL[k]),
                                                OptL[k],
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

}  /* Fix3_KiOption_t */



/*****  Fix3_Trigger_t  *******************************************************/
/*                                                            
*       Option price with exercise triggered by an index. It uses Arnon's
*       smoothing algorithm to american exercise premium only.
*/
int     Fix3_Trigger_t  
           (double      *Trigger,     /* (I/O) Trigger option prices     */
            double      *Under,       /* (I) Underlying swap prices      */
            double      *Index,       /* (I) Trigger values              */
            double      Notional,     /* (I) Notional of the option      */
            double      Strike,       /* (I) Strike of the option        */
            long        TriggerFlag,  /* (I) Trigger flag                */
            int         CoP,          /* (I) =1 for a call, -1 for a put */
            double      LowBarrier,   /* (I) Lower barrier               */
            double      HighBarrier,  /* (I) Higher barrier              */
            char        IoO,          /* (I) Trigger 'I'n or 'O'ut side  */
            char        Smoothing,    /* (I) Smoothing ('Y' or 'N')      */
            int         t,            /* (I) Current time point          */
            int         T,            /* (I) Last time point             */
            int         DCurve,       /* (I) Discount curve              */
            FIX3_DEV_DATA    *dev_data,    /* (I) Fix3_Dev data structure          */
            FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure         */
{

    double  *TriggerL;          /* Local slice pointers */
    double  *UnderL;
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  Intrinsic;          /* Intrinsic value of option   */
    double  IndexStep;          /* Diff between index values   */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Fix3_Dev (Trigger,
             t,
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE)
    {
        goto RETURN;
    }
    

    if (!TriggerFlag)                  /* Nothing to do */
    {
        return (SUCCESS);
    }

    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        TriggerL = Trigger + offset;
        UnderL   = Under   + offset;
        IndexL   = Index   + offset;

        if (Smoothing == 'N')
        {
            if (IoO == 'I')                                                
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] > LowBarrier * (1. - BARRIER_TOL)) 
                     && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        TriggerL[i] = 
                            MAX (TriggerL[i], Notional*CoP*(UnderL[i]-Strike));
                    }
                }  /* for i */  
            }           
            else                                                           
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] < LowBarrier * (1. + BARRIER_TOL)) 
                     || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        TriggerL[i] = 
                            MAX (TriggerL[i], Notional*CoP*(UnderL[i]-Strike));
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {   
                Intrinsic = Notional * CoP * (UnderL[i] - Strike);
                
                if(TriggerL[i] < Intrinsic)
                {
                    /* There is value to american exericse */
                    
                    IndexStep = Fix3_GetIndexStep (Index,
                                              1,                                                  
                                              i, 0, 0,
                                              t,
                                              tree_data);
                                                  
                    if (IoO == 'I')  
                    {                                               
                        if (Smooth_Step (&x,
                                         Intrinsic,  /* up=intrinsic   */
                                         TriggerL[i],/* down=option    */
                                         IndexL[i],
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        if (Smooth_Step (&(TriggerL[i]),
                                         TriggerL[i], /* up=option     */
                                         x,           /* down=smoothed */
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
                                         TriggerL[i],  /* up=option     */
                                         Intrinsic,    /* down=intrinsic*/
                                         IndexL[i],
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        if (Smooth_Step (&(TriggerL[i]),
                                         Intrinsic,       /* up=intrinsic */
                                         x,               /* down=smoothed*/
                                         IndexL[i],
                                         HighBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }  /* if then else */	
                } /* if there is value to american exercise */   
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
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    TriggerL = Trigger + offset;
                    UnderL   = Under   + offset;
                    IndexL   = Index   + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                         && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            TriggerL[j] = 
                              MAX(TriggerL[j],Notional*CoP*(UnderL[j]-Strike));
                        }
                    }  /* for j */
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    TriggerL = Trigger + offset;
                    UnderL   = Under   + offset;
                    IndexL   = Index   + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] < LowBarrier * (1. + BARRIER_TOL)) 
                         || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            TriggerL[j] = 
                              MAX(TriggerL[j],Notional*CoP*(UnderL[j]-Strike));
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                TriggerL = Trigger + offset;
                UnderL   = Under   + offset;
                IndexL   = Index   + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Intrinsic = Notional * CoP * (UnderL[j] - Strike);
                    
                    if(TriggerL[j] < Intrinsic)
                    {
                        /* There is value to american exericse */
                        
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  2,                                                  
                                                  i, j, 0,
                                                  t,
                                                  tree_data);
                                                  
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (&x,
                                             Intrinsic,
                                             TriggerL[j],
                                             IndexL[j],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            
                            if (Smooth_Step (&(TriggerL[j]),
                                             TriggerL[j],
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
                            if (Smooth_Step (&x,
                                             TriggerL[j],
                                             Intrinsic,
                                             IndexL[j],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            
                            if (Smooth_Step (&(TriggerL[j]),
                                             Intrinsic,
                                             x,
                                             IndexL[j],
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }  /* if then else */
                    } /* if american exercise has value */
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */
    else if (tree_data->NbFactor == 3)
    {
        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        TriggerL = Trigger+offset;
                        UnderL   = Under  +offset;
                        IndexL   = Index  +offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                          if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL)) 
                           && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                          {        
                            TriggerL[k] = 
                              MAX(TriggerL[k],Notional*CoP*(UnderL[k]-Strike));
                          }
                        }  /* for k */  
                    }  /* for j */  
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        TriggerL = Trigger + offset;
                        UnderL   = Under   + offset;
                        IndexL   = Index   + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                          if ((IndexL[k] < LowBarrier * (1. + BARRIER_TOL)) 
                           || (IndexL[k] > HighBarrier * (1. - BARRIER_TOL)))
                          {        
                            TriggerL[k] = 
                              MAX(TriggerL[k],Notional*CoP*(UnderL[k]-Strike));
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    TriggerL = Trigger + offset;
                    UnderL   = Under   + offset;
                    IndexL   = Index   + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)          
                    {
                        Intrinsic = Notional * CoP * (UnderL[k] - Strike);
                
                        if(TriggerL[k] < Intrinsic)
                        {
                            /* There is value to american exericse        */
                        
                            IndexStep = Fix3_GetIndexStep (Index,
                                                      3,                                                  
                                                      i, j, k,
                                                      t,
                                                      tree_data);
                                                  
                            if (IoO == 'I')
                            {                                               
                                if (Smooth_Step (&x,
                                                 Intrinsic,
                                                 TriggerL[k],
                                                 IndexL[k],
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step (&(TriggerL[k]),
                                                 TriggerL[k],
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
                                                 TriggerL[k],
                                                 Intrinsic,
                                                 IndexL[k],
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step (&(TriggerL[k]),
                                                 Intrinsic,
                                                 x,
                                                 IndexL[k],
                                                 HighBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }  /* if then else */
                        } /* if american exercise has value */
                    }  /* for k */  
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* TriggerOption_t */



/*****  Fix3_Top_t  ***********************************************************/
/*                                                            
 *       Utility function for the Top option with a trigger feature.  Must 
 *       be called  twice in  order to implement the callable and then the 
 *       puttable.It uses Arnon's smoothing algorithm to american exercise 
 *       premium only.
 *
 *       WARNING: This function  uses the CoP (call or put) flag to calc-
 *       ulate the price of 'long the call' or 'short  the put'.The call-
 *       er of  this function  must therefore be aware that, when pricing 
 *       a simple put, the value will always be <0.
 *
 */
int     Fix3_Top_t  
           (double      *TopOpt,      /* (I/O) Top option prices         */
            double      *Under,       /* (I) Underlying swap prices      */
            double      *Index,       /* (I) Trigger values              */
            double      Notional,     /* (I) Notional of the option      */
            double      Strike,       /* (I) Strike of the option        */
            long        ExerciseFlag, /* (I) Exercise flag               */
            int         CoP,          /* (I) Do the callable or puttable */
            double      LowBarrier,   /* (I) Lower barrier               */
            double      HighBarrier,  /* (I) Higher barrier              */
            char        IoO,          /* (I) Trigger 'I'n or 'O'ut side  */
            char        Smoothing,    /* (I) Smoothing ('Y' or 'N')      */
            int         t,            /* (I) Current time point          */
            int         T,            /* (I) Last time point             */
            int         DCurve,       /* (I) Discount curve              */
            FIX3_DEV_DATA    *dev_data,    /* (I) Fix3_Dev data structure          */
            FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure         */
{

    double  *TopOptL;           /* Local slice pointers */
    double  *UnderL;
    double  *IndexL;

    double  x;                  /* Intermediate smoothed value */
    double  Intrinsic;          /* Intrinsic value of option   */
    double  IndexStep;          /* Diff between index values   */

    int     Top1, Bottom1;      /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;    /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;  /* Tree limits (3rd dim)  */

    int     i, j, k;            /* Node indices           */
    int     offset;             /* Node offset            */
    int     status = FAILURE;   /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Fix3_Dev (TopOpt,
             t,
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE)
    {
        goto RETURN;
    }
    


    if (!ExerciseFlag)                  /* Nothing to do */
    {
        return (SUCCESS);
    }


    if (tree_data->NbFactor == 1)
    {
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        TopOptL = TopOpt + offset;
        IndexL  = Index  + offset;
        UnderL  = Under  + offset;

        if (Smoothing == 'N')
        {
            if (IoO == 'I')                                                
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] > LowBarrier * (1. - BARRIER_TOL)) 
                     && (IndexL[i] < HighBarrier * (1. + BARRIER_TOL)))
                    {        
                        TopOptL[i] = CoP * MAX(CoP*TopOptL[i], 
                                              Notional*CoP*(UnderL[i]-Strike));
                    }
                }  /* for i */  
            }           
            else                                                           
            {        
                for (i = Bottom1; i <= Top1; i ++)
                {
                    if ((IndexL[i] < LowBarrier * (1. + BARRIER_TOL)) 
                     || (IndexL[i] > HighBarrier * (1. - BARRIER_TOL)))
                    {        
                        TopOptL[i] = CoP * MAX(CoP*TopOptL[i], 
                                              Notional*CoP*(UnderL[i]-Strike));
                    }
                }  /* for i */  
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {   
                Intrinsic = Notional * CoP * (UnderL[i] - Strike);
                TopOptL[i] *= CoP;
                
                if(TopOptL[i] < Intrinsic)
                {
                    /* There is value to american exercise */
                    
                    IndexStep = Fix3_GetIndexStep (Index,
                                              1,                                                  
                                              i, 0, 0,
                                              t,
                                              tree_data);
                                                  
                    if (IoO == 'I')  
                    {                                               
                        if (Smooth_Step (&x,
                                         Intrinsic,     /* up=intrinsic   */
                                         TopOptL[i],    /* down=option    */
                                         IndexL[i],
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        if (Smooth_Step (&(TopOptL[i]),
                                         TopOptL[i],     /* up=option     */
                                         x,              /* down=smoothed */
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
                                         TopOptL[i],      /* up=option      */
                                         Intrinsic,       /* down=intrinsic */
                                         IndexL[i],
                                         LowBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                        
                        if (Smooth_Step (&(TopOptL[i]),
                                         Intrinsic,       /* up=intrinsic  */
                                         x,               /* down=smoothed */
                                         IndexL[i],
                                         HighBarrier,
                                         IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }  /* if then else */	
                } /* if there is value to american exercise */   

                TopOptL[i] *= CoP;

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
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    TopOptL = TopOpt + offset;
                    UnderL  = Under  + offset;
                    IndexL  = Index  + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] > LowBarrier * (1. - BARRIER_TOL)) 
                         && (IndexL[j] < HighBarrier * (1. + BARRIER_TOL)))
                        {        
                            TopOptL[j]=CoP*MAX (CoP * TopOptL[j], 
                                              Notional*CoP*(UnderL[j]-Strike));
                        }
                    }  /* for j */  
                }  /* for i */
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    TopOptL = TopOpt + offset;
                    UnderL  = Under  + offset;
                    IndexL  = Index  + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        if ((IndexL[j] < LowBarrier * (1. + BARRIER_TOL)) 
                         || (IndexL[j] > HighBarrier * (1. - BARRIER_TOL)))
                        {        
                            TopOptL[j]=CoP*MAX (CoP * TopOptL[j], 
                                              Notional*CoP*(UnderL[j]-Strike));
                        }
                    }  /* for j */  
                }  /* for i */
            }  /* if then else */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                TopOptL = TopOpt + offset;
                UnderL  = Under  + offset;
                IndexL  = Index  + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    Intrinsic = Notional * CoP * (UnderL[j] - Strike);
                    TopOptL[j] *= CoP;

                    if(TopOptL[j] < Intrinsic)
                    {
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  2,                                                  
                                                  i, j, 0,
                                                  t,
                                                  tree_data);
                                                  
                        if (IoO == 'I')
                        {                                               
                            if (Smooth_Step (&x,
                                             Intrinsic,
                                             TopOptL[j],
                                             IndexL[j],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            
                            if (Smooth_Step (&(TopOptL[j]),
                                             TopOptL[j],
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
                            if (Smooth_Step (&x,
                                             TopOptL[j],
                                             Intrinsic,
                                             IndexL[j],
                                             LowBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                            
                            if (Smooth_Step (&(TopOptL[j]),
                                             Intrinsic,
                                             x,
                                             IndexL[j],
                                             HighBarrier,
                                             IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }  /* if then else */
                    } /* if american exercise has value */

                    TopOptL[j] *= CoP;

                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */
    else if (tree_data->NbFactor == 3)
    {
        if (Smoothing == 'N')
        {
            if (IoO == 'I')
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        TopOptL = TopOpt + offset;
                        UnderL  = Under  + offset;
                        IndexL  = Index  + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] > LowBarrier * (1. - BARRIER_TOL)) 
                             && (IndexL[k] < HighBarrier * (1. + BARRIER_TOL)))
                            {        
                                TopOptL[k]=CoP*MAX(CoP*TopOptL[k], 
                                              Notional*CoP*(UnderL[k]-Strike));
                            }
                        }  /* for k */  
                    }  /* for j */  
            }           
            else
            {        
                for (i = Bottom1; i <= Top1; i ++)          
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        TopOptL = TopOpt + offset;
                        UnderL  = Under  + offset;
                        IndexL  = Index  + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            if ((IndexL[k] < LowBarrier * (1. + BARRIER_TOL)) 
                             || (IndexL[k] > HighBarrier * (1. - BARRIER_TOL)))
                            {        
                                TopOptL[k]=CoP*MAX(CoP*TopOptL[k], 
                                              Notional*CoP*(UnderL[k]-Strike));
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
                    offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                    TopOptL = TopOpt + offset;
                    UnderL  = Under  + offset;
                    IndexL  = Index  + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)          
                    {
                        Intrinsic = Notional * CoP * (UnderL[k] - Strike);
                        TopOptL[k] *= CoP;

                        if(TopOptL[k] < Intrinsic)
                        {
                            IndexStep = Fix3_GetIndexStep (Index,
                                                      3,                                                  
                                                      i, j, k,
                                                      t,
                                                      tree_data);
                                                  
                            if (IoO == 'I')
                            {                                               
                                if (Smooth_Step (&x,
                                                 Intrinsic,
                                                 TopOptL[k],
                                                 IndexL[k],
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step (&(TopOptL[k]),
                                                 TopOptL[k],
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
                                                 TopOptL[k],
                                                 Intrinsic,
                                                 IndexL[k],
                                                 LowBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                                
                                if (Smooth_Step (&(TopOptL[k]),
                                                 Intrinsic,
                                                 x,
                                                 IndexL[k],
                                                 HighBarrier,
                                                 IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }  /* if then else */
                        } /* if american exercise has value */

                        TopOptL[k] *= CoP;

                    }  /* for k */  
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* TopOption_t */
