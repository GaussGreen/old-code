/****************************************************************************/
/*      Chooser and survivor caps.             	                            */
/****************************************************************************/
/*      CHOOSER.c                                                           */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Fix3_Chooser_t  **********************************************************/
/*
*       Chooser cap price.
*/
int     Fix3_Chooser_t ( double    **Chooser,    /* (I/O) Chooser caps/floors     */
                    double     *Caplet,     /* (I) Current caplet/floorlet   */
                    long        CapletFlag, /* (I) Reset flag                */
                    int         NbCaplet,   /* (I) Maximum number of caplets */
                    int         t,          /* (I) Current time point        */
                    int         T,          /* (I) Last time point           */
                    int         DCurve,     /* (I) Discount curve            */
                    FIX3_DEV_DATA    *dev_data,	/* (I) Fix3_Dev data structure	     */
                    FIX3_TREE_DATA   *tree_data)	/* (I) Tree data structure 	     */
{

    double  *Chooser1L;                      /* Local slice pointers   */
    double  *Chooser2L;
    double  *CapletL;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */

    int     i, j, k;                        /* Node indices           */
    int     offset;                         /* Node offset            */
    int     l;                              /* Chooser cap index      */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* Discount all the state variables */
    for (l = 1; l <= NbCaplet; l++)
    {
        if (Fix3_Dev (   Chooser[l],
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
                    
        }  /* if */
    }  /* for l */


    if (CapletFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            for (l = NbCaplet; l > 0 ; l--) /* Need reverse order here */
            {
                offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

                Chooser1L = Chooser[l]   + offset;
                Chooser2L = Chooser[l-1] + offset;
                CapletL   = Caplet       + offset;
    
                for (i = Bottom1; i <= Top1; i ++)
                {
                    Chooser1L[i] = MAX (Chooser1L[i], Chooser2L[i] + CapletL[i]);
                }
            }  /* for l */
        }
        else if (tree_data->NbFactor == 2)
        {
            for (l = NbCaplet; l > 0 ; l--)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                    Chooser1L = Chooser[l]   + offset;
                    Chooser2L = Chooser[l-1] + offset;
                    CapletL   = Caplet       + offset;
    
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        Chooser1L[j] = MAX (Chooser1L[j], Chooser2L[j] + CapletL[j]);
                    }
                }  /* for i */
            }  /* for l */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (l = NbCaplet; l > 0 ; l--)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)	            
                    {
                        offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                        Chooser1L = Chooser[l]   + offset;
                        Chooser2L = Chooser[l-1] + offset;
                        CapletL   = Caplet       + offset;
    
                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                        {
                            Chooser1L[k] = MAX (Chooser1L[k], Chooser2L[k] + CapletL[k]);
                        }
                    }  /* for j */        
            }  /* for l */
        }  /* if then else */
    }  /* if */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Chooser_t */



/*****  Fix3_Survivor_t  *********************************************************/
/*
*       Survivor cap price with smoothing.
*/
int     Fix3_Survivor_t (double    **Survivor,   /* (I/O) Survivor caps/floors    */
                    double     *Caplet,     /* (I) Current caplet/floorlet   */
                    double     *Index,      /* (I) Floating index	         */
                    long        CapletFlag, /* (I) Reset flag                */
                    int         NbCaplet,   /* (I) Maximum number of caplets */
                    int         CoF,        /* (I) 1 for cap, -1 for floor 	 */
                    double      Strike,     /* (I) Strike 	                 */
                    char        Smoothing,	/* (I) Smoothing ('Y' or 'N')	 */
                    int         t,          /* (I) Current time point        */
                    int         T,          /* (I) Last time point           */
                    int         DCurve,     /* (I) Discount curve            */
                    FIX3_DEV_DATA    *dev_data,	/* (I) Fix3_Dev data structure	     */
                    FIX3_TREE_DATA   *tree_data)	/* (I) Tree data structure 	     */
{

    double  *Survivor1L;             /* Local slice pointers   */
    double  *Survivor2L;             /* Local slice pointers   */
    double  *CapletL;
    double  *IndexL;

    double  IndexStep;              /* Diff between index values */
            
    int     Top1, Bottom1;          /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;        /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;      /* Tree limits (3rd dim)  */

    int     i, j, k;                /* Node indices           */
    int     offset;                 /* Node offset            */
    int     l;                      /* Survivor cap index     */
    int     status = FAILURE;       /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    for (l = 1; l <= NbCaplet; l++)
    {
        if (Fix3_Dev (   Survivor[l],
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;
                    
        }  /* if */
    }  /* for l */


    if (CapletFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            if (Smoothing == 'N')
            {
                for (l = NbCaplet; l > 0 ; l--) /* Need reverse order here */
                {
                    offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

                    Survivor1L = Survivor[l]   + offset;
                    Survivor2L = Survivor[l-1] + offset;
                    CapletL    = Caplet        + offset;
    
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        if (CapletL[i] > ERROR)
                        {
                            Survivor1L[i] = Survivor2L[i] + CapletL[i];
                        }
                    }  /* for i */	
                }  /* for l */	
            }	        
            else
            {	        
                for (l = NbCaplet; l > 0 ; l--)
                {
                    offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

                    Survivor1L = Survivor[l]   + offset;
                    Survivor2L = Survivor[l-1] + offset;
                    CapletL    = Caplet        + offset;
                    IndexL     = Index         + offset;
    
                    for (i = Bottom1; i <= Top1; i ++)	        
                    {
                        IndexStep = Fix3_GetIndexStep (Index,
                                                  1,                                                  
                                                  i, 0, 0,
                                                  t,
                                                  tree_data);
                                                  
                        if (CoF == 1)
                        {
                            /* 'up' value = add exercised caplet */
                            /* 'down' value = keep current value */
                            if (Smooth_Step (	&(Survivor1L[i]),
                                                Survivor2L[i] + CapletL[i],	
                                                Survivor1L[i],                
                                                IndexL[i],
                                                Strike,
                                                IndexStep) == FAILURE)
                            {                           
                                goto RETURN;
                            }
                        }	
                        else    
                        {
                            /* 'up' value = keep current value     */
                            /* 'down' value = add exercised caplet */
                            if (Smooth_Step (	&(Survivor1L[i]),
                                                Survivor1L[i],                
                                                Survivor2L[i] + CapletL[i],	
                                                IndexL[i],
                                                Strike,
                                                IndexStep) == FAILURE)
                            {
                                goto RETURN;
                            }
                        }  /* if then else */
                    }  /* for i */	
                }  /* for l */
            }  /* if then else */	    
        }
        else if (tree_data->NbFactor == 2)
        {
            if (Smoothing == 'N')
            {
                for (l = NbCaplet; l > 0 ; l--)
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        Survivor1L = Survivor[l]   + offset;
                        Survivor2L = Survivor[l-1] + offset;
                        CapletL    = Caplet        + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            if (CapletL[j] > ERROR)
                            {
                                Survivor1L[j] = Survivor2L[j] + CapletL[j];
                            }
                        }  /* for j */	
                    }  /* for i */
                }  /* for l */	
            }	        
            else
            {	        
                for (l = NbCaplet; l > 0 ; l--)
                {
                    for (i = Bottom1; i <= Top1; i ++)	        
                    {
                        offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

                        Survivor1L = Survivor[l]   + offset;
                        Survivor2L = Survivor[l-1] + offset;
                        CapletL    = Caplet        + offset;
                        IndexL     = Index         + offset;
    
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            IndexStep = Fix3_GetIndexStep (Index,
                                                      2,                                                      
                                                      i, j, 0,
                                                      t,
                                                      tree_data);
                                                  
                            if (CoF == 1)
                            {
                                if (Smooth_Step (	&(Survivor1L[j]),
                                                    Survivor2L[j] + CapletL[j],
                                                    Survivor1L[j],
                                                    IndexL[j],
                                                    Strike,
                                                    IndexStep) == FAILURE)
                                {                           
                                    goto RETURN;
                                }
                            }	
                            else    
                            {
                                if (Smooth_Step (	&(Survivor1L[j]),
                                                    Survivor1L[j],
                                                    Survivor2L[j] + CapletL[j],
                                                    IndexL[j],
                                                    Strike,
                                                    IndexStep) == FAILURE)
                                {
                                    goto RETURN;
                                }
                            }  /* if then else */
                        }  /* for j */	
                    }  /* for i */	
                }  /* for l */
            }  /* if then else */	    
        }
        else if (tree_data->NbFactor == 3)
        {
            if (Smoothing == 'N')
            {
                for (l = NbCaplet; l > 0 ; l--)
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            Survivor1L = Survivor[l]   + offset;
                            Survivor2L = Survivor[l-1] + offset;
                            CapletL    = Caplet        + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                if (CapletL[k] > ERROR)
                                {
                                    Survivor1L[k] = Survivor2L[k] + CapletL[k];
                                }
                            }  /* for k */
                        }  /* for j */
                }  /* for l */	
            }	        
            else
            {	        
                for (l = NbCaplet; l > 0 ; l--)
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                            Survivor1L = Survivor[l]   + offset;
                            Survivor2L = Survivor[l-1] + offset;
                            CapletL    = Caplet        + offset;
                            IndexL     = Index         + offset;
    
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                            {
                                IndexStep = Fix3_GetIndexStep (Index,
                                                          3,
                                                          i, j, k,
                                                          t,
                                                          tree_data);
                                                  
                                if (CoF == 1)
                                {
                                    if (Smooth_Step (	&(Survivor1L[k]),
                                                        Survivor2L[k] + CapletL[k],
                                                        Survivor1L[k],
                                                        IndexL[k],
                                                        Strike,
                                                        IndexStep) == FAILURE)
                                    {                           
                                        goto RETURN;
                                    }
                                }	
                                else    
                                {
                                    if (Smooth_Step (	&(Survivor1L[k]),
                                                        Survivor1L[k],
                                                        Survivor2L[k] + CapletL[k],
                                                        IndexL[k],
                                                        Strike,
                                                        IndexStep) == FAILURE)
                                    {
                                        goto RETURN;
                                    }
                                }  /* if then else */
                            }  /* for k */	
                        }  /* for j */	
                    }  /* for i */	
                }  /* for l */
            }  /* if then else */	    
        }  /* if then else */
    }  /* if */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Survivor_t */
