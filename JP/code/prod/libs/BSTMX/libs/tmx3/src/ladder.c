/**************************************************************************/
/*      Payoff function for the ladder swap                               */
/**************************************************************************/
/*      ladder.c                                                          */
/**************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tmx123head.h"


/*****  LadderSwap_x  *****************************************************/
/*
*   Calculates the ladder swap price
*   no dev's, just add ladder and floating coupons
*/
int  LadderSwap_x
             (double    **Sticky,         /* (O) Prices for all states   */
              long        ResetFlagSt,    /* (I) Reset at timept         */
              long        ResetFlagFund,  /* (I) Reset at timept         */
              double     *Step,           /* (I) Step for ladder         */
              double     *Funding,        /* (I) Funding                 */
              double     *ZeroToPmtSt,    /* (I) PmtZero 4 this reset    */
              double      FloorSt,        /* (I) Floor spread            */
              double      CapSt,          /* (I) Cap spread              */
              double      OutsSt,         /* (I) Outstanding             */
              double      DcfSt,          /* (I) Day count fract         */
              char        SoZ,            /* (I) Swap/zero coupon        */
              char        CompSt,         /* (I) Simple/Compound pmt     */
              int         NbStates,       /* (I) Nb of states            */
              double     *CurrStates,     /* (I) Curr state levels at i-1*/
              double     *PrevStates,     /* (I) Prev state levels at i  */
              int         t,              /* (I) Current time point      */
              TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *StickyL[MAXNBSTATES+1];
    double  *FundingL;
    double  *ZeroToPmtStL;
    double  *StepL;

    /* state-variable variables */
    int     s;                        /* State variable index            */

    /* Payoff variables */
    double  tmpSticky[MAXNBSTATES];   /* Intermediate values of ladder   */
    double  PRate[MAXNBSTATES];

    double  X;
    int     isZeroSt, isSimpleSt;

    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */

    long    CurrDate = tree_data->TPDate[t];


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* STICKY LEG */


    isZeroSt   = (SoZ == 'Z');
    isSimpleSt = (CompSt == 'S');

    /* Add ladder payoff if it's a reset */
    if (ResetFlagSt)
    {
        if (PrevStates == NULL ||
            CurrStates == NULL)
        {
            DR_Error("LadderSwap_x: No previous or current states on %8d!", 
                     CurrDate);
            goto RETURN;
        }


        /* Note: outs is 1 for zero coupon */
        for (s=0; s<NbStates; s++)
        {
            PRate[s] =  ACC_FN(PrevStates[s],DcfSt,isSimpleSt);
        }

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            ZeroToPmtStL = ZeroToPmtSt + offset;
            StepL  = Step  + offset;


            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* ---------------------------- */
                /*   Update current value:      */
                /*   mulitply and add coupon    */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    tmpSticky[s] = StickyL[s][i];

                    if (isZeroSt) tmpSticky[s] = StickyL[s][i] * (1.+PRate[s]);

                    tmpSticky[s] += ZeroToPmtStL[i]*PRate[s]*OutsSt;
                }

                /* ---------------------------- */
                /*   Rearrange variables        */
                /* ---------------------------- */
                for (s = 0; s < NbStates; s++)
                {
                    /* ------------------------------------- */
                    /*   Find next coupon level              */
                    /* ------------------------------------- */
                    
                    X  = CurrStates[s];
                    X += StepL[i]; 
                    X  = COLLAR(X,CapSt,FloorSt);
        

                    /* ------------------------------------- */
                    /*   Payoff = interpolated Sticky price  */
                    /* ------------------------------------- */
                    if (DoubleQuadraticInterp(PrevStates, 
                                              tmpSticky, 
                                              NbStates,
                                              X, 
                                             &(StickyL[s][i])) == FAILURE)
                        goto RETURN;

                } /* for state idx */

            } /* for i */

        } /* if NbFactor == 1 */

        /************************  2 FACTORS   ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                StepL  = Step  + offset;
                ZeroToPmtStL = ZeroToPmtSt + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* ---------------------------- */
                    /*   Update current value:      */
                    /*   mulitply and add coupon    */
                    /* ---------------------------- */
                    for (s = 0; s < NbStates; s++)
                    {
                        tmpSticky[s] = StickyL[s][j];

                        if (isZeroSt) tmpSticky[s] = 
                                     StickyL[s][j] * (1.+PRate[s]);

                        tmpSticky[s] += ZeroToPmtStL[j]*PRate[s]*OutsSt;
                    }
                
                    for (s = 0; s < NbStates; s++)
                    {
                        /* ------------------------------------- */
                        /*   Find next coupon level              */
                        /* ------------------------------------- */
                    
                        X  = CurrStates[s];
                        X += StepL[j];
                        X  = COLLAR(X,CapSt,FloorSt);

                        /* ------------------------------------- */
                        /*   Payoff = interpolated Sticky price  */
                        /* ------------------------------------- */
        
                        if (DoubleQuadraticInterp(PrevStates, 
                                                  tmpSticky, 
                                                  NbStates,
                                                  X, 
                                                 &(StickyL[s][j])) == FAILURE)
                            goto RETURN;
    
                    } /* for state idx */
    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTORS   ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }
                    StepL  = Step  + offset;
                    ZeroToPmtStL = ZeroToPmtSt + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* ---------------------------- */
                        /*   Update current value:      */
                        /*   mulitply and add coupon    */
                        /* ---------------------------- */
                        for (s = 0; s < NbStates; s++)
                        {
                            tmpSticky[s] = StickyL[s][k];

                            if (isZeroSt) tmpSticky[s] = 
                                         StickyL[s][k] * (1.+PRate[s]);

                            tmpSticky[s] += ZeroToPmtStL[k]*PRate[s]*OutsSt;
                        }
        
                        for (s = 0; s < NbStates; s++)
                        {
                            /* ------------------------------------- */
                            /*   Find next coupon level              */
                            /* ------------------------------------- */
                    
                            X  = CurrStates[s];
                            X += StepL[k];
                            X  = COLLAR(X,CapSt,FloorSt);

                            /* ------------------------------------- */
                            /*   Payoff = interpolated Sticky price  */
                            /* ------------------------------------- */
            
                            DoubleQuadraticInterp(PrevStates, 
                                                  tmpSticky, 
                                                  NbStates,
                                                  X, 
                                                 &(StickyL[s][k]));
        
                        } /* for state idx */

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
        

    }  /* if ResetFlagSticky */

    /* FUNDING LEG */

    /* Add floating payoff if it's a reset BUT ONLY if it's ladder date */
    /* This restriction is achieved by setting flags in calc properly   */
    if (ResetFlagFund)
    {

        /************************  1 FACTOR    ****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            for (s = 0; s < NbStates ; s++)
            {
                StickyL[s] = Sticky[s] + offset;
            }
            FundingL = Funding + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                for (s = 0; s < NbStates ; s++)  
                {
                    StickyL[s][i] -= FundingL[i];
                }
                FundingL[i] = 0.;

            } /* for i */

        } 

        /************************  2 FACTOR    ****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                for (s = 0; s < NbStates ; s++)
                {
                    StickyL[s] = Sticky[s] + offset;
                }
                FundingL = Funding + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    for (s = 0; s < NbStates ; s++)  
                    {
                        StickyL[s][j] -= FundingL[j];
                    }
                    FundingL[j] = 0.;
                                    
                }  /* for j */
            }  /* for i */
        }

        /************************  3 FACTOR    ****************************/

        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    for (s = 0; s < NbStates ; s++)
                    {
                        StickyL[s] = Sticky[s] + offset;
                    }
                    FundingL = Funding + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        for (s=0; s<NbStates ; s++)  
                        {
                            StickyL[s][k] -= FundingL[k];
                        }
                        FundingL[k] = 0.;

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
 
    }  /* if ResetFlagFloat */

    status = SUCCESS;
    
RETURN:

    return (status);

}  /* LadderSwap_x */



/*****  LadderStep_x  ********************************************************/
/*                                                                           
*       Smooth 3-level step function.
*/
int     LadderStep_x (double      *Step,       /* (I/O) Smooth step          */
                      double      *Index,      /* (I) Index prices           */
                      double      LowBarrier,  /* (I) Lower barrier          */
                      double      HighBarrier, /* (I) Higher barrier         */
                      double      DownValue,   /* (I) Down value             */
                      double      MidValue,    /* (I) Mid value              */
                      double      UpValue,     /* (I) Up value               */
                      char        Smoothing,   /* (I) Smoothing ('Y' or 'N') */
                      int         t,           /* (I) Current time point     */
                      TREE_DATA   *tree_data)  /* (I) Structure of tree data */
{

    double  *StepL;             /* Local slice pointers */
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


    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        StepL = Step + offset;
        IndexL = Index + offset;

        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                if (IndexL[i] <= LowBarrier)
                {
                    StepL[i] = DownValue;
                }
                else 
                if (IndexL[i] <= HighBarrier)
                {        
                    StepL[i] = MidValue;
                }
                else
                {
                    StepL[i] = UpValue;
                }
            }  /* for i */  
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
                
                /* 'up' value = mid value    */
                /* 'down' value = down value */
                if (Smooth_Step (   &x,
                                    MidValue, 
                                    DownValue,
                                    IndexL[i],
                                    LowBarrier,
                                    IndexStep) == FAILURE)
                {
                    goto RETURN;
                }
                                
                /* 'up' value = up value         */
                /* 'down' value = smoothed value */
                if (Smooth_Step (   &(StepL[i]),
                                    UpValue,
                                    x,
                                    IndexL[i],
                                    HighBarrier,
                                    IndexStep) == FAILURE)
                {
                    goto RETURN;
                }
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 2)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);
  
                StepL = Step + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    if (IndexL[j] <= LowBarrier)
                    {
                        StepL[j] = DownValue;
                    }
                    else 
                    if (IndexL[j] <= HighBarrier)
                    {        
                        StepL[j] = MidValue;
                    }
                    else
                    {
                        StepL[j] = UpValue;
                    }
                }  /* for j */  
            }  /* for i */
        }           
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                StepL = Step + offset;
                IndexL = Index + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    IndexStep = GetIndexStep (Index,
                                              2,
                                              i, j, 0,
                                              t,
                                              tree_data);
                                                  
                    /* 'up' value = mid value    */
                    /* 'down' value = down value */
                    if (Smooth_Step (   &x,
                                        MidValue,
                                        DownValue,
                                        IndexL[j],
                                        LowBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                                
                    /* 'up' value = up value         */
                    /* 'down' value = smoothed value */
                    if (Smooth_Step (   &(StepL[j]),
                                        UpValue,
                                        x,
                                        IndexL[j],
                                        HighBarrier,
                                        IndexStep) == FAILURE)
                    {
                        goto RETURN;
                    }
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */   
    }
    else if (tree_data->NbFactor == 3)
    {
        if (Smoothing == 'N')
        {
            for (i = Bottom1; i <= Top1; i ++)          
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    StepL = Step + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        if (IndexL[k] <= LowBarrier)
                        {
                            StepL[k] = DownValue;
                        }
                        else 
                        if (IndexL[k] <= HighBarrier)
                        {        
                            StepL[k] = MidValue;
                        }
                        else
                        {
                            StepL[k] = UpValue;
                        }
                    }  /* for k */  
                }  /* for j */  
        }
        else
        {           
            for (i = Bottom1; i <= Top1; i ++)          
            {
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    StepL = Step + offset;
                    IndexL = Index + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        IndexStep = GetIndexStep (Index,
                                                  3,
                                                  i, j, k,
                                                  t,
                                                  tree_data);
                                                  
                        /* 'up' value = mid value    */
                        /* 'down' value = down value */
                        if (Smooth_Step (   &x,
                                            MidValue,
                                            DownValue,
                                            IndexL[k],
                                            LowBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                                    
                        /* 'up' value = up value         */
                        /* 'down' value = smoothed value */
                        if (Smooth_Step (   &(StepL[k]),
                                            UpValue,
                                            x,
                                            IndexL[k],
                                            HighBarrier,
                                            IndexStep) == FAILURE)
                        {
                            goto RETURN;
                        }
                    }  /* for k */  
                }  /* for j */  
            }  /* for i */  
        }  /* if then else */
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    return (status);

}  /* LadderStep_x */

