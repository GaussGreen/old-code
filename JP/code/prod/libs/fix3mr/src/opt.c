/****************************************************************************/
/*      Option price.                                                       */
/****************************************************************************/
/*      OPT.c                                                               */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/opt.c,v 1.6 1998/10/23 18:09:54 dfung Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Option_t  ***********************************************************/
/*
*       Standard option.
*/
int     Option_t (  double      *Opt,         /* (I/O) Option prices         */
                    double      *Under,       /* (I) Underlying prices       */
                    double      Notional,     /* (I) Notional                */
                    double      Strike,       /* (I) Strike                  */
                    long        ExerFlag,     /* (I) Exercise flag           */
                    int         CoP,          /* (I) =1 for call, -1 for put */
                    int         t,            /* (I) Current time point      */
                    int         T,            /* (I) Last time point         */
                    int         DCurve,       /* (I) Discount curve          */
                    DEV_DATA    *dev_data,    /* (I) Dev data structure      */
                    TREE_DATA   *tree_data)	  /* (I) Tree data structure     */
{

    double  *OptL;                          /* Local slice pointer    */
    double  *UnderL;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */

    int     i, j, k;                        /* Node indices           */
    int     offset;                         /* Node offset            */
    int     status = FAILURE;               /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    if (Dev (   Opt,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;
    }


    if (ExerFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            OptL   = Opt   + offset;
            UnderL = Under + offset;
    
            /* 
            *   If this is the last exercise date OptL has been initialized to
            *   zero so the exercise decision is Max (intrinsic value, 0).
            *   If this is not the last exercise date we have an American 
            *   exercise decision Max (intrinsic value, live option).
            */

            for (i = Bottom1; i <= Top1; i ++)
            {
                OptL[i] = MAX (Notional * CoP * (UnderL[i] - Strike), OptL[i]);
            }
        }
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                OptL   = Opt   + offset;
                UnderL = Under + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    OptL[j] = MAX (Notional * CoP*(UnderL[j]-Strike), OptL[j]);
                }
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    OptL   = Opt   + offset;
                    UnderL = Under + offset;
    
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        OptL[k]=MAX (Notional*CoP*(UnderL[k]-Strike), OptL[k]);
                    }
                }  /* for j */
        }  /* if then else */
    }  /* if */    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Option_t */


/*****  OptionPlus_t  ********************************************************/
/*
 *      Standard option plus the following stats:
 *      - prob of exericse
 *      - expected time to exercise
 *      - expected time to exercise squared (for variance calculation)
 *
 *      Only perform calculations for stats info if slice is not NULL
 *
 */
int OptionPlus_t (double     *Opt,          /* (I/O) Option prices         */
                  double     *ExerProb,     /* (I/O) exer prob             */
                  double     *ExerTime,     /* (I/O) time to exer          */
                  double     *ExerTimeSqr,  /* (I/O) time to exer squared  */
                  double     *Under,        /* (I) Underlying prices       */
                  double      Notional,     /* (I) Notional                */
                  double      Strike,       /* (I) Strike                  */
                  long        ExerFlag,     /* (I) Exercise flag           */
                  int         CoP,          /* (I) =1 for call, -1 for put */
                  int         SmoothingOn,  /* (I) TRUE = smoothing on     */
                  double     *AuxSlice,     /* (I) aux slice for smoothing */
                  int         t,            /* (I) Current time point      */
                  int         T,            /* (I) Last time point         */
                  int         DCurve,       /* (I) Discount curve          */
                  DEV_DATA   *dev_data,     /* (I) Dev data structure      */
                  TREE_DATA  *tree_data)	/* (I) Tree data structure     */
{
    double  *OptL;                          /* Local slice pointer         */
    double  *AuxSliceL;
    double  *ExerProbL    = NULL;
    double  *ExerTimeL    = NULL;
    double  *ExerTimeSqrL = NULL;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim)  */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)   */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)   */

    int     i, j, k;                        /* Node indices            */
    int     offset;                         /* Node offset             */
    int     status = FAILURE;               /* Error status            */

    double  CurrTimeInYrs;                  /* for ExerTime stats      */
    double  CurrTimeInYrsSqr;               /* for ExerTimeSqr stats   */
    double  Step;                           /* step size for smoothing */

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    if (Dev (Opt,
             t,
             T,
             DCurve,
             dev_data,
             tree_data) == FAILURE) goto RETURN;

    if (ExerProb != NULL)
    {
        if (Ev(ExerProb, t, T, dev_data, tree_data) == FAILURE) goto RETURN;
    }

    if (ExerTime != NULL)
    {
        if (Ev(ExerTime, t, T, dev_data, tree_data) == FAILURE) goto RETURN;
    }

    if (ExerTimeSqr != NULL)
    {
        if (Ev(ExerTimeSqr, t, T, dev_data, tree_data) == FAILURE) goto RETURN;
    }

    if (ExerFlag)
    {
        CurrTimeInYrs = Daysact(tree_data->TPDate[0],
                                tree_data->TPDate[t]) / 365.0;

        CurrTimeInYrsSqr = CurrTimeInYrs * CurrTimeInYrs;
        
        Step = 0.0;  /* initialise to smoothing off */

        /* Store "Notional * CoP * (Underlying - Strike) - Opt" in AuxSlice */

        if (LCombTwoSlices(AuxSlice,
                           Under,
                           (Notional * CoP),
                           Opt,  
                           -1.0, 
                           t, 
                           tree_data) == FAILURE) goto RETURN;

        if (AddScalar(AuxSlice,
                      (-1.0 * Notional * CoP * Strike),
                      t,
                      tree_data) == FAILURE) goto RETURN;


        /*************************** 1 Factor *****************************/

        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);

            OptL      = Opt      + offset;
            AuxSliceL = AuxSlice + offset;
            if (ExerProb    != NULL) ExerProbL    = ExerProb + offset;
            if (ExerTime    != NULL) ExerTimeL    = ExerTime + offset;
            if (ExerTimeSqr != NULL) ExerTimeSqrL = ExerTimeSqr + offset;
    
            for (i = Bottom1; i <= Top1; i ++)
            {
                if (SmoothingOn)
                {
                    Step = GetIndexStep(AuxSlice,1,i,0,0,t,tree_data);
                }

                if (ExerProb != NULL)
                    ExerProbL[i] = DrSmoothStep(1.0, 
                                                ExerProbL[i], 
                                                AuxSliceL[i],
                                                0.0, 
                                                Step);

                if (ExerTime != NULL)
                    ExerTimeL[i] = DrSmoothStep(CurrTimeInYrs,
                                                ExerTimeL[i],
                                                AuxSliceL[i],
                                                0.0,
                                                Step);

                if (ExerTimeSqr != NULL)
                    ExerTimeSqrL[i] = DrSmoothStep(CurrTimeInYrsSqr,
                                                   ExerTimeSqrL[i],
                                                   AuxSliceL[i],
                                                   0.0,
                                                   Step);

                OptL[i] = DrSmoothMax(AuxSliceL[i], 0.0, Step) + OptL[i];

            } /* for i */
        }

        /*************************** 2 Factor *****************************/

        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);

                OptL      = Opt      + offset;
                AuxSliceL = AuxSlice + offset;
                if (ExerProb    != NULL) ExerProbL    = ExerProb + offset;
                if (ExerTime    != NULL) ExerTimeL    = ExerTime + offset;
                if (ExerTimeSqr != NULL) ExerTimeSqrL = ExerTimeSqr + offset;
    
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    if (SmoothingOn)
                    {
                        Step = GetIndexStep(AuxSlice,2,i,j,0,t,tree_data);
                    }

                    if (ExerProb != NULL)
                        ExerProbL[j] = DrSmoothStep(1.0,
                                                    ExerProbL[j],
                                                    AuxSliceL[j], 
                                                    0.0, 
                                                    Step);

                    if (ExerTime != NULL)
                        ExerTimeL[j] = DrSmoothStep(CurrTimeInYrs,
                                                    ExerTimeL[j],
                                                    AuxSliceL[j], 
                                                    0.0, 
                                                    Step);

                    if (ExerTimeSqr != NULL)
                        ExerTimeSqrL[j] = DrSmoothStep(CurrTimeInYrsSqr,
                                                       ExerTimeSqrL[j],
                                                       AuxSliceL[j], 
                                                       0.0, 
                                                       Step);

                    OptL[j] = DrSmoothMax(AuxSliceL[j], 0.0, Step) + OptL[j];

                }
            }  /* for i */
        }

        /*************************** 3 Factor *****************************/
        
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);

                    OptL      = Opt      + offset;
                    AuxSliceL = AuxSlice + offset;
                    if (ExerProb    != NULL) ExerProbL    = ExerProb + offset;
                    if (ExerTime    != NULL) ExerTimeL    = ExerTime + offset;
                    if (ExerTimeSqr != NULL) 
                        ExerTimeSqrL = ExerTimeSqr + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        if (SmoothingOn)
                        {
                            Step = GetIndexStep(AuxSlice,3,i,j,k,t,tree_data);
                        }

                        if (ExerProb != NULL)
                            ExerProbL[k] = DrSmoothStep(1.0,
                                                        ExerProbL[k],
                                                        AuxSliceL[k],
                                                        0.0, 
                                                        Step);

                        if (ExerTime != NULL)
                            ExerTimeL[k] = DrSmoothStep(CurrTimeInYrs,
                                                        ExerTimeL[k],
                                                        AuxSliceL[k],
                                                        0.0, 
                                                        Step);

                        if (ExerTimeSqr != NULL)
                            ExerTimeSqrL[k] = DrSmoothStep(CurrTimeInYrsSqr,
                                                           ExerTimeSqrL[k],
                                                           AuxSliceL[k],
                                                           0.0, 
                                                           Step);

                        OptL[k] = DrSmoothMax(AuxSliceL[k], 0.0, Step) 
                                  + OptL[k];
                    }
                }  /* for j */
        }  /* if then else */
    }  /* if */    
   
    status = SUCCESS;

    RETURN:

    return (status);

}  /* OptionPlus_t */
