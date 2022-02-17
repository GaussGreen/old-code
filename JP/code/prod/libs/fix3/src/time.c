/****************************************************************************/
/*      Construct the time arrays of the tree (time steps, ...).            */
/****************************************************************************/
/*      TIME.c                                                              */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "fix123head.h"




/*****  Fix3_Time_Line  **********************************************************/
/*
*       Set up the time step array so that a time point falls on each 
*       critical date. Allocate memory for the tree arrays.
*       Each set of dates can be discrete (e.g.:coupon payments) or continuous
*       (e.g. Knock-out dates). For the latter we need to interpolate values
*       between input dates.
*
*       The time step for the timeline and jumpsize are set according to the
*       StepStyle:
*
*       'I'/'E': timeline is built with an increasing/equal timestep 
*                (a function of tree_data->Ppy) and jumpsize is calculated 
*                using the corresponding theoretical timestep (also a function
*                of tree_data->Ppy)
*           'J': timeline is built with an increasing timestep and jumpsize
*                is calculated using a modified timestep (a function of 
*                tree_data->JumpPpy)
*/
int     Fix3_Time_Line (   
                    long        ValueDate,  /* (I) Value date                */
                    int         NbDate,     /* (I) Number of dates           */
                    CRIT_DATE   *Date,      /* (I) Critical date list        */
                    char        StepStyle,  /* (I) Equal or increasing steps */
                    FIX3_TREE_DATA   *tree_data) /* (O) Structure of tree data    */
{

    double  x;

    int     *TPinIntval=NULL; /* Nb of time points in current interval    */
    int     Intval;           /* Interval between two critical dates      */
    int     NbStep;           /* Number of time steps in current interval */
    int     MaxStep=0;        /* Maximum length for a time step           */
    int     Step;             /* Length of current time step              */
    int     Step0;            /* Length of current time step after adjust */
    int     LastStep;         /* Length of last step of previous interval */
    int     Rest, Rest0;      /* Days left over in current interval       */
    int     NbDaysSoFar;
    int     NbDailyPts=0;     /* Number of daily pts to be inserted from
                                 value date on                            */

    int     i0;               /* Index of point crossing the boundary     */
    int     m;                /* Index for intermediate time points       */
    int     d;                /* Nb of pts from VD  to critical date      */
    int     UnderNDP = FALSE; /* Under number of daily points             */
    int     NbTP=0;           /* Total number of time points in the tree  */
    int     i, j, k=0, k1, k2, l;
    

    int     status = FAILURE; /* Error status = FAILURE initially         */


    if ( (StepStyle != 'I') && (StepStyle != 'E') && (StepStyle != 'J')  )
    {
        DR_Error ("Fix3_Time_Line: invalid step style!");
        goto FREE_MEM_AND_RETURN;
    }

    TPinIntval = (int *) DR_Array (INT, 0, NbDate+1);
    if (TPinIntval == NULL)
    {
        DR_Error ("Fix3_Time_Line: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }


    /* 
    *   Sort critical date list.
    */

    if (Sort_CritDate ( NbDate,
                        Date) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }
        
    if (Date[0].CritDate < ValueDate)
    {
        DR_Error("Fix3_Time_Line: first date %ld falls before value date %ld",
                 Date[0].CritDate,
                 ValueDate);
        goto FREE_MEM_AND_RETURN;
    }

   
    /* 
    *   Build time line i.e. put time points between critical dates.
    */

    /* Increasing time steps */
    if ( (StepStyle == 'I') || (StepStyle == 'J') )
    {
          /* Ppy parameter gives maximum time step in days */
        MaxStep = (int) floor (372. / tree_data->Ppy);   
        
        if (MaxStep < 1)
        {
            DR_Error ("Fix3_Time_Line: maximum step less than 1");
            goto FREE_MEM_AND_RETURN;
        }

        /* If the maximum distance between value date and last critical date
           is smaller than  desired NDP, insert daily points to the end   */
        if ((int)Daysact (ValueDate, Date[NbDate - 1].CritDate) < tree_data->NbDailyPts)
        { 
            NbDailyPts = (int)Daysact (ValueDate, Date[NbDate - 1].CritDate);
        }
        else
        {
            NbDailyPts = tree_data->NbDailyPts;
        }


        NbTP = 0;
        Step = 0;
        i0 = 0;

        /* 
        *	Calculate total number of time points.
        */
        for (i = 0; i < NbDate - 1; i++)                            
        {
            if ((int)Daysact (ValueDate, Date[i].CritDate) <= NbDailyPts)
            {
                Intval = 0; NbTP = NbDailyPts;
                if( (int)Daysact(ValueDate, Date[i + 1].CritDate) >  NbDailyPts)
                {
                    Intval = (int) Daysact (ValueDate, Date[i + 1].CritDate) 
                             - NbDailyPts ;
                    i0 = i;
                }      
            }
            else
            {
                Intval = (int) Daysact (Date[i].CritDate, Date[i+1].CritDate) ;
            }
                
            /* Equal dates: no time point required */
            if (Intval == 0)
                continue;
            
            NbStep = 0;
                
            do 
            {
                if (Step < MaxStep)
                    Step++;
                
                NbStep++;
                Intval -= Step;
                        
            }  while (Intval > 0);
        
            NbTP += NbStep;      

        }/*for i*/

        if (NbTP == 0)
        {
            DR_Error ("Fix3_Time_Line: no time points in the tree!");
            goto FREE_MEM_AND_RETURN;
        }

 
        tree_data->NbTP = NbTP;
       
        Fix3_Tree_Alloc (tree_data);
        
        /* 
        *   Place time points.
        */

        tree_data->TPDate[0] = Date[0].CritDate;
        
        k = 0;                                  
        LastStep = 0;
        
        for (i = 0; i < NbDate - 1; i++)                                
        {
            /* Record type of time point of current critical date */
            /* We have to do it before the "if" test below.       */
            for (j = 0; j < NBCRITDATE; j++)
            {
                if (Date[i].Type == j)
                {
                    tree_data->TPtype[j][k] = TRUE;       
                                                                                
                    tree_data->CritDate[j][k].CritDate    = Date[i].CritDate;
                    tree_data->CritDate[j][k].Type        = Date[i].Type;
                    tree_data->CritDate[j][k].Value[0]    = Date[i].Value[0];
                    tree_data->CritDate[j][k].Value[1]    = Date[i].Value[1];
                    tree_data->CritDate[j][k].Value[2]    = Date[i].Value[2];
                    tree_data->CritDate[j][k].Value[3]    = Date[i].Value[3];
                    tree_data->CritDate[j][k].Value[4]    = Date[i].Value[4];
                    tree_data->CritDate[j][k].SuppDate[0] = Date[i].SuppDate[0];
                    tree_data->CritDate[j][k].SuppDate[1] = Date[i].SuppDate[1];
                    tree_data->CritDate[j][k].SuppDate[2] = Date[i].SuppDate[2];
                }
            }  /* for j */
             
            d = (int)Daysact(ValueDate, Date[i].CritDate);

            /* If the number of daily time points has not been reached place 
               daily time points, otherwise increasing time steps */
            if((i0 == 0) && (i == 0) && 
                ((int)Daysact(ValueDate, Date[1].CritDate) > NbDailyPts) )
            {
                for (j = 0; j < NbDailyPts; j++)
                {
                    k++;
                    tree_data->TPDate[k] = Nxtday(tree_data->TPDate[k - 1], 1);
                }
               
                Intval = (int) Daysact (ValueDate, Date[i + 1].CritDate) 
                         - NbDailyPts;
                UnderNDP = FALSE; 
            }
            else if((i0 == 0) && 
                    (int)Daysact(ValueDate, Date[1].CritDate) <= NbDailyPts)
            {
                Intval = (int) Daysact(Date[i].CritDate, Date[i + 1].CritDate);
                UnderNDP = TRUE;
            }
            else if ((i0 > 0) && (d <= NbDailyPts))
            {
                Intval = 0;
                UnderNDP = FALSE;
                if (i < i0 )
                {
                    Intval = (int)Daysact(Date[i].CritDate,Date[i + 1].CritDate);
                    UnderNDP = TRUE;
                }
                else if(i == i0)
                {
                    Intval = 0;
                    UnderNDP = FALSE;
 
                    for ( j =  0; j  < (NbDailyPts - d); j++)
                    {
                        m = k + j + 1;
                        tree_data->TPDate[m] = Nxtday(tree_data->TPDate[m - 1],1);
                    }
               
                    k = k  + NbDailyPts  - d ;
                    Intval = (int) Daysact(ValueDate,Date[i + 1].CritDate)
                              -  NbDailyPts;   
                }   
            }
            else
            {
                Intval = (int) Daysact (Date[i].CritDate, Date[i+1].CritDate);
                UnderNDP = FALSE;
            }
        
                
            if (Intval == 0)
                continue;

            
            NbStep = 0;
        
            if (UnderNDP)
                Step = 0;
            else
                Step = LastStep;


            if(UnderNDP)
            {
                do 
                {
                   Step = 1;
                
                   NbStep++;
                   Intval -= Step;
                        
                }  while (Intval > 0);
           
            }
            else
            {
                do 
                {
                    if (Step < MaxStep)
                       Step++;
                
                   NbStep++;
                   Intval -= Step;
                        
                }  while (Intval > 0);
            }
                
            Rest  = -Intval;
            Rest0 = Rest / NbStep;
            Rest -= Rest0 * NbStep;
                        
            
            Step = LastStep;
            
            for (j = 0; j < NbStep ; j++)
            {
                if(UnderNDP)
                    Step = 1;
                else
                {
                    if(Step < MaxStep )
                       Step++;
                }

                Step0 = Step  - Rest0 - (j >= NbStep  - Rest);

            
                k++;
            
                tree_data->TPDate[k] = Nxtday (tree_data->TPDate[k-1], 
                                                            (long) Step0);
            }
                        
            if (tree_data->TPDate[k] != Date[i+1].CritDate)
            {
                DR_Error ("Fix3_Time_Line: time point does not fall on critical"
                            "date!");
                goto FREE_MEM_AND_RETURN;
            }

            if (UnderNDP)
                LastStep = 0;
            else
                LastStep = Step;

        }  /* for i */


        /* Record type of last date */
        for (j = 0; j < NBCRITDATE; j++)
        {                           
            if (Date[NbDate-1].Type == j)
            {
                tree_data->TPtype[j][k] = TRUE;       
                                                                                
                tree_data->CritDate[j][k].CritDate    = Date[NbDate-1].CritDate;
                tree_data->CritDate[j][k].Type        = Date[NbDate-1].Type;
                tree_data->CritDate[j][k].Value[0]    = Date[NbDate-1].Value[0];
                tree_data->CritDate[j][k].Value[1]    = Date[NbDate-1].Value[1];
                tree_data->CritDate[j][k].Value[2]    = Date[NbDate-1].Value[2];
                tree_data->CritDate[j][k].Value[3]    = Date[NbDate-1].Value[3];
                tree_data->CritDate[j][k].Value[4]    = Date[NbDate-1].Value[4];
                tree_data->CritDate[j][k].SuppDate[0] = Date[NbDate-1].SuppDate[0];
                tree_data->CritDate[j][k].SuppDate[1] = Date[NbDate-1].SuppDate[1];
                tree_data->CritDate[j][k].SuppDate[2] = Date[NbDate-1].SuppDate[2];
            }
        }  /* for j */
    }
    /* Equal time steps */
    else if (StepStyle == 'E')
    {
        NbTP = 0;

        for (i = 0; i < NbDate - 1; i++)
        {
            x=tree_data->Ppy*Daysact(Date[i].CritDate,Date[i+1].CritDate)/372.;
            TPinIntval[i] = (int) ceil (x);

            NbTP += TPinIntval[i];
        }
                         
        
        if (NbTP == 0)
        {
            DR_Error ("Fix3_Time_Line: no time points in the tree!");
            goto FREE_MEM_AND_RETURN;       
        }

        tree_data->NbTP = NbTP;
        Fix3_Tree_Alloc (tree_data);


        for (i = 0, k = 0; i < NbDate - 1; i++)
        {
            tree_data->TPDate[k] = Date[i].CritDate;

            for (j = 0; j < NBCRITDATE; j++)
            {
                if (Date[i].Type == j)
                {
                    tree_data->TPtype[j][k] = TRUE;       
                                                                                
                    tree_data->CritDate[j][k].CritDate    = Date[i].CritDate;
                    tree_data->CritDate[j][k].Type        = Date[i].Type;
                    tree_data->CritDate[j][k].Value[0]    = Date[i].Value[0];
                    tree_data->CritDate[j][k].Value[1]    = Date[i].Value[1];
                    tree_data->CritDate[j][k].Value[2]    = Date[i].Value[2];
                    tree_data->CritDate[j][k].Value[3]    = Date[i].Value[3];
                    tree_data->CritDate[j][k].Value[4]    = Date[i].Value[4];
                    tree_data->CritDate[j][k].SuppDate[0] = Date[i].SuppDate[0];
                    tree_data->CritDate[j][k].SuppDate[1] = Date[i].SuppDate[1];
                    tree_data->CritDate[j][k].SuppDate[2] = Date[i].SuppDate[2];
                }
            }  /* for j */
                
            if (!TPinIntval[i])
                continue;

            k++;

            /* Approximate length of interval */
            x = Daysact (Date[i].CritDate, Date[i+1].CritDate) 
                    / (double) TPinIntval[i];

            for (j = 1; j < TPinIntval[i]; j++)
            {
                tree_data->TPDate[k] = Nxtday(Date[i].CritDate,(long)(j*x+.5));

                k++;
            }
        }  /* for i */

        
        tree_data->TPDate[NbTP] = Date[NbDate - 1].CritDate;
        
        for (j = 0; j < NBCRITDATE; j++)
        {                           
            if (Date[NbDate-1].Type == j)
            {
                tree_data->TPtype[j][NbTP] = TRUE;       
                                                                                
                tree_data->CritDate[j][NbTP].CritDate    = Date[NbDate-1].CritDate;
                tree_data->CritDate[j][NbTP].Type        = Date[NbDate-1].Type;
                tree_data->CritDate[j][NbTP].Value[0]    = Date[NbDate-1].Value[0];
                tree_data->CritDate[j][NbTP].Value[1]    = Date[NbDate-1].Value[1];
                tree_data->CritDate[j][NbTP].Value[2]    = Date[NbDate-1].Value[2];
                tree_data->CritDate[j][NbTP].Value[3]    = Date[NbDate-1].Value[3];
                tree_data->CritDate[j][NbTP].Value[4]    = Date[NbDate-1].Value[4];
                tree_data->CritDate[j][NbTP].SuppDate[0] = Date[NbDate-1].SuppDate[0];
                tree_data->CritDate[j][NbTP].SuppDate[1] = Date[NbDate-1].SuppDate[1];
                tree_data->CritDate[j][NbTP].SuppDate[2] = Date[NbDate-1].SuppDate[2];

            }  /* if */
        }  /* for j */
    }  /* if then else */
   
       
        
    /* We did not count the number of time points correctly */
    if (k != NbTP)
    {
        DR_Error ("Fix3_Time_Line: problem in the construction of the time line: "
                    "k != NbTP");
        goto FREE_MEM_AND_RETURN;                           
    }


    /* 
    *   Calculate interval between time steps.
    */

    NbDaysSoFar = 0;

    for (k = 0; k < NbTP; k++)
    {
        int DeltaDays;

        DeltaDays = Daysact (tree_data->TPDate[k], tree_data->TPDate[k+1]);
        tree_data->Length[k] = DeltaDays / 365.;

        if (tree_data->Length[k] < ERROR)
        {
            DR_Error ("Fix3_Time_Line: problem in the construction of the time "
                        "line: Length = 0");
            goto FREE_MEM_AND_RETURN;      
        }

        /* If  Length is OK, then assign delta t for jumpsize. If on */
        /* 'I' it is the step size corresponding to the current time */
        /* and if 'E', then it is a constant.                        */
        if (StepStyle == 'I' )
        {
            if((k < NbDailyPts) && (NbDailyPts > 0))
                tree_data->LengthJ[k] = 1.0 / 365.;
            else
                tree_data->LengthJ[k] = MIN(0.5*(1.0+sqrt(1.0+8.0*NbDaysSoFar)),
                                        MaxStep)/365.;
        }
        else if (StepStyle == 'J')
        {
            int    MaxStepJ;
            double LengthJ;

            MaxStepJ = (int) floor (372. / tree_data->JumpPpy);

            LengthJ  = MIN(0.5*(1.0+sqrt(1.0+8.0*NbDaysSoFar)), MaxStepJ)/365.;

            if (( k < NbDailyPts) && (NbDailyPts > 0))
                tree_data->LengthJ[k] = 1.0/ 365.;
            else
                tree_data->LengthJ[k] = MAX(tree_data->Length[k], LengthJ);

        }
        else
        {
            tree_data->LengthJ[k] = 1.0/tree_data->Ppy;
        }
        NbDaysSoFar += DeltaDays;

    }  /* for k */

    /* We need 2 extra time points with arbitrary length for the calibration */
    tree_data->Length[NbTP] = tree_data->Length[NbTP - 1];
    tree_data->Length[-1]   = tree_data->Length[0];                  

    tree_data->LengthJ[NbTP] = tree_data->LengthJ[NbTP - 1];
    tree_data->LengthJ[-1]   = tree_data->LengthJ[0];

    tree_data->TPDate[NbTP+1] = Nxtday (tree_data->TPDate[NbTP], 
                                 (long) (tree_data->Length[NbTP] * 365. + .5));

    /* 
    *   Interpolate values between dates for American type critical dates.
    */

    for (i = 0; i < NBCRITDATE; i++)                    
    {                                                                       
        switch (tree_data->CritType[i])
        {
            /* Staircase profile */
            case 'S':
            {
                k1 = 0;
                while ((k1 <= NbTP) && (tree_data->TPtype[i][k1] != TRUE))
                k1++;

                if (k1 == NbTP + 1)
                    break;


                for (;;)
                {
                    k2 = k1 + 1;

                    while ((k2 <= NbTP) && (tree_data->TPtype[i][k2] != TRUE))
                    k2++;

                    if (k2 == NbTP + 1)
                        break;

                    for (k = k1+1; k < k2; k++)
                    {
                        tree_data->TPtype[i][k] = TRUE;

                        tree_data->CritDate[i][k].Value[0] = tree_data->CritDate[i][k1].Value[0];
                        tree_data->CritDate[i][k].Value[1] = tree_data->CritDate[i][k1].Value[1];
                        tree_data->CritDate[i][k].Value[2] = tree_data->CritDate[i][k1].Value[2];
                        tree_data->CritDate[i][k].Value[3] = tree_data->CritDate[i][k1].Value[3];
                        tree_data->CritDate[i][k].Value[4] = tree_data->CritDate[i][k1].Value[4];
                    }

                    k1 = k2;
        
                }  /* for */

                break;
            }        
            /* Linear interpolation */
            case 'L':
            {
                k1 = 0;
                while ((k1 <= NbTP) && (tree_data->TPtype[i][k1] != TRUE))
                k1++;

                if (k1 == NbTP + 1)
                    break;


                for (;;)
                {
                    k2 = k1 + 1;

                    while ((k2 <= NbTP) && (tree_data->TPtype[i][k2] != TRUE))  
                    k2++;

                    if (k2 == NbTP + 1)
                        break;

                    for (k = k1+1; k < k2; k++)
                    {
                        tree_data->TPtype[i][k] = TRUE;

                        for (l = 0; l < 5; l++)
                        {
                            dlinterp (  tree_data->TPDate[k], 
                                        &x,
                                        tree_data->TPDate[k1],
                                        tree_data->TPDate[k2],
                                        tree_data->CritDate[i][k1].Value[l],
                                        tree_data->CritDate[i][k2].Value[l]);

                            tree_data->CritDate[i][k].Value[l] = x;
                        }
                    }  /* for k */

                    k1 = k2;
        
                }  /* for */
        
                break;
            }       
            /* Discrete set: no need to interpolate */
            default:
            {
                break;
                                
            }        
        }  /* switch */
    }  /* for i */


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (TPinIntval, INT, 0, NbDate+1);

    return (status);

}  /* Fix3_Time_Line */


/*****  Fix3_Time_Line_Nmr  ******************************************************/
/*
*       Set up the time line including numeraire dates if necessary
*/   
int     Fix3_Time_Line_Nmr (   
                long        ValueDate,       /* (I) Value date                */
                int         *NbCritDate,     /* (I) Number of dates           */
                CRIT_DATE   **CritDate,      /* (I) Critical date list        */
                char        StepStyle,       /* (I) Equal or increasing steps */
                MKTVOL_DATA *mktvol_data,    /* (I) Volatility data           */
                FIX3_TREE_DATA   *tree_data) /* (O) Tree data                 */
{
    long LastProdDate = 0;
    int  status = FAILURE;
    int  i;

    /* Insert numeraire dates if dealing with numeraire model */
    if (mktvol_data->IsNmrModel)
    {
        for (i = 0; i < *NbCritDate; i++)
        {
            LastProdDate = MAX (LastProdDate, (*CritDate)[i].CritDate);
        }

        if (Nmr_Schedule (ValueDate,
                          LastProdDate,
                          mktvol_data) == FAILURE) goto RETURN;

        for (i = 0; i < mktvol_data->NbNmr; i++)
        {
            if (Add_To_DateList (NbCritDate,
                                 CritDate,
                                 mktvol_data->NmrDate[i],
                                 NMREVENT,
                                 i, 0, 0, 0, 0,
                                 0, 0, 0) == FAILURE) goto RETURN;
        }
    }

    if (Fix3_Time_Line (ValueDate,
                        *NbCritDate,
                        *CritDate,
                        StepStyle,
                        tree_data) == FAILURE) goto RETURN;


    status = SUCCESS;

RETURN:

    return (status);

}  /* Fix3_Time_Line_Nmr */


