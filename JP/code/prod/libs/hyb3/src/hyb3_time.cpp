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
#include "cupslib.h"

/*****  Hyb3_Time_Line  **********************************************************/
/*
*       Set up the time step array so that a node falls on each exercise
*       date, coupon date, benchmark volatility date, ... Set up the cash-
*       flow, exercise arrays, ... Allocate memory for the tree arrays.
*       There are 15 time arrays, of which: TimeArray[14]: Zero coupon resets.
*       Each set of dates can be discrete (e.g.:coupon payments) or continuous
*       (e.g. Knock-out dates). For the latter we need to interpolate values
*       between input dates.
*/
int     Hyb3_Time_Line (long      ValueDate,   /* (I) Value date             	             */
                   int       NbDate,	  /* (I) Number of dates in critical date list       */
                   CRIT_DATE *Date,	  /* (I) Critical date list	                     */
                   char      EoI,         /* (I) Equal or increasing time steps 	     */
                   HYB3_TREE_DATA *tree_data)  /* (O) Structure of tree data 	             */
{	                
    double
            x;
    int
            NbTP,                     /* Total number of nodes in the tree                   */
            i, j, k, k1, k2, l,       /*                                                     */
            Period,                   /* Length of current period between two critical dates */
            NbStep,                   /* Number of steps in current period                   */
            MaxStep=0,                /* Maximum length for a time step                      */
            Step,                     /* Length of current time step                         */
            Step0, 	              /* Length of current time step after adjustment        */
            LastStep,	              /* Length of last step of previous period              */
            Rest, Rest0,              /* Days left over in current period                    */
            status = FAILURE;	      /* Error status = FAILURE initially                    */
    int                               /*                                                     */
            NbDaysSoFar,              /*                                                     */
           *NodePerPeriod = NULL;     /* Number of Nodes between two critical dates          */

    NodePerPeriod = (int *)DR_Array(INT, 0, NbDate + 1);
    if (NodePerPeriod == NULL)
    {
        goto RETURN;
    }

    /* 
    *	Sort critical date list.
    */

    if (Sort_CritDate (	NbDate,
                        Date) == FAILURE)
    {
        goto RETURN;
       
    }  /* if */                                                
        

    if (Date[0].CritDate < ValueDate)	                                        /* Just in case */
    {
        DR_Error ("First date falls before value date! (Hyb3_Time_Line)");
        goto RETURN;	                        
       
    }  /* if */                                                


    /* 
    *	Build time line i.e. put nodes between critical dates.
    */

    if (EoI == 'I')                                                             /* Increasing time steps */
    {
        MaxStep = (int) floor (372. / tree_data->Ppy);	                        /* Maximum time step in days: linked to Ppy */
                    
        NbTP = 0;
        Step    = 0;

        /* 
        *	Calculate total number of time steps.
        */

        for (i = 0; i < NbDate - 1; i++)                            
        {
            Period = (int) Daysact (Date[i].CritDate, Date[i+1].CritDate);      /* Length of current period in days */
                
            if (Period == 0)                                                    /* If two consecutive dates are equal there are no nodes in the period */
                continue;
            
            NbStep = 0;
                
            do 
            {
                if (Step < MaxStep)	                                            /* Increase time step by one day if less than maximum time step */
                    Step++;
                
                NbStep++;	                                                    /* Add one step */
                Period -= Step;	                                                /* Subtract its length from the current period */
                        
            }  while (Period > 0);                	                            /* Until all current period is covered */
                
            NbTP += NbStep;	                                                /* Add the number of steps in current period to the total */
                                        
        }  /* for i */
                         
        
        if (NbTP == 0)	                                                    /* Just in case */
        {
            DR_Error ("No nodes in the tree! (Hyb3_Time_Line)");
            goto RETURN;	                        
       
        }  /* if */                                                


        tree_data->NbTP = NbTP;                                           /* Allocate memory for the arrays constituting the tree. We could not do it */
        Hyb3_Tree_Alloc (tree_data);                                                 /* before because we did not know the total number of nodes NbTP.        */

        
        /* 
        *	Place the nodes.
        */

        tree_data->TPDate[0] = Date[0].CritDate;                              /* First node is first date, i.e. value date */
        
        k = 0;	                                
        LastStep = 0;
        
        for (i = 0; i < NbDate - 1; i++)                                
        {
            for (j = 0; j < NBCRITDATE; j++)                                    /* Record the type of the node falling on the current critical date  */
            {	                                                                /* We have to do it before the "if" test below.	             */
                if (Date[i].Type == j)
                {
                    tree_data->TPType[j][k] = TRUE;       
                                                                                
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

                }  /* if */
            }  /* for j */
                
            Period = (int) Daysact (Date[i].CritDate, Date[i+1].CritDate);
                
            if (Period == 0)                                                    /* If two consecutive dates are equal there are no nodes in the period */
                continue;

            
            NbStep = 0;
            Step = LastStep;
                
            do 
            {
                if (Step < MaxStep)
                    Step++;
                
                NbStep++;
                Period -= Step;
                        
            }  while (Period > 0);                
                
            Rest  = -Period;
            Rest0 = Rest / NbStep;
            Rest -= Rest0 * NbStep;
                        
                
            Step = LastStep;
                        
            for (j = 0; j < NbStep; j++)
            {
                if (Step < MaxStep)
                    Step++;

                Step0 = Step - Rest0 - (j >= NbStep - Rest);                    /* Step0 = Step - Rest0 -1 while Rest > 0, Step - Rest0 otherwise */

                k++;
                tree_data->TPDate[k] = Nxtday (tree_data->TPDate[k-1], (long) Step0);   /* Nodes are falling on exact date */
                        
            }  /* for j */
                        
            if (tree_data->TPDate[k] != Date[i+1].CritDate)                   /* Last node should fall on end of current period */
            {
                DR_Error ("Node does not fall on critical date! (Hyb3_Time_Line)");
                goto RETURN;	                        
       
            }  /* if */                                                
                    
            LastStep = Step;                                                    /* Size of last time step of current period */
        
        }  /* for i */


        for (j = 0; j < NBCRITDATE; j++)                                        /* Record the type of the last date */
        {	                        
            if (Date[NbDate-1].Type == j)
            {
                tree_data->TPType[j][k] = TRUE;       
                                                                                
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

            }  /* if */
        }  /* for j */
    }
    else                                                                        /* Equal time steps */
    {
        for (i = 0, NbTP = 0; i < NbDate - 1; i++)                                   /* Calculate total number of periods and    */
        {	                                                                            /* number of periods between critical dates */
            x = tree_data->Ppy * Daysact (Date[i].CritDate, Date[i+1].CritDate) / 372.;	/* The number of periods between two critical dates  */
            NodePerPeriod[i] = (int) ceil (x);	                                        /* is the smallest  integer that gives a period      */
                                                                                        /* length smaller than 1/Ppy (in years). There are   */
            NbTP += NodePerPeriod[i];                                                /* 372 days in a year (1/12 year gives 31 day month) */

        }  /* for i */
                         
        
        if (NbTP == 0)	                                                    /* Just in case */
        {
            DR_Error ("No nodes in the tree! (Hyb3_Time_Line)");
            goto RETURN;	                        
       
        }  /* if */                                                

        tree_data->NbTP = NbTP;                                           /* Allocate memory for the arrays constituting the tree. We could not do it */
        Hyb3_Tree_Alloc (tree_data);                                                 /* before because we did not know the total number of nodes NbTP.        */


        for (i = 0, k = 0; i < NbDate - 1; i++)                                 /* Calculate length of each time step (period), stored in tree_data->Length */
        {
            tree_data->TPDate[k] = Date[i].CritDate;

            for (j = 0; j < NBCRITDATE; j++)                                    /* Record the type of the node falling on the current critical date  */
            {	                                                                /* We have to do it before the "if" test below.	             */
                if (Date[i].Type == j)
                {
                    tree_data->TPType[j][k] = TRUE;       
                                                                                
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

                }  /* if */
            }  /* for j */
                
            if (!NodePerPeriod[i])                                              /* If two consecutive dates are equal, NodePerPeriod[i] = 0. We ignore it */
                continue;

            k++;

            x = Daysact (Date[i].CritDate, Date[i+1].CritDate) / (double) NodePerPeriod[i];	/* Approximate length of periods between [i] and [i+1] */

            for (j = 1; j < NodePerPeriod[i]; j++)
            {
                tree_data->TPDate[k] = Nxtday (Date[i].CritDate, (long) (j * x + .5));    /* Nodes are falling on exact date */

                k++;

            }  /* for j */
        }  /* for i */

        
        tree_data->TPDate[NbTP]  = Date[NbDate - 1].CritDate;              /* This is the last node falling on the last date */
        
        for (j = 0; j < NBCRITDATE; j++)                                        /* Record the type of the last date */
        {	                        
            if (Date[NbDate-1].Type == j)
            {
                tree_data->TPType[j][NbTP] = TRUE;       
                                                                                
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
        
        
    if (k != NbTP)                                                           /* We did not count the number of periods correctly */
    {
        DR_Error ("Problem in the construction of the tree: k != NbTP (Hyb3_Time_Line)");
        goto RETURN;	                        
       
    }  /* if */                                                


    /* 
    *	Calculate time interval between time steps.
    */

    NbDaysSoFar = 0;
    for (k = 0; k < NbTP; k++)
    {
        int DeltaDays;

        DeltaDays = Daysact (tree_data->TPDate[k], tree_data->TPDate[k+1]);
        tree_data->Length[k] = DeltaDays / 365.;

        if (tree_data->Length[k] < ERROR)                                       /* If length=0, we made a mistake placing the nodes */
        {
            DR_Error ("Problem in the construction of the tree: Length = 0 (Hyb3_Time_Line)");
            goto RETURN;	                        
       
        }  /* if */

        /* If  Length is OK, then assign delta t for jumpsize. If on */
        /* 'I' it is the step size corresponding to the current time */
        /* and if 'E', then it is a constant.                        */
        if (EoI == 'I')
        {
            tree_data->LengthJ[k] = MIN(0.5*(1.0+sqrt(1.0+8.0*NbDaysSoFar)),
                                         MaxStep)/365.;
        }
        else
        {
            tree_data->LengthJ[k] = 1.0/tree_data->Ppy;
        }
        NbDaysSoFar += DeltaDays;

    }  /* for k */

    tree_data->Length[NbTP] = tree_data->Length[NbTP - 1];                /* We need two extra nodes for the calibration: their length is arbitrary */
    tree_data->Length[-1]      = tree_data->Length[0];            	    

    tree_data->LengthJ[NbTP] = tree_data->LengthJ[NbTP - 1];
    tree_data->LengthJ[-1]      = tree_data->LengthJ[0];

    tree_data->TPDate[NbTP+1] = Nxtday (tree_data->TPDate[NbTP], 
                                     (long) (tree_data->Length[NbTP] * 365. + .5));

    /* 
    *	Interpolate values between dates for American type critical dates.
    */

    for (i = 0; i < NBCRITDATE; i++)                    
    {                                                                       
        switch (tree_data->CritType[i])
        {
            case 'S':	                                                        /* Staircase profile */
            {
                k1 = 0;
                while ((k1 <= NbTP) && (tree_data->TPType[i][k1] != TRUE))	/* Search for first critical date of type [i] */
                k1++;

                if (k1 == NbTP + 1)	                                        /* No critical date of this type */
                    break;


                for (;;)
                {
                    k2 = k1 + 1;

                    while ((k2 <= NbTP) && (tree_data->TPType[i][k2] != TRUE))	/* Search for following critical date */
                    k2++;

                    if (k2 == NbTP + 1)
                        break;

                    for (k = k1+1; k < k2; k++)
                    {
                        tree_data->TPType[i][k] = TRUE;

                        tree_data->CritDate[i][k].Value[0] = tree_data->CritDate[i][k1].Value[0];
                        tree_data->CritDate[i][k].Value[1] = tree_data->CritDate[i][k1].Value[1];
                        tree_data->CritDate[i][k].Value[2] = tree_data->CritDate[i][k1].Value[2];
                        tree_data->CritDate[i][k].Value[3] = tree_data->CritDate[i][k1].Value[3];
                        tree_data->CritDate[i][k].Value[4] = tree_data->CritDate[i][k1].Value[4];

                    }  /* for k */

                    k1 = k2;
        
                }  /* for */

                break;
            }        
            case 'L':	                                                        /* Linear interpolation */
            {
                k1 = 0;
                while ((k1 <= NbTP) && (tree_data->TPType[i][k1] != TRUE))
                k1++;

                if (k1 == NbTP + 1)
                    break;


                for (;;)
                {
                    k2 = k1 + 1;

                    while ((k2 <= NbTP) && (tree_data->TPType[i][k2] != TRUE))	
                    k2++;

                    if (k2 == NbTP + 1)
                        break;

                    for (k = k1+1; k < k2; k++)
                    {
                        tree_data->TPType[i][k] = TRUE;

                        for (l = 0; l < 5; l++)
                        {
                            dlinterp (	tree_data->TPDate[k], 
                                        &x,
                                        tree_data->TPDate[k1],
                                        tree_data->TPDate[k2],
                                        tree_data->CritDate[i][k1].Value[l],
                                        tree_data->CritDate[i][k2].Value[l]);

                            tree_data->CritDate[i][k].Value[l] = x;

                        }  /* for l */
                    }  /* for k */

                    k1 = k2;
        
                }  /* for */
        
                break;
            }        
            default:     /* If the set of dates is discrete we don't have to interpolate values */
            {
                break;
                                
            }        
        }  /* switch */
    }  /* for i */


    status = SUCCESS;

    RETURN:   /* Passes control here on error if we have to free memory */

    Free_DR_Array(NodePerPeriod, INT, 0, NbDate + 1);
    
    return (status);

}  /* Hyb3_Time_Line */


