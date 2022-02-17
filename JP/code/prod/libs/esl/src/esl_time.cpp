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
#include "esl_time.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"
#include "esl_util.h"


/*****  Add_To_DateList  ******************************************************/
/**
*       Add a set of dates to the critical date list.
*/
int     Add_To_DateList (int        *NbCritDate  /** (I/O) Number of dates    */
                        ,CRIT_DATE  **CritDate 	 /** (I/O) Critical date list */
                        ,long       Date         /** (I) New date             */
                        ,int        Type         /** (I) New date type        */
                        ,double	    Value0       /** (I) New date values      */
                        ,double	    Value1       /** (I)                      */
                        ,double	    Value2       /** (I)                      */
                        ,double	    Value3       /** (I)                      */
                        ,double	    Value4       /** (I)                      */
                        ,long       SuppDate0    /** (I) Supplementary dates  */
                        ,long       SuppDate1    
                        ,long       SuppDate2)   
{

    int     status = FAILURE;   /* Error status = FAILURE initially */


    /* 
    *	Re-allocate memory.
    */

    *NbCritDate += 1;

    *CritDate =(CRIT_DATE *)DR_REALLOC(*CritDate,*NbCritDate*sizeof(CRIT_DATE));

    if (!(*CritDate))
    {
        DR_Error ("Add_To_DateList: allocation failure!\n");
        goto RETURN;
    }


    /* 
    *	Add the new date.
    */
	{
		CRIT_DATE *c = &((*CritDate)[*NbCritDate-1]);
	    c->CritDate    = Date;
	    c->Type        = Type;          
	    c->Value[0]    = Value0;
	    c->Value[1]    = Value1;
	    c->Value[2]    = Value2;
	    c->Value[3]    = Value3;
	    c->Value[4]    = Value4;
	    c->SuppDate[0] = SuppDate0;
	    c->SuppDate[1] = SuppDate1;
	    c->SuppDate[2] = SuppDate2;
	}

    status = SUCCESS;

    RETURN:

    return (status);

}  /* Add_To_DateList */


/*****  TimeInterp  *********************************************************/
/**
*       Utility routine. Eliminates dates falling before the value date and 
*       interpolate the corresponding value.
*       It assumes dates are entered in ascending order.
*/
int     TimeInterp (long    ValueDate  /** (I) Value date                     */
                   ,char    *Name      /** (I) Array name                     */
                   ,char    Type       /** (I) Type of interpolation          */
                   ,int     *NbDate    /** (I/O) Number of dates in the array */
                   ,long    *Date      /** (I/O) Dates in ascending order     */
                   ,double  *Curve0    /** (I/O) Set of associated curves     */
                   ,double  *Curve1 
                   ,double  *Curve2 
                   ,double  *Curve3 
                   ,double  *Curve4)
{

    double      *Curve[5];
    double      y;
    int         i, j, k;
    int         NbCurves;
    int         status = FAILURE;        /* Error status = FAILURE initially */
    char        Message[MAXBUFF];
                

    /* Nothing to do */
    if (*NbDate == 0)
    {
        return (SUCCESS);
    }
                
        
    if ((Type == 'L') && (*NbDate < 2))
    {
        DR_Error ("TimeInterp: only one date to do a linear interpolation!");
        goto RETURN;
    }


    /*
     *  Determine number of curves to interpolate.
     */

    if (Curve0 == NULL)
    {
        NbCurves = 0;
    }
    else
    {
        Curve[0] = Curve0;

        if (Curve1 == NULL)
        {
            NbCurves = 1;
        }
        else
        {
            Curve[1] = Curve1;

            if (Curve2 == NULL)
            {
                NbCurves = 2;
            }
            else
            {
                Curve[2] = Curve2;

                if (Curve3 == NULL)
                {
                    NbCurves = 3;
                }
                else
                {
                    Curve[3] = Curve3;

                    if (Curve4 == NULL)
                    {
                        NbCurves = 4;
                    }
                    else
                    {
                        Curve[4] = Curve4;

                        NbCurves = 5;
                    }
                }
            }
        }
    }

    
    /* 
     *  Find first date falling after maturity.
     */

    k = 0;
    while ((ValueDate > Date[k]) && (k < *NbDate - 1))                          
        k++;
        
        
    /* Nothing to do: all dates are falling after value date */
    if ((k == 0) && (ValueDate <= Date[0]))
    {
        return (SUCCESS);
    }
        
        
    /* There is no date left and we don't use stair interpolation */
    if ((ValueDate > Date[k]) && (Type != 'S'))
    {
        *NbDate = 0;
                
        return (SUCCESS);
    }
        

    /*
     *  Eliminate dates falling before value date and interpolate curves.
     */

    switch (Type)
    {
        /* Discrete set of dates */
        case 'D':
        {        
            *NbDate -= k;
        
            for (i = 0; i < *NbDate; i++)
            {
                Date[i] = Date[i+k];        
                                
                for (j = 0; j < NbCurves; j++)
                    Curve[j][i] = Curve[j][i+k];
            }
        
            break;
        }
        /* Staircase profile */
        /* If value date is falling before the kth date */
        /* we add one date equal to the value date      */
        case 'S':       
        {
            if (ValueDate < Date[k])                                            
                k--;                                                            

            *NbDate -= k;
                                
            for (i = 0; i < *NbDate; i++)                               
            {
                Date[i] = Date[i+k];        
                                
                for (j = 0; j < NbCurves; j++)
                    Curve[j][i] = Curve[j][i+k];
            }
                
            Date[0] = ValueDate;                    
                        
            break;
        } 
        /* Linear profile */
        case 'L':
        {
            if (ValueDate < Date[k])
            {        
                k--;
         
                for (j = 0; j < NbCurves; j++)
                {        
                    dlinterp (	ValueDate, &y,
                                Date[k], Date[k+1], 
                                Curve[j][k], Curve[j][k+1]);

                    Curve[j][k] = y;
                }
            }  /* if */
                        
            *NbDate -= k;
                                
            for (i = 0; i < *NbDate; i++)
            {
                Date[i] = Date[i+k];        
                                
                for (j = 0; j < NbCurves; j++)
                    Curve[j][i] = Curve[j][i+k];
            }
                
            Date[0] = ValueDate;
                
            break;
        }
        default:
        {
            sprintf (Message, "TimeInterp: unrecognised interpolation method "
                        "for %s!", Name);
            DR_Error (Message);
            goto RETURN;
        }
    }  /* switch */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* TimeInterp */



/*****  Sort_CritDate  ******************************************************/
/**
*       Sort routine of critical dates. Used in time.c.
*/
int     Sort_CritDate (	int	        n       /** Number of dates to sort */
                       ,CRIT_DATE	*Date   /** Dates to be sorted      */
		)
{

    CRIT_DATE   temp;                    /* Temporary storage */

    int     i, j;


    for (j = 1 ; j < n; j++)
    {
        temp.CritDate    = Date[j].CritDate;
        temp.Type        = Date[j].Type;
        temp.Value[0]    = Date[j].Value[0];
        temp.Value[1]    = Date[j].Value[1];
        temp.Value[2]    = Date[j].Value[2];
        temp.Value[3]    = Date[j].Value[3];
        temp.Value[4]    = Date[j].Value[4];
        temp.SuppDate[0] = Date[j].SuppDate[0];
        temp.SuppDate[1] = Date[j].SuppDate[1];
        temp.SuppDate[2] = Date[j].SuppDate[2];

        i = j-1;
        while (i >= 0 && Date[i].CritDate > temp.CritDate)
        {
            Date[i+1].CritDate    = Date[i].CritDate;
            Date[i+1].Type        = Date[i].Type;
            Date[i+1].Value[0]    = Date[i].Value[0];
            Date[i+1].Value[1]    = Date[i].Value[1];
            Date[i+1].Value[2]    = Date[i].Value[2];
            Date[i+1].Value[3]    = Date[i].Value[3];
            Date[i+1].Value[4]    = Date[i].Value[4];
            Date[i+1].SuppDate[0] = Date[i].SuppDate[0];
            Date[i+1].SuppDate[1] = Date[i].SuppDate[1];
            Date[i+1].SuppDate[2] = Date[i].SuppDate[2];

            i--;
        }

        Date[i+1].CritDate    = temp.CritDate;
        Date[i+1].Type        = temp.Type;
        Date[i+1].Value[0]    = temp.Value[0];
        Date[i+1].Value[1]    = temp.Value[1];
        Date[i+1].Value[2]    = temp.Value[2];
        Date[i+1].Value[3]    = temp.Value[3];
        Date[i+1].Value[4]    = temp.Value[4];
        Date[i+1].SuppDate[0] = temp.SuppDate[0];
        Date[i+1].SuppDate[1] = temp.SuppDate[1];
        Date[i+1].SuppDate[2] = temp.SuppDate[2];

    }  /* for j */


    return (SUCCESS);

}  /* Sort_CritDate */


/*****  DrExtendAmerSch  ****************************************************/
/**
 *       Utility routine to extend an input American date schedule to a full
 *       European style schedule using the 'Increasing time-step' method.
 *       It also does the following:
 *       - eliminates dates STRICTLY BEFORE value date
 *       - interpolates supplementary values (if available)
 *
 *       NOTE: the extension is performed on the NotifDates if they are given,
 *             the SchDates are then interpolated. Otherwise, the extension is
 *             done on the SchDates themselves.
 */

int DrExtendAmerSch(long      ValueDate    /** (I)   Value date          */
                   ,int       Ppy          /** (I)   ppy in tree         */
                   ,int      *NbSchDates   /** (I/O) Nb schedule dates   */
                   ,long    **SchDates     /** (I/O) Schedule dates      */
                   ,long    **NotifDates   /** (I/O) Notif dates         */
                   ,char      InterpStyle  /** (I)   Curve interp style  */
                   ,double  **Curve0       /** (I/O) Curve 0 values      */
                   ,double  **Curve1       /** (I/O) Curve 1 values      */
                   ,double  **Curve2       /** (I/O) Curve 2 values      */
                   ,double  **Curve3       /** (I/O) Curve 3 values      */
                   ,double  **Curve4       /** (I/O) Curve 4 values      */
		)
{
    int       status = FAILURE;

    int       i, j, k;
    long      newDate=0L;
    double    daysOfs;

    int       hasNotifDates; /* TRUE if there're notif dates for extension */
    int       isFwdStart;    /* TRUE if schedule is fwd starting           */
    int       Intval;        /* Interval between two sch dates             */
    int       NbStep;        /* Nb steps in current interval               */
    int       MaxStep=0;     /* Maximum length for a time step             */
    int       Step;          /* Length of current time step                */
    int       Step0;         /* Length of current time step after adjust   */
    int       LastStep;      /* Length of last step of previous interval   */
    int       Rest, Rest0;   /* Days left over in current interval         */
    int       NbTP;          /* Total nb of extended time points           */

    /* internal arrays */
    int       NbTmpDL      = 0;
    long     *TmpDL        = NULL;
    double   *NotifDaysOfs = NULL;
    double   *Curve[5]     = {NULL,NULL,NULL,NULL,NULL};
    long     *BaseDL;               /* the datelist we choose to extend    */
    long     *ExtendedDL   = NULL;  /* the extended datelist               */

    /* extended arrays for output */
    double   *ExtCurve[5] = {NULL,NULL,NULL,NULL,NULL};
    long     *ExtSchDL    = NULL;
    long     *ExtNotifDL  = NULL;
    int       NbExtSchDL  = 0;      /* size of all extended arrays         */

    /* basic checks */

    if ((NbSchDates  == NULL) ||
        (SchDates    == NULL)) goto FREE_MEM_AND_RETURN;
    if (*NbSchDates <= 0)
        return(SUCCESS);
    else
        if (*SchDates == NULL) goto FREE_MEM_AND_RETURN;

    if ((InterpStyle != 'L') && (InterpStyle != 'S'))
    {
        DR_Error("DrExtendAmerSch: InterpStyle is not 'L' or 'S'!");
        goto FREE_MEM_AND_RETURN;
    }
    if (*NbSchDates < 2)
    {
        DR_Error("DrExtendAmerSch: "
                 "cannot extend schedule with < 2 dates!");
        goto FREE_MEM_AND_RETURN;
    }
    for (k=1; k<*NbSchDates; k++)
    {
        if ((*SchDates)[k] <= (*SchDates)[k-1])
        {
            DR_Error("DrExtendAmerSch: "
                     "input schedule dates are not strictly increasing.");
            goto FREE_MEM_AND_RETURN;
        }
    }
    if ((*SchDates)[*NbSchDates-1] < ValueDate)
    {
        DR_Error("DrExtendAmerSch: all input dates are < Value Date");
        goto FREE_MEM_AND_RETURN;
    }

    hasNotifDates = ((NotifDates != NULL) && ((*NotifDates) != NULL));
    if (hasNotifDates)
    {
        for (k=0; k<*NbSchDates; k++)
        {
            if ((*NotifDates)[k] > (*SchDates)[k])
            {
                DR_Error("DrExtendAmerSch: "
                         "some notif dates are > schedule dates");
                goto FREE_MEM_AND_RETURN;
            }
        }
    }

    /* Prepare the curves */

    Curve[0] = (Curve0 != NULL) ? (*Curve0) : NULL;
    Curve[1] = (Curve1 != NULL) ? (*Curve1) : NULL;
    Curve[2] = (Curve2 != NULL) ? (*Curve2) : NULL;
    Curve[3] = (Curve3 != NULL) ? (*Curve3) : NULL;
    Curve[4] = (Curve4 != NULL) ? (*Curve4) : NULL;

    /********************************/
    /* extend the American schedule */
    /********************************/

    /* select the BaseDL to perform the extension */
    BaseDL = (hasNotifDates ? *NotifDates : *SchDates);

    /* prepare datelist containing value date and SchDates>value date */

    if (AddDateToList(&NbTmpDL, &TmpDL, ValueDate) == FAILURE) 
            goto FREE_MEM_AND_RETURN;

    for (k=0; k<*NbSchDates; k++)
    {
        if (BaseDL[k] > ValueDate)
        {
            if (AddDateToList(&NbTmpDL, &TmpDL, BaseDL[k]) == FAILURE) 
                    goto FREE_MEM_AND_RETURN;
        }
    }

    /* prepare extended datelist */
    isFwdStart = BaseDL[0] > ValueDate;

    if (isFwdStart)
    {
        if (AddDateToList(&NbExtSchDL, &ExtendedDL, BaseDL[0]) == FAILURE) 
                goto FREE_MEM_AND_RETURN;
    }
    else
    {
        if (AddDateToList(&NbExtSchDL, &ExtendedDL, ValueDate) == FAILURE) 
                goto FREE_MEM_AND_RETURN;
    }

    /* Ppy parameter gives maximum time step in days */
    MaxStep = (int)floor (372./Ppy);
    NbTP = 0;
    LastStep = 0;

    for (i=0; i<NbTmpDL-1; i++)
    {
        /* Calculate total number of time points */

        Intval = (int)Daysact(TmpDL[i], TmpDL[i+1]);        
        
        NbStep = 0;
        Step   = LastStep;
        do
        {
            if (Step < MaxStep) Step++;
            NbStep++;
            Intval -= Step;
            
        }while (Intval > 0);

        /* place the timepoints if appropriate */

        if (!isFwdStart || (i>0))
        {
            NbTP += NbStep;

            Rest  = -Intval;
            Rest0 = Rest / NbStep;
            Rest -= Rest0 * NbStep;
            
            Step = LastStep;
            for (j = 0; j < NbStep; j++)
            {
                if (Step < MaxStep) Step++;                
                Step0 = Step - Rest0 - (j >= NbStep - Rest);
                newDate = Nxtday(ExtendedDL[NbExtSchDL-1], (long)Step0);
                if (AddDateToList(&NbExtSchDL, &ExtendedDL, newDate) == FAILURE)
                        goto FREE_MEM_AND_RETURN;
            }

            if (newDate != TmpDL[i+1])
            {
                DR_Error ("DrExtendAmerSch: extension algorithm failed!");
                goto FREE_MEM_AND_RETURN;
            }

        }/* if (!isFwdStart || (i>0)) */

        LastStep = Step;
        
    }/* for i */

    if (NbExtSchDL != NbTP+1)
    {
        DR_Error("DrExtendAmerSch: fail to verify nb of extended dates!");
        goto FREE_MEM_AND_RETURN;
    }


    /**************************************/
    /* calculate the notif/schedule dates */
    /**************************************/

    if (hasNotifDates)
    {
        ExtNotifDL = ExtendedDL;
        ExtendedDL = NULL; /* reset to NULL for safe FREE */

        /* prepare the input days offset list */
        NotifDaysOfs = (double *)DR_Array(DOUBLE, 0, *NbSchDates-1);
        if (NotifDaysOfs == NULL) goto FREE_MEM_AND_RETURN;

        for (k=0; k<*NbSchDates; k++)
        {
            NotifDaysOfs[k] = (double)Daysact((*NotifDates)[k],
                                              (*SchDates)[k]);
        }

        /* interpolate the schedule dates from the extended notif dates */
        ExtSchDL = (long *)DR_Array(LONG, 0, NbExtSchDL-1);
        if (ExtSchDL == NULL) goto FREE_MEM_AND_RETURN;

        for (i=0; i<NbExtSchDL; i++)
        {
            /* find the nearest lower offset */
            k = GetDLOffset(*NbSchDates, *NotifDates, ExtNotifDL[i], CbkLOWER);
            if (k < 0)
            {
                DR_Error("DrExtendAmerSch: "
                         "unable to interpolate notif date!");
                goto FREE_MEM_AND_RETURN;
            }

            /* interpolate the days offset */
            if (ExtNotifDL[i] == (*NotifDates)[k])
            {
                daysOfs = NotifDaysOfs[k];
            }
            else
            {
                double tmpDaysOfs;

                if (k >= *NbSchDates-1)
                {
                    DR_Error("DrExtendAmerSch: "
                             "cannot extrapolate notif date!");
                    goto FREE_MEM_AND_RETURN;
                }
                dlinterp(ExtNotifDL[i],
                    &daysOfs,
                         (*NotifDates)[k], (*NotifDates)[k+1],
                         NotifDaysOfs[k],  NotifDaysOfs[k+1]);

                /* round up daysOfs if necessary */
                tmpDaysOfs = ceil(daysOfs);
                daysOfs = (fabs(tmpDaysOfs - daysOfs - 1.0) < TINY) ? \
                                (tmpDaysOfs - 1.0) :                  \
                                tmpDaysOfs;
            }

            /* calc the new notification date */
            ExtSchDL[i] = Nxtday(ExtNotifDL[i], (long)daysOfs);

            if (i>0)
            {
                if (ExtSchDL[i] <= ExtSchDL[i-1])
                {
                    /* This situation can occur when PPY is large, for example,
                     * 365, and the ExtNotifDL is close to one of the input
                     * notification date (NotifDates).  
                     *
                     *   ----------|------|----------------|--
                     *            Notif  Notif            Exe[i-1]
                     *            [i-1]   [i]             Exe[i]
                     *    
                     * Then we overwrite notification dates and exercise dates 
                     * as follows
                     *
                     *   ----------|------|----------------|--
                     *                  Notif[i-1]      Exe[i-1]
                     *                  Notif[i]        Exe[i]
                     *    
                     */
                    ExtSchDL[i] = ExtSchDL[i-1];
                    ExtNotifDL[i-1] = ExtNotifDL[i];
                }
            }
        }/* for i */
    }/* if (hasNotifDates) */
    else
    {
        ExtSchDL   = ExtendedDL;
        ExtNotifDL = NULL;
        ExtendedDL = NULL; /* reset to NULL for safe FREE */
    }

    
    /**********************************************/
    /* interpolate the supplementary curve values */
    /**********************************************/

    /* allocate memory for the extended supplementary lists */
    for (j=0; j<5; j++)
    {
        if (Curve[j] != NULL)
        {
            ExtCurve[j] = (double *)DR_Array(DOUBLE, 0, NbExtSchDL-1);
            if (ExtCurve[j] == NULL) goto FREE_MEM_AND_RETURN;
        }
    }

    for (i=0; i<NbExtSchDL; i++)
    {
        /* find the nearest lower offset */
        k = GetDLOffset(*NbSchDates, *SchDates, ExtSchDL[i], CbkLOWER);
        if (k < 0)
        {
            DR_Error("DrExtendAmerSch: unable to interpolate values!");
            goto FREE_MEM_AND_RETURN;
        }

        /* interpolate the values */
        if (ExtSchDL[i] == (*SchDates)[k])
        {
            for (j=0; j<5; j++) /* for each curve */
            {
                if (Curve[j] != NULL) ExtCurve[j][i] = Curve[j][k];
            }
        }
        else
        {
            if (k >= *NbSchDates-1)
            {
                DR_Error("DrExtendAmerSch: cannot extrapolate!");
                goto FREE_MEM_AND_RETURN;
            }

            for (j=0; j<5; j++) /* for each curve */
            {
                if (Curve[j] != NULL)
                {
                    switch (InterpStyle)                    
                    {
                    case 'S':
                        ExtCurve[j][i] = Curve[j][k];
                        break;
                    case 'L':
                        dlinterp(ExtSchDL[i],
                                 &(ExtCurve[j][i]),
                                 (*SchDates)[k], (*SchDates)[k+1],
                                 Curve[j][k],    Curve[j][k+1]);
                        break;
                    default:
                        goto FREE_MEM_AND_RETURN; /* should not get here */

                    }/* switch */
                }/* if not NULL */
            }/* for j */
        } /* else */
    } /* for i */

    /* transfer the arrays over for output */
    Free_DR_Array(*SchDates, LONG, 0, *NbSchDates-1);
    *SchDates = ExtSchDL;

    if (hasNotifDates)
    {
        Free_DR_Array(*NotifDates, LONG, 0, *NbSchDates-1);
        *NotifDates = ExtNotifDL;
    }

    for (j=0; j<5; j++) Free_DR_Array(Curve[j], DOUBLE, 0, *NbSchDates-1);
    if (Curve0 != NULL) *Curve0 = ExtCurve[0];
    if (Curve1 != NULL) *Curve1 = ExtCurve[1];
    if (Curve2 != NULL) *Curve2 = ExtCurve[2];
    if (Curve3 != NULL) *Curve3 = ExtCurve[3];
    if (Curve4 != NULL) *Curve4 = ExtCurve[4];

    *NbSchDates = NbExtSchDL;

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    Free_DR_Array(TmpDL,        LONG,   0, NbTmpDL-1);
    Free_DR_Array(NotifDaysOfs, DOUBLE, 0, *NbSchDates-1);

    if (status == FAILURE)
    {
        DR_Error("DrExtendAmerSch: Failed!");
        for (j=0; j<5; j++) Free_DR_Array(ExtCurve[j],DOUBLE,0,NbExtSchDL-1);
        Free_DR_Array(ExtSchDL,   LONG, 0, NbExtSchDL-1);
        Free_DR_Array(ExtNotifDL, LONG, 0, NbExtSchDL-1);
        Free_DR_Array(ExtendedDL, LONG, 0, NbExtSchDL-1);
    }

    return (status);

}  /* DrExtendAmerSch */

