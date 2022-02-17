#include <evtstat.h>
#include <fix123head.h>

#include <math.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#define  CAPACITY (sizeof(es->sched) / sizeof(es->sched[0]))
int
Fix3EventStatsScheduleInit(FIX3_EVENT_STATS_SCHEDULE* es, FIX3_TREE_DATA const* tree_data, FIX3_DEV_DATA const* dev_data)
{
    static char const* routine = "Fix3EventStatsScheduleInit";
    unsigned i;

    es->tree_data = tree_data;
    es->dev_data  = dev_data;

    for (i=0; i<CAPACITY; ++i)
    {
        es->sched[i].date      = 0;
        es->sched[i].time      =-1;
        es->sched[i].probSlice = NULL;
        es->sched[i].prob      =-1;
    }
    es->size = 0;

    es->probSlice = Fix3_Alloc_Slice(es->tree_data);
    es->timeSlice = Fix3_Alloc_Slice(es->tree_data);
    es->tsqrSlice = Fix3_Alloc_Slice(es->tree_data);

    es->fEventDate= 0;
    es->lEventDate= 0;
    es->lKnownDate= 0;
    es->valid     = 0;

    if (es->probSlice == NULL || es->timeSlice == NULL || es->tsqrSlice== NULL)
    {
        DR_Error("%s - slice allocation failed\n", routine);
        return FAILURE;
    }
    return SUCCESS;
}



void
Fix3EventStatsScheduleClear(FIX3_EVENT_STATS_SCHEDULE* es)
{
    unsigned i;
    for (i=0; i<CAPACITY; ++i)
        Fix3_Free_Slice(es->sched[i].probSlice, es->tree_data);

    Fix3_Free_Slice(es->probSlice, es->tree_data);
    Fix3_Free_Slice(es->timeSlice, es->tree_data);
    Fix3_Free_Slice(es->tsqrSlice, es->tree_data);
}


int 
Fix3EventStatsScheduleUpdate(FIX3_EVENT_STATS_SCHEDULE* es, double const* indicator, int add, int t)
{
    static char const* routine = "Fix3EventStatsScheduleUpdate";

    unsigned  i;

    IRDate curDate = es->tree_data->TPDate[t];

    /* move all slices present in the schedule from t+1 to t */
    for (i=0; i<es->size; ++i) {
        Fix3_Ev(es->sched[i].probSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);
    }

    Fix3_Ev(es->probSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);
    Fix3_Ev(es->timeSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);
    Fix3_Ev(es->tsqrSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);

    /* use last valid date to initialize time and time squared slices */
    if (es->lKnownDate == 0) {
        es->lKnownDate = curDate;
    }

    /* make exercise decision and update each probability slice */
    if (indicator)
    {
        double t2Event = Daysact(es->tree_data->TPDate[0], curDate) / 365.0;

        if (add)
        {
            struct Fix3EventStats stats;

            stats.date      = curDate;
            stats.time      = t2Event;
            stats.probSlice = Fix3_Alloc_Slice(es->tree_data);

            if (!stats.probSlice)
            {
                DR_Error("%s - slice allocation failed\n", routine);
                return FAILURE;
            }

            if (es->size == CAPACITY)
                Fix3_Free_Slice(es->sched[CAPACITY-1].probSlice, es->tree_data);

            /* inefficient and ugly, if concerned there is C++ version */
            memmove(&es->sched[1], &es->sched[0], (es->size < CAPACITY ? es->size : CAPACITY-1) * sizeof(es->sched[0]));
            es->sched[0] = stats; 

            if (es->size < CAPACITY)
                ++(es->size);
        }

        /* capture the last and first date when the event may occur */
        if (es->lEventDate == 0)
            es->lEventDate = curDate;

        es->fEventDate = curDate;
        es->valid = TRUE;

        /* based on the trigger indicator slice update probability slices: p = p | i */
        for (i=0; i<es->size; ++i) {
            Fix3_LogicalOr(es->sched[i].probSlice, es->sched[i].probSlice, indicator, t, es->tree_data); 
        }

        Fix3_LogicalOr(es->probSlice, es->probSlice, indicator, t, es->tree_data); 

        /* x = x * (1-p) + y * p */
        Fix3_Expectation(es->timeSlice, es->timeSlice, t2Event,         indicator, t, es->tree_data); 
        Fix3_Expectation(es->tsqrSlice, es->tsqrSlice, t2Event*t2Event, indicator, t, es->tree_data); 
    }

    if (t == 0)
    {
        unsigned offset = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);
        for (i=0; i<es->size; ++i)
            es->sched[i].prob = (es->sched[i].probSlice + offset)[0];
    }

    return SUCCESS;
}



double   
Fix3EventStatsScheduleGetProbability(FIX3_EVENT_STATS_SCHEDULE const* es)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeExp";

    unsigned offset = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);

    if (!es->valid)
        DR_Error("%s - does not have valid stats\n", routine);

    return (es->probSlice + offset)[0];
}



double   
Fix3EventStatsScheduleGetProbabilityAsOfDate(FIX3_EVENT_STATS_SCHEDULE const* es, IRDate date)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeExp";

    int i, idx=-1;

    if (!es->valid)
        DR_Error("%s - does not have valid stats\n", routine);

    if (date < es->fEventDate)
        return 0.0;

    for (i=0; i<(int)es->size; ++i)
    {
        if (es->sched[i].date <= date)
            idx = i;
        else
            break;
    }

    if (idx < 0)
    {
        DR_Error("%s - no event statistics calculated for %ld or earlier date\n", routine, date);
        return -1.0;
    }

    return es->sched[idx].prob;
}



double   
Fix3EventStatsScheduleGetTimeExp(FIX3_EVENT_STATS_SCHEDULE const* es)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeExp";
    unsigned offset = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);
    double   tte_exp= (es->timeSlice + offset)[0];

    if (!es->valid)
        DR_Error("%s - do not have valid stats\n", routine);

    if ((es->probSlice + offset)[0] > TINY) {
        tte_exp /= (es->probSlice + offset)[0];
    }
    else {
      if (fabs(tte_exp) > TINY)
        tte_exp = -1;  /* non-zero divided by 0 is not defined */
    }
    return tte_exp;
}



double   
Fix3EventStatsScheduleGetTimeStd(FIX3_EVENT_STATS_SCHEDULE const* es)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeStd";
    double   tte_exp = Fix3EventStatsScheduleGetTimeExp(es);

    unsigned offset  = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);
    double   tte_sqr = (es->tsqrSlice + offset)[0];
    double   tte_var = 0.0;

    if (!es->valid)
        DR_Error("%s - do not have valid stats\n", routine);

    if ((es->probSlice + offset)[0] > TINY) {
        tte_sqr /= (es->probSlice + offset)[0];
        tte_var = tte_sqr - pow(tte_exp, 2.0);
    }
    else {
        tte_sqr = 0.0;
    }

    if (tte_var < TINY)
        return 0.0;

    return sqrt(tte_var);
}

double   
Fix3EventStatsScheduleGetFugit(FIX3_EVENT_STATS_SCHEDULE const* es, int lastEvent)
{
    static char const* routine = "Fix3EventStatsScheduleGetFugit";

    /* here we have small problem - old code used unconditional expectation */
    double tte_exp = Fix3EventStatsScheduleGetTimeExp(es);
    double evt_prob= Fix3EventStatsScheduleGetProbability(es);
    double t2live  = Daysact(es->tree_data->TPDate[0], lastEvent ? es->lEventDate : es->lKnownDate) / 365.0;

    if (!es->valid)
        DR_Error("%s - do not have valid stats\n", routine);

    if (tte_exp < 0.0)
        return t2live;
    else
        /* tte_exp is conditional expectations - need to correct it */
        return tte_exp * evt_prob  + (1.0 - evt_prob) * t2live;
}

/*****  Fix3EventStatsInput  ********************************************************/
/**
*  	Read Event Statistic parameters and check validity of input.
*/
int     Fix3EventStatsInput (   
             FIX3_EVENT_STATS_INPUT *eventstat_data, /**< (O) Event Stats data              */
             char const             *FileName)       /**< (I) File name including extension */
{
    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    stream = fopen (FileName, "r");

    /*
     *  If there is no parameter file.
     */

    if (stream == NULL)
    {
         eventstat_data->OptCalcStats = 'N';
         eventstat_data->StatSchedStyle = 'I';
         eventstat_data->NbStatDates = 0;

	 for(i=0; i<MAXNBDATE;i++)
              eventstat_data->EventStatDate[i] = 0;

         status = SUCCESS;
	 goto RETURN;
    }

    if (FindAndSkipComLine(
                  stream, 
		  "Calc Opt Stats", 
		  "Fix3EventStatsInput", 
		  FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%c \n", 
                         &(eventstat_data->OptCalcStats));

    if (readerror != 1)
    {        
         sprintf (ErrorMsg, "Could not find Calc Opt Stats in file %s! (Fix3EventStatsInput)", FileName);
         DR_Error(ErrorMsg);
         goto RETURN;
    }

    eventstat_data->OptCalcStats = (char)toupper(eventstat_data->OptCalcStats);

    if (eventstat_data->OptCalcStats == 'N')
    {
         eventstat_data->StatSchedStyle = 'I';
         eventstat_data->NbStatDates = 0;

	 for(i=0; i<MAXNBDATE;i++)
              eventstat_data->EventStatDate[i] = 0;

         status = SUCCESS;
	 goto RETURN;
    }

    if (FindAndSkipComLine(
                  stream, 
		  "Calc Opt Sched Style", 
		  "Fix3EventStatsInput", 
		  FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%c \n", 
                         &(eventstat_data->StatSchedStyle));

    if (readerror != 1)
    {        
         sprintf (ErrorMsg, "Could not find Calc Sched Style in file %s! (Fix3EventStatsInput)", FileName);
         DR_Error(ErrorMsg);
         goto RETURN;
    }

    eventstat_data->StatSchedStyle = (char)toupper(eventstat_data->StatSchedStyle);

    if (FindAndSkipComLine(
                  stream, 
		  "Nb Stat Dates", 
		  "Fix3EventStatsInput", 
		  FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                         &(eventstat_data->NbStatDates));

    if (readerror != 1)
    {        
         sprintf (ErrorMsg, "Could not find Nb Stat Dates in file %s! (Fix3EventStatsInput)", FileName);
         DR_Error(ErrorMsg);
         goto RETURN;
    }

    if (eventstat_data->NbStatDates <= 0)
    {        
         sprintf (ErrorMsg, "Nb Stat Dates %d can not be <= 0! (Fix3EventStatsInput)", eventstat_data->NbStatDates);
         DR_Error(ErrorMsg);
         goto RETURN;
    }

    if (FindAndSkipComLine(
                  stream, 
		  "Stat Dates", 
		  "Fix3EventStatsInput", 
		  FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for(i=0; i<eventstat_data->NbStatDates; i++)
    {
        readerror = fscanf (stream, 
                            "%ld \n", 
                            &(eventstat_data->EventStatDate[i]));
        if (readerror != 1)
        {        
             sprintf (ErrorMsg, "Could not find Stat Dates nb %d in file %s! (Fix3EventStatsInput)", i+1, FileName);
             DR_Error(ErrorMsg);
             goto RETURN;
        }
        eventstat_data->EventStatDate[i] = IRDateFromYMDDate(eventstat_data->EventStatDate[i]);
    }

    if (   (eventstat_data->OptCalcStats != 'N')
        && (eventstat_data->OptCalcStats != 'Y'))
    {
        sprintf (ErrorMsg, "Opt Calc Stat Flag can only be Y or N, in file %s! (Fix3EventStatsInput)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (   (eventstat_data->StatSchedStyle != 'I')
        && (eventstat_data->StatSchedStyle != 'D')
        && (eventstat_data->StatSchedStyle != 'W')
        && (eventstat_data->StatSchedStyle != 'M')
        && (eventstat_data->StatSchedStyle != 'Q')
        && (eventstat_data->StatSchedStyle != 'S')
        && (eventstat_data->StatSchedStyle != 'A'))
    {
        sprintf (ErrorMsg, "Stat Sched Style can only be I,D,W,M,Q,S or A, in file %s! (Fix3EventStatsInput)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    for(i=0; i<eventstat_data->NbStatDates; i++)
    {
        if (Dateok(eventstat_data->EventStatDate[i]))
        {
            DR_Error ("Incorrect format for event stat date!");
            goto RETURN;
        }
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);
} /* Fix3EventStatsInput */


/**
*       Creat Event Statistic file
*/
int     Fix3EventStatsPrint (
             FIX3_EVENT_STATS_INPUT *eventstat_data, /**< (I) Event Stats data              */
             char const           *FileName)         /**< (I) File name including extension */
{

    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    stream = fopen (FileName, "w");

    if (stream == NULL)
    {
         sprintf (ErrorMsg, "Could not write Opt Stats file %s! (Fix3EventStatsPrint)", FileName);
         DR_Error(ErrorMsg);
         goto RETURN;
    }

    if (eventstat_data->OptCalcStats == 'N')
    {
        fprintf (stream,"# Exercise Stats Flag\n");
        fprintf (stream,"N\n");
        fprintf (stream,"# Exercise statistics schedule style - I,D,W,M,Q,S,A\n");
        fprintf (stream,"I\n");
        fprintf (stream,"# Exercise statistics number of dates\n");
        fprintf (stream,"0\n");
        fprintf (stream,"# Exercise statistics dates\n");
        fprintf (stream,"\n");
    }
    else
    {
        fprintf (stream,"# Exercise Stats Flag\n");
        fprintf (stream,"%c \n",eventstat_data->OptCalcStats);
        fprintf (stream,"# Exercise statistics schedule style - I,D,W,M,Q,S,A\n");
        fprintf (stream,"%c \n",eventstat_data->StatSchedStyle);
        fprintf (stream,"# Exercise statistics number of dates\n");
        fprintf (stream,"%d \n",eventstat_data->NbStatDates);
        fprintf (stream,"# Exercise statistics dates\n");

        for (i=0; i<eventstat_data->NbStatDates; i++)
            fprintf (stream,"%ld \n",eventstat_data->EventStatDate[i]);
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

} /* Fix3EventStatsPrint */

