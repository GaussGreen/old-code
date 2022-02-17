#include <evtstat.h>
#include <fix123head.h>

#include <math.h>
#include <string.h>
#include <assert.h>

#define  CAPACITY (sizeof(es->sched) / sizeof(es->sched[0]))
int
Fix3EventStatsScheduleInit(FIX3_EVENT_STATS_SCHEDULE* es, FIX3_TREE_DATA const* tree_data, FIX3_DEV_DATA const* dev_data, int smooth)
{
    static char const* routine = "Fix3EventStatsScheduleInit";
    unsigned i;

    es->tree_data = tree_data;
    es->dev_data  = dev_data;
    es->smooth    = smooth;

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
        DR_Error("%s - slice allocation failed", routine);
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
Fix3EventStatsScheduleUpdate(FIX3_EVENT_STATS_SCHEDULE* es, double const* trigIndicator, int add, int t)
{
    static char const* routine = "Fix3EventStatsScheduleUpdate";

    unsigned  i;

    ESL_DATE curDate = es->tree_data->TPDate[t];

    /* move all slices present in the schedule from t+1 to t */
    for (i=0; i<es->size; ++i)
        Fix3_Ev(es->sched[i].probSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);

    Fix3_Ev(es->probSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);
    Fix3_Ev(es->timeSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);
    Fix3_Ev(es->tsqrSlice, t, es->tree_data->NbTP, es->dev_data, es->tree_data);

    if (es->lKnownDate == 0)
        es->lKnownDate = curDate;

    /* make exercise decision and update each probability slice */
    if (trigIndicator)
    {
        double   t2Event = Daysact(es->tree_data->TPDate[0], curDate) / 365.0;

        if (add)
        {
            struct Fix3EventStats stats;

            stats.date      = curDate;
            stats.time      = t2Event;
            stats.probSlice = Fix3_Alloc_Slice(es->tree_data);

            if (!stats.probSlice)
            {
                DR_Error("%s - slice allocation failed", routine);
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

        /* based on the trigger indicator slice update probability slices */
        for (i=0; i<es->size; ++i)
            Fix3_SmoothStepUp(es->sched[i].probSlice, trigIndicator, 1.0, 0.0, es->smooth, t, es->tree_data); 

        Fix3_SmoothStepUp(es->probSlice, trigIndicator, 1.0,             0.0, es->smooth, t, es->tree_data); 
        Fix3_SmoothStepUp(es->timeSlice, trigIndicator, t2Event,         0.0, es->smooth, t, es->tree_data); 
        Fix3_SmoothStepUp(es->tsqrSlice, trigIndicator, t2Event*t2Event, 0.0, es->smooth, t, es->tree_data); 
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
        DR_Error("%s - does not have valid stats", routine);

    return (es->probSlice + offset)[0];
}



double   
Fix3EventStatsScheduleGetProbabilityAsOfDate(FIX3_EVENT_STATS_SCHEDULE const* es, ESL_DATE date)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeExp";

    int i, idx=-1;

    if (!es->valid)
        DR_Error("%s - does not have valid stats", routine);

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
        DR_Error("%s - no event statistics calculated for %d or earlier date", routine, date);
        return -1.0;
    }

    return es->sched[idx].prob;
}



double   
Fix3EventStatsScheduleGetTimeExp(FIX3_EVENT_STATS_SCHEDULE const* es, int  conditional)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeExp";
    unsigned offset = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);
    double   tte_exp= (es->timeSlice + offset)[0];

    if (!es->valid)
        DR_Error("%s - does not have valid stats", routine);

    if (conditional && (es->probSlice + offset)[0] > TINY)
        tte_exp /= (es->probSlice + offset)[0];

    return tte_exp;
}



double   
Fix3EventStatsScheduleGetTimeStd(FIX3_EVENT_STATS_SCHEDULE const* es, int  conditional)
{
    static char const* routine = "Fix3EventStatsScheduleGetTimeStd";
    double   tte_exp = Fix3EventStatsScheduleGetTimeExp(es, conditional);

    unsigned offset  = Fix3_Node_Offset(es->tree_data->NbFactor, 0, 0, 0, es->tree_data);
    double   tte_sqr = (es->tsqrSlice + offset)[0];
    double   tte_var;

    if (!es->valid)
        DR_Error("%s - does not have valid stats", routine);

    if (conditional && (es->probSlice + offset)[0] > TINY)
        tte_sqr /= (es->probSlice + offset)[0];

    tte_var = tte_sqr - pow(tte_exp, 2.0);

    if (tte_var < TINY)
        return 0.0;

    return sqrt(tte_var);
}



double   
Fix3EventStatsScheduleGetFugit(FIX3_EVENT_STATS_SCHEDULE const* es, int lastEvent)
{
    double tte_exp = Fix3EventStatsScheduleGetTimeExp(es, 0);
    double evt_prob= Fix3EventStatsScheduleGetProbability(es);
    double t2live  = Daysact(es->tree_data->TPDate[0], lastEvent ? es->lEventDate : es->lKnownDate) / 365.0;

    return tte_exp + (1.0 - evt_prob) * t2live;
}

