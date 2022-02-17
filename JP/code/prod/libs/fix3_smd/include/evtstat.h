#ifndef evtstat_dot_h
#define evtstat_dot_h

#include <esl.h>
#include <fix123.h>

struct Fix3EventStats
{
    ESL_DATE    date;      /**< event date                                   */
    double      time;      /**< time to the event date in years              */
    double*     probSlice; /**< event probability slice                      */
    double      prob;      /**< probability calculated at the end of rollback*/
};


/** Event (exercise, knoc-in/out) calculator. Calculates cumulative probability
 *  of the event occuring on given dates. From these it calculates expected 
 *  value and standard deviation of the event as well as the FUGIT number.
 */
typedef struct Fix3EventStatsSchedule
{
    struct Fix3EventStats      sched[128]; /**< event vector            */
    unsigned                   size;       /**< event vector size       */

    FIX3_TREE_DATA const*      tree_data;  /**< tree data               */
    FIX3_DEV_DATA  const*      dev_data;   /**< dev data                */
    int                        smooth;     /**< smoothing flag          */

    /* following are necessary to compute proper stats for long schedules */
    double*                    probSlice;  /**< total prob slice        */
    double*                    timeSlice;  /**< total time slice        */
    double*                    tsqrSlice;  /**< total time square slice */

    ESL_DATE                   fEventDate; /**< first event date        */
    ESL_DATE                   lEventDate; /**< last  event date        */
    ESL_DATE                   lKnownDate; /**< last  valuation date    */

    int                        valid;      /**< 1 if has valid stats    */

} FIX3_EVENT_STATS_SCHEDULE;


/** Event statistics constructor
 *
 *  @param  es      pointer to the event statistics object
 *  @param  tree    pointer to the tree object
 *  @param  dev     pointer to the DEV data object
 *
 *  @return FAILURE if unable to allocate statistic slices
 *
 *  NOTE: Must be called before first use.
 */
int
Fix3EventStatsScheduleInit(FIX3_EVENT_STATS_SCHEDULE* es, FIX3_TREE_DATA const* tree_data, FIX3_DEV_DATA const* dev_data, int smooth);



/** Event statistics destructor
 *
 *  @param  es      pointer to the event statistics object
 *
 *  NOTE: Must be called when finished using it. 
 */
void
Fix3EventStatsScheduleClear(FIX3_EVENT_STATS_SCHEDULE* es);



/** Update event statistics to the current timepoint. This method must be called 
 *  for every timepoint during tree rollback process. When the current timepoint
 *  is the event date, a valid indicator slice must be passed. The event is 
 *  triggered if the value of the trigger indicator slice is greater than 0. 
 *  If the add flag is not zero the new statistics for the current date is added
 *  to the schedule. If the trigger indicator slice pointer is null only the 
 *  statistics that are currently in the schedule are updated. Only when the tree 
 *  rollback reaches the beginning of the timeline the final event probabilities 
 *  are calculated.
 *
 *  @param  es              pointer to the event statistics object
 *  @param  trigIndicator   underlier/trigger slice
 *  @param  add             1 to add statistics for the current date
 *  @param  t               time index in the tree
 *
 *  @return FAILURE if unable to allocate new probability slice
 *
 *  NOTE:   add flag is ignored when indicator slice is null
 */
int 
Fix3EventStatsScheduleUpdate(FIX3_EVENT_STATS_SCHEDULE* es, double const* trigIndicator, int add, int t);


/** Get total probability of event occuring on the last date in the schedule.
 *
 *  @param  es          pointer to the event statistics object
 */
double   
Fix3EventStatsScheduleGetProbability(FIX3_EVENT_STATS_SCHEDULE const* es);

/** Get total probability of event occuring on or before given date.
 *
 *  @param  es          pointer to the event statistics object
 *  @param  date        date for which to report the probability
 */
double   
Fix3EventStatsScheduleGetProbabilityAsOfDate(FIX3_EVENT_STATS_SCHEDULE const* es, ESL_DATE date);

/** Get the expected value in years of the time to the event date. Expectation
 *  is calculated and as a sum of time to event date times the event probability 
 *  for every event date. If the conditional flag is true then the probability 
 *  is scaled by the total probability of the event.
 *
 *  @param  es          pointer to the event statistics object
 *  @param  conditional total probability scaling flag
 */
double   
Fix3EventStatsScheduleGetTimeExp(FIX3_EVENT_STATS_SCHEDULE const* es, int  conditional);



/** Get the standard deviation in years of the time to the event date. If 
 *  conditional flag is true then the probability used to calculate the variance 
 *  is scaled by the total probability of the event.
 *
 *  @param  es          pointer to the event statistics object
 *  @param  conditional total probability scaling flag
 */
double   
Fix3EventStatsScheduleGetTimeStd(FIX3_EVENT_STATS_SCHEDULE const* es, int  conditional);



/** Get fugit. Fugit number is calculated as the expected value of the time
 *  when the even may occur plus the time to the last event date multiplied by 
 *  the probability of the event not occuring before that date.
 *
 *  @param  es          pointer to the event statistics object
 *  @lastEvent          if 1 fugit is calculated to talst event date, else to last
 *                      known valuation date
 */
double   
Fix3EventStatsScheduleGetFugit(FIX3_EVENT_STATS_SCHEDULE const* es, int lastEvent);


#endif

