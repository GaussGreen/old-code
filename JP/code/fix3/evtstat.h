#ifndef evtstat_dot_h
#define evtstat_dot_h

#include <esl.h>
#include <fix123.h>

typedef struct Fix3InputsEventStats
{
    char        OptCalcStats;
    char        StatSchedStyle;
    int         NbStatDates;
    long        EventStatDate[MAXNBDATE];
} FIX3_EVENT_STATS_INPUT;

struct Fix3EventStats
{
    IRDate       date;      /**< event date                                   */
    double      time;      /**< time to the event date in years              */
    double*     probSlice; /**< event probability slice                      */
    double      prob;      /**< probability calculated at the end of rollback*/
};


/** Event (exercise, knoc-in/out) calculator. Calculates cumulative probability
 *  of the event occuring on or before given dates. From these it calculates 
 *  expected value and standard deviation of the event as well as the FUGIT 
 *  number.
 */
typedef struct Fix3EventStatsSchedule
{
    struct Fix3EventStats      sched[128]; /**< event vector            */
    unsigned                   size;       /**< event vector size       */

    FIX3_TREE_DATA const*      tree_data;  /**< tree data               */
    FIX3_DEV_DATA  const*      dev_data;   /**< dev data                */

    /* following are necessary to compute proper stats for long schedules */
    double*                    probSlice;  /**< total prob slice        */
    double*                    timeSlice;  /**< total time slice        */
    double*                    tsqrSlice;  /**< total time square slice */

    IRDate                      fEventDate; /**< first event date        */
    IRDate                      lEventDate; /**< last  event date        */
    IRDate                      lKnownDate; /**< last  valuation date    */

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
Fix3EventStatsScheduleInit(FIX3_EVENT_STATS_SCHEDULE* es, FIX3_TREE_DATA const* tree_data, FIX3_DEV_DATA const* dev_data);



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
 *  is the event date, a valid indicator slice must be passed. The indicator 
 *  slice is treated as the probability of the event occuring for a given node.
 *  The values must be either 0 or 1 but when smoothing is used, intermediate
 *  values are allowed.
 *  If the add flag is not zero the new statistics for the current date is added
 *  to the schedule. If the trigger indicator slice pointer is null only the 
 *  statistics that are currently in the schedule are updated. Only when the tree 
 *  rollback reaches the beginning of the timeline the final event probabilities 
 *  are calculated.
 *
 *  @param  es              pointer to the event statistics object
 *  @param  indicator       trigger probability slice
 *  @param  add             1 to add statistics for the current date
 *  @param  t               time index in the tree
 *
 *  @return FAILURE if unable to allocate new probability slice
 *
 *  NOTE:   add flag is ignored when indicator slice is null
 */
int 
Fix3EventStatsScheduleUpdate(FIX3_EVENT_STATS_SCHEDULE* es, double const* indicator, int add, int t);


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
Fix3EventStatsScheduleGetProbabilityAsOfDate(FIX3_EVENT_STATS_SCHEDULE const* es, IRDate date);

/** Get the expected value in years of the time to the event date conditional on
 *  the event occuring. 
 *
 *  @param  es          pointer to the event statistics object
 */
double   
Fix3EventStatsScheduleGetTimeExp(FIX3_EVENT_STATS_SCHEDULE const* es);



/** Get the standard deviation in years of the time to the event date conditional on
 *  the event occuring.
 *
 *  @param  es          pointer to the event statistics object
 */
double   
Fix3EventStatsScheduleGetTimeStd(FIX3_EVENT_STATS_SCHEDULE const* es);



/** Get fugit. Fugit number is calculated as the expected value of the time
 *  when the event may occur plus the time to the last event date multiplied by 
 *  the probability of the event not occuring before that date.
 *
 *  @param  es          pointer to the event statistics object
 *  @lastEvent          if 1 fugit is calculated to last event date, else to last
 *                      known valuation date
 */
double   
Fix3EventStatsScheduleGetFugit(FIX3_EVENT_STATS_SCHEDULE const* es, int lastEvent);

/**
*  	Read Event Statistic parameters and check validity of input.
*/
int     Fix3EventStatsInput (   
             FIX3_EVENT_STATS_INPUT *eventstat_data, /**< (O) Event Stats data              */
             char const             *FileName);      /**< (I) File name including extension */

/**
*       Creat Event Statistic file
*/
int     Fix3EventStatsPrint (
             FIX3_EVENT_STATS_INPUT *eventstat_data, /**< (I) Event Stats data              */
             char const              *FileName);     /**< (I) File name including extension */

#endif

