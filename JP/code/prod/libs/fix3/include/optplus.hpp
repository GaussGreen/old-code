#ifndef optplus_dot_hpp
#define optplus_dot_hpp

#include <vector>
#include <esl.h>
#include <fix123.h>

/** Option_Plus replacement
 */
class Fix3OptionPlus
{
    public:

    struct ExNode
    {
        IRDate                   m_date;     ///< exercise date
        double                  m_time;     ///< time to exercise
        double*                 m_prob;     ///< exercise probability slice
    };

    /** Option calculator constructor
     *  @param  optType option type - 1 for call -1 for put
     *  @param  smooth  true for payoff smoothing
     *  @param  tree    pointer to the tree object
     *  @param  dev     pointer to the DEV data object
     */
    Fix3OptionPlus(int optType, bool smooth, FIX3_TREE_DATA const* tree, FIX3_DEV_DATA const* dev);

    ~Fix3OptionPlus();

    /** Update calculator to the current timepoint.
     *  @param  opt         option price slice
     *  @param  under       underlier slice
     *  @param  curve       discount curve index
     *  @param  strike      strike value
     *  @param  exercise    exercise flag
     *  @param  stat        statistics calculation flag - add new stats for current date
     *  @param  t           time index in the tree
     */
    void update(double* opt, double const* under, int curve, double notional, double strike, 
            bool exercise, bool stat, int t);

    /// Return the current size of the exercise array
    size_t getSize() {
        return m_sched.size();
    }

    /// Get exercise date
    IRDate getExerciseDate(size_t i) {
        return m_sched[i].m_date;
    }

    /** Get exercise probability. It is a total cumulative probability of the exercise 
     *  on or before the last exercise date
     */
    double   getExerciseProbability();

    /** Get exercise probability. It is a cumulative probability of the exercise 
     *  on or before the date pointed to by the index
     */
    double   getExerciseProbability(size_t i);

    /** Get time to exercise expectation in years. Expected time to exercise is calculated as 
     *  sum of time to exercise times the exercise probability for every exercise date. If
     *  conditional flag is true then the probability is scaled by the total probability of the
     *  exercise.
     */
    double   getTimeToExerciseExp(bool conditional=true);

    /** Get time to exercise stdev in years. If conditional flag is true then the probability 
     *  used to calculate the variance is scaled by the total probability of the exercise.
     */
    double   getTimeToExerciseStd(bool conditional=true);

    /// Get fugit
    double   getFugit();

    private:

    std::vector<ExNode>         m_sched;    ///< exercise nodes along exercise schedule

    int                         m_optType;  ///< option type 1 call, -1 put
    bool                        m_smooth;   ///< smoothing flag

    FIX3_TREE_DATA const*       m_tree;     ///< tree data
    FIX3_DEV_DATA  const*       m_dev;      ///< dev data

    double*                     m_aux;      ///< temporary slice
};


#endif

