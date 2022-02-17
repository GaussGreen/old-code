#include <optplus.hpp>
#include <fix123oper.hpp>
#include <fix123head.h>
#include <esl_log.h>

#include <math.h>
#include <assert.h>

Fix3OptionPlus::Fix3OptionPlus(
        int optType, 
        bool smooth, 
        FIX3_TREE_DATA const* tree, 
        FIX3_DEV_DATA const* dev) 
: m_optType(optType), m_smooth(smooth), m_tree(tree), m_dev(dev)
{
    // allocate aux slice
    m_aux = Fix3_Alloc_Slice(tree);
}


Fix3OptionPlus::~Fix3OptionPlus()
{
    for (size_t i=0; i<m_sched.size(); ++i)
        Fix3_Free_Slice(m_sched[i].m_prob, m_tree);

    Fix3_Free_Slice(m_aux, m_tree);
}


void 
Fix3OptionPlus::update(double* opt, double const* under, int curve,
        double notional, double strike, bool exercise, bool stat, int t)
{
    // ------------------------------------------------------------------------
    // exercise decision function: lhs = notional * type * (under-strike) - opt
    // ------------------------------------------------------------------------
    struct ExIndicator : public BinaryOper
    {
        ExIndicator(double notional, double type, double strike) 
            : m_notional(notional), m_type(type), m_strike(strike) {}

        virtual void execute(double* lhs, double const* x, double const* y, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i)
                lhs[i] = m_notional * m_type * (x[i] - m_strike) - y[i];
        }
        double m_notional;
        double m_type;
        double m_strike;
    };

    // ------------------------------------------------------------------------
    // smooth step-up function lhs = smoothStep(up, lhs, x, barrier, step)
    // ------------------------------------------------------------------------
    struct StepUp : public UnarySmoothOper
    {
        StepUp(double up, double barrier, bool smooth) 
            : UnarySmoothOper(smooth), m_up(up), m_barrier(barrier) {}

        virtual void execute(double* lhs, double const* x, double const* step, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i)
                lhs[i] = DrSmoothStep(m_up, lhs[i], x[i], m_barrier, step[i]);
        }
        double m_up;
        double m_barrier;
    };

    // ------------------------------------------------------------------------
    // add exercise price to the option lhs = smoothMax(x, 0, step) + lhs
    // ------------------------------------------------------------------------
    struct AddEx : public UnarySmoothOper
    {
        AddEx(bool smooth) : UnarySmoothOper(smooth) {}
        virtual void execute(double* lhs, double const* x, double const* step, size_t cnt) const
        {
            for (size_t i=0; i<cnt; ++i)
                lhs[i] += DrSmoothMax(x[i], 0.0, step[i]);
        }
        double m_up;
        double m_barrier;
    };



    ESL_DATE curDate = m_tree->TPDate[t];

    // move option slice from t+1 to t
    Fix3_Dev(opt, t, m_tree->NbTP, curve, m_dev, m_tree);

    // move all slices present in the schedule from t+1 to t
    for (size_t i=0; i<m_sched.size(); ++i)
        Fix3_Ev(m_sched[i].m_prob, t, m_tree->NbTP, m_dev, m_tree);


    // make exercise decision and update each probability slice
    if (exercise)
    {
        // new probability slice needed on exercise
        if (stat)
        {
            ExNode node;

            node.m_date = curDate;
            node.m_time = Daysact(m_tree->TPDate[0], curDate) / 365.0;
            node.m_prob = Fix3_Alloc_Slice(m_tree);

            // bit inefficient but I can avoid the iterator
            m_sched.insert(m_sched.begin(), node);
        }

        // decide whether to exercise: aux = notional * type * (under-strike) - opt
        Fix3_SliceOper(m_tree, t, m_aux, under, opt, ExIndicator(notional, m_optType, strike)); 

        // add option exercise value: opt += smoothMax(aux, 0)
        Fix3_SliceOper(m_tree, t, opt, m_aux, AddEx(m_smooth)); 

        // based on the exercise decision update probability slices
        for (size_t i=0; i<m_sched.size(); ++i)
            Fix3_SliceOper(m_tree, t, m_sched[i].m_prob, m_aux, StepUp(1.0, 0.0, m_smooth)); 
    }
}



double   
Fix3OptionPlus::getExerciseProbability()
{
    static char const* routine = "Fix3OptionPlus::getExerciseProbability";

    if (!m_sched.size())
        throw EslException(routine) << " - empty statistics schedule";

    return getExerciseProbability(m_sched.size() - 1);
}

double   
Fix3OptionPlus::getExerciseProbability(size_t i)
{
    static char const* routine = "Fix3OptionPlus::getExerciseProbability";

    if (!m_sched.size())
        throw EslException(routine) << " - empty statistics schedule";

    size_t offset = Fix3_Node_Offset(m_tree->NbFactor, 0, 0, 0, m_tree);
    return (m_sched[i].m_prob + offset)[0];
}

double   
Fix3OptionPlus::getTimeToExerciseExp(bool conditional)
{
    static char const* routine = "Fix3OptionPlus::getTimeToExerciseExp";

    if (!m_sched.size())
        throw EslException(routine) << " - empty statistics schedule";

    if (getExerciseProbability() < TINY)
        return 999999.0;

    double prev_prob = 0.0;
    double curr_prob = 0.0;

    double tte = 0.0;
    for (size_t i=0; i<m_sched.size(); ++i)
    {
        prev_prob = curr_prob;
        curr_prob = getExerciseProbability(i);
        tte += m_sched[i].m_time * (curr_prob - prev_prob);
    }

    if (conditional)
    {
        tte /= getExerciseProbability();
    }

    return tte;
}

double   
Fix3OptionPlus::getTimeToExerciseStd(bool conditional)
{
    static char const* routine = "Fix3OptionPlus::getTimeToExerciseStd";

    if (!m_sched.size())
        throw EslException(routine) << " - empty statistics schedule";

    if (getExerciseProbability() < TINY)
        return sqrt(TINY);


    double prev_prob = 0.0;
    double curr_prob = 0.0;

    double cond_prob = 1.0;
    if (conditional)
    {
        cond_prob   = getExerciseProbability();
    }

    double tte_exp = getTimeToExerciseExp();
    double tte_var = 0.0;
    for (size_t i=0; i<m_sched.size(); ++i)
    {
        prev_prob = curr_prob;
        curr_prob = getExerciseProbability(i);
        tte_var += pow(m_sched[i].m_time - tte_exp, 2.0) * (curr_prob - prev_prob)/cond_prob;
    }

    return sqrt(tte_var);
}

double   
Fix3OptionPlus::getFugit()
{
    static char const* routine = "Fix3OptionPlus::getFugit";
    if (!m_sched.size())
        throw EslException(routine) << " - empty statistics schedule";


    // calculate FUGIT = E(timeToEx) + (1 - ProbOfEx) * TimeToLastEx

    size_t last = m_sched.size() - 1;

    double ex_prob = getExerciseProbability();
    double ex_time = m_sched[last].m_time;

    if (ex_prob < TINY)
        return ex_time;
    return getTimeToExerciseExp(false) + (1.0 - ex_prob) * ex_time;
}


