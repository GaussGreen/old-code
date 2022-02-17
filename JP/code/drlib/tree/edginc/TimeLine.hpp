//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TimeLine.hpp
//
//   Description : time line class
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#ifndef EDG_TIMELINE_HPP
#define EDG_TIMELINE_HPP

#include "edginc/Object.hpp"
#include "edginc/Array.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/TermStructure.hpp"
#include "edginc/TimeMetric.hpp"

DRLIB_BEGIN_NAMESPACE

#define MAX_TIME_STEP 10000 // maximum number of time steps (tree, MC or FD)
#define MIN_TREE_STEP 50 // minimum number of time steps for a tree or FD
#define MAX_SEG_STEP 10000
#define MIN_SEG_STEP 20


/////////////////////////////////////////////////////////
// CTimeLine class
/////////////////////////////////////////////////////////
/** Time line class generates a time line for given date intervals/segments.
    A time line can have different step densities for each segment. Spline or other
    interpolation will be used to connect density switching step when used in for a tree.
    When volatility term structure is supplied and use equal variance step
    is chosen time points will have equal variance assuming simple i.i.d. such as
    the usual log-normal model. */
class TREE_DLL CTimeLine
{
public:
    int             NumOfStep;    // num of pts = num of step + 1
    double          *TradeYrFrac; // trade year fraction between points
    double          *TradeTime;   // trade year fraction from start date to points
    CDateTimeArray  StepDates;    // Constructed date array

    // Usually there is just one segment so SegmentEnd = NumOfStep
    vector<int>     SegmentEnd;   // The end step of each segment. 
    vector<int>     SegmentEndNoCrit; // End step of segment before crit dates are inserted
    vector<int>     SegmentStart; // Start step of segment

    CTimeLine();
    virtual ~CTimeLine();

    void Clear();

    /** creates time points between [dates[0], dates[1]] and (dates[1], dates[2]], ... 
        returns number of steps */
    int CreateTimeLine(const DateTimeArray& dates, const vector<int>& density, 
                              TimeMetricConstSP metric, int stepsPerYear, double minGap, 
                              bool forMC, const CTermStructure& vol_square, 
                              DateTimeArray critDates, bool equalTime=true);


    int CreateTimeLineSimple(const DateTimeArray& dates, 
                              TimeMetricConstSP metric, int stepsPerYear,                              
                              DateTimeArray critDates, IntArray* isAddedSeg=0, int* minStepInAddedSeg = 0);


    // change tree step base
    void SetTreeStepBase(int stepBase);

protected:
    int             TreeStepBase;  // base step to be scaled for a tree if inputStep==0
            
    vector<int>     CritSteps; // time step corresponding to critical dates

    void CalcStepDate(const DateTime& startDate, const DateTime& endDate, 
                      const vector<double>& trade_t, int steps, TimeMetricConstSP metric, 
                      const DateTimeArray& dates);

    /** called by CreateTimeLine, inserting critical dates/time step after time line intialisation */ 
    void InsertCritDates(vector<double>& trade_t, double mat_trading, int step, 
                         TimeMetricConstSP metric, const DateTimeArray& critDates);
    /** is a timepoint the end of a segment */
    bool IsSegEnd(int j);
};

typedef refCountPtr<CTimeLine> TimeLineSP;

DRLIB_END_NAMESPACE
#endif
