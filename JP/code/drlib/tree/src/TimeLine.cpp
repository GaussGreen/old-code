//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TimeLine.cpp
//
//   Description : time line class
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// Constructor
CTimeLine::CTimeLine()
{
    TradeYrFrac = 0;
    TradeTime = 0;
    NumOfStep = 0;
    TreeStepBase = 200;
}

// Destructor
CTimeLine::~CTimeLine() {
    Clear();
}

void CTimeLine::Clear() {
    if (TradeTime != 0) {
        delete [] TradeTime;
        delete [] TradeYrFrac;
    }
    TradeYrFrac = 0;
    TradeTime = 0;
    NumOfStep = 0;
    SegmentEnd.clear();
    SegmentEndNoCrit.clear();
    SegmentStart.clear();
    StepDates.clear();
    CritSteps.clear();
}

/*======================================================================
*   Description:    Creates time points between [dates[0], dates[1]] and (dates[1], dates[2]], ...
*                   Returns number of steps created (there are NumOfStep+1 points)
*   Parameters:
*                   dates       : list of dates for creating time line intervals
*                   density     : density factor of num of steps in each interval
*                   metric      : TimeMetric used for year fractions etc.
*                   stepsPerYear: user input for total steps (if ==0 default mechanism used)
*                   minGap      : Now obsolete. All critical dates are inserted.
*                   forMC       : if time line is for MC use
*                   vol_square  : a vol^2 term structure for use to create equal variance step
*                   equalTime   : if equal time step is required
*
*   Return:         num of Effective steps
*
*======================================================================*/
int CTimeLine::CreateTimeLine(const DateTimeArray& dates, const vector<int>& density,
                              TimeMetricConstSP metric, int stepsPerYear, double minGap,
                              bool forMC, const CTermStructure& vol_square,
                              DateTimeArray critDates, bool equalTime)
{
    static string routine("CTimeLine::CreateTimeLine");
    try {

        int i, j, k;
        double dt, var_dt, mat_trading, vol2;
        vector<double> trade_t, var;

        vector<int> stepsInSeg;
        vector<double> seg_trading;
        vector<int> segmentStartNoCrit;
        double segStartTime;
        double varSegStart;

        double normFact;

        int totalSteps;

        int num_seg = dates.size() - 1;
        if (num_seg < 1)
            throw ModelException(routine, "must have at least two dates for generating a time line.");

        if (num_seg != (int) density.size())
            throw ModelException(routine, "num of density factor must be numOfDates -1.");

        Clear();

        // sort critical dates
        metric->SortDate(dates[0], dates[num_seg], true, critDates);

        mat_trading = metric->yearFrac(dates[0], dates[num_seg]);

        // force to do equal time steps if we have any zero vols or negative forward variances before maturity
        if (equalTime == false) {
            double varSum = 0.0;
            for (i=0; i<vol_square.size(); i++)
            {
                if ((vol_square.GetX(i)*vol_square.GetY(i) < varSum) ||
                    (Maths::isZero(vol_square.GetY(i))))                {
                    equalTime = true;
                    break;
                }
                varSum = vol_square.GetX(i)*vol_square.GetY(i);

                if (vol_square.GetX(i)>mat_trading)
                    break;
            }
        }

        // allow varying density segments
        SegmentEnd.resize(num_seg); 
        SegmentEndNoCrit.resize(num_seg);
        SegmentStart.resize(num_seg);

        stepsInSeg.resize(num_seg);
        seg_trading.resize(num_seg);
        segmentStartNoCrit.resize(num_seg);

        // compute a normalization factor for the density
        normFact = 0;
        for(i=0; i<num_seg; i++) {
            normFact += density[i]*metric->yearFrac(dates[i], dates[i+1])/mat_trading;
        }

        totalSteps = 0;
        for(i=0; i<num_seg; i++) {
            if (stepsPerYear==0 && forMC) // handle standard MC separately
            {
                stepsInSeg[i] = 1;
                totalSteps++;
                segmentStartNoCrit[i] = i;
                SegmentEndNoCrit[i] = i+1;
            } else {
                seg_trading[i] = metric->yearFrac(dates[i], dates[i+1]);
                if (stepsPerYear==0 && equalTime==false)
                {
                    // default steps in segment according to std dev over time period
                    vol2 = vol_square.InterpFlatFwd(i==0 ? 0 : seg_trading[i-1], seg_trading[i] + (i==0 ? 0 : seg_trading[i-1]) );
                    stepsInSeg[i] = (int)(TreeStepBase*sqrt(vol2*seg_trading[i]) * (density[i]/normFact));
                }
                else if (stepsPerYear<=0) // -1 => equal time daily step
                {
                    equalTime = true;
                    stepsInSeg[i] = (int)(dates[i+1].daysDiff(dates[i]) * (density[i]/normFact));
                }
                else
                {
                    stepsInSeg[i] = (int)(stepsPerYear*seg_trading[i] * (density[i]/normFact));
                }
                stepsInSeg[i] = Maths::max(stepsInSeg[i], MIN_SEG_STEP);
                stepsInSeg[i] = Maths::min(stepsInSeg[i], MAX_SEG_STEP);
                totalSteps += stepsInSeg[i];
                segmentStartNoCrit[i] = (i==0 ? 0 : SegmentEndNoCrit[i-1]);
                SegmentEndNoCrit[i] = segmentStartNoCrit[i] + stepsInSeg[i];
            }
        }

        // if there's only one segment put bounds on the number of steps
        if (num_seg == 1) {
            if (!forMC)
                stepsInSeg[0] = (int) Maths::max(stepsInSeg[0], MIN_TREE_STEP);
            // cap steps
            stepsInSeg[0] = (int) Maths::min(stepsInSeg[0], MAX_TIME_STEP);
            totalSteps = stepsInSeg[0];
            SegmentEndNoCrit[0] = stepsInSeg[0];
        }

        trade_t.resize(totalSteps + 1);
        var.resize(totalSteps + 1);

        trade_t[0] = 0.0;
        var[0] = 0.0;
        double dVar=0.0;

        for (k=0; k<num_seg; k++) {
            dt = seg_trading[k]/(double)stepsInSeg[k];
            segStartTime = trade_t[segmentStartNoCrit[k]];
            // set up time points
            if (equalTime) // use equal time step
            {
                for (i=segmentStartNoCrit[k]+1; i<=SegmentEndNoCrit[k]; i++)
                {
                    trade_t[i] = segStartTime + (i - segmentStartNoCrit[k]) * dt;
                }
            }
            else if (vol_square.size() == 1)
            {
                varSegStart = var[segmentStartNoCrit[k]];
                var_dt = dt * vol_square.InterpFlatFwd(segStartTime, segStartTime + seg_trading[k]);
                for (i=segmentStartNoCrit[k]+1; i<=SegmentEndNoCrit[k]; i++)
                {
                    trade_t[i] = segStartTime + (i - segmentStartNoCrit[k]) * dt;
                    var[i] = varSegStart + (i - segmentStartNoCrit[k]) * var_dt;
                }
            }
            else // term structures
            {// use constant var step instead of equal time step
                varSegStart = var[segmentStartNoCrit[k]];
                var_dt = dt * vol_square.InterpFlatFwd(segStartTime, segStartTime + seg_trading[k]);
                j = 0;
                for (i=segmentStartNoCrit[k]+1; i<=SegmentEndNoCrit[k]; i++)
                {
                    var[i] = varSegStart + (i - segmentStartNoCrit[k]) * var_dt;
                    if (j < vol_square.size()) // before last vol point
                    {
                        while (Maths::isNegative(vol_square.GetX(j)*vol_square.GetY(j) - var[i]))
                        {
                            j++;
                            if (j == vol_square.size()) break;
                        }
                    }
                    if (j == 0)
                        trade_t[i] = var[i]/vol_square.GetY(0);
                    else if (j == vol_square.size()) // beyond last vol point
                        trade_t[i] = vol_square.GetX(j-1)+(var[i]-vol_square.GetX(j-1)*vol_square.GetY(j-1))/vol_square.GetY(j-1);
                    else
                        trade_t[i] = vol_square.GetX(j-1)+(var[i]-vol_square.GetX(j-1)*vol_square.GetY(j-1))
                            *(vol_square.GetX(j)-vol_square.GetX(j-1))/(vol_square.GetX(j)
                                                                        *vol_square.GetY(j)-vol_square.GetX(j-1)*vol_square.GetY(j-1));

                    // Check that the variance interval between the successive times just calculated is within 1% of var_dt
                    if (Maths::isPositive(var_dt)){
                        dVar = trade_t[i]*vol_square.InterpFlatFwd(0,trade_t[i]) - trade_t[i-1]*vol_square.InterpFlatFwd(0,trade_t[i-1]);
                        if (!Maths::areEqualWithinTol(dVar, var_dt, 0.01*var_dt)) {
                            throw ModelException(routine, "Error constructing time line with equal variance steps");
                        }
                    }
                }
            }
        }
        // calculate dates for the steps
        CalcStepDate(dates[0], dates[num_seg], trade_t, totalSteps, metric, dates);

        // insert critical dates
        InsertCritDates(trade_t, mat_trading, totalSteps, metric, critDates);

        // calculate TradeYrFrac and CalTime array
        TradeYrFrac = new double[NumOfStep+1];  
        TradeTime = new double[NumOfStep+1];    
        TradeYrFrac[0] = TradeTime[0] = 0.0;
        for (i=1; i<=NumOfStep; i++)
        {
            TradeYrFrac[i] = metric->yearFrac(StepDates[i-1], StepDates[i]); // Length of time interval
            TradeTime[i] = TradeYrFrac[i] + TradeTime[i-1];                  // Time until StepDates[i]
        }

        // step back to find the segment ends
        int currSeg = SegmentEnd.size()-1;
        for (i=NumOfStep; i>0; i--) {
            if (StepDates[i].equals(dates[currSeg+1])) {
                SegmentEnd[currSeg] = i;
                currSeg--;
                if (currSeg < 0) {
                    break;
                }
            }
        }

        // set the SegmentStart array
        SegmentStart[0] = 0;
        for (i=1; i<(int)SegmentStart.size(); i++) {
            SegmentStart[i] = SegmentEnd[i-1];
        }

        // to do: EffectiveStep = MAX(for all segments, (segment end date - val date)/segment step size)
        // Should be done in tree as not needed for FD

        return NumOfStep;

    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

int CTimeLine::CreateTimeLineSimple(const DateTimeArray& dates, 
                              TimeMetricConstSP metric, int stepsPerYear,                              
                              DateTimeArray critDates, IntArray* isAddedSeg, int* minStepInAddedSeg)
{
    static string routine("CTimeLine::CreateTimeLineSimple");
    try {

        int i, k;
        double dt,  mat_trading;
        vector<double> trade_t, var;

        vector<int> stepsInSeg;
        vector<double> seg_trading;
        vector<int> segmentStartNoCrit;
        double segStartTime;

        int totalSteps;

        int num_seg = dates.size() - 1;
        if (num_seg < 1)
            throw ModelException(routine, "must have at least two dates for generating a time line.");

    //    if (num_seg != (int) density.size())
//            throw ModelException(routine, "num of density factor must be numOfDates -1.");

        Clear();

        // sort critical dates
        metric->SortDate(dates[0], dates[num_seg], true, critDates);

        mat_trading = metric->yearFrac(dates[0], dates[num_seg]);

        // allow varying density segments
        SegmentEnd.resize(num_seg); 
        SegmentEndNoCrit.resize(num_seg);
        SegmentStart.resize(num_seg);

        stepsInSeg.resize(num_seg);
        seg_trading.resize(num_seg);
        segmentStartNoCrit.resize(num_seg);



        totalSteps = 0;
        for(i=0; i<num_seg; i++) {
            seg_trading[i] = metric->yearFrac(dates[i], dates[i+1]);
            //we need to modify this line if we wants dt depends on var.
			stepsInSeg[i] = (int)(stepsPerYear*seg_trading[i]);


			//make sure step is no 0
            stepsInSeg[i] = Maths::max(stepsInSeg[i], 3);
            stepsInSeg[i] = Maths::min(stepsInSeg[i], MAX_SEG_STEP);

			//need to tes 3 and 15
			if (isAddedSeg ==0){
				//don't need to do special for seg
			}else{
				if ((*isAddedSeg)[0] != -1){
					double temp = 1.0/250.0;
					if ((seg_trading[i] <= temp) && ((*isAddedSeg)[i] == 0) ){ //is an added seg						
						stepsInSeg[i] = Maths::max(stepsInSeg[i], (*minStepInAddedSeg)); 
						//stepsInSeg[i] = Maths::max(stepsInSeg[i], 15); 
						stepsInSeg[i] = Maths::min(stepsInSeg[i], MAX_SEG_STEP);
					}
				}
			}

            totalSteps += stepsInSeg[i];
            segmentStartNoCrit[i] = (i==0 ? 0 : SegmentEndNoCrit[i-1]);
            SegmentEndNoCrit[i] = segmentStartNoCrit[i] + stepsInSeg[i];        
        }


        // if there's only one segment put bounds on the number of steps
        //don't think that we need it for FD
//        if (num_seg == 1) {
//            stepsInSeg[0] = (int) Maths::max(stepsInSeg[0], MIN_TREE_STEP);
//            // cap steps
//            stepsInSeg[0] = (int) Maths::min(stepsInSeg[0], MAX_TIME_STEP);
//            totalSteps = stepsInSeg[0];
//            SegmentEndNoCrit[0] = stepsInSeg[0];
//        }

        trade_t.resize(totalSteps + 1);
        var.resize(totalSteps + 1);

        trade_t[0] = 0.0;
        var[0] = 0.0;

        for (k=0; k<num_seg; k++) {
            dt = seg_trading[k]/(double)stepsInSeg[k];

            segStartTime = trade_t[segmentStartNoCrit[k]];
            // set up time points
            //for now, set dt as constant within seg.
            for (i=segmentStartNoCrit[k]+1; i<=SegmentEndNoCrit[k]; i++){
                trade_t[i] = segStartTime + (i - segmentStartNoCrit[k]) * dt;
            }
        }
        // calculate dates for the steps
        CalcStepDate(dates[0], dates[num_seg], trade_t, totalSteps, metric, dates);

        // insert critical dates
        InsertCritDates(trade_t, mat_trading, totalSteps, metric, critDates);

        // calculate TradeYrFrac and CalTime array
        TradeYrFrac = new double[NumOfStep+1];  
        TradeTime = new double[NumOfStep+1];    
        TradeYrFrac[0] = TradeTime[0] = 0.0;
        for (i=1; i<=NumOfStep; i++)
        {
            TradeYrFrac[i] = metric->yearFrac(StepDates[i-1], StepDates[i]); // Length of time interval
            TradeTime[i] = TradeYrFrac[i] + TradeTime[i-1];                  // Time until StepDates[i]
        }

        // step back to find the segment ends
        int currSeg = SegmentEnd.size()-1;
        for (i=NumOfStep; i>0; i--) {
            if (StepDates[i].equals(dates[currSeg+1])) {
                SegmentEnd[currSeg] = i;
                currSeg--;
                if (currSeg < 0) {
                    break;
                }
            }
        }

        // set the SegmentStart array
        SegmentStart[0] = 0;
        for (i=1; i<(int)SegmentStart.size(); i++) {
            SegmentStart[i] = SegmentEnd[i-1];
        }

        // to do: EffectiveStep = MAX(for all segments, (segment end date - val date)/segment step size)
        // Should be done in tree as not needed for FD

        return NumOfStep;

    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

////////////////////////////////////////////////////////
// CalcStepDate
////////////////////////////////////////////////////////
void CTimeLine::CalcStepDate(const DateTime& startDate, const DateTime& endDate,
                             const vector<double>& trade_t, int steps, TimeMetricConstSP metric,
                             const DateTimeArray& dates)
{
    static const string method = "CTimeLine::CalcStepDate";
    static const double YEAR_FRAC_TOL = 0.1/365.0;    // The tolerance used in verifying time line dates (0.1 days)
    try {
        int i;
        double yr;

        StepDates.resize(steps+1);
        StepDates[0] = startDate;

        int currSeg = 0;
        for (i=1; i<=steps; i++)
        {
            if (currSeg < (int)SegmentEndNoCrit.size()-1 && i == SegmentEndNoCrit[currSeg]) {
                StepDates[i] = dates[currSeg+1];
                currSeg++;
            } else {
                StepDates[i] = metric->impliedTime(startDate, trade_t[i], yr);

                if (!(Maths::areEqualWithinTol(trade_t[i], yr, YEAR_FRAC_TOL)) ||    // within absolute tolerance
                    !(Maths::areEqualWithinTol(trade_t[i], yr, 0.01*(trade_t[i]-trade_t[i-1])))) {    // within 1% of step size
                    // Will only happen if bug in impliedTime()
                    throw ModelException(method, "Error in date calculated by TimeMetric::impliedTime()");
                }
            }
        }
                    
        // Verify that computed end date of time line is equal to maturity date
        if (!(Maths::areEqualWithinTol(metric->yearFrac(StepDates[steps], endDate), 0.0, YEAR_FRAC_TOL))) {
            throw ModelException(method, "Calculated end date of time line (" + StepDates[steps].toString() + 
                                         ") is different to instrument maturity date (" +
                                         endDate.toString() + ").");
        }
        // assign end date to the exact date/time
        StepDates[steps] = endDate;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

////////////////////////////////////////////////////////
// InsertCritDates
////////////////////////////////////////////////////////
void CTimeLine::InsertCritDates(vector<double>& trade_t, 
                                double mat_trading, int step, TimeMetricConstSP metric,
                                const DateTimeArray& critDates)
{
    static const string method = "CTimeLine::InsertCritDates";
    try {

        int i, j, k;
        vector<double> crit_trade_time;

        crit_trade_time.resize(critDates.size());
        for(i=0; i<(int)critDates.size(); i++){
            crit_trade_time[i] = metric->yearFrac(i==0?
                                                  StepDates[0]: critDates[i-1],
                                                  critDates[i]);
            if (i > 0){
                crit_trade_time[i] += crit_trade_time[i-1];
            }
        }
        // insert critical dates
        k = 0;
        if (critDates.size() > 0)
        {
            j = 1;
            for (i=0; i<(int)critDates.size(); i++)
            {
                while (StepDates[j] < critDates[i]) j++; // locate nearby step

                // if the dates are different, look whether to insert the crit date
                if (StepDates[j] != critDates[i])
                {

                    // if they are really different ie Trading Time is more than zero, insert it
                    if (metric->yearFrac(critDates[i],StepDates[j]) > 0.0) {
                        StepDates.insert(StepDates.begin()+j, critDates[i]);
                        trade_t.insert(trade_t.begin()+j, crit_trade_time[i]);

                        k++; // num of critical points inserted
                    } else {

                    //they are different but Trading Time diff is zero (ie Monday EOD -> Tuesday SOD)
                    //keep the critical date rather than the tree one.
                        StepDates[j] = critDates[i];
                    }
                }
            }
        }
        NumOfStep = k+step; // Total num of steps in time line
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// change tree base step
void CTimeLine::SetTreeStepBase(int stepBase)
{
    TreeStepBase = stepBase;
}

bool CTimeLine::IsSegEnd(int j)
{
    int i;

    for (i=0; i<(int)SegmentEndNoCrit.size(); i++) {
        if (j == SegmentEndNoCrit[i]) {
            return true;
        }
    }
    return false;
}

DRLIB_END_NAMESPACE
