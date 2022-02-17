//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Tree1fLN.cpp
//
//   Description : one factor trinomial tree for log-normal process.
//
//   Author      : Ning Shen
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/LatticeProdEDR.hpp"

DRLIB_BEGIN_NAMESPACE

CTree1fLN::CTree1fLN():CTree1f(TYPE), 
                       volType(VolSurface::TYPE->getName()){}

/** get processed vol*/
// need to examine if this is general
void CTree1fLN::InitVol()
{
    //init volprodLN
    const IFDProductLN*  prodLN = dynamic_cast<const IFDProductLN*>(prod.get());

    if (prodLN){
        CVolRequestLNSP volRequests = prodLN->getVolInterp(0);    // [numAssets] 

        // interpolate the vol using our LN request
        //should retrieve vol from Underlier, not from mAsset since Underlier maybe changed!
        VolLN = CVolProcessedBSSP(Underlier->getProcessedVol(volRequests.get()));
    }else{
        //should retrieve vol from Underlier, not from mAsset since Underlier maybe changed!
        VolLN = CVolProcessedBSSP::dynamicCast(
            CVolProcessedSP(Underlier->getProcessedVol(latticeProd->GetLNRequest().get())));
    }        
    // get time metric
    timeMetric = VolLN->GetTimeMetric();
}

IModel::WantsRiskMapping CTree1fLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** calculate a term structure vol^2 or just one point at maturity */
void CTree1fLN::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
                           const DateTime& matDate, CTermStructure& v2)
{
    static const string method = "CTree1fLN::CalcV2Term";
    try {
        int i;
        DateTimeArray benchMark;
        DoubleArray vol_sq;
        // a simple but not most efficient way of creating vol term structure
        // access to vol term structure denied

        // write down the bench marks
        const int benchMarkSize = 19;
        string benchMarkString[benchMarkSize] = {"1D", "1W", "2W", "3W", "1M", "2M", "3M", "6M", "9M", 
                                                 "1Y", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "10Y", "20Y", "30Y"};
        // create bench mark dates
        DateTime valDatePM(valDate.getDate(), DateTime::END_OF_DAY_TIME);
        benchMark.resize(benchMarkSize*3-2);
        benchMark[0] = MaturityPeriod(benchMarkString[0]).toDate(valDatePM);
        for (i=1; i<benchMarkSize; i++)
        {
            benchMark[3*i-1] = MaturityPeriod(benchMarkString[i]).toDate(valDatePM);
            benchMark[3*i-2] = benchMark[3*i-1].rollDate(-1); // we won't be wrong by 1 day !
            benchMark[3*i] = benchMark[3*i-1].rollDate(1);
        }
        // use dates from start date onwards
        DateTimeArray benchMarkTrunc;
        for (i=0; i<benchMarkSize*3-2; i++)
        {
            if (benchMark[i] > startDate)
            {
                benchMarkTrunc.push_back(benchMark[i]);
            }
        }

        // retrieve variance for each date from the processed vol object
        vector<double> yrs;
        vol_sq.resize(benchMarkTrunc.size());
        yrs.resize(benchMarkTrunc.size());
        VolLN->CalcVar(startDate, benchMarkTrunc, 
                       CVolProcessedBS::fromFirst, vol_sq);

        for (i=0; i<vol_sq.size(); i++)
        {
            // for performance, calc trading time between dates
            yrs[i] = VolLN->calcTradingTime(i==0? 
                                            startDate: benchMarkTrunc[i-1],
                                            benchMarkTrunc[i]);
            if (i > 0){
                yrs[i] += yrs[i-1];
            }
            if (!Maths::isZero(yrs[i])){
                vol_sq[i] /= yrs[i]; // convert from variance to vol squared
            }
        }
        // populate it to a vol^2 term structure
        v2.Populate(startDate, vol_sq.size(), benchMarkTrunc.begin(),
                    yrs.begin(), vol_sq.begin());
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** set up variance array */
void CTree1fLN::PostSetup()
{
    Variance.resize(timeLine->NumOfStep);
    VolLN->CalcVar(timeLine->StepDates, CVolProcessedBS::fromFirst, Variance);
    Variance.insert(Variance.begin(), 0.0); // start with 0
}

/** calculate drift and variance in unit of fwd (ie. exp(v_dt)-1 for LN) for current step
    they must be either one per step or one per node (from -BotClip to TopClip)
    returns true if constant probabilities can be used for all nodes (fast roll) */
bool CTree1fLN::CalcStepDriftAndVar(const double* s, int start, int end, vector<double>& var_dt, vector<double>* drift)
{
    static const string method = "CTree1fLN::CalcStepDriftAndVar";
    try {
        ASSERT(CurrStep < timeLine->NumOfStep);

        // just one var per step
        var_dt.resize(1);
        var_dt[0] = exp(Variance[CurrStep+1] - Variance[CurrStep]) - 1.0;

        // fwd
        bool fastRoll = (RollDirection == -1);

        if (drift)
        {
            drift->resize(1);
            (*drift)[0] = StepForward[CurrStep+1]/StepForward[CurrStep];
        }

        // calculate prob for the step
        double mu = StepForward[CurrStep+1]/StepForward[CurrStep]; // step drift
        double m = Stock[1-CurrIdx][0]/Stock[CurrIdx][0]/mu - 1.0;
        double u = Stock[1-CurrIdx][1]/Stock[CurrIdx][0]/mu - 1.0;
        double d = Stock[1-CurrIdx][-1]/Stock[CurrIdx][0]/mu - 1.0;
        Pu = (var_dt[0] + d*m)/(u-d)/(u-m);
        Pd = (var_dt[0] + u*m)/(d-u)/(d-m);
        Pm = 1.0 - Pu - Pd;

        // protect rounding for ~0 vol interval
        if (var_dt[0]<1.0e-8)
        {
            if (Pu < 0.0 && Pu>-1.0e-10) Pu = 0.0;
            if (Pm < 0.0 && Pm>-1.0e-10) Pm = 0.0;
            if (Pd < 0.0 && Pd>-1.0e-10) Pd = 0.0;
        }
        if(Pu<0.0 || Pd<0.0 || Pm<0.0)
            fastRoll = false;     // -ve prob, needs per node recalc

        if (!fastRoll)
            PerNodeProb = true;

        return fastRoll;
        }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
/** calculate vol for current step for a set of spot levels
    returns number of vol calculated - one (flat for all node) or num */

int CTree1fLN::GetStepVol(int step, vector<double>& vol, const double*, int, int)
{
    if (vol.size() != 1)
        vol.resize(1);

    if (step >= timeLine->NumOfStep)
        step =timeLine->NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    
    if (timeLine->TradeYrFrac[step+1] != 0.) {
        vol[0] = (Variance[step+1]-Variance[step])/timeLine->TradeYrFrac[step+1];
        vol[0] = sqrt(vol[0]);
    } else { // if there is no trading time between steps return a one day vol to avoid a div by zero
        CDateTimeArray dateAndNextDay;
        dateAndNextDay.resize(2);
        dateAndNextDay[0] = timeLine->StepDates[step];
        dateAndNextDay[1] = dateAndNextDay[0].rollDate(1);
        CDoubleArray volTmp;
        volTmp.resize(1);
        VolLN->CalcVol(dateAndNextDay, CVolProcessedBS::forward, volTmp);
        vol[0] = volTmp[0];
    }
    return 1;
}

/** calculate node spacing for each time segment, usually just one segment */
void CTree1fLN::CalcNodeSpacing()
{
    static const string method = "CTree1fLN::CalcNodeSpacing";
    try {
        double v_dt;
        int startStep = 0;
        int startStepNoCrit = 0;
        int i, j, initStep;
        double root1, root2;

        int numOfSeg = timeLine->SegmentEnd.size();
        NodeSpace.resize(numOfSeg);
        for (i=0; i<numOfSeg; i++)
        {
            initStep = timeLine->SegmentEndNoCrit[i] - startStepNoCrit;
            v_dt = (Variance[timeLine->SegmentEnd[i]] - Variance[startStep])/initStep;
            if (Maths::isZero(v_dt)){
                // zero variance across segment so put NodeSpace to some
                // minimum. Any value > 0 will do
                NodeSpace[i] = 0.01; 
            } else {
            
                //NodeSpace[i] = sqrt(v_dt/(TreeAlpha/timeLine->PeakVarRatio)); // simple but does not work for high vol
                
                // search for max var in the segment, really for equal time step use
                // but keep it until equal variance step creation matches CalcVar() exactly
                for (j=startStep; j<timeLine->SegmentEnd[i]; j++)
                {
                    if (v_dt < Variance[j+1] - Variance[j])
                        v_dt = Variance[j+1] - Variance[j];
                }
                //     double m = exp(v_dt/2.); // this is for variance adjusted centre node
                double m =1.0;
                double alpha = 
                    (-exp(v_dt)*m*m + 1-2*TreeAlpha)/(TreeAlpha-1+m);
                SolveQuadratic(1.0, alpha, 1.0, &root1, &root2);
                ASSERT(root2>1.0);
                NodeSpace[i] = log(root2);
            }
            startStep = timeLine->SegmentEnd[i];
            startStepNoCrit = timeLine->SegmentEndNoCrit[i];
        }


        // calculate centre nodes for all steps
        for (i=0; i<=timeLine->NumOfStep; i++)
        {
            // variance adjusted centre nodes give better sampling for prob space 
            // and thus slightly better convergence than fwds centre nodes.
            // CentreNode[i] = StepForward[i]*exp(-Variance[i]/2.0);
 
            // but using fwds as centre nodes tend to range around usual boundaries (strikes etc.)
            // and may be more relvevant ?
            CentreNode[i] = StepForward[i];
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/**   tree node space set up */
void CTree1fLN::NodeSetupSpace(double* s, int bot, int top, 
                          double* sLast, int botLast, int topLast, 
                          double f0, double spacing, int step, double driftNext)
{
    static const string method = "CTree1fLN::NodeSetupSpace";
    try {
        // reset it to make sure default only calc one prob per step
        PerNodeProb = false;

        // spot tree is created with proportional spacing
        double spacingExp = exp(spacing);
        // starting node
        s[bot] = f0*exp(spacing*bot);
    for (int j = bot+1; j<=top; j++)
            s[j] = s[j-1]*spacingExp;

    // set floors and ceilings
    s[bot-2] = 0.0; // lowest floor set to zero
        if (step == timeLine->NumOfStep || (step == 0 && RollDirection==1))
        {
            int  n = Maths::min(step+1, Ndim);
            s[bot-1] = 0.1*s[-n-1];
        s[top+2] = 2*f0*exp(spacing*(Ndim+1)+Variance[timeLine->NumOfStep]/2.0);
        s[top+1] = s[n+1] +0.9*(s[n+3]-s[n+1]);
            if (step == 0)
            {
                double drift = StepForward[timeLine->NumOfStep]/StepForward[0];
                s[bot-1] /= drift;
            s[top+2] /= drift;
            s[top+1] /= drift;
            }
        }
        else
        {
            s[bot-1] = driftNext*sLast[botLast-1];
        s[top+2] = driftNext*sLast[topLast+2];
        s[top+1] = driftNext*sLast[topLast+1];
        }
    ASSERT(s[top+1] >= s[top]);
    ASSERT(s[top+2] >= s[top+1]);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** overwrite control of sameGrid from instrument, mainly for fwdStarting */
void CTree1fLN::controlSameGridFwdStart(const string& ccyTreatment)
{
    CTree1f::controlSameGridFwdStart(ccyTreatment);
    // if it's Quanto, it's not available to use SameGridVega
    if (ccyTreatment==CAsset::CCY_TREATMENT_PROTECTED){
        sameGridVega = false;
    }
}


/** make a Tree1fLN with most things safely defaulted */
CTree1fLN* CTree1fLN::make(const string& volType, bool useDivTreatment) {
    Tree1fLNSP tree(new Tree1fLN());

    tree->volType = volType;
    tree->StepsPerYear = 50;
    tree->DivAmountTreatment = useDivTreatment;
    tree->TreeAlpha = 0.7;

    return tree.release();
}

/** Override default createMDF in order to set the right MDF */
MarketDataFetcherSP CTree1fLN::createMDF() const{

    return MarketDataFetcherSP(new MarketDataFetcherLN(volType, useCcyBasis));
}

class Tree1fLNHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CTree1fLN, clazz);
        SUPERCLASS(CTree1f);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        EMPTY_SHELL_METHOD(defaultTree1fLN);
    }

    static IObject* defaultTree1fLN(){
        return new CTree1fLN();
    }
};

CClassConstSP const CTree1fLN::TYPE = CClass::registerClassLoadMethod(
    "Tree1fLN", typeid(CTree1fLN), Tree1fLNHelper::load);

DRLIB_END_NAMESPACE
