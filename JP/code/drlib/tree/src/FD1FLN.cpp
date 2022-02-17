//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FLN.cpp
//
//   Description : One factor finite difference base class for log-normal processes
//
//   Author      : André Segger
//
//   Date        : November 20, 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FLN.hpp"

DRLIB_BEGIN_NAMESPACE

FD1FLN::FD1FLN(): FD1F(TYPE), volType(VolSurface::TYPE->getName()) {
    // should really have a constructor in FD1F
    product = 0;
    
    useFwdGrid = true;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = true;

    varMethod = 0;

    gridType = 0;
}

    /** Less simple constructor */
FD1FLN::FD1FLN(int stepsPY, int stockSteps): 
    FD1F(TYPE), volType(VolSurface::TYPE->getName()) {
    product = 0;
    
    useFwdGrid = true;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    varMethod = 0;
    gridType = 0;

    this->stepsPerYear    = stepsPY;
    this->stockSteps      = stockSteps;
    this->TruncationStd   = 7.0;

}

FD1FLN::FD1FLN(CClassConstSP clazz): FD1F(clazz), volType(VolSurface::TYPE->getName())
{
    product = 0;

    useFwdGrid = true;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    varMethod = 0;
    gridType = 0;

    this->stepsPerYear    = -1;
    this->stockSteps      = -1;
    this->TruncationStd   = 7.0;
}


FD1FLN::~FD1FLN()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}

/** get processed vol */
// need to examine if this is general
void FD1FLN::InitVol()
{
    VolLN = CVolProcessedBSSP::dynamicCast(
        CVolProcessedSP(Underlier->getProcessedVol(product->GetLNRequest().get())));
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FLN::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
                        const DateTime& matDate, CTermStructure& v2)
{
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
    
    VolLN->CalcVol(startDate, benchMarkTrunc, CVolProcessedBS::fromFirst, vol_sq);
    
    for (i=0; i<vol_sq.size(); i++)
    {
        if (vol_sq[i] <= 0.0)
            throw ModelException("CTree1fLN::CalcV2Term", "received vol<=0 : cannot create tree with such vol");
        vol_sq[i] *= vol_sq[i];
        yrs[i] = VolLN->calcTradingTime(startDate, benchMarkTrunc[i]);
    }
    
    // populate it to a vol^2 term structure
    v2.Populate(startDate, vol_sq.size(), benchMarkTrunc.begin(),
                yrs.begin(), vol_sq.begin());
}

/** set up variance array */
void FD1FLN::PostSetup()
{
    Variance.resize(TimePts.NumOfStep);
    VolLN->CalcVar(TimePts.StepDates, CVolProcessedBS::fromFirst, Variance);
    Variance.insert(Variance.begin(), 0.0); // start with 0

    stockMaxSeg.resize(TimePts.SegmentEnd.size());
    stockMinSeg.resize(TimePts.SegmentEnd.size());
    strikeSeg.resize(TimePts.SegmentEnd.size());
    int i;
    for (i=0; i<int(stockMaxSeg.size()); i++) {
        stockMaxSeg[i] = -1.;
        stockMinSeg[i] = -1.;
        strikeSeg[i] = -1.;
    }
}

TimeMetricConstSP FD1FLN::GetTimeMetric() const
{
    if (!VolLN)
        throw ModelException("FD1FLN::GetTimeMetric", "Volatility not initialised.");
    
    return VolLN->GetTimeMetric();
}

int FD1FLN::GetStepVol(int step, vector<double>& vol, const double*, int, int end)
{
    if (int(vol.size()) != (end+1))
         vol.resize(end+1);

    if (step >= TimePts.NumOfStep)
        step =TimePts.NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    
    if (TimePts.TradeYrFrac[step+1] != 0.) {
        vol[0] = (Variance[step+1]-Variance[step])/TimePts.TradeYrFrac[step+1];
        vol[0] = ::sqrt(vol[0]);
    } else { // if there is no trading time between steps return a one day vol to avoid a div by zero
        CDateTimeArray dateAndNextDay;
        dateAndNextDay.resize(2);
        dateAndNextDay[0] = TimePts.StepDates[step];
        dateAndNextDay[1] = dateAndNextDay[0].rollDate(1);
        CDoubleArray volTmp;
        volTmp.resize(1);
        VolLN->CalcVol(dateAndNextDay, CVolProcessedBS::forward, volTmp);
        vol[0] = volTmp[0];
    }

	for(int i=1; i <= end; i++)
		vol[i] = vol[0];

    return 1;
}

CVolProcessedBSConstSP FD1FLN::getProcessedVol()
{
    return VolLN;
}


MarketDataFetcherSP FD1FLN::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}

class FD1FLNHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FLN, clazz);
        SUPERCLASS(FD1F);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        EMPTY_SHELL_METHOD(defaultFD1FLN);
    }
    
    static IObject* defaultFD1FLN(){
        return new FD1FLN();
    }
};


CClassConstSP const FD1FLN::TYPE = CClass::registerClassLoadMethod(
    "FD1FLN", typeid(FD1FLN), FD1FLNHelper::load);

DRLIB_END_NAMESPACE
