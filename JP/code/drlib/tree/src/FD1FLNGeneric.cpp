//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FLNGeneric.cpp
//
//   Description : One factor generic finite difference base class for log-normal processes
//
//   Author      : André Segger
//
//   Date        : 16 October 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FLNGeneric.hpp"

DRLIB_BEGIN_NAMESPACE

FD1FLNGeneric::FD1FLNGeneric(): 
    FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()) {
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
FD1FLNGeneric::FD1FLNGeneric(int stepsPY, int stockSteps): 
    FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()) {
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

FD1FLNGeneric::FD1FLNGeneric(CClassConstSP clazz): 
    FD1FGeneric(clazz), volType(VolSurface::TYPE->getName())
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


FD1FLNGeneric::~FD1FLNGeneric()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}


/** get processed vol */
// need to examine if this is general
void FD1FLNGeneric::InitVol()
{
    VolLN = CVolProcessedBSSP::dynamicCast(
        CVolProcessedSP(Underlier->getProcessedVol(product->GetLNRequest().get())));
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FLNGeneric::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
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
void FD1FLNGeneric::PostSetup()
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

TimeMetricConstSP FD1FLNGeneric::GetTimeMetric() const
{
    if (!VolLN)
        throw ModelException("FD1FLNGeneric::GetTimeMetric", "Volatility not initialised.");
    
    return VolLN->GetTimeMetric();
}

int FD1FLNGeneric::GetStepVol(int step, vector<double>& vol, const double*, int, int end)
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

CVolProcessedBSConstSP FD1FLNGeneric::getProcessedVol()
{
    return VolLN;
}

FDTermStructureSP FD1FLNGeneric::getDriftTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                              const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP driftTerm;
    // driftTerm = FDTermStructureSP(new FDTermStructure(0.0));

    if ( useFwdGrid ) {
        DoubleArraySP     termArray(new DoubleArray(end-start+1));
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = (irPert-divPert) * s[i];
        }
        driftTerm = FDTermStructureSP(new FDTermStructure(termArray));
    } else {
        DoubleArraySP     termArray(new DoubleArray(end-start+1));
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = ((inIr[step]+irPert)-(inDivy[step]+divPert)) * s[i];

        }
        driftTerm = FDTermStructureSP(new FDTermStructure(termArray));
    }
    return driftTerm;
}

FDTermStructureSP FD1FLNGeneric::getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                                  const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP diffusionTerm;
    vector<double>    vol;
    DoubleArraySP     termArray(new DoubleArray(end-start+1));

    GetStepVol(step-1, vol, s, start, end);

    for (int i=0; i<termArray->size() ; ++i) {
        (*termArray)[i] = .5 * vol[i] * vol[i] * s[i] * s[i];
    }

    diffusionTerm = FDTermStructureSP(new FDTermStructure(termArray));
    return diffusionTerm;
}

FDTermStructureSP FD1FLNGeneric::getCouponTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                               const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP couponTerm;
    couponTerm = FDTermStructureSP(new FDTermStructure(0.0));
    return couponTerm;
}

FDTermStructureSP FD1FLNGeneric::getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                                 const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP discountTerm;
    if ( useFwdGrid ) {
        discountTerm = FDTermStructureSP(new FDTermStructure(irPert));
    } else {
        discountTerm = FDTermStructureSP(new FDTermStructure(inIr[step]+irPert));
    }
    return discountTerm;
}

MarketDataFetcherSP FD1FLNGeneric::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}


class FD1FLNGenericHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FLNGeneric, clazz);
        SUPERCLASS(FD1FGeneric);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        EMPTY_SHELL_METHOD(defaultFD1FLNGeneric);
    }
    
    static IObject* defaultFD1FLNGeneric(){
        return new FD1FLNGeneric();
    }
};


CClassConstSP const FD1FLNGeneric::TYPE = CClass::registerClassLoadMethod(
    "FD1FLNGeneric", typeid(FD1FLNGeneric), FD1FLNGenericHelper::load);

DRLIB_END_NAMESPACE
