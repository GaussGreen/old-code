//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FLNResettable.cpp
//
//   Description : One factor finite difference base class for log-normal processes (Resettable)
//
//   Author      : André Segger
//
//   Date        : 05 June 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FLNResettable.hpp"

DRLIB_BEGIN_NAMESPACE

FD1FLNResettable::FD1FLNResettable(): FD1FLN(TYPE), num3DSteps(50) {
     // empty
}

/** Constructor */
FD1FLNResettable::FD1FLNResettable(int stepsPY, int stockSteps): 
    FD1FLN(TYPE), num3DSteps(50) {
    // empty
}

FD1FLNResettable::~FD1FLNResettable()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}


/** get processed vol */
void FD1FLNResettable::InitVol()
{
    FD1FLN::InitVol();
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FLNResettable::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
                                  const DateTime& matDate, CTermStructure& v2)
{
    FD1FLN::CalcV2Term(valDate, startDate, matDate, v2);
}

/** set up variance array */
void FD1FLNResettable::PostSetup()
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

TimeMetricConstSP FD1FLNResettable::GetTimeMetric() const
{
    if (!VolLN)
        throw ModelException("FD1FLNResettable::GetTimeMetric", "Volatility not initialised.");
    
    return VolLN->GetTimeMetric();
}

int FD1FLNResettable::GetStepVol(int step, vector<double>& vol, const double*, int, int end)
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

CVolProcessedBSConstSP FD1FLNResettable::getProcessedVol()
{
    return VolLN;
}

class FD1FLNResettableHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FLNResettable, clazz);
        SUPERCLASS(FD1FLN);
        FIELD(num3DSteps,          "number of steps in strike dimension");
        FIELD_MAKE_OPTIONAL(num3DSteps);
        EMPTY_SHELL_METHOD(defaultFD1FLNResettable);
    }
    
    static IObject* defaultFD1FLNResettable(){
        return new FD1FLNResettable();
    }
};


CClassConstSP const FD1FLNResettable::TYPE = CClass::registerClassLoadMethod(
    "FD1FLNResettable", typeid(FD1FLNResettable), FD1FLNResettableHelper::load);

DRLIB_END_NAMESPACE
