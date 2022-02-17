//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FDDE.cpp
//
//   Description : One factor finite difference base class for DDE
//
//   Author      : Qing Hou
//
//   Date        : Nov 21, 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FDDE.hpp"
#include "edginc/AssetDDE.hpp"

DRLIB_BEGIN_NAMESPACE

FD1FDDE::FD1FDDE(): FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()), 
                    ddeType("RG"), calibParams(new DDECalib()), defProbStep(0)
{
    // should really have a constructor in FD1F
    product = 0;
    
    useFwdGrid = false;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = true;

    varMethod = 0;

    // default exponential grid. the sinh grid is problematic 
    // for high vol/long horizon case, or low strike case
    gridType = 2;
}

    /** Less simple constructor */
FD1FDDE::FD1FDDE(int stepsPY, int stockSteps): 
    FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()), 
    ddeType("RG"), calibParams(new DDECalib()), defProbStep(0)
{
    product = 0;
    
    useFwdGrid = false;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    varMethod = 0;
    gridType = 0;

    this->stepsPerYear    = stepsPY;
    this->stockSteps      = stockSteps;
    this->TruncationStd   = 7.0;

}

FD1FDDE::FD1FDDE(CClassConstSP clazz): FD1FGeneric(clazz),
                                       volType(VolSurface::TYPE->getName()), 
                                       ddeType("RG"),
                                       calibParams(new DDECalib()), defProbStep(0)
{
    product = 0;

    useFwdGrid = false;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    varMethod = 0;
    gridType = 0;

    this->stepsPerYear    = -1;
    this->stockSteps      = -1;
    this->TruncationStd   = 7.0;
}


FD1FDDE::~FD1FDDE()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}


/** get processed vol and prepare the asset for DDE calculations */
// need to examine if this is general
void FD1FDDE::InitVol()
{
	if( !AssetDDE::TYPE->isInstance(Underlier.get()) )
		throw(ModelException("FD1FDDE::InitVol", "Underlier must be of type AssetDDE"));
	const AssetDDE *asset = (dynamic_cast<const AssetDDE *> (Underlier.get()));

	const IDDEInitiator *initiator = dynamic_cast<const IDDEInitiator *>(product);
	if( !initiator )
		throw(ModelException("FD1FDDE::InitVol", "Product does not support DDEInitiator"));
	
	if( !ddeModule )
	{
		ddeModule = DDEModule::createModule(ddeType, calibParams.get(), asset, initiator);
	}
	else
		ddeModule->update(asset, initiator);

    VolLN = CVolProcessedBSSP::dynamicCast(
        CVolProcessedSP(ddeModule->getProcessedVol(product->GetLNRequest().get())));
		
	// get credit spreads function
	csFunc = SpreadEquityFuncConstSP(ddeModule->getSpreadFunc());
}

const CleanSpreadCurve *FD1FDDE::getNonEqSpreads() const
{
	return ddeModule->getNonEqSpreads();
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FDDE::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
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
void FD1FDDE::PostSetup()
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

void FD1FDDE::recordOutputRequests(Control* control, CResults* results)
{
	ddeModule->recordOutputRequests(control, results);
}


TimeMetricConstSP FD1FDDE::GetTimeMetric() const
{
    if (!VolLN)
        throw ModelException("FD1FDDE::GetTimeMetric", "Volatility not initialised.");
    
    return VolLN->GetTimeMetric();
}

int FD1FDDE::GetStepVol(int step, vector<double>& vol, const double*, int, int end)
{
    if (int(vol.size()) != (end+1))
         vol.resize(end+1);

    if (step >= TimePts.NumOfStep)
        step =TimePts.NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    
    if (TimePts.TradeYrFrac[step+1] != 0.) {
        vol[0] = (Variance[step+1]-Variance[step])/TimePts.TradeYrFrac[step+1];
        vol[0] = sqrt(vol[0]);
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

CVolProcessedBSConstSP FD1FDDE::getProcessedVol()
{
    return VolLN;
}

FDTermStructureSP FD1FDDE::getDriftTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                              const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP driftTerm;

    if ( useFwdGrid ) {
        throw ModelException("FD1FDDE::getDriftTerm", "UseFwdGrid=TRUE not allowed.");
    } else {
        DoubleArraySP termArray(new DoubleArray(end-start+1));
		DoubleArray	csArray(end - start + 1);
		double dt = times[step];
		csFunc->getSpreadCC(TimePts.StepDates[step-1], TimePts.StepDates[step], end - start + 1, &s[start], &csArray.front(), dt);
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = ((inIr[step]+irPert)-(inDivy[step]+divPert) + csArray[i]) * s[i];
        }
        driftTerm = FDTermStructureSP(new FDTermStructure(termArray));
    }
    return driftTerm;
}

FDTermStructureSP FD1FDDE::getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGrid,
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

FDTermStructureSP FD1FDDE::getCouponTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                               const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP couponTerm;
    couponTerm = FDTermStructureSP(new FDTermStructure(0.0));
    return couponTerm;
}

FDTermStructureSP FD1FDDE::getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                                 const double &irPert, const double &divPert, const bool doEquityLayer)
{
    FDTermStructureSP discountTerm;
    if ( useFwdGrid ) {
        throw ModelException("FD1FDDE::getDiscountTerm", "UseFwdGrid=TRUE not allowed.");
    } else {
		// may want to somehow cache the credit spreads
        DoubleArraySP	termArray(new DoubleArray(end-start+1));
		DoubleArray		csArray(end - start + 1);
		double dt = times[step];
		csFunc->getSpreadCC(TimePts.StepDates[step-1], TimePts.StepDates[step], end - start + 1, &s[start], &csArray.front(), dt);
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = inIr[step]+irPert + csArray[i];

        }
        discountTerm = FDTermStructureSP(new FDTermStructure(termArray));
    }
    return discountTerm;
}

// initialize the default prob for the next step
void FD1FDDE::setDefaultProb(int step, double* defProbs, int start, int end)
{
    ASSERT( step > 0 );
    double irFwd = 1./DiscountCurve->pv(TimePts.StepDates[step-1], TimePts.StepDates[step]);
    for(int i=start; i<=end; i++) defProbs[i] = irFwd;
    defProbStep = step;
}

void FD1FDDE::adjustPriceByDefPO(int nextStep, double nextDefPayoff, double *defProbs, double *price, int start, int end)
{
    ASSERT( nextStep == defProbStep );
    double irPv = DiscountCurve->pv(TimePts.StepDates[nextStep-1], TimePts.StepDates[nextStep]);
    for(int i=start; i<=end; i++)
        price[i] += ( 1 - defProbs[i] ) * nextDefPayoff * irPv;
}

MarketDataFetcherSP FD1FDDE::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}


class FD1FDDEHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FDDE, clazz);
        SUPERCLASS(FD1FGeneric);
        IMPLEMENTS(CAsset::ICanHaveDDE);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(ddeType, "Type of DDE model: RG, LN1F, LN2F");
        FIELD_MAKE_OPTIONAL(ddeType);
        FIELD(calibParams,    "Calibration method and param");
		FIELD_MAKE_OPTIONAL(calibParams);
		FIELD(ddeModule, "");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(ddeModule);
        EMPTY_SHELL_METHOD(defaultFD1FDDE);
    }
    
    static IObject* defaultFD1FDDE(){
        return new FD1FDDE();
    }
};


CClassConstSP const FD1FDDE::TYPE = CClass::registerClassLoadMethod(
    "FD1FDDE", typeid(FD1FDDE), FD1FDDEHelper::load);

DRLIB_END_NAMESPACE
