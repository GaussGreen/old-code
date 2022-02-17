//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1FLV.cpp
//
//   Description : One factor finite difference base class for local vol processes
//
//   Author      : Oleg Divinskiy
//
//   Date        : April 23 2002
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FLV.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/VolSpline.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"


DRLIB_BEGIN_NAMESPACE

FD1FLV::FD1FLV(): FD1F(TYPE), volType(VolSurface::TYPE->getName()) {
    // should really have a constructor in FD1F
    product = 0;
    
    useFwdGrid = false;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    gridType = 0;
}

FD1FLV::~FD1FLV()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}

/** get processed vol */
// need to examine if this is general
void FD1FLV::InitVol()
{
    DateTime startDate = product->GetFwdStartDateLV();

    LocVolRequestSP volRequest(new LocVolRequest(startDate,product->GetFwdStartLV(),false,true,false, 0.005,0.001,0.01));

    CAssetConstSP plainAsset = Underlier;

    if (StruckEquity::TYPE->isInstance(Underlier) && !(AssetUtil::isBasket(Underlier)))
    {
        throw ModelException("FD1FLV::InitVol", "Struck assets are not supported yet");
        /*
        plainAsset = StruckEquityConstSP::dynamicCast(Underlier)->getPlainAsset();
        fxAsset = StruckEquityConstSP::dynamicCast(Underlier)->getFX();
        eqFXCorr = StruckEquityConstSP::dynamicCast(Underlier)->getCorrelation();

        if (!fxAsset){
            throw ModelException("CTree1fLV::InitVol", "NULL fx asset");
        }
        if (!eqFXCorr){
            throw ModelException("CTree1fLV::InitVol", "NULL correlation");
        }
        // create an atm vol request for the fx side of things
        CVolRequestSP atmInterp(new ATMVolRequest());
        CVolProcessedSP interpFXVol(fxAsset->getProcessedVol(atmInterp.get()));
        // cast to the type of vol we need
        volFXBS = CVolProcessedBSSP::dynamicCast(interpFXVol);
        */
    }
    else if (ProtEquity::TYPE->isInstance(Underlier) && !(AssetUtil::isBasket(Underlier)))
    {
        throw ModelException("FD1FLV::InitVol", "Currency protected assets are not supported yet");
        /*
        plainAsset = ProtEquityConstSP::dynamicCast(Underlier)->getPlainAsset();
        FXVolBaseWrapper fxVol = ProtEquityConstSP::dynamicCast(Underlier)->getFXVol();
        eqFXCorr = ProtEquityConstSP::dynamicCast(Underlier)->getCorrelation();
        // interpolate the fx vols atm
        ATMVolRequestSP fxVolRequest(new ATMVolRequest());
        // interpolate the vol
        CVolProcessedSP  fxVoltmp(fxVol->getProcessedVol(fxVolRequest.get(), NULL));
        // cast to the type of vol we're expecting
        volFXBS = CVolProcessedBSSP::dynamicCast(fxVoltmp);
        Underlier = plainAsset; // assign it to plain asset for the tree
        */
    }

    VolLV = CVolProcessedDVFSP::dynamicCast(
        CVolProcessedSP(plainAsset->getProcessedVol(volRequest.get())));
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FLV::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
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
    
    // *** this is only for setting up time line, using ATM vol should be ok ?
    double strike = Underlier->getSpot();

    for (i=0; i<vol_sq.size(); i++)
    {
         vol_sq[i] = VolLV->computeImpVol(benchMarkTrunc[i], strike);
         if (vol_sq[i] <= 0.0)
             throw ModelException("FD1FLV::CalcV2Term", "received vol<=0 : cannot create tree with such vol");
         vol_sq[i] *= vol_sq[i];
         yrs[i] = VolLV->calcTradingTime(valDate, benchMarkTrunc[i]); // need review for fwd start !!!
    }

    // populate it to a vol^2 term structure
    v2.Populate(startDate, vol_sq.size(), benchMarkTrunc.begin(), 
                yrs.begin(), vol_sq.begin());
}


/** set up variance array */
void FD1FLV::PostSetup()
{
    Variance.resize(TimePts.NumOfStep+1);
    // *** requiring ATM vol should be ok, Variance is used for setting up nodes only
    double strike = Underlier->getSpot();

    DateTime startDate = product->GetFwdStartDateLV();
    strike = Underlier->fwdValue(startDate);

    Variance[0] = 0.0; // start with 0
	int i;
    for (i=1; i<=TimePts.NumOfStep; i++)
    {
        Variance[i] = VolLV->computeImpVol(TimePts.StepDates[i], strike);
        if (Variance[i] <= 0.0)
            throw ModelException("CTree1fLV::PostSetup", "received vol<=0 : cannot create tree with such vol");
        Variance[i] *= Variance[i] * TimePts.TradeTime[i];
    }

    stockMaxSeg.resize(TimePts.SegmentEnd.size());
    stockMinSeg.resize(TimePts.SegmentEnd.size());
    strikeSeg.resize(TimePts.SegmentEnd.size());

    for (i=0; i<(int)stockMaxSeg.size(); i++) {
        stockMaxSeg[i] = -1.;
        stockMinSeg[i] = -1.;
        strikeSeg[i] = -1.;
    }
}

TimeMetricConstSP FD1FLV::GetTimeMetric() const
{
    if (!VolLV)
        throw ModelException("FD1FLV::GetTimeMetric", "Volatility not initialised.");
    
    return VolLV->GetTimeMetric();
}

int FD1FLV::GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end)
{

	if ((int)vol.size() != (end-start+1))
        vol.resize(end-start+1);

    if (step >= TimePts.NumOfStep)
        step =TimePts.NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;

    DateTimeArray t(2);
    t[0] = TimePts.StepDates[step];
    // go from now to next timepoint or to tomorrow whichever is further.
    // Local vol is not very accurate if you ask for it over too short an interval.
    if (TimePts.StepDates[step].rollDate(1) > TimePts.StepDates[step+1]) {
        t[1] = TimePts.StepDates[step].rollDate(1);
    } else {
        t[1] = TimePts.StepDates[step+1];
    }

    /*
    // flattern <=0 spots here
    double* s = const_cast<double*>(s_inp);
    for (int n_zero=start; n_zero<=end; n_zero++)
    {
         if (s[n_zero] < FD_MIN)
             s[n_zero] = FD_MIN;
    }

    VolLV->CalcLocVol(t,(s+start),(end - start + 1), vol);

    return (end-start+1);
    */

    if(false /* Payoff->getCcyTreatment()=="S" */ )
    {// struck
        /*
        vector<double> stmp;
        stmp.resize(end-start+1);
        double* s = stmp.begin();
        for (int n_zero=0; n_zero<=end - start; n_zero++)
        {	
            s[n_zero] = s_inp[n_zero + start];
            s[n_zero] = s[n_zero]/fxAsset->fwdValue(t[0]);
            if (s[n_zero] < FD_MIN)
                s[n_zero] = FD_MIN;
        }
        CSliceDouble spots(s, end - start + 1);
        CSliceDouble locVol(vol.begin(), end - start + 1);
        VolLV->CalcLocVol(&spots, t, &locVol);

        double vol_fx = volFXBS->CalcVol(t[0], t[1]);
        double corr = eqFXCorr->getCorrelation();   
        for (int i=0; i<=end - start; i++)
            vol[i] = sqrt(vol[i]*vol[i] + vol_fx*vol_fx + 2.0*corr*vol[i]*vol_fx);
        */
    }
    else
    {
        double* s = const_cast<double*>(s_inp);
        for (int n_zero=start; n_zero<=end; n_zero++)
        {			   
            if (s[n_zero] < FD_MIN)
                s[n_zero] = FD_MIN;
        }
        CSliceDouble spots(s+start, end - start + 1);
        CSliceDouble locVol(&vol[0], end - start + 1);
        VolLV->CalcLocVol(&spots, t, &locVol);

    }

    return (end-start+1);
}

MarketDataFetcherSP FD1FLV::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherLNSpline(volType));
}

class FD1FLVHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FLV, clazz);
        SUPERCLASS(FD1F);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        EMPTY_SHELL_METHOD(defaultFD1FLV);
    }
    
    static IObject* defaultFD1FLV(){
        return new FD1FLV();
    }
};


CClassConstSP const FD1FLV::TYPE = CClass::registerClassLoadMethod(
    "FD1FLV", typeid(FD1FLV), FD1FLVHelper::load);

DRLIB_END_NAMESPACE
