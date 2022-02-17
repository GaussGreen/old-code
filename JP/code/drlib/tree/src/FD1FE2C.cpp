//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FE2C.cpp
//
//   Description : One factor finite difference base class for E2C dynamics
//
//   Author      : André Segger
//
//   Date        : 16 October 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FE2C.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/FirmAsset.hpp"
#include "edginc/BadDayConventionFactory.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/DayCountConventionFactory.hpp"

DRLIB_BEGIN_NAMESPACE


FD1FE2C::FD1FE2C(): FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()) {
    // should really have a constructor in FD1F
    product = 0;
    
    useFwdGrid = true;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = true;

    varMethod = 0;

    gridType = 0;

    riskyEquityGrowth     = false;

    useDivsForAssetDrift = false;
    additionalAssetDrift = 0.0;
}

    /** Less simple constructor */
FD1FE2C::FD1FE2C(int stepsPY, int stockSteps): 
    FD1FGeneric(TYPE), volType(VolSurface::TYPE->getName()) {
    product = 0;
    
    useFwdGrid            = true;
    DEBUG_UseCtrlVar      = true;
    DEBUG_SameGridDelta   = true;
    DEBUG_SameGridVegaRho = false;

    varMethod             = 0;
    gridType              = 0;

    useDivsForAssetDrift  = false;
    additionalAssetDrift  = 0.0;
    defaultBarrier        = 0.0;

    riskyEquityGrowth     = false;

    this->stepsPerYear    = stepsPY;
    this->stockSteps      = stockSteps;
    this->TruncationStd   = 7.0;
}

FD1FE2C::FD1FE2C(CClassConstSP clazz): 
    FD1FGeneric(clazz), volType(VolSurface::TYPE->getName())
{
    product               = 0;

    useFwdGrid            = true;
    DEBUG_UseCtrlVar      = true;
    DEBUG_SameGridDelta   = true;
    DEBUG_SameGridVegaRho = false;

    varMethod             = 0;
    gridType              = 0;

    useDivsForAssetDrift  = false;
    additionalAssetDrift  = 0.0;
    defaultBarrier        = 0.0;

    riskyEquityGrowth     = false;

    this->stepsPerYear    = -1;
    this->stockSteps      = -1;
    this->TruncationStd   = 7.0;
}


FD1FE2C::~FD1FE2C()
{
    Clear();

    if (product) {
        delete product;
        product = 0;
    }
}

/** get processed vol */
// need to examine if this is general
void FD1FE2C::InitVol()
{
    // asset vol - the asset vol must be a single strike surface, so any vol request with the correct
    // maturity date and fwd start flag will be ok, but it's not strictly ok
    VolLN = CVolProcessedBSSP::dynamicCast(
        CVolProcessedSP(Underlier->getProcessedVol(product->GetLNRequest().get())));

    // equity vol
    FirmAsset* firmAsset = 0;
    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
    } else {
        throw ModelException("FD1FE2C::getDiscountTerm",
                "Underlying credit asset must be a FirmAsset");
    }

    CAssetWrapper equityAsset = firmAsset->getEquityAsset();

    VolLNEquity = CVolProcessedBSSP::dynamicCast(
        CVolProcessedSP(equityAsset->getProcessedVol(product->GetLNRequest().get())));
}

/** calculate a term structure vol^2 or just one point at maturity */
/* Hard coded vol benchmarks? Very bad. */
void FD1FE2C::CalcV2Term(const DateTime& valDate,const DateTime& startDate,
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
void FD1FE2C::PostSetup()
{
    Variance.resize(TimePts.NumOfStep);
    VolLN->CalcVar(TimePts.StepDates, CVolProcessedBS::fromFirst, Variance);
    Variance.insert(Variance.begin(), 0.0); // start with 0

    equityVariance.resize(TimePts.NumOfStep);
    VolLNEquity->CalcVar(TimePts.StepDates, CVolProcessedBS::fromFirst, equityVariance);
    equityVariance.insert(equityVariance.begin(), 0.0); // start with 0

    stockMaxSeg.resize(TimePts.SegmentEnd.size());
    stockMinSeg.resize(TimePts.SegmentEnd.size());
    strikeSeg.resize(TimePts.SegmentEnd.size());
    int i;
    for (i=0; i<int(stockMaxSeg.size()); i++) {
        stockMaxSeg[i] = -1.;
        stockMinSeg[i] = -1.;
        strikeSeg[i] = -1.;
    }

    equityMinSeg.resize(TimePts.SegmentEnd.size());
    equityMaxSeg.resize(TimePts.SegmentEnd.size());
    equityStrikeSeg.resize(TimePts.SegmentEnd.size());
    for (i=0; i<int(TimePts.SegmentEnd.size()); i++) {
        equityMinSeg[i]    = -1.;
        equityMaxSeg[i]    = -1.;
        equityStrikeSeg[i] = -1.;
    }

    // bootstrap the liquidity spread curve if required
    try {

        FirmAsset* firmAsset = 0;
        if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
            CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
            firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
        } else {
            throw ModelException("FD1FE2C::PostSetup",
                    "Underlying credit asset must be a FirmAsset");
        }

        defaultBarrier = firmAsset->getDefaultBarrier();

        CDSParSpreadsWrapper  liquiditySpreads = firmAsset->getLiquiditySpreadCurve();
        if (!liquiditySpreads.getName().empty())  {
            BadDayConventionSP bdc(BadDayConventionFactory::make("None"));
            DayCountConventionSP dcc(DayCountConventionFactory::make("Actual/360"));

            YieldCurveSP discount(copy(riskFreeCurve.get()));

            liquiditySpreadRates = CDSHelper::CParSpreadDefaultRatesSP( new 
                    CDSHelper::CParSpreadDefaultRates(
                        liquiditySpreads.getSP(), 
                        discount, 
                        TimePts.StepDates[0], 
                        TimePts.StepDates[0].rollDate(1), 
                        bdc.get(),
                        dcc));
        }

        if (riskyEquityGrowth) {
            equitySpreadRates = CDSHelper::CFirmAssetDefaultRatesSP( new CDSHelper::CFirmAssetDefaultRates(firmAsset, *liquiditySpreadRates.get()));
        }

    } catch (exception &e) {
        throw ModelException(&e, "CredDefSwaption::basePrice: Failed to calculate clean liquidity spread curve");
    }


}

TimeMetricConstSP FD1FE2C::GetTimeMetric() const
{
    if (!VolLN)
        throw ModelException("FD1FLN::GetTimeMetric", "Volatility not initialised.");
    
    return VolLN->GetTimeMetric();
}

int FD1FE2C::GetStepVol(int step, vector<double>& vol, const double*, int, int end)
{
    if (int(vol.size()) != (end+1))
         vol.resize(end+1);

    if (step >= TimePts.NumOfStep)
        step =TimePts.NumOfStep - 1; // returns the same vol for maturity step !!!
    if (step < 0)
        step = 0;
    
    if (TimePts.TradeYrFrac[step+1] != 0.) {
        vol[0] = (equityVariance[step+1]-equityVariance[step])/TimePts.TradeYrFrac[step+1];
        vol[0] = sqrt(vol[0]);
    } else { // if there is no trading time between steps return a one day vol to avoid a div by zero
        CDateTimeArray dateAndNextDay;
        dateAndNextDay.resize(2);
        dateAndNextDay[0] = TimePts.StepDates[step];
        dateAndNextDay[1] = dateAndNextDay[0].rollDate(1);
        CDoubleArray volTmp;
        volTmp.resize(1);
        VolLNEquity->CalcVol(dateAndNextDay, CVolProcessedBS::forward, volTmp);
        vol[0] = volTmp[0];
    }

	for(int i=1; i <= end; i++)
		vol[i] = vol[0];

    return 1;
}


int FD1FE2C::GetAssetStepVol(int step, vector<double>& vol, const double*, int, int end)
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

CVolProcessedBSConstSP FD1FE2C::getProcessedVol()
{
    return VolLN;
}

void FD1FE2C::SetGridData() {
    static const string method = "FD1D::SetGridData";
    try {
        Actual365F fwdRateDCC;
        times.resize(TimePts.NumOfStep+1);

        // the forward space grid stuff is not yet implemented ... 
        if (useFwdGrid == true) {
            throw ModelException("FD1FE2C::SetGridData",
                                 "The forward space grid is not currently supported for E2C");
        }

        int i;
        for (i=0; i<TimePts.NumOfStep+1; i++) {
            times[i] = TimePts.TradeYrFrac[i];
        }

        inForwards.resize(TimePts.NumOfStep+1);
        Underlier->fwdValue(TimePts.StepDates, inForwards);

        if (hasEquityLayer()) {
            FirmAsset* firmAsset = 0;
            if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
                CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
                firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
            } else {
                throw ModelException("FD1FE2C::getDiscountTerm",
                        "Underlying credit asset must be a FirmAsset");
            }

            inFwdEQ.resize(TimePts.NumOfStep+1);
            CAssetWrapper equityAsset = firmAsset->getEquityAsset();
            equityAsset->fwdValue(TimePts.StepDates, inFwdEQ);
        }

        if (useFwdGrid == true) {
            inPVs.resize(TimePts.NumOfStep+1);

            for (i=1; i<TimePts.NumOfStep+1; i++) {
                inPVs[i] = DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
            }
        } else {

            inIr.resize(TimePts.NumOfStep+1);
            inDivy.resize(TimePts.NumOfStep+1);

            inIrEQ.resize(TimePts.NumOfStep+1);
            inDivyEQ.resize(TimePts.NumOfStep+1);
            inPVRiskFree.resize(TimePts.NumOfStep+1);
    

            inPVs.resize(TimePts.NumOfStep+1);
            inPVs[0] = 1.0;
            for (i=1; i<TimePts.NumOfStep+1; i++) {
                inPVs[i]        = DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
                inPVRiskFree[i] = riskFreeCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
            }

            for (i=1; i<TimePts.NumOfStep+1; i++) {

                if (times[i] != 0.) {
                    // remember to change for riskFreeCurve here
                    // inIr[i] = -log(riskFreeCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]))/TimePts.TradeYrFrac[i];
                    // inIr[i]   = riskFreeCurve->fwd(TimePts.StepDates[i-1], TimePts.StepDates[i].rollDate(1),
                    //                               &fwdRateDCC, CompoundBasis::CONTINUOUS);
                    inIr[i] = -log(inPVRiskFree[i])/TimePts.TradeYrFrac[i];



                    inDivy[i] = inIr[i] - log(inForwards[i]/inForwards[i-1])
                        /TimePts.TradeYrFrac[i];
                } else { // If the time is zero DEV routine will simply map accross. 
                    // inIr is the ir forward factor
                    inIr[i] = 1./riskFreeCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
                    // inDivy is the ir forward factor divided by the equity fwd factor 
                    inDivy[i] = inForwards[i-1]/inForwards[i]*inIr[i];
                }

                if ( hasEquityLayer() ) {
                    if (times[i] != 0.) {
                        inIrEQ[i] = -log(inPVs[i])/TimePts.TradeYrFrac[i];

                        inDivyEQ[i] = inIrEQ[i] - log(inFwdEQ[i]/inFwdEQ[i-1])
                            /TimePts.TradeYrFrac[i];
                    } else { // If the time is zero DEV routine will simply map accross. 
                        // inIr is the ir forward factor
                        inIrEQ[i] = 1./DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
                        // inDivy is the ir forward factor divided by the equity fwd factor 
                        inDivyEQ[i] = inFwdEQ[i-1]/inFwdEQ[i]*inIrEQ[i];
                    }
                }
            }
        }


        inSegEnd.resize(TimePts.SegmentEnd.size());
        for (i=0; i<int(TimePts.SegmentEnd.size()); i++) {
            inSegEnd[i] = TimePts.SegmentEnd[i];
        }            
        
        inSegMax.resize(TimePts.SegmentEnd.size());
        inSegMin.resize(TimePts.SegmentEnd.size());
        inSegStrike.resize(TimePts.SegmentEnd.size());

        // set-up data for the equity grid ... pretty simplistic at the moment, but works, could do with some
        // fine tuning
        FirmAsset* firmAsset = 0;
        if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
            CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
            firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
        } else {
            throw ModelException("FD1FE2C::getDiscountTerm",
                    "Underlying credit asset must be a FirmAsset");
        }

        CAssetWrapper equityAsset = firmAsset->getEquityAsset();
        double equityForward  = equityAsset->fwdValue(TimePts.StepDates[TimePts.NumOfStep]);
        double recoveryVariance = firmAsset->getLambda() * firmAsset->getLambda();
        for (i=0; i<int(TimePts.SegmentEnd.size()); i++) {
            if (useFwdGrid == false) {
                double truncSegMax = inForwards[TimePts.NumOfStep]*exp(TruncationStd*sqrt( Variance[TimePts.NumOfStep])+recoveryVariance);
                double truncSegMin = inForwards[TimePts.NumOfStep]*exp(-TruncationStd*sqrt(Variance[TimePts.NumOfStep])+recoveryVariance);
                if ( !Maths::isPositive(stockMaxSeg[i]) ) {
                    inSegMax[i] = truncSegMax;
                } else {
                    inSegMax[i] = Maths::min(stockMaxSeg[i], truncSegMax);
                }
                if ( !Maths::isPositive(stockMinSeg[i])) {
                    inSegMin[i] = truncSegMin;
                } else {
                    inSegMin[i] = stockMinSeg[i];
                    // inSegMin[i] = Maths::min(stockMinSeg[i], truncSegMin);;
                }
                
                // never let the max be less than 1.01 times the min
                inSegMax[i] = Maths::max(inSegMax[i], inSegMin[i]*1.01);

                if ( !Maths::isPositive(strikeSeg[i])) {
                    inSegStrike[i] = inForwards[inSegEnd[i]];
                } else {
                    inSegStrike[i] = strikeSeg[i];
                }

                if (hasEquityLayer()) {
                    FirmAsset* firmAsset = 0;
                    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
                        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
                        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
                    } else {
                        throw ModelException("FD1FE2C::getDiscountTerm",
                                "Underlying credit asset must be a FirmAsset");
                    }

                    equityInSegMin.resize(TimePts.SegmentEnd.size());
                    equityInSegMax.resize(TimePts.SegmentEnd.size());
                    equityInSegStrike.resize(TimePts.SegmentEnd.size());

                    double eqTruncSegMax = inFwdEQ[TimePts.NumOfStep]*exp(TruncationStd*sqrt(equityVariance[TimePts.NumOfStep]));
                    double eqTruncSegMin = inFwdEQ[TimePts.NumOfStep]*exp(-TruncationStd*sqrt(equityVariance[TimePts.NumOfStep]));
                    if ( !Maths::isPositive(equityMaxSeg[i]) ) {
                        equityInSegMax[i] = eqTruncSegMax;
                    } else {
                        equityInSegMax[i] = Maths::min(equityMaxSeg[i], truncSegMax);
                    }
                    if ( !Maths::isPositive(equityMinSeg[i])) {
                        equityInSegMin[i] = eqTruncSegMin;
                    } else {
                        equityInSegMin[i] = Maths::max(equityMinSeg[i], truncSegMin);;
                    }
                
                    // never let the max be less than 1.01 times the min
                    equityInSegMax[i] = Maths::max(equityInSegMax[i], equityInSegMin[i]*1.01);

                    if ( !Maths::isPositive(equityStrikeSeg[i])) {
                        equityInSegStrike[i] = inFwdEQ[inSegEnd[i]];
                    } else {
                        equityInSegStrike[i] = strikeSeg[i];
                    }
                }
            } else {
                double pvFact = 
                    DiscountCurve->pv(TimePts.StepDates[TimePts.SegmentStart[i]], TimePts.StepDates[TimePts.SegmentEnd[i]]);
        
                double truncSegMax = exp(TruncationStd*sqrt(Variance[TimePts.NumOfStep]));
                double truncSegMin = exp(-TruncationStd*sqrt(Variance[TimePts.NumOfStep]));

                if ( !Maths::isPositive(stockMaxSeg[i])) {
                    inSegMax[i] = truncSegMax;
                 //   topMult = exp(genericFDModel->TruncationStd*sqrt(genericFDModel->Variance[genericFDModel->TimePts.SegmentEnd[currSeg]]));
                } else {
                    inSegMax[i] = Maths::min(stockMaxSeg[i]/inForwards[TimePts.SegmentEnd[i]]/pvFact, truncSegMax);
                }

                if ( !Maths::isPositive(stockMinSeg[i])) {
                   inSegMin[i] = truncSegMin;
                 //   botMult = exp(-genericFDModel->TruncationStd*sqrt(genericFDModel->Variance[genericFDModel->TimePts.SegmentEnd[currSeg]]));
                } else {
                    inSegMin[i] = Maths::max(stockMinSeg[i]/inForwards[TimePts.SegmentStart[i]]/pvFact, truncSegMin);
                }

                // never let the max be less than 1.01 times the min
                inSegMax[i] = Maths::max(inSegMax[i], inSegMin[i]*1.01);

                if ( !Maths::isPositive(strikeSeg[i])) {
                    inSegStrike[i] = exp(-Variance[inSegEnd[i]]/4.);
                } else {
                    inSegStrike[i] = strikeSeg[i]/inForwards[inSegEnd[i]];
                }
            }
        }
    } catch (exception& e){
        throw ModelException(&e, method);
    }
}

void FD1FE2C::preProcessGrid(const double* s, int start, int end)
{
    double jumpSpread   = 0.0;
    double equitySpread = 0.0;
    int    numStep      = end-start+1;
    int    step,j;
    if ( riskyEquityGrowth ) {

        FirmAsset* firmAsset = 0;
        if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
            CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
            firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
        }

        equitySpreads =     CDoubleMatrix(TimePts.NumOfStep+1, numStep);
        for (step=1; step<TimePts.NumOfStep+1; step++) {
            jumpSpread = 0.0;
            if ( riskyEquityGrowth && !!liquiditySpreadRates) {
                jumpSpread = liquiditySpreadRates->calcDefaultPV(TimePts.StepDates[step-1], TimePts.StepDates[step]);
            }
            for (j=0; j<numStep; j++) {
                equitySpread = jumpSpread * firmAsset->CalcNoDefaultProb(TimePts.StepDates[step], TimePts.StepDates[step+1], s[j]);
                equitySpreads[step][j] = -log(equitySpread) / TimePts.TradeYrFrac[step];
            }
        }
    }
}


FDTermStructureSP FD1FE2C::getDriftTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                        const double &irPert, const double &divPert, const bool doEquityLayer)
{
    static FDTermStructureSP driftTerm;
    double equitySpread;

    FirmAsset* firmAsset = 0;
    if ( doEquityLayer && riskyEquityGrowth && FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
    }

    // calculate liquidity spread rate
    double jumpSpread = 0.0;
    if ( riskyEquityGrowth && !!liquiditySpreadRates) {
        jumpSpread = liquiditySpreadRates->calcDefaultPV(TimePts.StepDates[step-1], TimePts.StepDates[step]);
    }

    if ( useFwdGrid ) {
        throw ModelException("FD1FE2C::getDriftTerm",
                "The forward space grid is not currently supported for E2C");
    } else {
        
        DoubleArraySP  termArray;
        if (!driftTerm || driftTerm->termStructure->size() != end-start+1) {
            termArray = DoubleArraySP(new DoubleArray(end-start+1));
            driftTerm = FDTermStructureSP(new FDTermStructure(termArray));
        } else {
            termArray = driftTerm->termStructure;
        }

        if ( doEquityLayer ) {

            // log-normal equity diffusion
            if ( useFwdGrid ) {
                for (int i=0; i<termArray->size() ; ++i) {
                    (*termArray)[i] = (irPert-divPert) * s[i];
                }
            } else {
                for (int i=0; i<termArray->size() ; ++i) {
                    if ( riskyEquityGrowth ) {
                        equitySpread = equitySpreads[step][i];
                    } else {
                        equitySpread = 0.0;
                    }

                    (*termArray)[i] = ((equitySpread+inIrEQ[step]+irPert)-(inDivyEQ[step]+divPert)) * s[i];

                }
            }
        } else {
            // log-normal asset diffusion
            for (int i=0; i<termArray->size() ; ++i) {
                if (useDivsForAssetDrift) {
                   if ( useFwdGrid ) {
                      (*termArray)[i] = - (s[i] - defaultBarrier) * additionalAssetDrift;
                   } else {
                      (*termArray)[i] = - (s[i] - defaultBarrier) * (additionalAssetDrift + inDivyEQ[step] + divPert);
                   }
                } else {
                    if ( useFwdGrid ) {
                        (*termArray)[i] = - (s[i] - defaultBarrier) * additionalAssetDrift;
                    } else {
                       (*termArray)[i] = - (s[i] - defaultBarrier) * additionalAssetDrift;
                    }
                }
            }
        }

    }
    return driftTerm;
}

FDTermStructureSP FD1FE2C::getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                            const double &irPert, const double &divPert, const bool doEquityLayer)
{
    static FDTermStructureSP diffusionTerm;
    vector<double>    vol;

    DoubleArraySP  termArray;
    if (!diffusionTerm || diffusionTerm->termStructure->size() != end-start+1) {
        termArray     = DoubleArraySP(new DoubleArray(end-start+1));
        diffusionTerm = FDTermStructureSP(new FDTermStructure(termArray));
    } else {
        termArray = diffusionTerm->termStructure;
    }

    // need to have asset and equity specific arrays here
    if ( doEquityLayer ) {
        GetStepVol(step-1, vol, s, start, end);
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = .5 * vol[i] * vol[i] * s[i] * s[i];
        }
        diffusionTerm = FDTermStructureSP(new FDTermStructure(termArray));
    } else {
        GetAssetStepVol(step-1, vol, s, start, end);
        for (int i=0; i<termArray->size() ; ++i) {
            (*termArray)[i] = .5 * vol[i] * vol[i] * s[i] * s[i];
        }
        // diffusionTerm = FDTermStructureSP(new FDTermStructure(termArray));
    }
    return diffusionTerm;
}

FDTermStructureSP FD1FE2C::getCouponTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                         const double &irPert, const double &divPert, const bool doEquityLayer)
{
    static FDTermStructureSP couponTerm;

    DoubleArraySP  termArray;
    if (!couponTerm || couponTerm->termStructure->size() != end-start+1) {
        termArray  = DoubleArraySP(new DoubleArray(end-start+1));
        couponTerm = FDTermStructureSP(new FDTermStructure(termArray));
    } else {
        termArray = couponTerm->termStructure;
    }

    if ( doEquityLayer ) {
        couponTerm->term      = 0.0;
        couponTerm->dimension = 0;
    } else {
        FirmAsset* firmAsset = 0;
        if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
            CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
            firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

            if (!firmAsset) {
                throw ModelException("FD1FE2C::getDriftTerm",
                        "Underlying credit asset must be a FirmAsset");
            }
        } else {
            throw ModelException("FD1FE2C::getDriftTerm",
                    "Underlying credit asset must be a FirmAsset");
        }

        // calculate liquidity spread rate
        double defRate = 0.0;
        if ( !!liquiditySpreadRates) {
            DoubleArraySP     termArray(new DoubleArray(end-start+1));

            double liquidityJumpProb = liquiditySpreadRates->calcDefaultPV(
                    TimePts.StepDates[step-1], TimePts.StepDates[step]);

            defRate = -log(liquidityJumpProb) / TimePts.TradeYrFrac[step];

            // calculate default payment
            double defaultPayment = product->getCoupon(step, s, start, end);

            for (int i=0; i<termArray->size() ; ++i) {
                (*termArray)[i] = defRate * defaultPayment;
            }
            couponTerm->dimension = 1;
        } else {
            couponTerm->term      = 0.0;
            couponTerm->dimension = 0;
        }
    }

    return couponTerm;
}

FDTermStructureSP FD1FE2C::getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGrid,
                                           const double &irPert, const double &divPert, const bool doEquityLayer)
{
    static FDTermStructureSP discountTerm;

    double         defRate              = 0.0;
    double         equityDiscountSpread = 0.0;
    DoubleArraySP  termArray;


    // only discount at credit spread for the bond floor ... there might be a problem for exchangeables here
    FirmAsset* firmAsset = 0;
    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

        if (!firmAsset) {
            throw ModelException("FD1FE2C::getDiscountTerm",
                    "Underlying credit asset must be a FirmAsset");
        }
    } else {
        throw ModelException("FD1FE2C::getDiscountTerm",
                "Underlying credit asset must be a FirmAsset");
    }

    if ( !!liquiditySpreadRates) {
        double liquidityJumpProb = liquiditySpreadRates->calcDefaultPV(
                TimePts.StepDates[step-1], TimePts.StepDates[step]);

        defRate = -log(liquidityJumpProb) / TimePts.TradeYrFrac[step];
    }

    if (!discountTerm || discountTerm->termStructure->size() != end-start+1) {
        termArray     = DoubleArraySP(new DoubleArray(end-start+1));
        discountTerm = FDTermStructureSP(new FDTermStructure(termArray));
    } else {
        termArray = discountTerm->termStructure;
    }

    if ( !doEquityLayer) {
        for (int i=0; i<termArray->size() ; ++i) {
            if ( useFwdGrid ) {
                (*termArray)[i] = defRate + irPert;
            } else {
                (*termArray)[i] = defRate + inIr[step] + irPert;
            }
        }
    } else {
        for (int i=0; i<termArray->size() ; ++i) {
            if ( useFwdGrid ) {
                (*termArray)[i] = defRate + irPert;
            } else {
                equityDiscountSpread = 0.0;
                if ( riskyEquityGrowth ) {
                    equityDiscountSpread = equitySpreads[step][i];
                }
                (*termArray)[i] = defRate + inIrEQ[step] + irPert;
            }
        }
    }

    // discountTerm = FDTermStructureSP(new FDTermStructure(termArray));
    return discountTerm;
}

bool FD1FE2C::doLambdaAdjust()
{
    return true;
}

double FD1FE2C::getLambda()
{
    FirmAsset* firmAsset = 0;
    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);
    } else {
        throw ModelException("FD1FE2C::getDiscountTerm",
                "Underlying credit asset must be a FirmAsset");
    }

    return firmAsset->getLambda();
}

double FD1FE2C::getLambdaAdjustedSpot(const double spot)
{
    double lambda = getLambda();
    return (spot * exp(lambda*lambda));
}

bool FD1FE2C::isE2C()
{
    return true;
}

bool FD1FE2C::hasEquityLayer()
{
    return product->hasEquityLayer();
}

void FD1FE2C::mapToEquitySpace(double* stockArray, int numnStockSteps)
{
    FirmAsset* firmAsset = 0;
    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

        if (!firmAsset) {
            throw ModelException("FD1FE2C::getDiscountTerm",
                    "Underlying credit asset must be a FirmAsset");
        }
    } else {
        throw ModelException("FD1FE2C::getDiscountTerm",
                "Underlying credit asset must be a FirmAsset");
    }

    double lambda         = firmAsset->getLambda();
    double fxRate         = firmAsset->getFXRate();
    double expLambda      = exp(lambda*lambda);
    double defaultBarrier = firmAsset->getDefaultBarrier();
    for(int i=0 ; i<numnStockSteps ; ++i) {
        stockArray[i] = (stockArray[i] / fxRate + defaultBarrier) * expLambda;
    }
}

double FD1FE2C::getFXRate()
{
    FirmAsset* firmAsset = 0;
    if ( FirmAsset::TYPE->isInstance(Underlier.get())) {
        CAsset*    ncAsset   = const_cast<CAsset*>(Underlier.get());
        firmAsset = dynamic_cast<FirmAsset*>(ncAsset);

        if (!firmAsset) {
            throw ModelException("FD1FE2C::getDiscountTerm",
                    "Underlying credit asset must be a FirmAsset");
        }
    } else {
        throw ModelException("FD1FE2C::getDiscountTerm",
                "Underlying credit asset must be a FirmAsset");
    }

    return firmAsset->getFXRate();
}


MarketDataFetcherSP FD1FE2C::createMDF() const {
    return MarketDataFetcherSP(new MarketDataFetcherLN(volType));
}


class FD1FE2CHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FE2C, clazz);
        SUPERCLASS(FD1FGeneric);
        FIELD(volType, "Type of vol to use");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(useDivsForAssetDrift, "whether to drift the asset with the equity dividends");
        FIELD_MAKE_OPTIONAL(useDivsForAssetDrift);
        FIELD(additionalAssetDrift, "additional explicit drift for the asset");
        FIELD_MAKE_OPTIONAL(additionalAssetDrift);
        FIELD(riskyEquityGrowth,    "whether to grow the equity at the risky local spread");
        FIELD_MAKE_OPTIONAL(riskyEquityGrowth);

        EMPTY_SHELL_METHOD(defaultFD1FE2C);
    }
    
    static IObject* defaultFD1FE2C(){
        return new FD1FE2C();
    }
};


CClassConstSP const FD1FE2C::TYPE = CClass::registerClassLoadMethod(
    "FD1FE2C", typeid(FD1FE2C), FD1FE2CHelper::load);

DRLIB_END_NAMESPACE
