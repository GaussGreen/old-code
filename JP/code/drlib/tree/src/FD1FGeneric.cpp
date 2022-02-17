//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FGeneric.cpp
//
//   Description : One factor generic finite difference base class
//
//   Author      : Andr� Segger
//
//   Date        : 15 October 2003
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD1FGeneric.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Actual365F.hpp"
#include "edginc/CrossGamma.hpp"
#include "edginc/FXCrossGamma.hpp"
#include "edginc/BasketDelta.hpp"
#include "edginc/DeltaDDE.hpp"
#include "edginc/MuParallel.hpp"
#include "edginc/MuPointwise.hpp"
#include "edginc/MuSpecial.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"
#include "edginc/RollingTheta.hpp"
#include "edginc/DeltaNextDay.hpp"
#include "edginc/FXDelta.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DeltaSurface.hpp"

DRLIB_BEGIN_NAMESPACE

static bool rebuildRequest(const SensitivitySP sens)
{
    const Sensitivity* sensPtr = sens.get();
    bool rebuild = Theta::TYPE->isInstance(sensPtr)
                || MuParallel::TYPE->isInstance(sensPtr)
                || MuPointwise::TYPE->isInstance(sensPtr)
                || MuSpecial::TYPE->isInstance(sensPtr)
                || FXDelta::TYPE->isInstance(sensPtr)
                || DeltaNextDay::TYPE->isInstance(sensPtr)
                || DeltaSurface::TYPE->isInstance(sensPtr)
                || RollingTheta::TYPE->isInstance(sensPtr);

    return rebuild;
}

bool FD1FGeneric::SameGrid(const CControl* control)
{
    bool same = false;
    
    if (!control || control->isPricing()) { 
        same = false;
    } else {
        SensitivitySP sens = control->getCurrentSensitivity();

        if (rebuildRequest(sens) == true) {
            same = false;
        } else {
            Sensitivity* sensPtr = sens.get();
            if (Delta::TYPE->isInstance(sensPtr)       ||
                BasketDelta::TYPE->isInstance(sensPtr) ||
                DeltaDDE::TYPE->isInstance(sensPtr)    ||
                CrossGamma::TYPE->isInstance(sensPtr)  ||
                FXCrossGamma::TYPE->isInstance(sensPtr) ) 
            {
                same = DEBUG_SameGridDelta;
            } else {
                same = DEBUG_SameGridVegaRho;
            }
        }
    }

    return same;
}


FD1FGeneric::FD1FGeneric(CClassConstSP clazz): CModel(clazz){
    product = 0;
    
    useFwdGrid = false;
    DEBUG_UseCtrlVar = true;
    DEBUG_SameGridDelta = true;
    DEBUG_SameGridVegaRho = false;

    varMethod = 0;

    gridType = 2;

    stockStepsToUse = -1;
    stepsPerYearToUse = -1;
}

FD1FGeneric::~FD1FGeneric()
{
    Clear();

    if (product != 0)
        delete product;
}


/** clean up */
void FD1FGeneric::Clear()
{
    TimePts.Clear();
}

/** calculate single price and store result in CResult */
void FD1FGeneric::Price(CInstrument*  instrument, 
                    CControl*     control, 
                    CResults*     results){
    if (!IIntoProduct::TYPE->isInstance(instrument)){
        throw ModelException("FD1FGeneric::Price", "Instrument of type "+
                             instrument->getClass()->getName() +
                             " does not support FD1FGeneric::IIntoProduct");
    }

    FD1FGenericSP modelToUse;

    try {
        // check and price a dead instrument
        if (instrument->priceDeadInstrument(control, results))
            return; // done for a dead instrument

        if (!control || control->isPricing()) {
            modelToUse = FD1FGenericSP::attachToRef(this);
        } else if (SameGrid(control) == false) {
            modelToUse = FD1FGenericSP(dynamic_cast<FD1FGeneric*>(this->clone()));
        } else {
            modelToUse = FD1FGenericSP::attachToRef(this);
        }

        // cast to FD1FGeneric::IIntoProduct
        IIntoProduct& intoProd = dynamic_cast<IIntoProduct&>(*instrument);
        // create the product
        if ( !(!modelToUse) &&modelToUse->product != 0)
            delete modelToUse->product;
        modelToUse->product = intoProd.createProduct(modelToUse.get());
        // Get the stock and discount curves
        modelToUse->Underlier = modelToUse->product->GetAssetRef();
        modelToUse->DiscountCurve = modelToUse->product->GetDiscCurveRef();
        modelToUse->riskFreeCurve = modelToUse->product->GetRiskFreeCurve();

        modelToUse->InitVol();

        modelToUse->product->InitGridSize();

        if (SameGrid(control) == false) {
            modelToUse->product->InitFD(control);
        }
        modelToUse->PostSetup();

        modelToUse->product->InitProdFD();

        modelToUse->fdEngine.initModel(modelToUse);
 
        if (SameGrid(control) == false) {
            // invoke the pricing

            modelToUse->SetGridData();

            modelToUse->fdEngine.init(modelToUse->TimePts.NumOfStep, 
                                      modelToUse->stockStepsToUse, 
                                      &*modelToUse->times.begin(), 
                                      &*modelToUse->inForwards.begin(), 
                                      &*modelToUse->inPVs.begin(),
                                      modelToUse->TimePts.SegmentEnd.size(), 
                                      &*modelToUse->inSegEnd.begin(), 
                                      &*modelToUse->inSegMax.begin(), 
                                      &*modelToUse->inSegMin.begin(),
                                      &*modelToUse->inSegStrike.begin(), 
                                      modelToUse->useFwdGrid, 
                                      &*modelToUse->inIr.begin(), 
                                      &*modelToUse->inDivy.begin(), 
                                      modelToUse->numPriceArrays, 
                                      modelToUse->gridType);
        
            modelToUse->PriceEnd.resize(modelToUse->numPriceArrays);

            modelToUse->divPert.resize(modelToUse->TimePts.NumOfStep+1);
            modelToUse->irPert.resize(modelToUse->TimePts.NumOfStep+1);

            int i;
            for (i = 1; i<=modelToUse->TimePts.NumOfStep; i++){
                modelToUse->divPert[i] = 0.;
                modelToUse->irPert[i] = 0.;
            }

            modelToUse->fdEngine.loop(
                modelToUse->product, 
                modelToUse->Underlier->fwdValue(modelToUse->
                                                TimePts.StepDates[0]), 
                &*modelToUse->PriceEnd.begin(), 
                &*modelToUse->divPert.begin(),
                &*modelToUse->irPert.begin());
            
        } else {
            // same grid. Don't bother with modelToUse as it's definitely the base one
            int i;
            CDoubleArray newForwards;

            newForwards.resize(TimePts.NumOfStep+1);

            Underlier->fwdValue(TimePts.StepDates, newForwards);

            irPert.resize(TimePts.NumOfStep+1);
            double newPV;
            divPert.resize(TimePts.NumOfStep+1);
            for (i = 1; i<=TimePts.NumOfStep; i++){
                newPV = DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
                
                if (times[i] != 0.) {
                    irPert[i] = log(inPVs[i]/newPV)/times[i];
                    divPert[i] = log((inForwards[i]/inForwards[i-1]*inPVs[i])/(newForwards[i]/newForwards[i-1]*newPV))/times[i];
                } else { // times[i] == 0.
                    if (useFwdGrid == true) {
                        irPert[i] = inPVs[i]/newPV;
                        divPert[i] = (inForwards[i]/inForwards[i-1]*inPVs[i])/(newForwards[i]/newForwards[i-1]*newPV);
                    } else {
                        irPert[i] = 1./newPV - inIr[i];
                        divPert[i] = newForwards[i-1]/newForwards[i]/newPV - inDivy[i];
                    }
                }
            } 

            PriceEnd.resize(numPriceArrays);
    
            fdEngine.loop(product, 
                          Underlier->fwdValue(TimePts.StepDates[0]),
                          &*PriceEnd.begin(), &*divPert.begin(),
                          &*irPert.begin());
        }


        double price;
        
        // discount factor for forward starting effect if tree does not start today
        double discFactor = modelToUse->DiscountCurve->pv(modelToUse->TimePts.StepDates[0]);
        
        price = modelToUse->product->RefinePrice(modelToUse->PriceEnd[0],discFactor,modelToUse->DEBUG_UseCtrlVar);

        // record price and additional outputs
        results->storePrice(price, modelToUse->DiscountCurve->getCcy());
        modelToUse->product->recordOutputRequests(control, results, price);

		// record output from model itself. Used for DDE output
        modelToUse->recordOutputRequests(control, results);
        
    } catch (exception& e){
        if ( !(!modelToUse) && modelToUse->product) {
            delete modelToUse->product;
            modelToUse->product = 0;
        }
        throw ModelException(&e, "FD1FGeneric::Price");
    }


    delete modelToUse->product;
    modelToUse->product = 0;
}

void FD1FGeneric::SetGridData(){
    static const string method = "FD1D::SetGridData";
    try
    {
        Actual365F fwdRateDCC;
        times.resize(TimePts.NumOfStep+1);

        int i;
        for (i=0; i<TimePts.NumOfStep+1; i++) {
            times[i] = TimePts.TradeYrFrac[i];
        }

        inForwards.resize(TimePts.NumOfStep+1);
        Underlier->fwdValue(TimePts.StepDates, inForwards);       
        
        if (useFwdGrid == true) {
            inPVs.resize(TimePts.NumOfStep+1);

            for (i=1; i<TimePts.NumOfStep+1; i++) {
                inPVs[i] = DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
            }
        } else {

            inIr.resize(TimePts.NumOfStep+1);
            inDivy.resize(TimePts.NumOfStep+1);

            inPVs.resize(TimePts.NumOfStep+1);
            for (i=1; i<TimePts.NumOfStep+1; i++) {
                inPVs[i] = DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
            }

            for (i=1; i<TimePts.NumOfStep+1; i++) {

                if (times[i] != 0.) {
                    inIr[i] = -log(inPVs[i])/TimePts.TradeYrFrac[i];

                    inDivy[i] = inIr[i] - log(inForwards[i]/inForwards[i-1])
                        /TimePts.TradeYrFrac[i];
                } else { // If the time is zero DEV routine will simply map accross. 
                    // inIr is the ir forward factor
                    inIr[i] = 1./DiscountCurve->pv(TimePts.StepDates[i-1], TimePts.StepDates[i]);
                    // inDivy is the ir forward factor divided by the equity fwd factor 
                    inDivy[i] = inForwards[i-1]/inForwards[i]*inIr[i];
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
        
        double truncSegMax = exp(TruncationStd*sqrt(Variance[TimePts.NumOfStep]));
        double truncSegMin = exp(-TruncationStd*sqrt(Variance[TimePts.NumOfStep]));
        
        if (false == useFwdGrid) {
            truncSegMax *= inForwards[TimePts.NumOfStep];
            truncSegMin *= inForwards[TimePts.NumOfStep];
        }
        
        for (i=0; i<int(TimePts.SegmentEnd.size()); i++) {

            if ( !Maths::isPositive(stockMaxSeg[i]) ) {
                inSegMax[i] = truncSegMax;
            } else {
                inSegMax[i] = Maths::min(stockMaxSeg[i], truncSegMax);
            }
            if ( !Maths::isPositive(stockMinSeg[i])) {
                inSegMin[i] = truncSegMin;
            } else {
                inSegMin[i] = Maths::max(stockMinSeg[i], truncSegMin);;
            }
                
            // never let the max be less than 1.01 times the min
            inSegMax[i] = Maths::max(inSegMax[i], inSegMin[i]*1.01);

            if (useFwdGrid == false) {
                if ( !Maths::isPositive(strikeSeg[i])) {
                    inSegStrike[i] = inForwards[inSegEnd[i]];
                } else {
                    inSegStrike[i] = strikeSeg[i];
                }
            }
            else {
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

void FD1FGeneric::Setup(const DateTime& valDate, const DateTimeArray& segDates, 
               const vector<int>& density, const DateTimeArray* critDatesIn,
               double minGap, bool equalTime, int numOfPrice){

    static const string method = "FD1D::Setup";
    try
    {
        Clear();
        
        TimeMetricConstSP metric = GetTimeMetric();   
        
        DateTimeArray critDates;

        // assign user supplied critical dates
        if (critDatesIn)
        {
            for (int i=0; i<critDatesIn->size(); i++)
	            critDates.push_back((*critDatesIn)[i]);
        }

        CTermStructure v_term; 
        CalcV2Term(valDate, segDates[0], segDates[segDates.size()-1], v_term);
        // int StepsPerYear = (int)(timeSteps/(double)(segDates[segDates.size()-1].getDate() - segDates[0].getDate())*365.);
        // int stepsPerYearToUse = (stepsPerYear==-1)?100:stepsPerYear;
        // create time line
        TimePts.CreateTimeLine(segDates, density, metric, stepsPerYearToUse, minGap, false, v_term, critDates, equalTime);
        
        numPriceArrays = numOfPrice;

    } catch (exception& e){
        throw ModelException(&e, method);
    }
   
}

void FD1FGeneric::preProcessGrid(const double* s, int start, int end)
{}

IModel::WantsRiskMapping FD1FGeneric::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

class FD1FGenericHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
		clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FD1FGeneric, clazz);
        SUPERCLASS(CModel);
        FIELD(stepsPerYear, "number of steps per year for T");
        FIELD(stockSteps, "number of steps for S");
        FIELD(TruncationStd, "number of standard deviations to build stock grid");
        FIELD(useFwdGrid, "build grid along stock (false) or along forward (true)");
        FIELD_MAKE_OPTIONAL(useFwdGrid);
        FIELD(varMethod, "0 for straight var, 1 for exp(var)-1");
        FIELD_MAKE_OPTIONAL(varMethod);
        FIELD(DEBUG_UseCtrlVar, "For DR use only : true(default)=use control variate when possible");
        FIELD_MAKE_OPTIONAL(DEBUG_UseCtrlVar);
        FIELD(DEBUG_SameGridDelta, "For DR use only : true(default)=same grid tweak for delta/gamma");
        FIELD_MAKE_OPTIONAL(DEBUG_SameGridDelta);
        FIELD(DEBUG_SameGridVegaRho, "For DR use only : false(default)=new grid tweak for vega, rho, etc.");
        FIELD_MAKE_OPTIONAL(DEBUG_SameGridVegaRho);
        FIELD(gridType, "0 for sin hyperbolic, 1 for linear, 2 for exponential");
        FIELD_MAKE_OPTIONAL(gridType);
        FIELD(stepsPerYearToUse, "transient");
        FIELD_MAKE_TRANSIENT(stepsPerYearToUse);
        FIELD(stockStepsToUse, "transient");
        FIELD_MAKE_TRANSIENT(stockStepsToUse);
    }
};


CClassConstSP const FD1FGeneric::TYPE = CClass::registerClassLoadMethod(
    "FD1FGeneric", typeid(FD1FGeneric), FD1FGenericHelper::load);

CClassConstSP const FD1FGeneric::IIntoProduct::TYPE =
CClass::registerInterfaceLoadMethod("FD1FGeneric::IIntoProduct",
                                    typeid(FD1FGeneric::IIntoProduct), 0);


DRLIB_END_NAMESPACE
