//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RangeAccrue.cpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/RangeAccrue.hpp"

DRLIB_BEGIN_NAMESPACE

// Make RangeAccrue Class
// 1.  set valueDate
// 2.  all range coupon is adjusted to be present value.
// 3.  return RangeAccrue Class.
RangeAccrue::RangeAccrue(const DateTime    baseDate, 
                         const YieldCurve* discount,
                         RangeAccrueMaker* maker): maker(maker),
                                                   pvAmounts(maker->amounts.size(), 0.) {
    valueDate = baseDate;       // set valueDate
    remainedValue = 0.0;        // initialize the remained value.
    for (int iStep = 0; iStep<maker->amounts.size(); iStep++){
        if (baseDate < maker->paymentDates[iStep])
            pvAmounts[iStep] = discount->pv(baseDate,maker->paymentDates[iStep])*maker->amounts[iStep];
        else
            pvAmounts[iStep] = 0.0;  // already paid coupon are not included.
    }
}; 

// Set bool array which tells the time point is accrue or not.
// the size of isMonitorStep must be same to engine's time step (MC simulation steps or Tree steps).
void RangeAccrue::setIsMonitorStep(const DateTimeArray& timeStepDate){
    static const string method = "RangeAccrue::setIsMonitorStep";
    try{
        int i=0, j=0;
        int numSteps = timeStepDate.size();
        int numMon  = maker->monitoringDates.size();
        isMonitorStep.resize(numSteps);
        for (i=0; i<numSteps; i++){
            if (maker->monitoringDates[j] == timeStepDate[i]){
                isMonitorStep[i] = true;
                if (j++ >= numMon)
                    break;      // no more range monitoring date
            }
            else
                isMonitorStep[i] = false;
        }
        // validate all range monitoring date are found on timeStepDate
        if ( j != numMon)
            throw ModelException(method,"Not all rangeAccrue monitoringDate are found in simulation time points");
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
}

// for MC product, calculate the range coupon value.
double RangeAccrue::getValue(const IPathGenerator*  pathGen,            //(I) simulated path
                             const int              endStep,            //(I) calculate up to endStep
                             IAggregateSP&          assetBasketForRange,//(O) each asset performance at each step
                             SimpleDoubleArray&     assetCompsForRange  //(O) basket defenition for range accrue
                             ){
    static const string method = "RangeAccrue::getValue";
    try{
        double value =remainedValue;        // start from the not yet paid coupon.
        int iStep, iAsset;
        int nbAssets=pathGen->NbSimAssets();

        if (endStep > isMonitorStep.size()){
            throw ModelException(method,"the endStep = " + Format::toString(endStep) + "is over the array size "
                                        + Format::toString(isMonitorStep.size()) + "of isMonitorStep" );
        }
        // add the range coupon before endStep
        for(iStep = pathGen->begin(0); iStep<endStep; iStep++){
            if (isMonitorStep[iStep]){
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetCompsForRange[iAsset] = pathGen->Path(iAsset, 0)[iStep]
                                                /pathGen->refLevel(iAsset, 0);
                }
                value += getValue(assetBasketForRange->aggregate(), iStep);
            }
        }
        return value;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }

}
/*
double RangeAccrue::getValue(const DoubleArray& levels, const int endStep) {
    // not tested, but for Tree case.
    double value = remainedValue;
    // add the range coupon before endStep
    for (int iStep=0;iStep<endStep;iStep++) {
        value += getValue(levels[iStep], iStep);
    }
    return value;
};*/

double RangeAccrue::getValue(const double level, const int iStep, const bool isPVed) {
    static const string method = "RangeAccrue::getValue";
    try{
        double value;
        if (maker->lowLevels[iStep] < level && level < maker->highLevels[iStep]){
            value = maker->isInside ? (isPVed ? pvAmounts[iStep] : maker->amounts[iStep]) : 0.0;
        }
        else{
            value = maker->isInside ? 0.0 : (isPVed ? pvAmounts[iStep] : maker->amounts[iStep]);
        }
        return value;
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
};

// make the known cash flow.  
// calculate the already monitored coupon to be paid in future!!
void RangeAccrue::makeKnownCashFlow(const IPathGenerator*  pathGen,            //(I) simulated path
                             const int              endStep,            //(I) calculate up to endStep
                             const bool             hasFuture,          //(I) if hasFuture = false, no need to calc Ramained Value
                             IAggregateSP&          assetBasketForRange,//(O) each asset performance at each step
                             SimpleDoubleArray&     assetCompsForRange  //(O) basket defenition for range accrue
                             ){         
    static const string method = "RangeAccrue::makeKnownCashFlow";
    try{
        double value =0.0;    
        int iStep, iAsset, j=0;
        int nbAssets=pathGen->NbSimAssets();
        
        // initialization
        knownCFL = CashFlowArray(0);
        
        // add the range coupon before endStep
        for(iStep = pathGen->begin(0); iStep<endStep; iStep++){
            if (isMonitorStep[iStep]){
                for(iAsset=0; iAsset<nbAssets; iAsset++) {
                    assetCompsForRange[iAsset] = pathGen->Path(iAsset, 0)[iStep]
                                                /pathGen->refLevel(iAsset, 0);
                }
                value = getValue(assetBasketForRange->aggregate(), iStep, false);   // no need to PVed value.
                knownCFL.push_back(CashFlow(maker->paymentDates[j], value));
                if(maker->paymentDates[j] > valueDate)
                    remainedValue += getValue(assetBasketForRange->aggregate(), iStep);   // set remainedValue. Need to be pved.
                j++;
            }
        }
        
        // if option terminated, pathGen->begin(0) returns 0.
        // So it will calculate all remained coupons at getValues.
        // thus, we need to reset the remainedValue, otherwise
        // it will double count.
        if (!hasFuture)
            remainedValue = 0.0;        
        CashFlow::aggregate(knownCFL);
    }
    catch (exception& e)
    {            
        throw ModelException(e, method);
    }
};

DateTimeArray RangeAccrueMaker::getMonitorDates(){
    return monitoringDates;
};

// validation
void RangeAccrueMaker::validatePop2Object() {
    static const string method = "RangeAccrueMaker::validatePop2Object";
    try {
        // not wmpty data
        if (monitoringDates.size() ==0 ) { 
            throw ModelException(method, "there is no monitoring dates ");
        }
        
        // the input arrays have all the same size
        if (monitoringDates.size() != paymentDates.size() ) { 
            throw ModelException(method,
                                 "number of monitoring dates (" + 
                                 Format::toString(monitoringDates.size()) + 
                                 ") must be the same as number of payment dates (" +
                                 Format::toString(paymentDates.size()) + ").");
        }
        
        if (monitoringDates.size() != amounts.size() ) { 
            throw ModelException(method,
                                 "number of monitoring dates (" + 
                                 Format::toString(monitoringDates.size()) + 
                                 ") must be the same as number of amounts (" +
                                 Format::toString(amounts.size()) + ").");
        }
        
        if (monitoringDates.size() != lowLevels.size() ) { 
            throw ModelException(method,
                                 "number of monitoring dates (" + 
                                 Format::toString(monitoringDates.size()) + 
                                 ") must be the same as number of low levels (" +
                                 Format::toString(lowLevels.size()) + ").");
        }
        
        if (monitoringDates.size() != highLevels.size() ) { 
            throw ModelException(method,
                                 "number of monitoring dates (" + 
                                 Format::toString(monitoringDates.size()) + 
                                 ") must be the same as number of high levels (" +
                                 Format::toString(highLevels.size()) + ").");
        }
        
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

class RangeAccrueMakerHelper {
public:
    


    static void load(CClassSP& clazz){
        REGISTER(RangeAccrueMaker, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultRangeAccrueMaker);
        FIELD(isInside, "is Inside monitoring?");
        FIELD(monitoringDates, "list of range monitoring date.");
        FIELD(paymentDates, "list of payment date of range accrued coupons");
        FIELD(amounts, "list of range accrue coupon amounts");
        FIELD(lowLevels, "lower level or range");
        FIELD(highLevels, "higher level or range");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Range Accrue Maker");
    }
    static IObject* defaultRangeAccrueMaker(){
        return new RangeAccrueMaker();
    }
};

CClassConstSP const RangeAccrueMaker::TYPE = CClass::registerClassLoadMethod(
    "RangeAccrueMaker", typeid(RangeAccrueMaker), RangeAccrueMakerHelper::load);
bool  RangeAccrueLoad() {
    return (RangeAccrueMaker::TYPE != 0);
   }


RangeAccrue* RangeAccrueMaker::getDiscounted(const DateTime baseDate, 
                                             const YieldCurve* discount){
    return new RangeAccrue(baseDate,
                           discount,
                           this);
}

DRLIB_END_NAMESPACE

    

