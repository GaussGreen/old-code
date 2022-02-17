//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DblBarrier.cpp
//
//   Description   double barrier instrument. 
//                 A KO barrier takes priority over KI barrier.
//                 If both are KI barriers, one KI counts as in.
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/Barrier.hpp"
#include "edginc/DblBarrier.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/InstrumentUtil.hpp"
#include "edginc/Tree1fLV.hpp"
#include "edginc/Tree1fLN.hpp"
#include "edginc/PhysicalDelivery.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/FD1DLV.hpp"

DRLIB_BEGIN_NAMESPACE

// some helper functions
// validate Barrier
void validateBar(ScheduleSP BarSched, ScheduleSP RebateSched){
    static const string method("validateBar");
    try {
        if (BarSched->getInterp() == Schedule::INTERP_NONE){
            if (BarSched->isSameSchedule(RebateSched->getDates(), true) == false)
                throw ModelException(method, "Barrier and Rebate should be same schedule for Interp='N' case");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** add barrier dates to critical date list */
void CDblBarrier::addCritBarDates(const ScheduleSP bar, const DateTime& valueDate, const DateTime& matDate,
                            CTree1f* tree1f, DateTimeArray& critDates)
{
    // place critical dates
    if (bar->getInterp() != Schedule::INTERP_NONE)
    {// add monitor start and monitor end to critical dates list if needed
        if (bar->firstDate() > valueDate)
            critDates.push_back(bar->firstDate());
        if (bar->lastDate() < matDate)
            critDates.push_back(bar->lastDate());
    }
    else if (tree1f)
    {
        const DateTimeArray& barDates = bar->getDates();
        if (tree1f->GetStepsPerYear() == 0)
            tree1f->SetStepsPerYear(-1); // default to daily step
        else if ( (tree1f->GetStepsPerYear() < 250 && tree1f->GetStepsPerYear() > 0)
                  ||(bar->getInterp() == Schedule::INTERP_NONE))
        {// input is less than one step per day ory on Date Case, add steps around barrier dates
            for (int i=0; i<bar->length(); i++)
            {
                if (barDates[i] > valueDate && barDates[i] <= matDate)
                {
                    // 2 days gap seem to be best for convergence
                    critDates.push_back(barDates[i].rollDate(-2)); 
                    critDates.push_back(barDates[i]);
                    // more test needed for 2 days gap following
                    // critDates.push_back(inst->BarDates[i].rollDate(2)); 
/*                    for (int j = 0; j<=5; j++)
                      {
                      DateTime tmpDate(barDates[i].getDate(),5+j*1440);
                      critDates.push_back(tmpDate);
                      critDates.push_back(tmpDate.rollDate(-1));
                      }*/
                }
            }
        }
    }
    else
    {
        const DateTimeArray& barDates = bar->getDates();
        for (int i=0; i<bar->length(); i++)
        {
            if (barDates[i] > valueDate && barDates[i] <= matDate)
            {
                // 1 days gap seem to be best for convergence
                critDates.push_back(barDates[i].rollDate(-1)); 
                critDates.push_back(barDates[i]);
            }
        }
    }
}

// set up the "can I KO flag" per timepoint for a barrier 
// static to avoid changing the product in header file
static void setKOCondition(vector<bool>&   stepCanKO,
                           const DateTimeArray& stepDates,
                           CAssetConstSP   Underlier) {
    static const string method("setKOCondition");
    try {
        // all points are potentially KO 
        // the simple rule is: a timepoint is KO if there is a chance to KO 
        // any time after the previous point
        if (AssetUtil::getHoliday(Underlier.get())->isBusinessDay(stepDates[0])) {
            stepCanKO[0] = true;
        }
        else {
            stepCanKO[0] = false;
        }

        // loop through steps, maturity step cannot KO
        for (int iStep = 1; iStep < (int)stepDates.size()-1; iStep++)
        {
            if (AssetUtil::hasTradingTime(Underlier,
                                          stepDates[iStep-1],
                                          stepDates[iStep]))
            {
                stepCanKO[iStep] = true;
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// tweak the spot and price
static double tweakSpotAndPrice(double spotTweak, const IModel* model, const CInstrument* inst) {
    static string method = "tweakSpotAndPrice";

    IModelSP modelCopy(copy(model));  
    CInstrumentSP tweaked(copy(inst));
    TweakGroupSP tweakGroup(new TweakGroup(tweaked,modelCopy));
    SensMgr sensMgr(tweakGroup);

    // shifth by 1% of the barrier level and
    SpotLevelSP shift(new SpotLevel(spotTweak));

/* removed as now it finds two assets - one in Model and one in Instrument
    OutputNameArrayConstSP names = sensMgr.allNames(shift.get());

    // only one underlying allowed
    if (names->size()>1) {
        throw ModelException(method, "only one underlying allowed");
    }
*/

    CControlSP ctrl(Control::makeFromFlags("", 0.0));

    shift->findAndShift(tweakGroup, OutputNameConstSP());
    ResultsSP tweakResults(tweakGroup->getModel()->Run(tweakGroup->getInstrument(), ctrl.get()));

    // returns the difference of prices
    return tweakResults->retrievePrice();
}

////// end of helper functions

// helpers
class CDblBarrierHelper {
public:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(CDblBarrier, clazz);
        SUPERCLASS(Generic1Factor);
        EMPTY_SHELL_METHOD(defaultDblBarrier);
        // same as vanilla
        IMPLEMENTS(FD1F::IIntoProduct);
        IMPLEMENTS(FD1FGeneric::IIntoProduct);
        IMPLEMENTS(FDModel::IIntoProduct);
        IMPLEMENTS(ISensitiveStrikes);
        //IMPLEMENTS(Theta::Shift);
        IMPLEMENTS(LastSensDate);
        IMPLEMENTS(LegalTerms::Shift);
        FIELD(PayoffMode, "CALL, PUT, BINARY, FORWARD");
        FIELD(exerciseSchedule,        "Exercise Schedule");
        FIELD(canExerciseEarly, "Can option be exercised early");
        FIELD(spotAtMaturity,   "underlying spot level when exercised");
        FIELD_MAKE_OPTIONAL(spotAtMaturity);
        FIELD(isExercised, "Indicates whether option has been exercised");
        FIELD_MAKE_OPTIONAL(isExercised);
        FIELD(dateExercised,        "Date on which option has been exercised");
        FIELD_MAKE_OPTIONAL(dateExercised);
        // barrier parts
        FIELD(UpperBarrier, "upper barrier schedule");
        FIELD_MAKE_OPTIONAL(UpperBarrier);
        FIELD(LowerBarrier, "lower barrier schedule");
        FIELD_MAKE_OPTIONAL(LowerBarrier);
        FIELD(UpperBarType, "upper barrier type: NA, KI or KO");
        FIELD(LowerBarType, "lower barrier type: NA, KI or KO");
        FIELD(IntraDayMonitor, "true for hi/low continuous monitoring, false if once a day");
        FIELD(UpperRebate, "upper rebate schedule");
        FIELD_MAKE_OPTIONAL(UpperRebate);
        FIELD(LowerRebate, "lower rebate schedule");
        FIELD_MAKE_OPTIONAL(LowerRebate);
        FIELD(RebateAtMat, "true if rebate paid at maturity, false if paid at touch");
        FIELD(RebateNotScaled, "true if rebate is not scaled for fwd start option, false=it is scaled");

        FIELD(upperEcoBarrier, "upper economic barrier");
        FIELD_MAKE_OPTIONAL(upperEcoBarrier);
        FIELD(lowerEcoBarrier, "lower economic barrier");
        FIELD_MAKE_OPTIONAL(lowerEcoBarrier);

        FIELD(UpperBarBreached, "true if upper barrier breached");
        FIELD(UpperBarBreachDate, "date upper barrier breached");
        FIELD_MAKE_OPTIONAL(UpperBarBreachDate);
        FIELD(LowerBarBreached, "true if lower barrier breached");
        FIELD(LowerBarBreachDate, "date lower barrier breached");
        FIELD_MAKE_OPTIONAL(LowerBarBreachDate);
        FIELD(BarrierDependence, "barrier dependence: KI_KEEP_KO, KI_CANCEL_KO");
        FIELD_MAKE_OPTIONAL(BarrierDependence);
        // As barrier adjustment is applied whent IntraDayMonitor = false, 
        // the Dependence is also active when IntraDayMonitor = false.
        FIELD(MonitoringDependence, "When IntraDayMonitor=false, it's applied to : BOTH(deafault), UPPER, LOWER.");
        FIELD_MAKE_OPTIONAL(MonitoringDependence);

        // used internally
        FIELD(isCall, "Is it a call option");
        FIELD_MAKE_TRANSIENT(isCall);

        FIELD(DEBUG_num_segments, "For DR use only. Number of time line segments, default=1.");
        FIELD_MAKE_OPTIONAL(DEBUG_num_segments);

        FIELD(DEBUG_SmoothBarrierWidth, "For DR use Only.  0(default) is no use.  double value are used as SmoothWidth, e.g. 0.05 (5%)");
        FIELD_MAKE_OPTIONAL(DEBUG_SmoothBarrierWidth);
    }

    static IObject* defaultDblBarrier(){
        return new CDblBarrier();
    }
};


CClassConstSP const CDblBarrier::TYPE = CClass::registerClassLoadMethod(
    "DblBarrier", typeid(CDblBarrier), CDblBarrierHelper::load);

bool  CDblBarrierLoad() {
    return (CDblBarrier::TYPE != 0);
}



bool CDblBarrier::avoidVegaMatrix(const IModel* model)
{
    if (CTree1fLV::TYPE->isInstance(model)) {
        return true; // do pointwise instead
    }
    return false;
}

// for reflection 
CDblBarrier::CDblBarrier(CClassConstSP clazz): Generic1Factor(clazz), isCall(false) {
    canExerciseEarly = false;
    LowerBarrier = ScheduleSP(   );
    UpperBarrier = ScheduleSP(   );
    LowerRebate = ScheduleSP(   );
    UpperRebate = ScheduleSP(   );
    BarrierDependence = "KI_KEEP_KO";
    MonitoringDependence = "BOTH";
    DEBUG_num_segments = 1;
    isExercised = false;
    DEBUG_SmoothBarrierWidth = -1.0;
    DEBUG_num_segments       = 1;

    spotAtMaturity           = 0.0;
    isExercised              = false;
    LowerBarBreachedDynamic         = false;
    UpperBarBreachedDynamic         = false;
}


// constructor
CDblBarrier::CDblBarrier(): Generic1Factor(TYPE), isCall(false)
{
    canExerciseEarly = false;
    LowerBarrier = ScheduleSP(   );
    UpperBarrier = ScheduleSP(   );
    LowerRebate = ScheduleSP(   );
    UpperRebate = ScheduleSP(   );
    BarrierDependence = "KI_KEEP_KO";
    MonitoringDependence = "BOTH";
    DEBUG_num_segments = 1;
    isExercised = false;
    DEBUG_SmoothBarrierWidth = -1.0;

    spotAtMaturity           = 0.0;
    isExercised              = false;
    LowerBarBreachedDynamic         = false;
    UpperBarBreachedDynamic         = false;
}

CDblBarrier::CDblBarrier(const DateTime&        valueDate, 
                         const DateTime&        startDate,
                         bool                   fwdStarting,
                         bool                   oneContract,
                         double                 notional,
                         double                 initialSpot,
                         InstrumentSettlementSP instSettle,
                         InstrumentSettlementSP premiumSettle,
                         YieldCurveWrapper      discount,
                         CAssetWrapper          asset,
                         const string&          ccyTreatment,
                         const string&          PayoffMode,
                         ScheduleSP             exerciseSchedule,
                         bool                   canExerciseEarly,
                         ScheduleSP             upperBarrier,
                         const string&          upperBarType,    
                         ScheduleSP             upperRebate,
                         ScheduleSP             lowerBarrier,
                         const string&          lowerBarType,
                         ScheduleSP             lowerRebate,
                         bool                   intraDayMonitor,
                         bool                   rebateAtMat,
                         bool                   rebateNotScaled,
                         const string&          barrierDependence,
                         const string&          monitoringDependence)
    : Generic1Factor(TYPE), isCall(false)
{
    this->valueDate         = valueDate;
    this->startDate         = startDate;
    this->fwdStarting       = fwdStarting;
    this->oneContract       = oneContract;
    this->notional          = notional;
    this->initialSpot       = initialSpot;
    this->instSettle        = instSettle;
    this->premiumSettle     = premiumSettle;
    this->discount          = YieldCurveWrapper(copy(discount.get()));
    this->asset             = CAssetWrapper(copy(asset.get()));
    this->ccyTreatment      = ccyTreatment;
    this->PayoffMode        = PayoffMode;
    this->exerciseSchedule  = exerciseSchedule;
    this->canExerciseEarly  = canExerciseEarly;
    this->UpperBarrier      = upperBarrier;
    this->UpperBarType      = upperBarType;
    this->UpperRebate       = upperRebate;
    this->LowerBarrier      = lowerBarrier;
    this->LowerBarType      = lowerBarType;
    this->LowerRebate       = lowerRebate;

    this->IntraDayMonitor   = intraDayMonitor;
    this->RebateAtMat       = rebateAtMat;
    this->RebateNotScaled   = rebateNotScaled;
    this->BarrierDependence = barrierDependence;
    this->MonitoringDependence = monitoringDependence;
    
    DEBUG_num_segments       = 1;
    DEBUG_SmoothBarrierWidth = -1.0;
    spotAtMaturity           = 0.0;
    isExercised              = false;
    UpperBarBreached                = false;
    LowerBarBreached                = false;
    LowerBarBreachedDynamic         = false;
    UpperBarBreachedDynamic         = false;
}

CDblBarrier1fProd::CDblBarrier1fProd(const CDblBarrier* instr)
: NumOfPriceArray(3), hasDefPayoff(false), inst(instr)
{}


#define USE_SEGMENTS
/** initialise FD - allow product customisation */
void CDblBarrier1fProd::InitFD(CControl* control)
{
#ifndef USE_SEGMENTS
    // add barrier dates to critical dates)
    if (inst->UpperBarType == "KO" || inst->UpperBarType == "KI")
        AddCritBarDates(inst->UpperBarrier);
    if (inst->LowerBarType == "KO" || inst->LowerBarType == "KI")
        AddCritBarDates(inst->LowerBarrier);
    // set base step used only for default input (StepsPerYear ==0)
    //   model1F->TimePts.SetTreeStepBase(500);

    // ************ below is almost the same as in Vanilla  ***************
    /** customize tree parameters here and set up the tree */
    DateTimeArray segDates;
    segDates.resize(2); // to try using more segments
    // this needs change if fwd start tree has to start today !!!
    if (inst->fwdStarting && inst->startDate>inst->valueDate)
    {
        segDates[0] = inst->startDate;  // to do: this needs to be checked !!!
    }
    else
        segDates[0] = inst->valueDate;

    segDates[1] = inst->exerciseSchedule->lastDate();
    vector<int> density;
    density.resize(1);
    density[0] = 1;
    double minGap = 0; // insert all points
    bool useEqualTime = false; // false = equal variance steps
 
#else
        
    double minGap = -1; // default value
    int i;
        
    DateTimeArray upBar, downBar;
    if (inst->UpperBarrier.get() && inst->UpperBarType != "NA") {
        upBar = inst->UpperBarrier->getDates();
    }
    if (inst->LowerBarrier.get() && inst->LowerBarType != "NA") {
        downBar = inst->LowerBarrier->getDates();
    }
    DateTimeArray allDates = DateTime::merge(
        upBar, 
        downBar);           


    DateTimeArray segDates;
    segDates.resize(1);
    if (inst->fwdStarting && inst->startDate>inst->valueDate)
    {
        segDates[0] = inst->startDate;  // to do: this needs to be checked !!!
    }
    else
        segDates[0] = inst->valueDate;
    
    for (i=0; i<allDates.size(); i++) {
        if (allDates[i] > inst->valueDate && allDates[i] < inst->exerciseSchedule->lastDate()) {
            segDates.push_back(allDates[i]);
        }
    }
        
    segDates.push_back(inst->exerciseSchedule->lastDate());

    vector<int> density;
    density.resize(segDates.size()-1);
    for (i=0; i<segDates.size()-1; i++) {
        density[i] = 1;
    }
    bool useEqualTime = true;
 
#endif
    model1F->TimePts.SetTreeStepBase(500);
    // all exercise dates are copied to critical dates
    DateTimeArray critDates = inst->exerciseSchedule->getDates();
    // remove exercise date from crit date
    critDates.erase(critDates.end()-1);
    // add div event dates if needed
    EventAssetMove divEvent;
    DateTimeArraySP divCritDates;
    if (inst->canExerciseEarly)
    {// American exercise or dollar div interp treatment
        // create div events, this call is very expensive for basket with lots of div dates
        // should only need once for pricing call but need to think how to store/copy for tweaks
        int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates

        if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                     segDates[0], 
                                     segDates[1],
                                     numDivs,
                                     divEvent)){
            // calculate critical dates
            divCritDates = divEvent.getCritDate(numDivs, inst->isCall);
            for (int i=0; i<divCritDates->size(); i++)
                critDates.push_back((*divCritDates)[i]);
        }
    }

 

    // call tree step set up routine
    if( fdModel )
        fdModel->Setup(inst->valueDate, segDates, density, &critDates, 
                   minGap, useEqualTime, NumOfPriceArray);
    else
        genericFDModel->Setup(inst->valueDate, segDates, density, &critDates, 
                   minGap, useEqualTime, NumOfPriceArray);
}

double CDblBarrier1fProd::scalePremium(const double& fairValue)
{
    double fwdAtStart = 0.0;
    if (inst->fwdStarting)
    {
        fwdAtStart = inst->asset->fwdValue(inst->startDate);
    }

    double scalingFactor = InstrumentUtil::scalePremium(
        inst->oneContract,
        inst->fwdStarting,
        inst->notional,
        fwdAtStart,
        inst->initialSpot);

    return fairValue * scalingFactor;
}

/** initialise product specific data */
void CDblBarrier1fProd::InitProd()
{
    static const string method("CDblBarrier1fProd::InitProd");
    try {
        int i;
        const CTimeLine* timeline = &model1F->TimePts;
        int numStep = timeline->NumOfStep;
        double fwdAtStart = 1.0, scaleFactor = 1.0;

        // set discrete adjustment factor. Estimate trading days per year using one year from now.
        // use 364 because inclusive counting
        DateTime endDate(inst->valueDate.getDate()+364, inst->valueDate.getTime()); 
        VolDaysPerYear = 252.0; // model1F->GetTimeMetric()->volDays(inst->valueDate, endDate);

        stepStrike.resize(numStep+1);
        stepCanExercise.resize(numStep+1);
        // ask tree to decide first about steps that can exercise
        AssetUtil::setStepExercise(stepCanExercise,
                                 model1F->TimePts.StepDates,
                                 inst->exerciseSchedule,
                                 inst->canExerciseEarly,
                                 inst->asset.getSP());

        bool canInterpExSched = (inst->exerciseSchedule->length() > 1);
        // get spot at start if needed
        if (inst->fwdStarting && timeline->StepDates[0]>inst->valueDate)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
            if (!inst->RebateNotScaled)
                scaleFactor = fwdAtStart; 
        }

        // get strikes
        stepStrike[numStep] = fwdAtStart*inst->exerciseSchedule->lastValue();

        for (i=0; i<numStep; i++)
        {
            if (stepCanExercise[i])
            {
                if (canInterpExSched)
                {
                    stepStrike[i] = fwdAtStart*inst->exerciseSchedule->interpolate(timeline->StepDates[i]);
                }
                else
                {
                    stepStrike[i] = stepStrike[numStep]; // strike level is constant for american with no schedule
                }
            }
        }

        // barrier parts
        if (inst->UpperBarType == "KO")
            UType = KO;
        else if (inst->UpperBarType == "KI")
            UType = KI;
        else
            UType = NA;

        if (inst->LowerBarType == "KO")
            LType = KO;
        else if (inst->LowerBarType == "KI")
            LType = KI;
        else
            LType = NA;

        if (inst->BarrierDependence == "KI_KEEP_KO")
            BarDepend = KI_KEEP_KO;
        else if (inst->BarrierDependence == "KI_CANCEL_KO")
            BarDepend = KI_CANCEL_KO;
        else if (inst->BarrierDependence == "ONCE_TOUCH")
            BarDepend = ONCE_TOUCH;
        else if (inst->BarrierDependence == "TWO_TOUCH")
            BarDepend = TWO_TOUCH;
        else
            throw ModelException("CDblBarrier1fProd::InitProd", inst->BarrierDependence + "unkown Barrier dependence");

        if (inst->MonitoringDependence == "BOTH")
            MonDepend = BOTH;
        else if (inst->MonitoringDependence == "UPPER")
            MonDepend = UPPER;
        else if (inst->MonitoringDependence == "LOWER")
            MonDepend = LOWER;
        else
            throw ModelException("CDblBarrier1fProd::InitProd", inst->MonitoringDependence + "unkown Barrier dependence");

        if (inst->PayoffMode == "CALL")
            PayMode = CALL;
        else if (inst->PayoffMode == "PUT")
            PayMode = PUT;
        else if (inst->PayoffMode == "BINARY")
            PayMode = BINARY;
        else if (inst->PayoffMode == "FORWARD")
            PayMode = FORWARD;
        else
            throw ModelException("CDblBarrier1fProd::InitProd", inst->PayoffMode + "unkown payoff mode");

        UBarInput.resize(numStep+1);
        LBarInput.resize(numStep+1);
        UBarRebate.resize(numStep+1);
        LBarRebate.resize(numStep+1);

        // init barriers to far away levels
        double spot = inst->asset->getSpot();
        for(i=0; i<=numStep; i++)
        {
            UBarInput[i] = spot/FP_MIN;
            LBarInput[i] = - 2.0*FP_MIN;
        }
        // init rebate to 0
        for(i=0; i<=numStep; i++)
            UBarRebate[i] = LBarRebate[i] = 0.0;

        const DateTime& valDate= inst->valueDate;
//    const DateTime& matDate= inst->exerciseSchedule->lastDate();
        string uInterp;
        string lInterp;
//    const DateTime& matDate= inst->exerciseSchedule->lastDate();
        if (inst->UpperBarType != "NA")
            uInterp =  inst->UpperBarrier->getInterp();
        if (inst->LowerBarType != "NA")
            lInterp =  inst->LowerBarrier->getInterp();
        // for continuous KI check if it's in already
        if (inst->IntraDayMonitor)
        {
            if (inst->UpperBarType != "NA")
                if (valDate>= inst->UpperBarrier->firstDate() &&  valDate<= inst->UpperBarrier->lastDate() &&
                    (uInterp!= Schedule::INTERP_NONE && spot>=inst->UpperBarrier->interpolate(valDate) ||
                     uInterp== Schedule::INTERP_NONE && inst->UpperBarrier->coversDateRange(valDate, valDate, true)))
                {
                    if (UType == KI)
                    {
                        UType = NA; // make it not applicable
                        if ((BarDepend == KI_CANCEL_KO && LType == KO) || (BarDepend == ONCE_TOUCH && LType == KI))
                            LType = NA;
                    }
                }
            if (inst->LowerBarType != "NA")
                if (valDate>= inst->LowerBarrier->firstDate() &&  valDate<= inst->LowerBarrier->lastDate() &&
                    (lInterp!= Schedule::INTERP_NONE && spot<=inst->LowerBarrier->interpolate(valDate) ||
                     lInterp== Schedule::INTERP_NONE && inst->LowerBarrier->coversDateRange(valDate, valDate, true)))
                {
                    if (LType == KI)
                    {
                        LType = NA; // make it not applicable
                        if ((BarDepend == KI_CANCEL_KO && UType == KO) || (BarDepend == ONCE_TOUCH && UType == KI))
                            UType = NA;
                    }
                }
        }

        // set up barrier and rebate levels at each step
        if (UType != NA)
        {
            Barrier::setStepLevel(inst->UpperBarrier, fwdAtStart, timeline->StepDates, UBarInput);
            if (!!inst->UpperRebate)
            {
                Barrier::setStepLevel(inst->UpperRebate, scaleFactor, timeline->StepDates, UBarRebate);
                if (inst->RebateAtMat)
                {
                    for (i=0; i<= numStep; i++)
                    {
                        if(!Maths::isZero(UBarRebate[i]))
                            UBarRebate[i] *= inst->discount->pv(timeline->StepDates[i],
                                                                    inst->instSettle->settles(inst->exerciseSchedule->lastDate(), inst->asset.get()));
                    }
                }
                /*  Need to be considered about Rebate settlment PV adjsutment
                else{
                    for (i=0; i<= numStep; i++)
                    {
                        if(!Maths::isZero(UBarRebate[i])){
                            if (CashSettleDate::TYPE->isInstance(inst->instSettle.get())) {
                                // It's really wired case, that rebate pay at hit, but settlement is fixed.
                                // It may better to disallow this case, but it's too late.  
                                // Many live position could have this case.
                                // Those case are just calculated as pay immidiatelly (T+0 days)!
                            }
                            else                            
                                UBarRebate[i] *= inst->discount->pv(timeline->StepDates[i],
                                                                    inst->instSettle->settles(timeline->StepDates[i], inst->asset.get()));
                        }
                    }
                }*/
            }
        }
        if (LType != NA)
        {
            Barrier::setStepLevel(inst->LowerBarrier, fwdAtStart, timeline->StepDates, LBarInput);
            if (!!inst->LowerRebate)
            {
                Barrier::setStepLevel(inst->LowerRebate, scaleFactor, timeline->StepDates, LBarRebate);
                if (inst->RebateAtMat)
                {
                    for (i=0; i<= numStep; i++)
                    {
                        if(!Maths::isZero(LBarRebate[i]))
                            LBarRebate[i] *= inst->discount->pv(timeline->StepDates[i],
                                                                    inst->instSettle->settles(inst->exerciseSchedule->lastDate(), inst->asset.get()));
                    }
                }
                /*  Need to be considered about Rebate settlment PV adjsutment
                else{
                    for (i=0; i<= model->getLastStep(); i++)
                    {
                        if(!Maths::isZero(LBarRebate[i])){
                            if (CashSettleDate::TYPE->isInstance(inst->instSettle.get())) {
                                // It's really wired case, that rebate pay at hit, but settlement is fixed.
                                // It may better to disallow this case, but it's too late.  
                                // Many live position could have this case.
                                // Those case are just calculated as pay immidiatelly (T+0 days)!
                            }
                            else
                                LBarRebate[i] *= inst->discount->pv(timeline->StepDates[i],
                                                                    inst->instSettle->settles(timeline->StepDates[i], inst->asset.get()));
                        }
                    }
                }*/
            }
        }

        // disallow non-zero Rabate for KI barrier.
        if (UType == KI){
            for (i=0; i<= numStep; i++)
            {
                if(!Maths::isZero(UBarRebate[i]))
                    throw ModelException(method, "Upper Barrier is KI. non-zero rebate for KI is disallowed.");
            }                    

        }
        if (LType == KI){
            for (i=0; i<= numStep; i++)
            {
                if(!Maths::isZero(LBarRebate[i]))
                    throw ModelException(method, "Lower Barrier is KI. non-zero rebate for KI is disallowed.");
            }                    

        }


        // now we need to set up an indication of whether or not we can KO at a step
        // UBarInput/LBarInput are set up with the "right" barrier levels so we can
        // interrogate them at each step (i.e. KO on dates only has zero/infinite
        // barriers on dates other than barrier dates) - HOWEVER we want to forbid
        // KO on holidays/weekends (c.f. exercise), so set up a vector of flags
        // indicating if we can check the KO condition
        stepCanKOUp.resize(numStep+1);
        stepCanKODown.resize(numStep+1);

        if (UType == NA) {
            for(i = 0; i <= numStep; i++) {        
                stepCanKOUp[i] = false;
            }
        }
        else {     
            setKOCondition(stepCanKOUp, timeline->StepDates, GetAssetRef());
        }
     
        if (LType == NA) {
            for(i = 0; i <= numStep; i++) {        
                stepCanKODown[i] = false;
            }
        }
        else {       
            setKOCondition(stepCanKODown, timeline->StepDates, GetAssetRef());
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calculate barriers for fd  */
void CDblBarrier1fProd::preCalcFD(int step, int idx, int pStart, int pEnd)
{
    // set barrier, make adjustment for once a day monitoring if needed
    // to do: forward start adjustment
    UpperBar = UBarInput[step];
    LowerBar = LBarInput[step];

    // taking care of discrete monitoring, treat as daily monitor
    if (!inst->IntraDayMonitor)
    {
        vector<double> vol;
        // adjust barrier if needed
        if (UType != NA)
        {
            model1F->GetStepVol(step, vol, &UpperBar, 0, 0); // get vol at barrier
            Barrier::BarrierAdjustment(vol[0], true, UpperBar);
        }
        if (LType != NA)
        {
            model1F->GetStepVol(step, vol, &LowerBar, 1, 1); // get vol at barrier
            Barrier::BarrierAdjustment(vol[0], false, LowerBar);
        }
    }

    int NumOfStep = fdModel?fdModel->TimePts.NumOfStep:genericFDModel->TimePts.NumOfStep;
    double *upBarrierNew = fdModel?fdModel->fdEngine.upBarrierNew:genericFDModel->fdEngine.upBarrierNew;
    double *upPayoutNew = fdModel?fdModel->fdEngine.upPayoutNew:genericFDModel->fdEngine.upPayoutNew;
    double *upPayoutDeltaNew = fdModel?fdModel->fdEngine.upPayoutDeltaNew:genericFDModel->fdEngine.upPayoutDeltaNew;
    double *upBarrierOld = fdModel?fdModel->fdEngine.upBarrierOld:genericFDModel->fdEngine.upBarrierOld;
    double *upPayoutOld = fdModel?fdModel->fdEngine.upPayoutOld:genericFDModel->fdEngine.upPayoutOld;
    double *upPayoutDeltaOld = fdModel?fdModel->fdEngine.upPayoutDeltaOld:genericFDModel->fdEngine.upPayoutDeltaOld;
    double *downBarrierNew = fdModel?fdModel->fdEngine.downBarrierNew:genericFDModel->fdEngine.downBarrierNew;
    double *downPayoutNew = fdModel?fdModel->fdEngine.downPayoutNew:genericFDModel->fdEngine.downPayoutNew;
    double *downPayoutDeltaNew = fdModel?fdModel->fdEngine.downPayoutDeltaNew:genericFDModel->fdEngine.downPayoutDeltaNew;
    double *downBarrierOld = fdModel?fdModel->fdEngine.downBarrierOld:genericFDModel->fdEngine.downBarrierOld;
    double *downPayoutOld = fdModel?fdModel->fdEngine.downPayoutOld:genericFDModel->fdEngine.downPayoutOld;
    double *downPayoutDeltaOld = fdModel?fdModel->fdEngine.downPayoutDeltaOld:genericFDModel->fdEngine.downPayoutDeltaOld;

    for(int i=pStart; i<=pEnd; i++)
    {
        upBarrierNew[i] = -1;
        upPayoutNew[i] = 0.;
        upPayoutDeltaNew[i] = 0.;
        upBarrierOld[i] = -1;
        upPayoutOld[i] = 0.;
        upPayoutDeltaOld[i] = 0.;
        
        downBarrierNew[i] = -1;
        downPayoutNew[i] = 0.;
        downPayoutDeltaNew[i] = 0.;
        downBarrierOld[i] = -1;
        downPayoutOld[i] = 0.;
        downPayoutDeltaOld[i] = 0.;           
    }
    
    if (step < NumOfStep) {
        if (UType == KO) {
            upBarrierNew[0] = UpperBar;
            upPayoutNew[0] = UBarRebate[step];
            upPayoutDeltaNew[0] = 0.;
        } 
        
        if (LType == KO) {
            downBarrierNew[0] = LowerBar;
            downPayoutNew[0] = LBarRebate[step];
            downPayoutDeltaNew[0] = 0.;
        }

        if (UType == KI || LType == KI)
        {
            if (UType == KO) {
                upBarrierNew[1] = UpperBar;
                upPayoutNew[1] = UBarRebate[step];
                upPayoutDeltaNew[1] = 0.;
            } 

            if (LType == KO) {
                downBarrierNew[1] = LowerBar;
                downPayoutNew[1] = LBarRebate[step];
                downPayoutDeltaNew[1] = 0.;
            } 
        }
    }
}

/** check if knocked out already */
bool CDblBarrier1fProd::CheckDeadInstr()
{
    UpIsOut = inst->UpperBarBreachedDynamic && UType == KO;
    DownIsOut = inst->LowerBarBreachedDynamic && LType == KO;

    if(inst->UpperBarBreachedDynamic && UType == KI)
        UType = NA;

    if(inst->LowerBarBreachedDynamic  && LType == KI)
        LType = NA;

    if (inst->IntraDayMonitor)
    {
        double s = inst->asset->getSpot();
        if (UType == KO && s > UBarInput[0]*(1.0-FP_MIN))
            UpIsOut = true;
        if(LType == KO && s < LBarInput[0]*(1.0+FP_MIN))
            DownIsOut = true;
    }

    bool isDead = (((UpIsOut || DownIsOut) && BarDepend != TWO_TOUCH)
                   || ((UpIsOut && DownIsOut) && BarDepend == TWO_TOUCH));

    return isDead;
}

/** calculate dead value */
void CDblBarrier1fProd::CalcDeadInstr(int step, int bot, int top, double * const * price)
{
    int i;

    bool hasKnockIn = (UType == KI  || LType == KI );

    // discount rebate if needed
    if ( hasKnockIn ) {
        if (UType != NA && !!inst->UpperRebate)
        {
            if (!Maths::isZero(UBarRebate[step])) {
                for (i=-bot; i<=top; i++)
                    price[2][i] = UBarRebate[step]; // index [2] is upper barrier option
            }
        }
        if (LType != NA && !!inst->LowerRebate)
        {
            if (!Maths::isZero(LBarRebate[step])) {
                for (i=-bot; i<=top; i++)
                    price[1][i] = LBarRebate[step]; // index [1] is lower barrier option
            }
        }

        // to do: need to decide which rebate to pay for NO_CANCEL on double KO ?
        if (UpIsOut)
        {
            for (i=-bot; i<=top; i++)
                price[0][i] = price[2][i];
        }
        else if (DownIsOut)
        {
            for (i=-bot; i<=top; i++)
                price[0][i] = price[1][i];
        }
    } else {
        // to do: need to decide which rebate to pay for NO_CANCEL on double KO ?
        if (UpIsOut)
        {
            if (!!inst->UpperRebate) {
                if (!Maths::isZero(UBarRebate[step])) {
                    for (i=-bot; i<=top; i++)
                        price[0][i] = UBarRebate[step];
                }
            }
        } else if (DownIsOut) {
            if (LType != NA && !!inst->LowerRebate) {
                if (!Maths::isZero(LBarRebate[step])) {
                    for (i=-bot; i<=top; i++)
                        price[0][i] = LBarRebate[step];
                }
            }
        }
    }
}

/** product payoff method at steps earlier than maturity */
void CDblBarrier1fProd::PayoffAtMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                    double * const * price)
{
    int j; 

    if (CheckDeadInstr())
    {   
        CalcDeadInstr(step, bot, top, price);
    }

    double settlementPV = inst->instSettle->pvAdjust(inst->exerciseSchedule->lastDate(),
                                                         inst->discount.get(), 
                                                         inst->asset.get());
 
    if (UType == KI || LType == KI)
    {
        for (j=-bot; j<=top; j++)
        {
            if (PayMode == BINARY)
                (price[2])[j] = settlementPV; // binary
            else
                (price[2])[j] = settlementPV*GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD);

            if (UType == KO && s[j] > UpperBar * (1.0-FP_MIN)) // dead on upper barrier
                (price[0])[j] = (price[1])[j] = UBarRebate[step];
            else if(LType == KO && s[j] < LowerBar * (1.0+FP_MIN)) // dead on lower barrier
                (price[0])[j] = (price[1])[j] = LBarRebate[step];
            else
            {
                (price[1])[j] = (price[2])[j];

                if( (UType == KI && s[j] > UpperBar * (1.0-FP_MIN))
                    || (LType == KI && s[j] < LowerBar * (1.0+FP_MIN))
                    || (LType != KI && UType != KI) )
                    (price[0])[j] = (price[1])[j];
                else
                {
                    if (UType == KI)
                        (price[0])[j] = UBarRebate[step]; // failed to KI upper bar
                    else if (LType == KI)
                        (price[0])[j] = LBarRebate[step]; // failed to KI lower bar
                    else
                        throw ModelException("CDblBarrier1fProd::PayoffAtMatD", "There is a bug in here.  Contact DR!!");
                    //validate
                    if (UType == KI && LType == KI && UBarRebate[step] != LBarRebate[step])
                        throw ModelException("CDblBarrier1fProd::PayoffAtMat", "UpperRebate and LowerRebate should be same for Double KI case.");
                }
            }
        }

    }
    else
    {
        for (j=-bot; j<=top; j++)
        {
            if (UType == KO && s[j] > UpperBar * (1.0-FP_MIN)) // dead on upper barrier
                (price[0])[j] = UBarRebate[step];
            else if(LType == KO && s[j] < LowerBar * (1.0+FP_MIN)) // dead on lower barrier
                (price[0])[j] = LBarRebate[step];
            else
            {
                if (PayMode == BINARY)
                    (price[0])[j] = settlementPV; // binary
                else
                    (price[0])[j] = settlementPV*GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD);
            }
        }

    }

    // DDE part: compute default probability
    if( hasDefPayoff )
    {
        FD1FDDE *fdDDE = dynamic_cast<FD1FDDE *>(genericFDModel);
        fdDDE->setDefaultProb(step, price[pEnd], -bot, top);
    }
}

/** product payoff method at steps earlier than maturity */
void CDblBarrier1fProd::PayoffBeforeMat(const double * s, int step, int bot, int top, int pStart, int pEnd,
                                        double * const * price)
{
    if ( (((UpIsOut || DownIsOut) && BarDepend != TWO_TOUCH)
          || ((UpIsOut && DownIsOut) && BarDepend == TWO_TOUCH)) )
    {
        CalcDeadInstr(step, bot, top, price);
        return;
    }

    // DDE part. need to adjust price slices before any exercise/rebate decisions for current slice
    // we should really try to shift most of the burden here into FD1FDDE so that product simply need
    // to preprocess the default payoffs for each step and pass it to FD1FDDE one time
    if( hasDefPayoff )
    {
        FD1FDDE *fdDDE = dynamic_cast<FD1FDDE *>(genericFDModel);
        double *defProbSlice = price[pEnd];
        
        double defStrikeVal = 0.0, defValue = 0.0;
        if (PayMode != CALL)
        {
            // for american, find next exercisable point to get the def payoff
            // may want to preprocess this calculation to save performance
            int strkIdx=step+1;
            while( strkIdx<model1F->TimePts.NumOfStep && !stepCanExercise[strkIdx] ) strkIdx++;
            double pv = inst->discount->pv(model1F->TimePts.StepDates[step + 1], 
                                               model1F->TimePts.StepDates[strkIdx]);
            if (PayMode==PUT)
                defStrikeVal = stepStrike[strkIdx] * pv;
            if (PayMode==FORWARD)
                defStrikeVal = -stepStrike[strkIdx] * pv;
            else // (PayMode==BINARY)
                defStrikeVal = pv;
        }

        if( UType == KI || LType == KI )
        {
            // this part is not implemented yet. we only state out the logic to be implemented

            // adjust KIKO case. we know, upon default, upper barrier is irrelevant, so only combinations are
            // LType = KI, find next KI step, then find relevant strike for that step or closest to that step,
            //             calc defStrikeVal2 (diff from above defStrikeVal which is for the strike val right 
            //             after current step, not the next KI step)
            // LType = NA, if UType==KI, no mat payoff, no need to adjust, 
            //             if UType==KO, we know mat payoff is here to stay, ie. defStrikeVal
            // LType = KO, find next KO step's rebate
            if( !Maths::isZero( defValue ) )
                fdDDE->adjustPriceByDefPO(step+1, defValue, defProbSlice, (price[0]), -bot, top);

            // adjust KO price only ((price[1]))
            // LType = KI, this price slice will not be KO anymore, ie. defStrikeVal
            // LType = NA, must be UType==KI, this slice not relevant, no need to adjust
            // LType = KO, find next KO step's rebate
            if( !Maths::isZero( defValue ) )
                fdDDE->adjustPriceByDefPO(step+1, defValue, defProbSlice, (price[1]), -bot, top);

            // adjust vanilla price ((price[2]))
            if( !Maths::isZero( defStrikeVal ) )
                fdDDE->adjustPriceByDefPO(step+1, defStrikeVal, defProbSlice, (price[2]), -bot, top);

        } else {
            // adjust (price[0]) by default rebate and put-strike/binary/fwd-strike val
            // LType = NA, must have UType==KO, we know mat payoff is here to stay, ie. defStrikeVal
            //             THIS IS NOT IMPLEMENTED SINCE WE CONCENTRATE ON SINGLE BARRIER DOWNOUT 
            // LType = KO, find next KO step's rebate
            //
            // this lookup and discounting calculation should really be preprocessed!
            int rebateIdx=step+1;
            while( rebateIdx<model1F->TimePts.NumOfStep &&
                !(stepCanKODown[rebateIdx] && 0 < LBarInput[rebateIdx]*(1.0+FP_MIN)) ) rebateIdx++;
            if( rebateIdx<=model1F->TimePts.NumOfStep )
            {
                double pv = inst->discount->pv(model1F->TimePts.StepDates[step + 1], 
                                                   model1F->TimePts.StepDates[rebateIdx]);
                defValue = pv * LBarRebate[rebateIdx];
                fdDDE->adjustPriceByDefPO(step+1,defValue, defProbSlice, (price[0]), -bot, top);
            }
        }
        
        // prep for the next time step
        if( step ) fdDDE->setDefaultProb(step, defProbSlice, -bot, top);
    }

    int j;

    double settlementPV = inst->instSettle->pvAdjust(model1F->TimePts.StepDates[step],
                                                         inst->discount.get(), 
                                                         inst->asset.get());
    double callput = settlementPV*(inst->isCall?1.0 : -1.0);

    // for KI condition below
    bool canKO = stepCanKOUp[step] || stepCanKODown[step];
    bool onlyKI = (UType == KI && LType == KI);

    if (UType == KI || LType == KI)
    {
        for ( j=-bot; j<=top; j++)
        {
            if (stepCanExercise[step]) // exercise KO
            {
                double intrinsic = GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD);
                (price[1])[j] = Maths::max(settlementPV*intrinsic,(price[1])[j]);
            }

            if (UType == KO && s[j] > UpperBar * (1.0-FP_MIN) && stepCanKOUp[step]) {
                // dead on upper barrier
                (price[0])[j] = (price[1])[j] = UBarRebate[step];
            }
            else if(LType == KO && s[j] < LowerBar * (1.0+FP_MIN) && stepCanKODown[step]) {
                // dead on lower barrier
                (price[0])[j] = (price[1])[j] = LBarRebate[step];
            }
            else if (s[j] > UpperBar * (1.0-FP_MIN) || s[j] < LowerBar * (1.0+FP_MIN)) {
                // KI already
                if ((canKO && BarDepend == KI_CANCEL_KO) || onlyKI ) {
                    (price[0])[j] = (price[2])[j]; // becomes vanilla
                }
                else {
                    (price[0])[j] = (price[1])[j]; // becomes KO (if (price[1]) is KO otherwise vanilla)
                }
            } 

        }
                
    }
    else
    {
        for (j=-bot; j<=top; j++)
        {
            if (UType == KO && s[j] > UpperBar * (1.0-FP_MIN) && stepCanKOUp[step]) // dead on upper barrier
                (price[0])[j] = UBarRebate[step];
            else if(LType == KO && s[j] < LowerBar * (1.0+FP_MIN) && stepCanKODown[step]) // dead on lower barrier
                (price[0])[j] = LBarRebate[step];
            else if (stepCanExercise[step]) // early exercise
                (price[0])[j] = Maths::max(callput*(s[j] - stepStrike[step]),(price[0])[j]);
        }

    }

}

CAssetConstSP CDblBarrier1fProd::GetAssetRef()
{
    return CAssetConstSP::dynamicCast((IObjectConstSP)inst->asset.getSP());
}

bool CDblBarrier1fProd::GetFwdStartLV() {
    return inst->fwdStarting;
}

DateTime CDblBarrier1fProd::GetFwdStartDateLV() {
    return inst->fwdStarting ? inst->startDate : inst->valueDate;
}


YieldCurveConstSP CDblBarrier1fProd::GetDiscCurveRef()  {
    return YieldCurveConstSP::dynamicCast(
        (IObjectConstSP)inst->discount.getSP());
}


string CDblBarrier1fProd::getCcyTreatment()
{
    return inst->ccyTreatment;
}

/** create a fd payoff product */
FD1F::IProduct* CDblBarrier::createProduct(FD1F* model) const
{
    
    CDblBarrier1fProd* fdProd = new CDblBarrier1fProd(this);
    fdProd->fdModel = model;
    fdProd->model1F = model;
    if ( !(UpperBarType == "KI" || LowerBarType == "KI") )
        fdProd->NumOfPriceArray = 1;

    return fdProd;
}

/** create a fd payoff product */
FD1FGeneric::IProduct* CDblBarrier::createProduct(FD1FGeneric* model) const
{
    static const string method = "CDblBarrier::createProduct";

    CDblBarrier1fProd* fdProd = new CDblBarrier1fProd(this);
    if ( !(UpperBarType == "KI" || LowerBarType == "KI") )
        fdProd->NumOfPriceArray = 1;

    if(FD1FDDE::TYPE->isInstance(model))
    {
        // don't support fwd starting
        if( fwdStarting )
            throw ModelException(method, "Do not support fwdStarting under DDE");
        if( ccyTreatment == CAsset::CCY_TREATMENT_STRUCK || ccyTreatment == CAsset::CCY_TREATMENT_PROTECTED )
            throw ModelException(method, "Do not support ccy feature under DDE");
        if( LowerBarType != "KO" || (UpperBarType == "KO" || UpperBarType == "KI") )
            throw ModelException(method, "only simple down and out barriers are currently available for dde");

        // may need to add slice for default prob, (p[3]) for default probability
        if( PayoffMode != "CALL" || LowerBarType == "KO" )
        {
            fdProd->hasDefPayoff = true;
            fdProd->NumOfPriceArray++;
        }
    }
    else
        throw ModelException("DblBarrier::createProduct", "Do not support non-DDE FD1FGeneric");

    fdProd->genericFDModel = model;
    fdProd->model1F = model;

    return fdProd;
}


// is the economic barrier equivalent to the risk barrier?
static void equivalentBarrier(const Schedule* risk, const Schedule* eco) {
    static const string method("CDblBarrier::equivalentBarrier");
    if (risk->getInterp() != eco->getInterp()) {
        throw ModelException(method, "risk & economic barrier interp "
                             "styles differ");
    }

    // what about dates? be a bit draconian and make dates match
    const DateTimeArray &riskDates = risk->getDateArray();
    const DateTimeArray &ecoDates  = eco->getDateArray();
    if (riskDates.size() != ecoDates.size()) {
        throw ModelException(method, "risk & economic barriers are "
                             "different lengths");
    }

    for (int i = 0; i < riskDates.size(); i++) {
        if (riskDates[i] != ecoDates[i]) {
            throw ModelException(method, "risk & economic barriers have "
                                 "different dates - risk (" + 
                                 riskDates[i].toString() + 
                                 ") vs economic (" + 
                                 ecoDates[i].toString() + ")");
        }
    }
}

void CDblBarrier::Validate()
{
    static const string method = "CDblBarrier::Validate";
    // just check the things that aren't/cannot be checked in 
    // validatePop2Object
    if (!asset){
        throw ModelException(method, "Asset is null");
    }
    if (!discount){
        throw ModelException(method, "Discount YC is null");
    }

    // validate against struck until thoroughly tested and agreed methodology
    if (ccyTreatment == Asset::CCY_TREATMENT_STRUCK) {
        throw ModelException(method, "ccy struck type not supported.");
    }

    AssetUtil::assetCrossValidate(asset.get(),
                                  fwdStarting,
                                  startDate,
                                  valueDate,
                                  discount,
                                  this);  
        
    if (PayoffMode == "PUT")
        isCall = false;
    else
        isCall = true;

    // can not be physical settle if Binary type
    // disable this validation for now since Tracer.cpp is using this!
    //if ( PayoffMode == "BINARY" && instSettle->isPhysical() )
    //{
    //    throw ModelException(method,
    //                         "can not physical settle in BINARY payoff mode");
    //}
    if (!(PayoffMode == "PUT" || PayoffMode == "CALL") && canExerciseEarly == true)
    {
        throw ModelException(method,
                             "can only exercise early in CALL or PUT payoff mode");
    }
    if (!(UpperBarType == "KO" || UpperBarType == "KI" || UpperBarType == "NA"))
    {
        throw ModelException(method,
                             "UpperBarType has illegal string.  It should be KO, KI or NA");
    }
    if ((UpperBarType == "KO" || UpperBarType == "KI") && (!UpperBarrier || UpperBarrier->length()==0))
    {
        throw ModelException(method,  "UpperBarrier not provided.");
    }

    if (!(LowerBarType == "KO" || LowerBarType == "KI" || LowerBarType == "NA"))
    {
        throw ModelException(method,
                             "LowerBarType has illegal string.  It should be KO, KI or NA");
    }
    if ((LowerBarType == "KO" || LowerBarType == "KI") && (!LowerBarrier || LowerBarrier->length()==0))
    {
        throw ModelException(method, "LowerBarrier  not provided.");
    }

    if (!!UpperRebate)
    {
        if (UpperRebate->length() == 0)
            UpperRebate = ScheduleSP(   );
    }

    if (!!LowerRebate)
    {
        if (LowerRebate->length() == 0)
            LowerRebate = ScheduleSP(   );
    }

    // Check the schedule doesn't contains holiday, for the case barrier interp = "N".
    if(!!UpperBarrier && UpperBarrier->getInterp() == "N" && UpperBarType != "NA"){
        DateTimeArray dates = UpperBarrier->getDates();
        for (int i = 0; i<dates .size(); i++){
            if (AssetUtil::getHoliday(asset.get())->isBusinessDay(dates[i])==false) {
                throw ModelException(method, "You can not set holiday when Interp = 'N' case.  ("
                    +Format::toString(i+1) + "-th dates in Upper Barrier is holiday).");
            }
        }
        if (!!UpperRebate && UpperRebate->getInterp() == "N")
            validateBar(UpperBarrier, UpperRebate);
    }
    if(!!LowerBarrier && LowerBarrier->getInterp() == "N" && LowerBarType != "NA"){
        DateTimeArray dates = LowerBarrier->getDates();
        for (int i = 0; i<dates.size(); i++){
            if (AssetUtil::getHoliday(asset.get())->isBusinessDay(dates[i])==false) {
                throw ModelException(method, "You can not set holiday when Interp = 'N' case.  ("
                    +Format::toString(i+1) + "-th dates in Lower Barrier is holiday).");
            }
        }
        if (!!LowerRebate && LowerRebate->getInterp() == "N")
            validateBar(LowerBarrier, LowerRebate);
    }

    // if I have economic barriers (used for reporting only, not pricing),
    // make sure they're least date consistent with the risk barriers
    if (UpperBarType != "NA" && upperEcoBarrier.get()) {
        try {
            equivalentBarrier(UpperBarrier.get(), upperEcoBarrier.get());
        }
        catch (exception& e) {
            throw ModelException(e, method, 
                                 "upper risk and economic barriers differ");
        }
    }
    if (LowerBarType != "NA" && lowerEcoBarrier.get()) {
        try {
            equivalentBarrier(LowerBarrier.get(), lowerEcoBarrier.get());
        }
        catch (exception& e) {
            throw ModelException(e, method, 
                                 "lower risk and economic barriers differ");
        }
    }
    // to do : double KI must have identical rebate schedule

}

/** returns a vol request for log-normal vol */
CVolRequestConstSP CDblBarrier1fProd::GetLNRequest()
{
    // get strike and maturity date from instrument
    DateTime        matDate = inst->exerciseSchedule->lastDate();

    double volStrike  = inst->exerciseSchedule->lastValue();

    DateTime imntStartDate = inst->fwdStarting? 
        inst->startDate: inst->valueDate;

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, imntStartDate, 
                                   matDate, inst->fwdStarting));
    return volRequest;
}

// initiate GetMarket 
void CDblBarrier::GetMarket(const IModel*         model, 
                            const CMarketDataSP    market)
{
    market->GetReferenceDate(valueDate);

    if (asset.usingCache() || !Tree1fLN::TYPE->isInstance(model))
    {// should always come in - just to cope with old regression convertion
        CAsset::getAssetMarketData(model, market.get(), ccyTreatment, 
                                   discount, asset);
    }

    discount.getData(model, market);
    instSettle->getMarket(model, market.get());
    if (premiumSettle.get())
    {
        premiumSettle->getMarket(model, market.get());
    }
}

/** returns the current value date */
DateTime CDblBarrier::getValueDate() const {
   return valueDate;
}

/** price a dead instrument until settlement - exercised, expired, knocked out etc.
returns true if it is dead (and priced), false if it is not dead */
bool CDblBarrier::priceDeadInstrument(CControl* control, CResults* results) const
{
    static string method = "CDblBarrier::priceDeadInstrument";
    try {
        double    strike        = 0.0;
        double    value         = 0.0;
        bool      foundExerDate = false;
        DateTime  exerDate;

        // reset dynamic bools to instrument level values; 
        UpperBarBreachedDynamic = UpperBarBreached;
        LowerBarBreachedDynamic = LowerBarBreached;

        bool monitorUpNow = false;
        bool monitorDownNow = false;
        if (!IntraDayMonitor)
        {
            if ( UpperBarType == "KO" || UpperBarType == "KI" )
            {
                monitorUpNow = (valueDate.getTime() == UpperBarrier->getDates()[0].getTime());
            }
            if ( LowerBarType == "KO" || LowerBarType == "KI" )
            {
                monitorDownNow = (valueDate.getTime() == LowerBarrier->getDates()[0].getTime());
            }
        }    
        double s = asset->getSpot(); 
        if ((IntraDayMonitor || monitorUpNow) && !fwdStarting)
        {
            // check breach of barriers on value date
            if (   (UpperBarType == "KO" || UpperBarType == "KI") 
                   && !UpperBarBreachedDynamic
                   && valueDate.daysDiff(UpperBarrier->getDates()[0]) >=0  // allow same date
                   && UpperBarrier->lastDate().daysDiff(valueDate)  >= 0) // barrier started and not finished.
            {
                if(UpperBarrier->getInterp() == "N")
                {
                    for(int i = 0; i < UpperBarrier->length(); i++)
                    {
                        if (UpperBarrier->getDates()[i].equals(valueDate,false))
                        {
                            if (s > UpperBarrier->interpolate(valueDate)*(1.0-FP_MIN))
                            {      
                                if (control && control->isPricing())
                                {  // set instrument level status flags only during pricing run
                                    UpperBarBreached = true; 
                                }

                                UpperBarBreachedDynamic = true;
                                UpperBarBreachDate = valueDate;
                            }
                            break;
                        }
                    }
                }
                else
                {
                    if (valueDate >= UpperBarrier->firstDate() &&
                        s > UpperBarrier->interpolate(valueDate)*(1.0-FP_MIN))
                    {                           
                        if (control && control->isPricing())
                        {  // set instrument level status flags only during pricing run
                            UpperBarBreached = true; 
                        }
                    
                        UpperBarBreachedDynamic = true;
                        UpperBarBreachDate = valueDate;
                    }
                }
            }
        }
        if ((IntraDayMonitor || monitorDownNow) && !fwdStarting)
        {
            if ( (LowerBarType == "KO" || LowerBarType == "KI") 
                 && !LowerBarBreachedDynamic
                 && valueDate.daysDiff(LowerBarrier->getDates()[0]) >= 0 //allow same date
                 && LowerBarrier->lastDate().daysDiff(valueDate) >= 0) // barrier started and not finished.
            {
                if(LowerBarrier->getInterp() == "N")
                {
                    for(int i = 0; i < LowerBarrier->length(); i++)
                    {
                        if (LowerBarrier->getDates()[i].equals(valueDate,false))
                        {
                            if (s < LowerBarrier->interpolate(valueDate)*(1.0+FP_MIN))
                            {                    
                                if (control && control->isPricing())
                                {  // set instrument level status flags only during pricing run
                                    LowerBarBreached = true; 
                                }
                            
                                LowerBarBreachedDynamic = true;
                                LowerBarBreachDate = valueDate;
                            }
                            break;
                        }
                    }
                }
                else
                {
                    if (valueDate >= LowerBarrier->firstDate() && 
                        s < LowerBarrier->interpolate(valueDate)*(1.0+FP_MIN))
                    {
                        if (control && control->isPricing())
                        {  // set instrument level status flags only during pricing run
                            LowerBarBreached = true; 
                        }
                    
                        LowerBarBreachedDynamic = true;
                        LowerBarBreachDate = valueDate;
                    }
                }
            }
        }  

        if (isExercised && ((LowerBarBreachedDynamic && !LowerBarBreached && LowerBarType == "KO") ||
                            (UpperBarBreachedDynamic && !UpperBarBreached && UpperBarType == "KO"))) 
        {
            throw ModelException(method, 
                                 "Can not breach KO barrier and exercise.");
        }

        DateTime matDate = exerciseSchedule->lastDate();

        bool lowerOut = LowerBarBreachedDynamic  && LowerBarType == "KO";
        bool upperOut = UpperBarBreachedDynamic  && UpperBarType == "KO";

        bool isOut = ((lowerOut || upperOut) && BarrierDependence != "TWO_TOUCH")
            || (lowerOut && upperOut);

        bool isDead = valueDate >= matDate || (isExercised && canExerciseEarly) || isOut;
        if (!isDead)
            return false; // not dead yet

        DateTime settlementDate = instSettle->settles(matDate, asset.get());

        if (valueDate >= settlementDate)
        {// settled already
            results->storePrice(0.0, discount->getCcy());
            addOutputRequests(control, results, 0.0, 0.0);
            return true;
        }
    
        // sort out ko case first
        if (isOut) // two touch case has problem choosing rebate - use lower for now
        {
            if (!RebateAtMat)
                settlementDate = instSettle->settles(
                    lowerOut ? LowerBarBreachDate : UpperBarBreachDate,
                    asset.get());

            if (settlementDate >= valueDate) // should really be > but then can't record cash flow if no delay
            {
                if (lowerOut && !!LowerRebate)
                    value = LowerRebate->interpolate(LowerBarBreachDate);
                else if (upperOut && !!UpperRebate)
                    value = UpperRebate->interpolate(UpperBarBreachDate);
            }
        }
        // we may simply consider isExercised flag
        // but here we give the instrinsic value at maturity regardless
        else if (!isExercised && valueDate >= matDate)
        {// maturity instrinsic value
            if((UpperBarType == "KI" && LowerBarType == "KI") && (!LowerBarBreachedDynamic &&
                                                                  !UpperBarBreachedDynamic))
            {
                if (!!LowerRebate)
                    value = LowerRebate->interpolate(matDate);
                else if (!!UpperRebate)
                    value = UpperRebate->interpolate(matDate);
            }
            else if(UpperBarType == "KI" && LowerBarType == "KO")
            {
                if(LowerBarBreachedDynamic)
                    value = LowerRebate->interpolate(LowerBarBreachDate);
                else if(!UpperBarBreachedDynamic)
                    value = UpperRebate->interpolate(matDate);
            }
            else if(UpperBarType == "KO" && LowerBarType == "KI")
            {
                if(UpperBarBreachedDynamic)
                    value = UpperRebate->interpolate(UpperBarBreachDate);
                else if(!LowerBarBreachedDynamic)
                    value = LowerRebate->interpolate(matDate);
            }
            else if(UpperBarType == "KI" && !UpperBarBreachedDynamic){
                value = UpperRebate->interpolate(matDate);
            }
            else if(LowerBarType == "KI" && !LowerBarBreachedDynamic){
                value = LowerRebate->interpolate(matDate);
            }
            else
            { 
                strike  = exerciseSchedule->lastValue();
                exerDate = exerciseSchedule->lastDate();
                foundExerDate = true;
            }
        }
        else
        {// early exercised
            // check that exercise date is not in the Future
            if ( dateExercised.isGreater(valueDate))
            {
                throw ModelException(method,
                                     "Option exercised on " + 
                                     dateExercised.toString() + ". " +
                                     "This date is after the current value date (" + 
                                     valueDate.toString());
            }

            // check whether exercise date is valid and find corresponding strike
            if ( !canExerciseEarly ) {
                // multi-european case
                DateTimeArray exerciseDates = exerciseSchedule->getDates();

                for (int i=0 ; i<exerciseDates.size() ; ++i)
                {
                    if ( dateExercised.equals(exerciseDates[i])) {
                        foundExerDate = true;
                        exerDate      = dateExercised;
                        strike        = exerciseSchedule->interpolate(
                            dateExercised);
                        break;
                    }
                    else 
                    {
                        throw ModelException(method, 
                                             "Option has been exercised on " + 
                                             dateExercised.toString() + ". " +
                                             "This date could not be found in the exercise schedule.");
                    }
                }
            } else {
                if (  dateExercised.isGreaterOrEqual(
                    exerciseSchedule->firstDate()) &&
                      !dateExercised.isGreater(
                          exerciseSchedule->lastDate())  ) {
                    // american case - note: not checking for weekend yet
                    foundExerDate = true;
                    exerDate      = dateExercised;
                    strike        = exerciseSchedule->interpolate(
                        dateExercised);
                }
                else 
                {
                    throw ModelException(method, 
                                         "Option has been exercised on " + 
                                         dateExercised.toString() + " but valid " +
                                         "exercise range is from " + 
                                         exerciseSchedule->firstDate().toString() +
                                         " to " +
                                         exerciseSchedule->lastDate().toString());
                }
            }
        }

        if (foundExerDate)
        {
            if( PayoffMode == "BINARY" )
                value = 1.0;
            else
                value = GetIntrinsic(spotAtMaturity,
                                     strike,
                                     isCall, 
                                     PayoffMode != "FORWARD");
  
            settlementDate = instSettle->settles(exerDate, asset.get());
        }

        double scalingFactor = InstrumentUtil::scalePremium(
            oneContract,
            false,
            notional,
            0.0,
            initialSpot);

        value *= scalingFactor;

        // write known cash flows before pv'ing
        if (!Maths::isZero(value))
        {
            if (valueDate >= matDate || 
                (isOut && 
                 (valueDate >= UpperBarBreachDate 
                  || valueDate >= LowerBarBreachDate )))
            {
                OutputRequest* request = 0;
                // record known cash flow for processing
                bool isCash = isOut || !instSettle->isPhysical();
                if ( isCash && (request = control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS))){
                    CashFlow flow(settlementDate, value);
                    CashFlowArray flowArr(1, flow);
                    OutputRequestUtil::recordKnownCashflows(control, 
                                                            results, 
                                                            discount->getCcy(), 
                                                            &flowArr);
                }
                else if( !isCash && (request = control->requestsOutput(OutputRequest::PHYSICAL_DELIVERY)))
                {
                    double shares = isCall?scalingFactor:(-scalingFactor);
                    PhysicalDelivery::recordPhysicalDelivery(shares,
                                                             strike,
                                                             matDate,
                                                             asset.get(),
                                                             control, 
                                                             results);
                }
                if( !isCash ) value = 0.0; // drop the physical delivery value from price
            }
        }

        // pv from settlement to today
        value *= discount->pv(valueDate, settlementDate);
        // store results
        results->storePrice(value, discount->getCcy());
        addOutputRequests(control, results, value, 0.0);

        return true;
    }
    catch (exception &e) {
        throw ModelException(e, method);
    }   
}

/** Rolls the value date and sets initial spot if rolling over start date */
bool CDblBarrier::sensShift(Theta* shift)
{    
    DateTime newDate = shift->rollDate(valueDate);

    DateTime matDate = exerciseSchedule->lastDate();

    if ( ( newDate >= matDate && valueDate < matDate ) ||
         ( valueDate == matDate && Maths::isZero(spotAtMaturity)))
        spotAtMaturity = asset->getThetaSpotOnDate(shift, matDate);

    if (fwdStarting && newDate.isGreaterOrEqual(startDate) &&
        startDate.isGreaterOrEqual(valueDate))
    {
        fwdStarting = false;
        initialSpot = asset->getThetaSpotOnDate(shift, startDate);
        exerciseSchedule->scale(initialSpot);
        if (!!UpperBarrier)
            UpperBarrier->scale(initialSpot);
        if (!!LowerBarrier)
            LowerBarrier->scale(initialSpot);
        if (!RebateNotScaled)
        {
            if (!!UpperRebate)
                UpperRebate->scale(initialSpot);
            if (!!LowerRebate)
                LowerRebate->scale(initialSpot);
        }
    }

    // roll today 
    valueDate = newDate;
    
    return true;
}

/** when to stop tweaking */
DateTime CDblBarrier::endDate(const Sensitivity* sensControl) const {
    DateTime matDate = exerciseSchedule->lastDate();
    DateTime instEnd  = instSettle->settles(matDate, asset.get());
    DateTime assetEnd = asset->settleDate(matDate);
    DateTime end = instEnd.isGreater(assetEnd) ? instEnd : assetEnd;
    return end;
}

/** for ITaxableInst::Basic */
const DateTime CDblBarrier::getFinalPaymentDate() const {
    // can't support anything which may have early cash flows
    if (canExerciseEarly || isExercised) {
        throw ModelException("CDblBarrier::getFinalPaymentDate",
                             "Tax is not supported for early exercisable barriers");
    }
    DateTime matDate = exerciseSchedule->lastDate();
    return instSettle->settles(matDate, asset.get());
}

bool CDblBarrier::sensShift(LegalTerms* shift) {
    // 1. upperEcoBarrier -> UpperBarrier
    // 2. lowerEcoBarrier -> LowerBarrier
    // 3. Anything else?
    if (UpperBarType != "NA" && upperEcoBarrier.get()) {
        UpperBarrier = upperEcoBarrier;
    }
    if (LowerBarType != "NA" && lowerEcoBarrier.get()) {
        LowerBarrier = lowerEcoBarrier;
    }
    return true; // continue shifting
}

void CDblBarrier::addOutputRequests(Control* control,
                                 Results* results,
                                 const double& fairValue,
                                 const double& indVol) const
{
    // take care of additional outputs
    if (control && control->isPricing()) {
        DateTime matDate = exerciseSchedule->lastDate();
        // DELAY_PRICE
        InstrumentUtil::delayPriceHelper(control,
                                         results,
                                         fairValue,
                                         valueDate,
                                         discount.get(),
                                         asset.get(),
                                         premiumSettle.get());
        // IND_VOL
        InstrumentUtil::recordIndicativeVol(control,results,indVol);
        // FWD_AT_MAT
        try {
            InstrumentUtil::recordFwdAtMat(control,
                                           results,
                                           matDate,
                                           valueDate,
                                           asset.get());
        }
        catch(exception&) {
            // continue if fwd failed - this hapens now for quanto asset with CEVj vol
        }

        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request && !fwdStarting) {
            // report barrier levels over a date range
            DateTime upperDate = BarrierLevel::barrierWindow(valueDate);

            BarrierLevelArraySP levels(new BarrierLevelArray(0));

            if (LowerBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = lowerEcoBarrier.get() ? lowerEcoBarrier.get(): 
                    LowerBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(false,(*subset)[i].date,(*subset)[i].amount,IntraDayMonitor);
                    levels->push_back(bl);
                }
            }

            if (UpperBarType != "NA") {
                // use economic barrier (if it exists)
                Schedule* s = upperEcoBarrier.get() ? upperEcoBarrier.get(): 
                    UpperBarrier.get();
                CashFlowArraySP subset(s->subset(valueDate, upperDate));
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(true,(*subset)[i].date,(*subset)[i].amount,IntraDayMonitor);
                    levels->push_back(bl);
                }
            }

            if (!levels->empty()) {
                OutputRequestUtil::recordBarrierLevels(control,
                                                       results,
                                                       asset->getTrueName(),
                                                       levels.get());
            }
        }            
    }
}

void CDblBarrier1fProd::recordOutputRequests(Control* control, Results* results, double fairValue)
{
    // take care of additional outputs
    if (control && control->isPricing()) {
        DateTime       matDate = inst->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        try {
            if ( matDate.isGreater(inst->valueDate) )
            {

                DateTime imntStartDate = inst->fwdStarting? 
                    inst->startDate: inst->valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest = GetLNRequest();

                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("CDblBarrier1fProd::recordOutputRequests", 
                                         "No Black Scholes Vol");
                }

                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            else
            {
                indVol = 0.0;
            }
        }
        catch(exception&)
        {// continue if indicative vol fail - non BS vol
            indVol = 0.0;
        }

        inst->addOutputRequests(control,
                                    results,
                                    fairValue,
                                    indVol);

        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_DISCONTINUITY_RISK);
        
        // only available if the instrument has started
        // and only one barrier
        if (request) {
            if (                
                // not foward starting
                (!inst->fwdStarting) &&
                // at lease one NA barrier
                (inst->UpperBarType=="NA" ||  inst->LowerBarType=="NA") &&
                // no more than one NA barrier
                (inst->UpperBarType!="NA" ||  inst->LowerBarType!="NA")) {
                
                // select the barrier, rebate and spot tweak
                ScheduleSP barrier;
                ScheduleSP rebate;
                double spotTweak;
                DateTime today = inst->valueDate;
                
                if (inst->UpperBarType!="NA") {
                    // upper barrier
                    barrier = ScheduleSP(inst->UpperBarrier);
                    rebate = ScheduleSP(inst->UpperRebate);
                    spotTweak = 0.99;
                }
                else {
                    //down barrier
                    barrier = ScheduleSP(inst->LowerBarrier);
                    rebate = ScheduleSP(inst->LowerRebate);
                    spotTweak = 1.01;
                }
                
                // the barrier and rebate may be bend 
                // we use today barrier level
                // already started so absolute level
               
                
                // returns the difference of prices
                try {
                    double barrierLevel = barrier->interpolate(today);
                    // calculate the difference of price when the spot is 1% out and 1%in the barrier
                    double tweakPrice = tweakSpotAndPrice(spotTweak*barrierLevel, fdModel, inst)
                                      - tweakSpotAndPrice((2.0-spotTweak)*barrierLevel, fdModel, inst);
                    results->storeRequestResult(request,tweakPrice);
                }
                catch (exception&){
                    // nothing to do
                }
            }
        }    
    }
}


/** returns all strikes on the vol surface to which 
    this instrument is sensitive */
DoubleArraySP CDblBarrier::getSensitiveStrikes(OutputNameConstSP outputName,
                                               const IModel*      model)
{
    DoubleArraySP sensStrikes = DoubleArraySP(new DoubleArray(0));

    if (avoidVegaMatrix(model)) {
        throw ModelException("CDblBarrier::getSensitiveStrikes", 
                             "VEGA_MATRIX is not valid for this instrument");
    }

    // get start date for vol interpolation
    DateTime imntStartDate = fwdStarting ? startDate:valueDate;
    // get last exercise date in exercise schedule
    DateTime maturityDate = exerciseSchedule->lastDate();
    // get last strike in exercise schedule
    double   strike       = exerciseSchedule->lastValue();

    // create a vol request object to be passed on
    LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(strike, 
                                                                   imntStartDate, 
                                                                   maturityDate,
                                                                   fwdStarting));

    SensitiveStrikeDescriptor sensStrikeDesc;
    sensStrikeDesc.forwardOnly = false;

    asset->getSensitiveStrikes(volRequest.get(), outputName, 
                               sensStrikeDesc, sensStrikes);

    return sensStrikes;
}

///////////////////////////////
//    Smoothing Price
// Need to be more Tested
void CDblBarrier1fProd::SmoothPrice(const double* s, int step, int bot, int top, int pStart, int pEnd, const vector< double * > & p)
{
    int index_Upper = -bot, index_Lower = -bot;
    double upSideValue, downSideValue;
    double boundaryWidth = inst->DEBUG_SmoothBarrierWidth;

    if (LType != NA && LowerBar > s[-bot])    
    {   //// search barrier level  ////
        for (int iPrice = pStart; iPrice <= pEnd; iPrice++)
        {
            for (int j=-bot; j<top; j++)
            {
                if (s[j] < LowerBar * (1.0+FP_MIN) && s[j+1] > LowerBar * (1.0-FP_MIN)) 
                {
                    index_Lower = j;
                    break;
                }
            }
            //// Shit Price on Smoothed function ////
            if (index_Lower > -bot && index_Lower < top)
            {
                upSideValue   = (p[iPrice])[index_Lower+1];
                downSideValue = (p[iPrice])[index_Lower];

                (p[iPrice])[index_Lower]   = SmoothValue(upSideValue, LowerBar * (1.0+boundaryWidth),
                                                           downSideValue, LowerBar * (1.0 - boundaryWidth),
                                                           s[index_Lower]);
                (p[iPrice])[index_Lower+1] = SmoothValue(upSideValue, LowerBar * (1.0+boundaryWidth),
                                                           downSideValue, LowerBar * (1.0-boundaryWidth),
                                                           s[index_Lower+1]);
            }
        }
    }

    if (UType != NA && UpperBar < s[top])
    {   //// search barrier level  ////
        for (int iPrice = pStart; iPrice <= pEnd; iPrice++)
        {
            for (int j=index_Lower; j<top; j++)
            {
                if (UType != NA && s[j] < UpperBar * (1.0+FP_MIN) && s[j+1] > UpperBar * (1.0-FP_MIN)) 
                {
                    index_Upper = j;
                    break;
                }
            }
            //// Shit Price on Smoothed function ////
            if (index_Upper > -bot && index_Upper < top)
            {
                upSideValue   = (p[iPrice])[index_Upper+1];
                downSideValue = (p[iPrice])[index_Upper];
                (p[iPrice])[index_Upper]   = SmoothValue(upSideValue, UpperBar * (1.0+boundaryWidth),
                                                           downSideValue, UpperBar * (1.0-boundaryWidth),
                                                           s[index_Upper]);
                (p[iPrice])[index_Upper+1] = SmoothValue(upSideValue, UpperBar * (1.0+boundaryWidth),
                                                           downSideValue, UpperBar * (1.0-boundaryWidth),
                                                           s[index_Upper+1]);
            }
        }
    }

}

/** Implementation of DDEInitiator interface, built on FD1FGeneric
 */
DateTime CDblBarrier1fProd::maxMaturity() const { return inst->exerciseSchedule->lastDate(); }

void CDblBarrier1fProd::sensitiveDates(  DateTimeArray  &dates) const
{
    dates = inst->exerciseSchedule->getDates();
}

void CDblBarrier1fProd::sensitiveStrikes(  const DateTimeArray     dates,
                                         DoubleArray             &strikes,   // same dimension as dates
                                         bool                    &strikeIsPct) const
{
    if( dates.size() != strikes.size() )
        throw ModelException("CDblBarrier1fProd::sensitiveStrikes", "Dates and strikes dimension mismatch");
    
    /************ doesn't support fwd start ***************/
    if ( inst->fwdStarting ) 
        throw ModelException("CDblBarrier1fProd::sensitiveStrikes", "Forward starting option not supported under DDE");
    
    double strike = inst->exerciseSchedule->lastValue();
    for(int i=dates.size()-1; i>=0; i--)
        strikes[i] = strike;
    
    strikeIsPct = false;
}


//*************** new FD products *******************
DblBarrierFDProd::DblBarrierFDProd(const CDblBarrier* inst, FDModel* m) :
    LatticeProdEDRIns(m,0,0),
    inst(inst)
{
    static const string method = "CDblBarrierFDProd::CDblBarrierFDProd";

    // first: set discount curve
    if( tree1f )
        tree1f->setDiscountCurve( inst->discount.getSP() );

    // second: create spot payoff
    payoffIndex = model->createProduct( IProdCreatorSP( new
        IndexSpecEQ( inst->asset.getName(), inst->asset, inst->ccyTreatment ) ) );

    // barrier parts
    if (inst->UpperBarType == "KO")
        UType = KO;
    else if (inst->UpperBarType == "KI")
        UType = KI;
    else
        UType = NA;

    if (inst->LowerBarType == "KO")
        LType = KO;
    else if (inst->LowerBarType == "KI")
        LType = KI;
    else
        LType = NA;

    if (inst->BarrierDependence == "KI_KEEP_KO")
        BarDepend = KI_KEEP_KO;
    else if (inst->BarrierDependence == "KI_CANCEL_KO")
        BarDepend = KI_CANCEL_KO;
    else if (inst->BarrierDependence == "ONCE_TOUCH")
        BarDepend = ONCE_TOUCH;
    else if (inst->BarrierDependence == "TWO_TOUCH")
        BarDepend = TWO_TOUCH;
    else
        throw ModelException(method, inst->BarrierDependence + " - unkown Barrier dependence");

    if (inst->MonitoringDependence == "BOTH")
        MonDepend = BOTH;
    else if (inst->MonitoringDependence == "UPPER")
        MonDepend = UPPER;
    else if (inst->MonitoringDependence == "LOWER")
        MonDepend = LOWER;
    else
        throw ModelException(method, inst->MonitoringDependence + " - unkown Barrier dependence");

    if (inst->PayoffMode == "CALL")
        PayMode = CALL;
    else if (inst->PayoffMode == "PUT")
        PayMode = PUT;
    else if (inst->PayoffMode == "BINARY")
        PayMode = BINARY;
    else if (inst->PayoffMode == "FORWARD")
        PayMode = FORWARD;
    else
        throw ModelException(method, inst->PayoffMode + " - unkown payoff mode");

    // set number of prices
    if( UType == KI || LType == KI )
    {
        if( ( UType == KO || LType == KO ) && BarDepend == KI_CANCEL_KO )
            numPrices = 3;
        else
            numPrices = 2;
    }
    else
        numPrices = 1;

    // set number of inserted nodes
    if( tree1f )
    {
        // default to NODE_INSERTION smoothing
        if( tree1f->GetSmoothMethod() == CTree1f::DEFAULT )
            tree1f->SetSmoothMethod( CTree1f::NODE_INSERTION );

        if( tree1f->GetSmoothMethod() != CTree1f::NODE_INSERTION )
            numIns = 0;
        else
            numIns = 2;
    }
    else if( fd1dRet )
        numIns = 2 * numPrices;
}

/** returns a vol request for log-normal vol */
CVolRequestConstSP DblBarrierFDProd::GetLNRequest() const
{
    // get strike and maturity date from instrument
    DateTime matDate = inst->exerciseSchedule->lastDate();
    double volStrike = inst->exerciseSchedule->lastValue();

    CVolRequestConstSP volRequest(
        new LinearStrikeVolRequest(volStrike, getStartDate(),
                                   matDate, inst->fwdStarting));
    return volRequest;
}

/** initialise tree1f - allow product customisation */
void DblBarrierFDProd::init(CControl* control) const
{
    static const string method = "DblBarrierFDProd::init()";
    try
    {
        int i;

        if( ! fd1dRet )
        {
            if( tree1f )
            {
                // default to DollarDivTreament = false !!! to be removed after EAS/EIS can handle it correctly !!!
                tree1f->SetDivAmountTreatment(false);

                tree1f->NumOfPrice = numPrices;
                tree1f->NumOfInsertNode = numIns;

                // this needs change if fwd start tree has to start today !!!
                if (inst->fwdStarting)
                    tree1f->controlSameGridFwdStart(inst->ccyTreatment);                
            }

            // ************ below is similar to Vanilla  ***************
            /** customize tree parameters here and set up the tree */
            DateTimeArray segDates;
            segDates.resize(inst->DEBUG_num_segments+1); // to try using more segments

            segDates[0] = getStartDate();

            DateTime mat = inst->exerciseSchedule->lastDate();
            IntArray density (inst->DEBUG_num_segments);
            // ********** testing simple density ratios, 3:2:1 ...need to change this after testing !!!
            double t = segDates[0].yearFrac(mat);
            for (i=0; i<inst->DEBUG_num_segments; i++)
            {
                //density[i] = 1;
                density[i] = inst->DEBUG_num_segments-i;
                segDates[i+1] = segDates[i].rollDate((int)(365*t/inst->DEBUG_num_segments));
            }
            // make sure the last date is spot on
            segDates[inst->DEBUG_num_segments] = mat;
                                                     
            // all exercise dates are copied to critical dates
            DateTimeArray critDates = inst->exerciseSchedule->getDates();

            // add div event dates if needed
            EventAssetMove divEvent;
            DateTimeArraySP divCritDates;
            if (inst->canExerciseEarly)
            {// American exercise or dollar div interp treatment
                // create div events, this call is very expensive for basket with lots of div dates
                // should only need once for pricing call but need to think how to store/copy for tweaks
                int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates

                if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                             segDates[0], 
                                             segDates[1],
                                             numDivs,
                                            divEvent)){
                    // calculate critical dates
                    divCritDates = divEvent.getCritDate(numDivs, inst->isCall);
                }
            }

            // add barrier dates to critical dates)
            if (UType == KO || UType == KI)
                CDblBarrier::addCritBarDates(inst->UpperBarrier, inst->getValueDate(), mat, tree1f, critDates);
            if (LType == KO || LType == KI)
                CDblBarrier::addCritBarDates(inst->LowerBarrier, inst->getValueDate(), mat, tree1f, critDates);

            // add critical dates
            model->addCritDates( critDates );
            if( divCritDates.get() )
                model->addCritDates( *divCritDates );

            // prepare timeline set up
            model->initSegments( segDates, density );
        }
        else
        {
            //need to add for forward starting, 
            
            //set up segment info.
            //set bar dates as seg dates
            IntArray isAddedSeg ;
            DateTimeArray segDates;
            DateTime t0;
            TimeMetricConstSP metric = model->getTimeMetric();
            const DateTime& matDate= inst->exerciseSchedule->lastDate();

            if (inst->fwdStarting && inst->startDate>inst->valueDate)
            {
                t0 = inst->startDate;  // to do: this needs to be checked !!!
            }
            else{
                t0 = inst->valueDate;
            }

            vector<ScheduleSP> bar;
            vector<string> barType;
            bar.resize(2);
            barType.resize(2);

            bar[0] = inst->LowerBarrier;
            bar[1] = inst->UpperBarrier;
            barType[0] = inst->LowerBarType;
            barType[1] = inst->UpperBarType;

            FDUtils::SetSegDatesGen(t0,
                    matDate,
                    metric,
                    bar, 
                    barType, 
                    segDates, &isAddedSeg);

            //end with seg

        //-----------------------------------------------------------------------------

            // all exercise dates are copied to critical dates    
            DateTimeArray critDates = inst->exerciseSchedule->getDates();
            // remove exercise date from crit date
            critDates.erase(critDates.end()-1);
            // add div event dates if needed
            EventAssetMove divEvent;
            DateTimeArraySP divCritDates;

        /*
            if (inst->canExerciseEarly || tree1f->DEBUG_DollarDivMethod > 0)
            {// American exercise or dollar div interp treatment
                // create div events, this call is very expensive for basket with lots of div dates
                // should only need once for pricing call but need to think how to store/copy for tweaks
                int numDivs = (int)(4*segDates[0].yearFrac(segDates[1]))+1; // 4 divs per year selected as critical dates

                if (AssetUtil::getJumpEvents(inst->asset.get(), 
                                            segDates[0], 
                                            segDates[1],
                                            numDivs,
                                            divEvent)){
                    // calculate critical dates
                    divCritDates = divEvent.getCritDate(numDivs, inst->isCall);
                }
            }
        */

            //using var grid
            DoubleArray outCritSpacePts;

            outCritSpacePts.resize(1);

            //get roughly the barriers level (average) from products 

            //general case: allow more than 2 barriers
            FDUtils::setCritSpacePtsAllGen(t0,
                                        inst->exerciseSchedule->lastDate(),
                                        bar,
                                        barType,
                                        outCritSpacePts);


            fd1dRet->setNumOfInsertNode( numIns );
            fd1dRet->setNbOfProd(1);
            fd1dRet->setMaxNumOfValue( numPrices );

            // add critical dates
            model->addCritDates( critDates );
            if( divCritDates.get() )
                model->addCritDates( *divCritDates );

            IntArray density (1, 1);

            // prepare model set up   
            //dblBarrier need extra data, call fd1dRet->addInitData to fill iniDataExtra
            //otherwise, just call model->addInitData to fill initData.
            fd1dRet->initSegments(segDates, density, outCritSpacePts, &isAddedSeg);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** initialising and setting product variables
    this is called per pricing call before tree sweep call (after InitTree)  */
void DblBarrierFDProd::initProd()
{
    static const string method = "DblBarrierFDProd::initProd()";
    try {
        int i;
        int lastStep = model->getLastStep();
        double fwdAtStart = 1.0, scaleFactor = 1.0;

        // to customize for fd

        stepStrike.resize(lastStep+1);
        stepCanExercise.resize(lastStep+1);
        // ask tree to decide first about steps that can exercise
        AssetUtil::setStepExercise(stepCanExercise,
                                 model->getDates(),
                                 inst->exerciseSchedule,
                                 inst->canExerciseEarly,
                                 inst->asset.getSP());

        bool canInterpExSched = (inst->exerciseSchedule->length() > 1);
        // get spot at start if needed
        if (inst->fwdStarting && model->getDate(0)>inst->valueDate)
        {
            fwdAtStart = inst->asset->fwdValue(inst->startDate);
            if (!inst->RebateNotScaled)
                scaleFactor = fwdAtStart; 
        }

        // scale adjust of gammaThreshold
        if (tree1f){
            if (tree1f->gammaNodesInterval>1){
                // on tree, price is always treated as one contract in DblBarrier                
                double refStrike = inst->fwdStarting ? fwdAtStart : inst->initialSpot;
                tree1f->SetGammaThresholdScaled(refStrike, refStrike);
            }
        }

        // get strikes
        stepStrike[lastStep] = fwdAtStart*inst->exerciseSchedule->lastValue();

        for (i=0; i<lastStep; i++)
        {
            if (stepCanExercise[i])
            {
                if (canInterpExSched)
                {
                    stepStrike[i] = fwdAtStart*inst->exerciseSchedule->interpolate(model->getDate(i));
                }
                else
                {
                    stepStrike[i] = stepStrike[lastStep]; // strike level is constant for american with no schedule
                }
            }
        }

        UBarInput.resize(lastStep+1);
        LBarInput.resize(lastStep+1);
        UBarRebate.resize(lastStep+1);
        LBarRebate.resize(lastStep+1);

        // init barriers to far away levels
        double spot = inst->asset->getSpot();
        for(i=0; i<=lastStep; i++)
        {
            UBarInput[i] = spot/FP_MIN;
            LBarInput[i] = - 2.0*FP_MIN;
        }
        // init rebate to 0
        for(i=0; i<=lastStep; i++)
            UBarRebate[i] = LBarRebate[i] = 0.0;

        const DateTime& valDate= inst->valueDate;
//    const DateTime& matDate= inst->exerciseSchedule->lastDate();
        string uInterp;
        string lInterp;
//    const DateTime& matDate= inst->exerciseSchedule->lastDate();
        if (UType != NA)
            uInterp =  inst->UpperBarrier->getInterp();
        if (LType != NA)
            lInterp =  inst->LowerBarrier->getInterp();
    
        // check KI barrier is breached or not.
        if(inst->UpperBarBreachedDynamic && UType == KI)
            UType = NA;
        if(inst->LowerBarBreachedDynamic  && LType == KI)
            LType = NA;
        // for continuous KI check if it's in already
        if (inst->IntraDayMonitor)
        {
            if (UType != NA)
                if (valDate>= inst->UpperBarrier->firstDate() &&  valDate<= inst->UpperBarrier->lastDate() &&
                    (uInterp!= Schedule::INTERP_NONE && spot>=inst->UpperBarrier->interpolate(valDate) ||
                     uInterp== Schedule::INTERP_NONE && inst->UpperBarrier->coversDateRange(valDate, valDate, true)))
                {
                    if (UType == KI)
                    {
                        UType = NA; // make it not applicable
                        if ((BarDepend == KI_CANCEL_KO && LType == KO) || (BarDepend == ONCE_TOUCH && LType == KI))
                            LType = NA;
                    }
                }
            if (LType != NA)
                if (valDate>= inst->LowerBarrier->firstDate() &&  valDate<= inst->LowerBarrier->lastDate() &&
                    (lInterp!= Schedule::INTERP_NONE && spot<=inst->LowerBarrier->interpolate(valDate) ||
                     lInterp== Schedule::INTERP_NONE && inst->LowerBarrier->coversDateRange(valDate, valDate, true)))
                {
                    if (LType == KI)
                    {
                        LType = NA; // make it not applicable
                        if ((BarDepend == KI_CANCEL_KO && UType == KO) || (BarDepend == ONCE_TOUCH && UType == KI))
                            UType = NA;
                    }
                }
        }

        // set up barrier and rebate levels at each step
        if (UType != NA)
        {
            Barrier::setStepLevel(inst->UpperBarrier, fwdAtStart, model->getDates(), UBarInput);
            if (!!inst->UpperRebate)
            {
                Barrier::setStepLevel(inst->UpperRebate, scaleFactor, model->getDates(), UBarRebate);
                if (inst->RebateAtMat)
                {
                    for (i=0; i<= lastStep; i++)
                    {
                        if(!Maths::isZero(UBarRebate[i]))
                            UBarRebate[i] *= inst->discount->pv(model->getDate(i),
                                             inst->instSettle->settles(inst->exerciseSchedule->lastDate(), inst->asset.get()));
                    }
                }
                /*  Need to be considered about Rebate settlment PV adjsutment
                else{
                    for (i=0; i<= lastStep; i++)
                    {
                        if(!Maths::isZero(UBarRebate[i])){
                            if (CashSettleDate::TYPE->isInstance(inst->instSettle.get())) {
                                // It's really wired case, that rebate pay at hit, but settlement is fixed.
                                // It may better to disallow this case, but it's too late.  
                                // Many live position could have this case.
                                // Those case are just calculated as pay immidiatelly (T+0 days)!
                            }
                            else                            
                                UBarRebate[i] *= inst->discount->pv(model->getDate(i),
                                                                    inst->instSettle->settles(model->getDate(i), inst->asset.get()));
                        }
                    }
                }*/
            }
        }
        if (LType != NA)
        {
            Barrier::setStepLevel(inst->LowerBarrier, fwdAtStart, model->getDates(), LBarInput);
            if (!!inst->LowerRebate)
            {
                Barrier::setStepLevel(inst->LowerRebate, scaleFactor, model->getDates(), LBarRebate);
                if (inst->RebateAtMat)
                {
                    for (i=0; i<= model->getLastStep(); i++)
                    {
                        if(!Maths::isZero(LBarRebate[i]))
                            LBarRebate[i] *= inst->discount->pv(model->getDate(i),
                                                                    inst->instSettle->settles(inst->exerciseSchedule->lastDate(), inst->asset.get()));
                    }
                }
                /*  Need to be considered about Rebate settlment PV adjsutment
                else{
                    for (i=0; i<= model->getLastStep(); i++)
                    {
                        if(!Maths::isZero(LBarRebate[i])){
                            if (CashSettleDate::TYPE->isInstance(inst->instSettle.get())) {
                                // It's really wired case, that rebate pay at hit, but settlement is fixed.
                                // It may better to disallow this case, but it's too late.  
                                // Many live position could have this case.
                                // Those case are just calculated as pay immidiatelly (T+0 days)!
                            }
                            else
                                LBarRebate[i] *= inst->discount->pv(model->getDate(i),
                                                                    inst->instSettle->settles(model->getDate(i), inst->asset.get()));
                        }
                    }
                }*/
            }
        }

        // disallow non-zero Rabate for KI barrier.
        if (UType == KI){
            for (i=0; i<= lastStep; i++)
            {
                if(!Maths::isZero(UBarRebate[i]))
                    throw ModelException(method, "Upper Barrier is KI. non-zero rebate for KI is disallowed.");
            }                    

        }
        if (LType == KI){
            for (i=0; i<= lastStep; i++)
            {
                if(!Maths::isZero(LBarRebate[i]))
                    throw ModelException(method, "Lower Barrier is KI. non-zero rebate for KI is disallowed.");
            }                    

        }


        // now we need to set up an indication of whether or not we can KO at a step
        // UBarInput/LBarInput are set up with the "right" barrier levels so we can
        // interrogate them at each step (i.e. KO on dates only has zero/infinite
        // barriers on dates other than barrier dates) - HOWEVER we want to forbid
        // KO on holidays/weekends (c.f. exercise), so set up a vector of flags
        // indicating if we can check the KO condition
        stepCanKOUp.resize(lastStep+1);
        stepCanKODown.resize(lastStep+1);

        if (UType == NA) {
            for(i = 0; i <= lastStep; i++) {        
                stepCanKOUp[i] = false;
            }
        }
        else {     
            setKOCondition(stepCanKOUp, model->getDates(), inst->asset.getSP());
        }
     
        if (LType == NA) {
            for(i = 0; i <= lastStep; i++) {        
                stepCanKODown[i] = false;
            }
        }
        else {       
            setKOCondition(stepCanKODown, model->getDates(), inst->asset.getSP());
        }

        //extra init prod data
        if( fd1dRet )
        {
            if (UType == KI || LType == KI){

                if (LType != NA){ //Ko, KI
                    fd1dRet->barrier[0]->needSpecialFDDownBar = true; //KIKO
                    if (LType == KO){
                        fd1dRet->barrier[1]->needSpecialFDDownBar = true; //KIKO
                    }
                }

                if (UType != NA){//KO, KI
                    fd1dRet->barrier[0]->needSpecialFDUpBar = true; //KIKO
                    if(UType == KO){
                        fd1dRet->barrier[1]->needSpecialFDUpBar = true; //KIKO
                    }
                }
            }else{        
                if (UType == KO){
                    fd1dRet->barrier[0]->needSpecialFDUpBar = true; 
                }

                if(LType == KO){
                    fd1dRet->barrier[0]->needSpecialFDDownBar = true; 
                }
            }
        }

        prepDeadInstr();

        initSlices( numPrices, inst->discount->getName() );
        initInsertNode();

        // if tree1f engine then use original way (for performance) otherwise slice operators
        if( tree1f )
        {
            prod_BWD_T = &DblBarrierFDProd::prod_BWD_T_orig;
            prod_BWD = &DblBarrierFDProd::prod_BWD_orig;
        }
        else
        {
            prod_BWD_T = &DblBarrierFDProd::prod_BWD_T_oper;
            prod_BWD = &DblBarrierFDProd::prod_BWD_oper;
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** calculate barriers and place barriers at inserted node if needed */
void DblBarrierFDProd::preCalc(int step)
{
    // set barrier, make adjustment for once a day monitoring if needed
    // to do: forward start adjustment
    UpperBar = UBarInput[step];
    LowerBar = LBarInput[step];

    // taking care of discrete monitoring, treat as daily monitor
    if( tree1f )
    {
        // taking care of discrete monitoring, treat as daily monitor
        if (!inst->IntraDayMonitor)
        {
            vector<double> vol;
            // adjust barrier if needed
            if (UType != NA && MonDepend != LOWER)
            {
                tree1f->GetStepVol(step, vol, &UpperBar, 0, 0); // get vol at barrier
                Barrier::BarrierAdjustment(vol[0], true, UpperBar);
            }
            if (LType != NA && MonDepend != UPPER)
            {
                tree1f->GetStepVol(step, vol, &LowerBar, 0, 0); // get vol at barrier
                Barrier::BarrierAdjustment(vol[0], false, LowerBar);
            }
        }

        if (tree1f->GetSmoothMethod() == CTree1f::NODE_INSERTION) 
        { //assuming Lower Barrier is index 0 and Upper uses 1 !!!
            bool isActiveInsNode = (LType != NA || UType != NA);
            int idx = tree1f->getSliceIndex( step );
            tree1f->SetInsertNode(idx, 0, LowerBar, (int)LType, isActiveInsNode); //lower barrier
            tree1f->SetInsertNode(idx, 1, UpperBar, (int)UType, isActiveInsNode); // upper barrier
        }
    }
    else if( fd1dRet )
    {
        //to change, 
        //for FD1D
        // set barrier, make adjustment for once a day monitoring if needed
        // to do: forward start adjustment

        // taking care of discrete monitoring, treat as daily monitor
        if (!inst->IntraDayMonitor)
        {
            vector<double> vol;
            // adjust barrier if needed
            if (UType != NA)
            {
                fd1dRet->GetStepVol(step, vol, &UpperBar, 0, 0); // get vol at barrier
                Barrier::BarrierAdjustment(vol[0], true, UpperBar);
            }
            if (LType != NA)
            {
                fd1dRet->GetStepVol(step, vol, &LowerBar, 0, 0); // get vol at barrier
                Barrier::BarrierAdjustment(vol[0], false, LowerBar);
            }
        }

        if (UType == KI || LType == KI){

            if (UType == KI)   {
                fd1dRet->getUpBarrier( 0, UpperBar);

                if (LType == KO){ 
                    //KOKI, use p[1] =KO
                    fd1dRet->getDownBarrier( 0, LowerBar);                
                    fd1dRet->getDownBarrier( 1, LowerBar);
                }
            }

            if (LType == KI)   {
                fd1dRet->getDownBarrier( 0, LowerBar);

                if (UType == KO){ 
                    //KOKI, use p[1] =KO
                    fd1dRet->getUpBarrier( 0, UpperBar);                
                    fd1dRet->getUpBarrier( 1, UpperBar);
                }
            }
        }else{
            if (UType != NA) {
                fd1dRet->getUpBarrier( 0, UpperBar);
            } 

            if (LType != NA) {
                fd1dRet->getDownBarrier( 0, LowerBar);
            }
        }
    }
    else
    {
        const FD1DLV * modelLV = dynamic_cast< const FD1DLV * >( model );
        if( modelLV && ! inst->IntraDayMonitor )
        {
            vector< double > vol;
            // adjust barrier if needed
            if( UType != NA && MonDepend != LOWER )
            {
                modelLV->GetStepVol( step, vol, &UpperBar, 0, 0 ); // get vol at barrier
                Barrier::BarrierAdjustment( vol[0], true, UpperBar );
            }
            if( LType != NA && MonDepend != UPPER )
            {
                modelLV->GetStepVol( step, vol, &LowerBar, 0, 0 ); // get vol at barrier
                Barrier::BarrierAdjustment( vol[0], false, LowerBar );
            }
        }

        const TreeSlice & s = payoffIndex->getValue( step );
        if( UType != NA )
        {
            if( UType == KO )
            {
                model->addCriticalLevel(
                    step, s, UpperBar, *slices[ 0 ], FDModel::LEVEL_BARRIER_UP, UBarRebate[step] );

                if( LType == KI )
                {
                    model->addCriticalLevel(
                        step, s, UpperBar, *slices[ 1 ], FDModel::LEVEL_BARRIER_UP, UBarRebate[step] );
                }
            }
            else
            {
                model->addCriticalLevel(
                    step, s, UpperBar, *slices[ 0 ], FDModel::LEVEL_BARRIER_UP );
            }
        }
        if( LType != NA )
        {
            if( LType == KO )
            {
                model->addCriticalLevel(
                    step, s, LowerBar, *slices[ 0 ], FDModel::LEVEL_BARRIER_DOWN, LBarRebate[step] );

                if( UType == KI )
                {
                    model->addCriticalLevel(
                        step, s, LowerBar, *slices[ 1 ], FDModel::LEVEL_BARRIER_DOWN, LBarRebate[step] );
                }
            }
            else
            {
                model->addCriticalLevel(
                    step, s, LowerBar, *slices[ 0 ], FDModel::LEVEL_BARRIER_DOWN );
            }
        }
    }
}

/** isInitValue == true, payoff at T for backward or value at t=0 for fwd induction
    isInitValue == false, payoff boundary condition, for KO, early exercise etc. */
void DblBarrierFDProd::update(int& step, FDProduct::UpdateType type)
{
    // we assume just need one und level for spot here
    const TreeSlice & spot = payoffIndex->getValue( step );

    if( type == FDProduct::BWD_T )
    {
        ( this->*prod_BWD_T )( step, spot, slices );

        if( tree1f && numIns )
            ( this->*prod_BWD_T )( step, *insNodes, *insPrices );
    }
    else if( type == FDProduct::BWD )
    {
        ( this->*prod_BWD )( step, spot, slices );

        if( tree1f && numIns )
            ( this->*prod_BWD )( step, *insNodes, *insPrices );
    }
    else if( type == BWD_NODE_INSERTION )
    {
        ( this->*prod_BWD )( step, *insNodes, *insPrices );
    }
    else
    {
        // to do fwd induction
    }
}

/** check if knocked out already */
void DblBarrierFDProd::prepDeadInstr()
{
    UpIsOut = inst->UpperBarBreachedDynamic && UType == KO;
    DownIsOut = inst->LowerBarBreachedDynamic && LType == KO;

    if( inst->IntraDayMonitor )
    {
        double s = inst->asset->getSpot();
        if( UType == KO && s > UBarInput[0]*(1.0-FP_MIN) )
            UpIsOut = true;
        if( LType == KO && s < LBarInput[0]*(1.0+FP_MIN) )
            DownIsOut = true;
    }

    isDead =
        ( UpIsOut || DownIsOut ) && BarDepend != TWO_TOUCH || 
        ( UpIsOut && DownIsOut ) && BarDepend == TWO_TOUCH;
}

/** calculate dead value */
void DblBarrierFDProd::calcDeadInstr(
    int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
{
    // discount rebate if needed
    if( UType != NA && inst->UpperRebate.get() )
    {
        if( ! Maths::isZero( UBarRebate[step] ) )
            *price[2] = UBarRebate[step]; // index [2] is upper barrier option
    }
    if( LType != NA && inst->LowerRebate.get() )
    {
        if( ! Maths::isZero( LBarRebate[step] ) )
            *price[1] = LBarRebate[step]; // index [1] is lower barrier option
    }

    // to do: need to decide which rebate to pay for NO_CANCEL on double KO ?
    if( UpIsOut )
        *price[0] = *price[2];
    else if( DownIsOut )
        *price[0] = *price[1];
}

/** product payoff method at steps earlier than maturity */
void DblBarrierFDProd::prod_BWD_T_orig(
    int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
{
    static const string method = "CDblBarrier::prod_BWD_T";    

    int bot, top;
    spot.getCalcRange( bot, top );
    double * s = spot.getValues();
    const vector< double * > & p = getValues( price );

    if( isDead )
    {
        calcDeadInstr( step, spot, price );
        // shrink node array to be calculated, Move to step 0 will be even more efficient
        if( tree1f )
            tree1f->ResetNodeBoundary(step, tree1f->GetSliceIdx(), 2, 2);

        return;
    }

    double settlementPV = inst->instSettle->pvAdjust(
        inst->exerciseSchedule->lastDate(),
        inst->discount.get(),
        inst->asset.get() );

    if( UType == KI || LType == KI )
    {
        //validate
        if( UType == KI && LType == KI && UBarRebate[step] != LBarRebate[step] )
            throw ModelException( method, "UpperRebate and LowerRebate should be same for Double KI case." );

        if( ( UType == KO || LType == KO ) && BarDepend == KI_CANCEL_KO )
        {
            // KIKO array [0], KO only array [1], Vanilla array [2]

            for( int j = bot; j <= top; ++j )
            {
                p[2][j] = settlementPV * ( PayMode == BINARY ?
                    1. : GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD) );

                if( UType == KO && s[j] > UpperBar * (1.0-FP_MIN) )
                    p[0][j] = p[1][j] = UBarRebate[step]; // dead on upper barrier
                else if( LType == KO && s[j] < LowerBar * (1.0+FP_MIN) )
                    p[0][j] = p[1][j] = LBarRebate[step]; // dead on lower barrier
                else
                {
                    p[1][j] = p[2][j];

                    if( s[j] > UpperBar * (1.0-FP_MIN) || s[j] < LowerBar * (1.0+FP_MIN) )
                        p[0][j] = p[1][j];
                    else
                    {
                        if( UType == KI )
                            p[0][j] = UBarRebate[step]; // failed to KI upper bar
                        else // if( LType == KI )
                            p[0][j] = LBarRebate[step]; // failed to KI lower bar
                    }      
                }
            }
        }
        else
        {
            // KIKO array [0], KO only array [1]

            for( int j = bot; j <= top; ++j )
            {
                if( UType == KO && s[j] > UpperBar * (1.0-FP_MIN) )
                    p[0][j] = p[1][j] = UBarRebate[step]; // dead on upper barrier
                else if( LType == KO && s[j] < LowerBar * (1.0+FP_MIN) )
                    p[0][j] = p[1][j] = LBarRebate[step]; // dead on lower barrier
                else
                {
                    p[1][j] = settlementPV * ( PayMode == BINARY ?
                        1. : GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD) );

                    if( s[j] > UpperBar * (1.0-FP_MIN) || s[j] < LowerBar * (1.0+FP_MIN) )
                        p[0][j] = p[1][j];
                    else
                    {
                        if( UType == KI )
                            p[0][j] = UBarRebate[step]; // failed to KI upper bar
                        else // if( LType == KI )
                            p[0][j] = LBarRebate[step]; // failed to KI lower bar
                    }      
                }
            }
        }
    }
    else
    {
        // KO array [0]

        for( int j = bot; j <= top; ++j )
        {
            if( UType == KO && s[j] > UpperBar * (1.0-FP_MIN) )
                p[0][j] = UBarRebate[step]; // dead on upper barrier
            else if( LType == KO && s[j] < LowerBar * (1.0+FP_MIN) )
                p[0][j] = LBarRebate[step]; // dead on lower barrier
            else
            {
                p[0][j] = settlementPV * ( PayMode == BINARY ?
                    1. : GetIntrinsic(s[j], stepStrike[step], inst->isCall, PayMode!=FORWARD) );
            }
        }
    }
}
void DblBarrierFDProd::prod_BWD_T_oper(
    int step, const TreeSlice & spot, const vector< TreeSliceSP > & price )
{
    static const string method = "CDblBarrier::prod_BWD_T";    

    if( isDead )
    {
        calcDeadInstr( step, spot, price );
        // shrink node array to be calculated, Move to step 0 will be even more efficient
        if( tree1f )
            tree1f->ResetNodeBoundary(step, tree1f->GetSliceIdx(), 2, 2);

        return;
    }

    double settlementPV = inst->instSettle->pvAdjust(
        inst->exerciseSchedule->lastDate(),
        inst->discount.get(),
        inst->asset.get() );

    double callput = inst->isCall ? settlementPV : -settlementPV;

    #define payoff()                                                                               \
        IF( PayMode == BINARY )                                                                    \
            settlementPV                                                                           \
        ELSE                                                                                       \
            IF( PayMode == FORWARD )                                                               \
                callput * ( spot - stepStrike[step] )                                              \
            ELSE                                                                                   \
                smax( callput * ( spot - stepStrike[step] ), 0. )                                  \
            ENDIF                                                                                  \
        ENDIF                                                                                      \

    #define ko_payoff( x )                                                                         \
        IF( UType == KO && UpperBar * (1.0-FP_MIN) < spot )                                        \
            UBarRebate[step]                                                                       \
        ELSE                                                                                       \
            IF( LType == KO && spot < LowerBar * (1.0+FP_MIN) )                                    \
                LBarRebate[step]                                                                   \
            ELSE                                                                                   \
                x                                                                                  \
            ENDIF                                                                                  \
        ENDIF                                                                                      \

    #define ki_payoff()                                                                            \
        IF( spot < LowerBar * (1.0+FP_MIN) || UpperBar * (1.0-FP_MIN) < spot )                     \
            *price[1]                                                                              \
        ELSE                                                                                       \
            IF( UType == KI )                                                                      \
                UBarRebate[step]                                                                   \
            ELSE                                                                                   \
                LBarRebate[step]                                                                   \
            ENDIF                                                                                  \
        ENDIF                                                                                      \

    if( UType == KI || LType == KI )
    {
        //validate
        if( UType == KI && LType == KI && UBarRebate[step] != LBarRebate[step] )
            throw ModelException( method, "UpperRebate and LowerRebate should be same for Double KI case." );

        if( ( UType == KO || LType == KO ) && BarDepend == KI_CANCEL_KO )
        {
            // KIKO array [0], KO only array [1], Vanilla array [2]

            //*price[2] = payoff;
            if( PayMode == BINARY )
                *price[2] = settlementPV;
            else if( PayMode == FORWARD )
                *price[2] = callput * ( spot - stepStrike[step] );
            else
            {
                *price[2] = smax( 0., callput * ( spot - stepStrike[step] ) );
            }

            *price[1] = ko_payoff( *price[2] );
            *price[0] = ko_payoff( ki_payoff() );
        }
        else
        {
            // KIKO array [0], KO only array [1]

            *price[1] = ko_payoff( payoff() );
            *price[0] = ko_payoff( ki_payoff() );
        }
    }
    else
    {
        // KO array [0]

        *price[0] = ko_payoff( payoff() );
    }

    #undef payoff
    #undef ko_payoff
    #undef ki_payoff
}

/** product payoff method at steps earlier than maturity */
void DblBarrierFDProd::prod_BWD_orig(
    int step, const TreeSlice & spot, const vector< TreeSliceSP > & price)
{
    static const string method = "DblBarrierFDProd::prodBCDS";
    try
    {
        int bot, top;
        spot.getCalcRange( bot, top );
        double * s = spot.getValues();
        const vector< double * > & p = getValues( price );

        if( ( UpIsOut || DownIsOut ) && BarDepend != TWO_TOUCH || 
            ( UpIsOut && DownIsOut ) && BarDepend == TWO_TOUCH )
        {
            calcDeadInstr( step, spot, price );
            return;
        }

        // ******* penultimate soothing not done ***************
        // forbid KO/KI on holidays/weekends via stepCanKOUp/stepCanKODown

        double settlementPV = inst->instSettle->pvAdjust(
            model->getDate(step),
            inst->discount.get(),
            inst->asset.get() );

        double callput = inst->isCall ? settlementPV : -settlementPV;

        if( UType == KI || LType == KI )
        {
            // for KI condition below
            bool canKO = stepCanKOUp[step] || stepCanKODown[step];

            if( ( UType == KO || LType == KO ) && BarDepend == KI_CANCEL_KO )
            {
                // KIKO array [0], KO only array [1], Vanilla array [2]

                for( int j = bot; j <= top; ++j )
                {
                    if( stepCanExercise[step] ) // American
                    {
                        p[2][j] = Maths::max(callput*(s[j] - stepStrike[step]),p[2][j]);
                        p[1][j] = Maths::max(callput*(s[j] - stepStrike[step]),p[1][j]);
                    }

                    if( UType == KO && stepCanKOUp[step] && s[j] > UpperBar * (1.0-FP_MIN) )
                        p[0][j] = p[1][j] = UBarRebate[step]; // dead on upper barrier
                    else if( LType == KO && stepCanKODown[step] && s[j] < LowerBar * (1.0+FP_MIN) )
                        p[0][j] = p[1][j] = LBarRebate[step]; // dead on lower barrier
                    else if( canKO && ( s[j] > UpperBar * (1.0-FP_MIN) || s[j] < LowerBar * (1.0+FP_MIN) ) )
                        p[0][j] = p[2][j]; // KI already
                }
            }
            else
            {
                // KIKO array [0], KO only array [1]

                for( int j = bot; j <= top; ++j )
                {
                    if( stepCanExercise[step] ) // American
                        p[1][j] = Maths::max(callput*(s[j] - stepStrike[step]),p[1][j]);

                    if( UType == KO && stepCanKOUp[step] && s[j] > UpperBar * (1.0-FP_MIN) )
                        p[0][j] = p[1][j] = UBarRebate[step]; // dead on upper barrier
                    else if( LType == KO && stepCanKODown[step] && s[j] < LowerBar * (1.0+FP_MIN) )
                        p[0][j] = p[1][j] = LBarRebate[step]; // dead on lower barrier
                    else if( canKO && ( s[j] > UpperBar * (1.0-FP_MIN) || s[j] < LowerBar * (1.0+FP_MIN) ) )
                        p[0][j] = p[1][j]; // KI already
                }
            }
        }
        else
        {
            // KO array [0]

            for( int j = bot; j <= top; ++j )
            {
                if( UType == KO && stepCanKOUp[step] && s[j] > UpperBar * (1.0-FP_MIN) )
                    p[0][j] = UBarRebate[step]; // dead on upper barrier
                else if( LType == KO && stepCanKODown[step] && s[j] < LowerBar * (1.0+FP_MIN) )
                    p[0][j] = LBarRebate[step]; // dead on lower barrier
                else if( stepCanExercise[step] ) // American
                    p[0][j] = Maths::max( callput*( s[j] - stepStrike[step] ), p[0][j] );
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}
void DblBarrierFDProd::prod_BWD_oper(
    int step, const TreeSlice & spot, const vector< TreeSliceSP > & price)
{
    static const string method = "DblBarrierFDProd::prodBCDS";
    try
    {
        if( ( UpIsOut || DownIsOut ) && BarDepend != TWO_TOUCH || 
            ( UpIsOut && DownIsOut ) && BarDepend == TWO_TOUCH )
        {
            calcDeadInstr( step, spot, price );
            return;
        }

        // ******* penultimate soothing not done ***************
        // forbid KO/KI on holidays/weekends via stepCanKOUp/stepCanKODown

        double settlementPV = inst->instSettle->pvAdjust(
            model->getDate(step),
            inst->discount.get(),
            inst->asset.get() );

        double callput = inst->isCall ? settlementPV : -settlementPV;

        #define payoff( x )                                                                        \
            IF( (bool)stepCanExercise[step] )                                                      \
                smax( callput * ( spot - stepStrike[step] ), x )                                   \
            ELSE                                                                                   \
                x                                                                                  \
            ENDIF                                                                                  \

        #define ko_payoff( x )                                                                     \
            IF( UType == KO && stepCanKOUp[step] && UpperBar * (1.0-FP_MIN) < spot )               \
                UBarRebate[step]                                                                   \
            ELSE                                                                                   \
                IF(LType == KO && stepCanKODown[step] && spot < LowerBar * (1.0+FP_MIN) )          \
                    LBarRebate[step]                                                               \
                ELSE                                                                               \
                    x                                                                              \
                ENDIF                                                                              \
            ENDIF                                                                                  \

        #define ki_payoff( x, y )                                                                  \
            IF( canKO && ( spot < LowerBar * (1.0+FP_MIN) || UpperBar * (1.0-FP_MIN) < spot ) )    \
                x                                                                                  \
            ELSE                                                                                   \
                y                                                                                  \
            ENDIF                                                                                  \

        if( UType == KI || LType == KI )
        {
            // for KI condition below
            bool canKO = stepCanKOUp[step] || stepCanKODown[step];

            if( ( UType == KO || LType == KO ) && BarDepend == KI_CANCEL_KO )
            {
                // KIKO array [0], KO only array [1], Vanilla array [2]

                // *price[2] = payoff( *price[2] );
                if( stepCanExercise[step] ) // American
                    *price[2] = smax( callput * ( spot - stepStrike[step] ), *price[2] );

                *price[1] = ko_payoff( payoff( *price[1] ) );
                *price[0] = ko_payoff( ki_payoff( *price[2], *price[0] ) );
            }
            else
            {
                // KIKO array [0], KO only array [1]

                *price[1] = ko_payoff( payoff( *price[1] ) );
                *price[0] = ko_payoff( ki_payoff( *price[1], *price[0] ) );
            }
        }
        else
        {
            // KO array [0]

            *price[0] = ko_payoff( payoff( *price[0] ) );
        }

        #undef payoff
        #undef ko_payoff
        #undef ki_payoff
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** premium scaling */
double DblBarrierFDProd::scalePremium(const double& fairValue, YieldCurveConstSP disc)
{
    double fwdAtStart = 0.0;
    double fwdStartDF = 1.0;
    if (inst->fwdStarting)
    {
        fwdAtStart = inst->asset->fwdValue(inst->startDate);
        fwdStartDF = disc->pv(inst->valueDate, model->getDate(0));
    }

    double scalingFactor = InstrumentUtil::scalePremium(
                                    inst->oneContract,
                                    inst->fwdStarting,
                                    inst->notional,
                                    fwdAtStart,
                                    inst->initialSpot);

    return (fairValue*scalingFactor*fwdStartDF);
}

//output results 
void DblBarrierFDProd::recordOutput(Control* control, YieldCurveConstSP disc, Results* results)
{
    // get prices at t=0
    // save price
    double price = scalePremium(model->getPrice0( *slices[0] ), disc);
    results->storePrice(price, disc->getCcy());

    if (control && control->isPricing()) {
        DateTime       matDate = inst->exerciseSchedule->lastDate();
        double         indVol;
        // calculate indicative vol
        try {
            if ( matDate.isGreater(inst->valueDate) )
            {

                DateTime imntStartDate = inst->fwdStarting? 
                    inst->startDate: inst->valueDate;

                // get vol request
                CVolRequestConstSP lnVolRequest = GetLNRequest();

                // interpolate the vol
                CVolProcessedSP  vol(inst->asset->getProcessedVol(lnVolRequest.get()));
                // cast to the type of vol we're expecting
                CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
                // this should never happen if our get market data has worked properly
                if (!vol){
                    throw ModelException("DblBarrierFDProd::recordOutput", 
                                         "No Black Scholes Vol");
                }

                // calculate the indicative vol
                indVol = volBS->CalcVol(imntStartDate, matDate);
            }
            else
            {
                indVol = 0.0;
            }
        }
        catch(exception&)
        {// continue if indicative vol fail - non BS vol
            indVol = 0.0;
        }

        inst->addOutputRequests(control,
                                    results,
                                    price,
                                    indVol);

        OutputRequest* request = control->requestsOutput(OutputRequest::BARRIER_DISCONTINUITY_RISK);
        
        // only available if the instrument has started
        // and only one barrier
        if (request) {
            if (                
                // not foward starting
                (!inst->fwdStarting) &&
                // at least one NA barrier
                (UType==NA || LType==NA) &&
                // no more than one NA barrier
                (UType!=NA || LType!=NA)) {
                
                // select the barrier, rebate and spot tweak
                ScheduleSP barrier;
                ScheduleSP rebate;
                double spotTweak;
                DateTime today = inst->valueDate;
                
                if (UType!=NA) {
                    // upper barrier
                    barrier = ScheduleSP(inst->UpperBarrier);
                    rebate = ScheduleSP(inst->UpperRebate);
                    spotTweak = 0.99;
                }
                else {
                    //down barrier
                    barrier = ScheduleSP(inst->LowerBarrier);
                    rebate = ScheduleSP(inst->LowerRebate);
                    spotTweak = 1.01;
                }
                
                // the barrier and rebate may be bend 
                // we use today barrier level
                // already started so absolute level
               
                
                // returns the difference of prices
                try {
                    double barrierLevel = barrier->interpolate(today);
                    // calculate the difference of price when the spot is 1% out and 1%in the barrier
                    double tweakPrice = tweakSpotAndPrice(spotTweak*barrierLevel, model, inst) 
                                        - tweakSpotAndPrice((2.0-spotTweak)*barrierLevel, model, inst);
                    results->storeRequestResult(request,tweakPrice);
                }
                catch (exception&){
                    // nothing to do
                }
            }
        }    
    }
}

/** create a fd payoff product */
FDProductSP CDblBarrier::createProduct(FDModel* model) const
{
    return FDProductSP( new DblBarrierFDProd(this, model) );
}

/// ------------  end of new FD product -----------

DRLIB_END_NAMESPACE
