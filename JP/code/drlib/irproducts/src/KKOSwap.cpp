//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKOSwap.hpp
//
//   Description : one touch KO swap component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KKOSwap.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/BarrierLevel.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

/****************************** KKOSwapTree ********************************/

class KKOSwapTree : public FDProduct {

public:
    /******************** methods ********************/
    KKOSwapTree(const KKOSwap* inst, FDModel* model);

    virtual DateTime getStartDate(void) const { return model->getDate(0);}

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);

    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;

protected:
    const KKOSwap*  inst;
    FDProductSP     underlying;
    FDProductSP     obs;

    TreeSliceSP value;
};

void KKOSwapTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    if (inst->recordOutputName)
        recordSliceToOutputName(ctrl, results, model, 
            inst->isTopComponent(), inst->outputName, "", getValue(0, model->getDate(0)));

    OutputRequest* request;
    request = ctrl->requestsOutput(OutputRequest::BARRIER_LEVEL);

    if (!request) {
        inst->recordExtraOutput(ctrl, results);
        return;
    }

    DateTime valueDate = model->getValueDate();
    DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
    upperDate = valueDate.rollDate(100000); // to ease debug
 
    Schedule* s = inst->barrier.get();
    bool isUp = inst->isUpOut;

        
    CashFlowArraySP subset(s->subset(valueDate, upperDate));
    if (!subset->empty()) {
        BarrierLevelArray levels;
        for (int i = 0; i < subset->size(); i++) {
            BarrierLevel bl(isUp, (*subset)[i].date, (*subset)[i].amount);
            levels.push_back(bl);
        }
        OutputRequestUtil::recordBarrierLevels(
            ctrl, results, inst->outputName, &levels);
    }
     

    inst->recordExtraOutput(ctrl, results);
}

KKOSwapTree::KKOSwapTree(const KKOSwap* inst, FDModel* model) :
    FDProduct(model), inst(inst) {

    try {
        underlying = model->createProduct(inst->underlier);
        underlying->addModelResetDates(inst->obsDates);

        obs = model->createProduct(inst->obsUnd);
        obs->addModelResetDates(inst->obsDates);
    }
    catch (exception& e){
        string errMsg = "Error creating dependent product of KKOSwap";
        throw ModelException(e, "KKOSwapTree::KKOSwapTree", errMsg);
    }
}

/** initializing and setting product variables */
void KKOSwapTree::initProd(void){
    const string curve = inst->discount->getName();

    value = model->createSlice(curve);
    value->name = inst->outputName+"_value";
    *value = 0.0;
    startDEV(value);
}

void KKOSwapTree::update(int& step, UpdateType type){
    DateTime stepDate = model->getDate(step);

    
    const TreeSlice &cashFlow = underlying->getCashFlow(step);
    if (!cashFlow.isZero()) 
        *value += cashFlow;

    if (! inst->isObsDate(stepDate) )
        return;

    double koLevel = inst->getKOLevel(stepDate);
    const TreeSlice &obsSlice = obs->getValue(step, stepDate);
    if (inst->isUpOut) 
        *value = cond(obsSlice < koLevel, *value, 0.0);
    else
        *value = cond(obsSlice < koLevel, 0.0, *value);
}

const TreeSlice & KKOSwapTree::getValue(int step, DateTime eventDate) const {
    if (eventDate != model->getDate(step)) {
        throw ModelException(__FUNCTION__, "Cannot be valued at a date != currentDate");
    }
    return *value;
}

/*********************************** KKOSwap ***********************************/
bool KKOSwap::isObsDate(const DateTime& date) const{
    return RatesUtils::happensNow(obsDates, date);
}

void KKOSwap::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        if (obsDates.empty()) {
            KComponentSP comp = KComponentSP::dynamicCast(underlier);
            comp->getCfDates(obsDates);
            DateTime::doSortUniq(obsDates);
        }
        if (obsUnd.get()) {
    	    obsUnd->addResetDates(obsDates);
            obsUnd->setup(model, market);
        }
        if (underlier.get()) {
            underlier->addResetDates(obsDates);
            underlier->setup(model, market);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KKOSwap::setup "+outputName);
    }
}

// get ko level
double KKOSwap::getKOLevel(const DateTime& date) const{
    // !!! temporary to make test work !!!
    double koLevel;
    if (barrier->lastDate() < date)
        koLevel = barrier->lastValue();
    else if (barrier->firstDate() > date)
        koLevel = barrier->firstValue();
    else
        koLevel = barrier->interpolate(date);

    return koLevel;
}

FDProductSP KKOSwap::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KKOSwapTree(this, model));
}

void KKOSwap::load(CClassSP& clazz) {

    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KKOSwap, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(underlier, "underlying swap leg"); FIELD_MAKE_OPTIONAL(underlier);
    FIELD(obsUnd, "observation index for KO"); FIELD_MAKE_OPTIONAL(obsUnd);
    FIELD(obsDates, "KO observation dates. If not provided, will take the underlying's payment dates.");
    FIELD_MAKE_OPTIONAL(obsDates);
    FIELD(isUpOut, "true= up and out. false = down and out"); 
    FIELD_MAKE_OPTIONAL(isUpOut);
    FIELD(barrier, "barrier schedule"); 
    FIELD_MAKE_OPTIONAL(barrier);
    Addin::registerConstructor(Addin::UTILITIES, KKOSwap::TYPE);
}

CClassConstSP const KKOSwap::TYPE = CClass::registerClassLoadMethod(
    "KKOSwap", typeid(KKOSwap), KKOSwap::load);

/******************************/
// for type linking
bool KKOSwapLoad(void){
    return (KKOSwap::TYPE != 0);
}

DRLIB_END_NAMESPACE

