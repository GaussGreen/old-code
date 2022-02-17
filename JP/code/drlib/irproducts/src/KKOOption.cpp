//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KKOOption.cpp
//
//   Description : KO option component - KKOOption, a subclass of KOption
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KKOOption.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/RatesUtils.hpp"
#include "edginc/OutputRequestUtil.hpp"

DRLIB_BEGIN_NAMESPACE

class KKOOptionTree : public KOptionTree {
public:
    /************************ methods ************************/
    KKOOptionTree(const KKOOptionConstSP &inst, FDModel* model);

    //virtual void init(Control*) const;
    virtual void initProd();
    virtual void update(int& step, UpdateType type);
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);

protected:
    /************************ variables ************************/
    KKOOptionConstSP    inst;
    FDProductSP         obs;
    TreeSliceSP         KIValue;
};

/******************************* KKOOptionTree ********************************/

KKOOptionTree::KKOOptionTree(const KKOOptionConstSP &inst, FDModel* model) :
    KOptionTree(inst, model), inst(inst)
{
    try {
        obs = model->createProduct(inst->obsUnd);
        obs->addModelResetDates(inst->obsDates);
    }
    catch (exception& e){
        string errMsg = "Error creating dependent product of KKOOption";
        throw ModelException(e, "KKOOptionTree::KKOOptionTree", errMsg);
    }
}

void KKOOptionTree::recordOutput(Control* ctrl, YieldCurveConstSP yc, Results* results) {
    string ccy = inst->discount->getCcy();

    if (inst->recordOutputName)
        recordSliceToOutputName(ctrl, results, model, 
        inst->isTopComponent(), inst->outputName, ccy, getValue(0, model->getDate(0)));

    try {
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
        
    } 
    catch (exception& e) {
        throw ModelException(e, "KKOOptionTree::recordOutput, failed to report barriers for outputName " +
                             getOutputName());
    }
}

/** initialising and setting product variables */
void KKOOptionTree::initProd(void){

    KOptionTree::initProd();

    if (inst->isKnockIn) {
        const string curve = inst->discount->getName();
    
        KIValue = model->createSlice(curve);
        KIValue->name = inst->outputName + "_KIValue";
        *KIValue = 0.;
        startDEV(KIValue);
    }
}

/** update products slice at each step */
void KKOOptionTree::update(int& step, UpdateType type) {
    try {
        KOptionTree::update(step, type);
        
        DateTime stepDate = model->getDate(step);

        if (! inst->isObsDate(stepDate) )
            return;

        double koLevel = inst->getKOLevel(stepDate);

        if (inst->isKnockIn) {
            if (inst->isUpOut)
                *KIValue = cond(obs->getValue(step, stepDate) < koLevel, *KIValue, *optionPrice);
            else
                *KIValue = cond(obs->getValue(step, stepDate) < koLevel, *optionPrice, *KIValue);
        }
        else {
            if (inst->isUpOut) 
                *optionPrice = cond(obs->getValue(step, stepDate) < koLevel, *optionPrice, 0.0);
            else
                *optionPrice = cond(obs->getValue(step, stepDate) < koLevel, 0.0, *optionPrice);
        }

    } 
    catch (exception& e) {
        string errMsg = "KKOOptionTree::update at date " + model->getDate(step).toString();
        throw ModelException(e, errMsg);
    }
}

const TreeSlice & KKOOptionTree::getValue(int step, DateTime eventDate) const {
    // only support callable and puttable if timeStep = 0 at this stage
    // so it is a simple callable = underlying - option , puttable = underlying + option

    if (eventDate != model->getDate(step)) {
        throw ModelException(__FUNCTION__, "Cannot be valued at a date != currentDate");
    }

    double los;
    if (inst->longOrShort == KOption::LongOrShort::OPTION_LONG)
        los = 1.0;
    else
        los = -1.0;
    
    if (inst->isKnockIn) 
        *getValuePrice = *KIValue * los;
    else 
        *getValuePrice = *optionPrice * los;

    return *getValuePrice;
        
}

/********************************** KKOOption **********************************/

void KKOOption::setup(const IModel* model, const MarketData* market) {

    try {
        KOption::setup(model, market);

        if (obsUnd.get()) {
            obsUnd->addResetDates(obsDates);
            obsUnd->setup(model, market);
        }

        /***** validity checks & calc transient */
    }
    catch (exception& e) {
        throw ModelException(e, "KKOOption::setup");
    }
}

/*********************************** KKOSwap ***********************************/
bool KKOOption::isObsDate(const DateTime& date) const{
    return RatesUtils::happensNow(obsDates, date);
}


// get ko level
double KKOOption::getKOLevel(const DateTime& date) const{
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


/** implement FDModel product interface */
FDProductSP KKOOption::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KKOOptionTree(KKOOptionConstSP(this), model));
}


void KKOOption::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Knock In/Out Option component");
    REGISTER(KKOOption, clazz)
    SUPERCLASS(KOption);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(obsUnd, "observation index for KoOption");
    FIELD(obsDates, "KO observation dates");
    FIELD(isUpOut, "true= up and out (default). false = down and out"); 
    FIELD_MAKE_OPTIONAL(isUpOut);
    FIELD(isKnockIn, "Knock In/Out selection. default is KnockOut"); 
    FIELD_MAKE_OPTIONAL(isKnockIn);
    FIELD(barrier, "barrier schedule");
    Addin::registerConstructor(Addin::UTILITIES, KKOOption::TYPE);
}

CClassConstSP const KKOOption::TYPE = CClass::registerClassLoadMethod(
    "KKOOption", typeid(KKOOption), KKOOption::load);


bool KKOOptionLoad() {return KKOOption::TYPE !=0;}

DRLIB_END_NAMESPACE
