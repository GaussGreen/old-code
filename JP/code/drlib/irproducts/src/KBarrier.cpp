//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Description : KnockOut component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KBarrier.hpp"
#include "edginc/UtilFuncs.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/TreeSliceOperExt.hpp"
#include "edginc/Addin.hpp"

#define  BARRIER_TOL   1E-4   /* Error tolerance for barriers */

DRLIB_BEGIN_NAMESPACE

/*************/

class KBarrierTree : public FDProduct {
public:

    /************************ variables ************************/
    const KBarrier* inst;    
    FDProductSP     und;
    /************************ methods ************************/

    KBarrierTree(const KBarrier* inst, FDModel* model);

    DateTime getStartDate(void) const {return model->getValueDate();}
    virtual string getOutputName(void) const { return inst->outputName; }

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);
    void addModelResetDates(const DateTimeArray &modelResetDates, 
                            const DateTimeArray &lateResetDates);
    virtual void init(Control*) const {}
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual void initProd(void);

private:
    TreeSliceSP slice;
};

/******************************** KBarrierTree ********************************/

KBarrierTree::KBarrierTree(const KBarrier* inst, FDModel* model) 
: FDProduct(model), inst(inst)
{   
    try {
        und = model->createProduct(inst->und);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KBarrierTree::addModelResetDates(const DateTimeArray &modelResetDates, 
                                      const DateTimeArray &lateResetDates)
{
    try {
        und->addModelResetDates(modelResetDates, lateResetDates);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/** initialising and setting product variables */
void KBarrierTree::initProd(void) {
    try {
        slice = model->createSlice();
        slice->name = inst->outputName + "_KBarrier";
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

struct oper_smooth_barrier : SliceOper< double >
{
    static const char* symbol(){return "oper_ko_smooth_barrier";}
    static Type apply( double value, double level, double smoothStep, double upVal, double downVal) 
    {        
        return SmoothValue( upVal,   level + smoothStep, 
                            downVal, level - smoothStep,
                            value);
    }
};

struct oper_smooth_double_barrier : SliceOper< double >
{
    static const char* symbol(){return "oper_ko_smooth_double_barrier";}
    static Type apply( double value, double levelLo, double levelHi, double smoothStep, double upVal, double downVal) 
    {        
        double x = SmoothValue( upVal,   levelLo + smoothStep, 
                                downVal, levelLo - smoothStep,
                                value);

        return SmoothValue( downVal, levelHi + smoothStep,
                            x,       levelHi - smoothStep,
                            value);
    }
};

const TreeSlice & KBarrierTree::getValue(int step, DateTime eventDate) const { 
    try {
        if (historicalValueAvailable(*inst, *slice, eventDate))
            return *slice;

        double upVal, downVal;
        if (inst->upIsTrue) {
            upVal = 1.;
            downVal = 0.;
        }
        else {
            upVal = 0.;
            downVal = 1.;
        }
        double level = inst->levels->interpolate(eventDate);
        TreeSlice const& undSl = und->getValue(step, eventDate);
        double tol = BARRIER_TOL;

        if (inst->levelsHi.get()) {
            double levelHi = inst->levelsHi->interpolate(eventDate);
            if (inst->smoothing == KBarrier::Smoothing::STEP) {
                TreeSliceSP smoothStep = undSl.calcSmoothStep();
                *slice = Slice6Expr<oper_smooth_double_barrier, const TreeSlice, double, double, const TreeSlice, double, double>
                    (undSl, level, levelHi, *smoothStep, upVal, downVal);
                return *slice;
            }
            else {
                if (inst->onBarrierIsUp) {
                    *slice = cond(undSl <= (level - tol), 
                                  downVal, 
                                  cond(undSl < (levelHi + tol), upVal, downVal));
                }
                else *slice = cond(undSl < (level + tol), 
                                   downVal, 
                                   cond(undSl <= (levelHi - tol), upVal, downVal));
            }
        }
        else {
            if (inst->smoothing == KBarrier::Smoothing::STEP) {
                TreeSliceSP smoothStep = undSl.calcSmoothStep();
                *slice = Slice5Expr<oper_smooth_barrier, const TreeSlice, double, const TreeSlice, double, double>
                    (undSl, level, *smoothStep, upVal, downVal);
                return *slice;
            }
            else if (inst->onBarrierIsUp) {
                *slice = cond(undSl <= (level - tol), downVal, upVal);
            }
            else *slice = cond(undSl < (level + tol), downVal, upVal);
        }
        return *slice; 
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KBarrierTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) 
{
    try {
        OutputRequest* request
            = ctrl->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (!request) {
            return;
        }
        DateTime valueDate = model->getValueDate();
        DateTime upperDate = BarrierLevel::barrierWindow(valueDate);
        upperDate = valueDate.rollDate(100000); // to ease debug

        CashFlowArraySP subset(inst->levels->subset(valueDate, upperDate));
        bool isUp = !inst->upIsTrue;
        BarrierLevelArray levels;
        for (int b=0;;++b) {
            if (!subset->empty()) {
                for (int i = 0; i < subset->size(); i++) {
                    BarrierLevel bl(isUp, (*subset)[i].date, (*subset)[i].amount);
                    levels.push_back(bl);
                }
            }
            if (b==1)
                break;
            if (!inst->levelsHi) 
                break;
            subset.reset(inst->levelsHi->subset(valueDate, upperDate));
            isUp = !isUp;
        }
        if (levels.size()) {
            OutputRequestUtil::recordBarrierLevels(
                ctrl, results, inst->outputName, &levels);
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}
/********************************** KBarrier **********************************/

void KBarrier::addResetDates(const DateTimeArray &resetDates) {
    try {
        und->addResetDates(resetDates);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

// the IGetKnownValue interface for retrieving a single sample
double KBarrier::getValue(DateTime date, CashflowInfo &cfi) const
{
    try {
        CashflowInfo cfiLoc; // local CashflowInfo, in case we want to cancel the computation
        double undVal = und->getValue(date, cfiLoc);
        double level = levels->interpolate(date);
        double result = 1.;
        double upVal, downVal;
        bool onIsDown = !onBarrierIsUp;
        if (upIsTrue) {
            upVal = 1.;
            downVal = 0.;
        }
        else {
            upVal = 0.;
            downVal = 1.;
        }
        for (int b=0;;++b) {
            result *= (undVal < level || (undVal==level && onIsDown) ? downVal : upVal);
            if (b==1)
                break;
            if (!levelsHi)
                break;
            level = levelsHi->interpolate(date);
            onIsDown = !onIsDown;
            swap(downVal, upVal);
        }

        if (cfiLoc.amountType != CashflowInfo::AmountType::KNOWN) {
            // we need to know the exact value, an estimation is not good enough.
            cfi.updateAmountType(CashflowInfo::AmountType::UNKNOWN);
            return 0.;
        }
        cfi.merge(cfiLoc, outputName);
        return result;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

void KBarrier::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        und->setup(model, market);
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

/** implement FDModel product interface */
FDProductSP KBarrier::createProduct(FDModel * model) const {    
    if (!setupCalled) throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KBarrierTree(this, model));
}

void KBarrier::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("KBarrier composite component");
    REGISTER(KBarrier, clazz)
    SUPERCLASS(KComponent);
    IMPLEMENTS(IIndicCreator);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(KBarrier::defaultConstructor);
    FIELD(und, "Underlying index to observe")
    FIELD(levels, "");
    FIELD(levelsHi, "If this is provided, it is mean to be the high barrier.");
    FIELD_MAKE_OPTIONAL(levelsHi);
    FIELD(onBarrierIsUp, "");
    FIELD_MAKE_OPTIONAL(onBarrierIsUp);
    FIELD(upIsTrue, "");
    FIELD_MAKE_OPTIONAL(upIsTrue);
    FIELD(smoothing, "");
    FIELD_MAKE_OPTIONAL(smoothing);
    Addin::registerConstructor(Addin::UTILITIES, KBarrier::TYPE);
}
   
CClassConstSP const KBarrier::TYPE = CClass::registerClassLoadMethod(
    "KBarrier", typeid(KBarrier), KBarrier::load);

START_PUBLIC_ENUM_DEFINITION(KBarrier::Smoothing::Enum, "");
ENUM_VALUE_AND_NAME(KBarrier::Smoothing::NONE, "NONE", "");
ENUM_VALUE_AND_NAME(KBarrier::Smoothing::STEP, "STEP", "");
END_ENUM_DEFINITION(KBarrier::Smoothing::Enum);

/**********************************/
bool KBarrierLoad(void) {
    return (KBarrier::TYPE !=0);
}

DRLIB_END_NAMESPACE
