//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Description : Simple boolean expression (!)A (&& (!)B)
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KBoolExpr.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

class KBoolExprTree : public FDProduct {
public:
    KBoolExprConstSP     inst;
    FDProductSP arg1, arg2;

    /******************** methods ********************/
    KBoolExprTree(const KBoolExprConstSP &inst, FDModel* model);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void addModelResetDates(const DateTimeArray &modelResetDates, 
                                    const DateTimeArray &lateResetDates);
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {}
    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;

private:
    TreeSliceSP slice;
};

/****************************** KBoolExprTree ********************************/

KBoolExprTree::KBoolExprTree(const KBoolExprConstSP &inst, FDModel* model) 
: FDProduct(model), inst(inst) 
{
    try {
        arg1 = model->createProduct(inst->arg1);
        if (inst->arg2.get()) {
            arg2 = model->createProduct(inst->arg2);
        }
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

void KBoolExprTree::addModelResetDates(const DateTimeArray &modelResetDates, 
                                       const DateTimeArray &lateResetDates)
{
    arg1->addModelResetDates(modelResetDates, lateResetDates);
    if (arg2.get()) {
        arg2->addModelResetDates(modelResetDates, lateResetDates);
    }
}

/** initializing and setting product variables */
void KBoolExprTree::initProd(void) {
    slice = model->createSlice();
    slice->name = inst->outputName+"_slice";
}

const TreeSlice & KBoolExprTree::getValue(int step, DateTime eventDate) const {
    try {
        if (historicalValueAvailable(*inst, *slice, eventDate))
            return *slice;

        TreeSlice const& val1 = arg1->getValue(step, eventDate);
        bool not1 = inst->not1;
        if (!arg2) {
            if (not1)
                *slice = 1. - val1;
            else *slice = val1;
        }
        else {
            TreeSlice const& val2 = arg2->getValue(step, eventDate);
            bool not2 = inst->not2;
            if (not1 && not2)
                *slice = (1. - val1) * (1. - val2);
            else if (!not1 && not2) 
                *slice = val1 * (1. - val2);
            else if (not1 && !not2) 
                *slice = (1. - val1) * val2;
            else if (!not1 && !not2) 
                *slice = val1 * val2;
        }
        return *slice;
    }
    catch (exception& e){
        throw makeException(e, __FUNCTION__);
    }
}

/*********************************** KBoolExpr ***********************************/

void KBoolExpr::addResetDates(const DateTimeArray &resetDatesP) {
    arg1->addResetDates(resetDatesP);
    if (arg2.get()) {
        arg2->addResetDates(resetDatesP);
    }
}

void KBoolExpr::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        arg1->setup(model, market);
        if (arg2.get()) {
            arg2->setup(model, market);
        }
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

double KBoolExpr::getValue(DateTime date, CashflowInfo &cfi) const {
    try {
        CashflowInfo cfiLoc;
        double val;
        double val1 = arg1->getValue(date, cfiLoc);
        val = (not1 ? (1. - val1) : val1);
        if (arg2.get()) {
            double val2 = arg2->getValue(date, cfiLoc);
            val *= (not2 ? (1. - val2) : val2);
        }
        cfi.merge(cfiLoc, outputName);
        return val;
    }
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);
    }
}

FDProductSP KBoolExpr::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KBoolExprTree(KBoolExprConstSP(this), model));
}

void KBoolExpr::load(CClassSP& clazz) {

    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KBoolExpr, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(IIndicCreator);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(arg1, "");
    FIELD(arg2, "");
    FIELD_MAKE_OPTIONAL(arg2);
    FIELD(not1, "Apply a not operator to the first argument"); 
    FIELD_MAKE_OPTIONAL(not1);
    FIELD(not2, "Apply a not operator to the first argument"); 
    FIELD_MAKE_OPTIONAL(not2);
    Addin::registerConstructor(Addin::UTILITIES, KBoolExpr::TYPE);
}

CClassConstSP const KBoolExpr::TYPE = CClass::registerClassLoadMethod(
    "KBoolExpr", typeid(KBoolExpr), KBoolExpr::load);

/******************************/

bool KBoolExprLoad(void){
    return (KBoolExpr::TYPE != 0);
}

DRLIB_END_NAMESPACE

