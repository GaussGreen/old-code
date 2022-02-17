    //----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KSum.cpp
//
//   Description : sum component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KSum.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

class KSumTree : public FDProduct {
public:
    KSumConstSP     inst;
    FDProductArray  underlying;
    FDProductArray  fxProds;

    /******************** methods ********************/
    KSumTree(const KSumConstSP &inst, FDModel* model);
    virtual DateTime getStartDate(void) const { return model->getDate(0);}
    virtual string getOutputName(void) const { return inst->outputName; }
    virtual void addModelResetDates(
        const DateTimeArray &modelResetDates, 
        const DateTimeArray &lateResetDates);

    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
        if (inst->recordOutputName)
            recordSliceToOutputName(ctrl, results, model, 
                inst->isTopComponent(), inst->outputName, "", getValue(0, model->getDate(0)));
        inst->recordExtraOutput(ctrl, results);
    }

    virtual void init(Control*) const {}
    virtual void initProd(void);
    virtual void update(int& step, UpdateType type) {}
    virtual const TreeSlice & getValue(int step, DateTime eventDate) const;
    virtual const TreeSlice & getCashFlow(int step) const;

private:
    TreeSliceSP slice;
    TreeSliceSP cashFlowSlice;
};

/****************************** KSumTree ********************************/

KSumTree::KSumTree(const KSumConstSP &inst, FDModel* model) :
    FDProduct(model), inst(inst) {

    int i;
    try {
        int s=inst->listK.size();
        underlying.resize(s);
        fxProds.resize(s);
        for (i=0; i<s; ++i) {
            if (Maths::isZero(inst->weights[i])) 
                continue;
            underlying[i] = model->createProduct(inst->listK[i]);
            if (inst->fx[i].get())
                fxProds[i] = model->createProduct(inst->fx[i]);
        }
    }
    catch (exception& e){
        string errMsg = "Error creating dependent product of KSum (listK array index = " +
                        Format::toString(i) + "), object outputName = " + inst->outputName;
        throw ModelException(e, "KSumTree::KSumTree", errMsg);
    }
}

void KSumTree::addModelResetDates(
    const DateTimeArray &modelResetDates, 
    const DateTimeArray &lateResetDates)
{
    for (size_t i=0; i < underlying.size(); ++i) {
        if (!Maths::isZero(inst->weights[i])) {
            underlying[i]->addModelResetDates(modelResetDates, lateResetDates);
            if (inst->fx[i].get()) {
                fxProds[i]->addModelResetDates(modelResetDates, lateResetDates);
            }
        }
    }
}

/** initializing and setting product variables */
void KSumTree::initProd(void){
    // nothing to do as slice allocated when first required
    slice = model->createSlice();
    slice->name = inst->outputName+"_KSum_slice";

    cashFlowSlice = model->createSlice();
    cashFlowSlice->name = inst->outputName+"_KSum_cashFlowSlice";
}

const TreeSlice & KSumTree::getValue(int step, DateTime eventDate) const {
    try {
        if (historicalValueAvailable(*inst, *slice, eventDate))
            return *slice;

        int i;
        int nbUnderlying = underlying.size();

        DateTime currentDate = model->getDate(step);

        if (inst->constSchedule.get() && 
            currentDate.isGreaterOrEqual(inst->constSchedule->firstDate()) )
            *slice = inst->constant + inst->constSchedule->interpolate(currentDate);
        else
            *slice = inst->constant;
        
        if ( nbUnderlying==0 ) 
            return *slice;

        TreeSliceSP tmp = slice->clone(false);
        for (i=0; i<nbUnderlying; ++i) {
            double weight = inst->weights[i];

            if (Maths::isZero(weight)) 
                continue;

            *tmp = underlying[i]->getValue(step, eventDate);
            *slice += weight * (*tmp);
        }
        return *slice;
    }
    catch (exception& e) {
        throw ModelException(e, "KSumTree::getValue, timeStep = " + Format::toString(step) +
                             ", eventDate = " + eventDate.toString() + ", outputName = " + 
                             getOutputName());
    }
}

const TreeSlice & KSumTree::getCashFlow(int step) const {
    try {        
        *cashFlowSlice = 0.;
        TreeSliceSP tmp = cashFlowSlice->clone(false);
        int nbUnderlying = underlying.size();

        for (int i=0; i<nbUnderlying; ++i) {
            double weight = inst->weights[i];

            if (Maths::isZero(weight)) 
                continue;

            *tmp = underlying[i]->getCashFlow(step);
            if (fxProds[i].get()) {
                *tmp = *tmp * fxProds[i]->getValue(step, model->getDate(step));
            }
            *cashFlowSlice += weight * *tmp;
        }
        return *cashFlowSlice;
    }
    catch (exception& e) {
        throw ModelException(e, "KSumTree::getCashFlow, " + inst->outputName);
    }
}

/*********************************** KSum ***********************************/

void KSum::add(IProdCreatorSP el, double weight) {
    listK.push_back(el);
    weights.push_back(weight);
}

void KSum::addResetDates(const DateTimeArray &resetDatesP) {
    resetDates.insert(resetDates.end(), resetDatesP.begin(), resetDatesP.end());
}

void KSum::validatePop2Object(void) {
    try {
        if (listK.size() != weights.size())
            throw ModelException("Number of underlying elements in sum component (=" + 
                Format::toString(listK.size()) +
                ") does not equal number of supplied weights in sum component (=" + 
                Format::toString(weights.size()) + ")");
    }
    catch (exception& e) {
        throw ModelException(e, "KSum::validatePop2Object "+outputName);
    }
}

void KSum::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);
        validatePop2Object();

		// build IndexSpecFX IndexSpecs
		fx.resize(listK.size());
		string disc = discountYieldCurveName();

		for (int i=0; i<listK.size(); ++i) {
            if (Maths::isZero(weights[i])) 
                continue;
			KComponent* kComp = dynamic_cast<KComponent*>(listK[i].get());
            if (kComp) {
                string discComp = kComp->discountYieldCurveName();

                if (discComp.size() && disc.size()
                    && discComp != disc) {

				    fx[i].reset(new IndexSpecFX(discComp, disc));    
                    fx[i]->addResetDates(resetDates);
                    fx[i]->setup(model, market);
			    }
            }
            listK[i]->addResetDates(resetDates);
            listK[i]->setup(model, market);
        }
    }
    catch (exception& e) {
        throw ModelException(e, "KSum::setup "+outputName);
    }
}

double KSum::getValue(DateTime date, CashflowInfo &cfi) const {
    try {
        int i, s=listK.size();
        double val = constant;
        CashflowInfo cfiLoc;
        cfiLoc.updateAmountType(CashflowInfo::AmountType::KNOWN);
        if (constant!=0.) {
            cfiLoc.push("constant", constant);
        }
        for (i=0; i<s; ++i) {
            if (Maths::isZero(weights[i])) 
                continue;

            double fxValue = 1.;
            if (fx[i].get()) {
                fxValue = fx[i]->getValue(date, cfiLoc);
            }
            val += weights[i] * fxValue * listK[i]->getValue(date, cfiLoc);
            cfiLoc.push("weights["+Format::toString(i)+"]", weights[i]);
        }
        cfi.merge(cfiLoc, outputName);
        return val;
    }
    catch (exception& e) {
        throw ModelException(e, "KSum::getPastValue");
    }
}

KComponent::CashflowNodeSP KSum::reportCashflowsTree(bool amountsNeeded) const {
    CashflowNodeSP cfn(new CashflowNode);
    cfn->comp = KComponentConstSP(this);
    int i, s=listK.size();
    for (i=0; i<s; ++i) {
        if (!Maths::isZero(weights[i])) {
            KComponent *und = dynamic_cast<KComponent*>(listK[i].get());
            if (und) {
                CashflowNodeSP node = und->reportCashflowsTree(amountsNeeded);
                if (node.get()) {
                    node->setMultiplier(weights[i]);
                    cfn->underlyings.push_back(node);
                }
            }
        }
    }
    return cfn;
}

FDProductSP KSum::createProduct(FDModel * model) const {
    if (!setupCalled) 
        throw ModelException(getClass()->getName()+" \""+outputName+"\" not setup");
    return FDProductSP(new KSumTree(KSumConstSP(this), model));
}

void KSum::load(CClassSP& clazz) {

    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setDescription("Wrapper style Simple IR swap interface component");
    REGISTER(KSum, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(listK, "product list");
    FIELD_MAKE_OPTIONAL(listK);
    FIELD(weights, "weight for each product");
    FIELD_MAKE_OPTIONAL(weights);
    FIELD(constant, "constant to add"); 
    FIELD_MAKE_OPTIONAL(constant);
    FIELD(constSchedule, "Schedule of constants to add"); 
    FIELD_MAKE_OPTIONAL(constSchedule);
    FIELD(fx, ""); FIELD_MAKE_TRANSIENT(fx);
    FIELD(resetDates,""); FIELD_MAKE_TRANSIENT(resetDates);
    Addin::registerConstructor(Addin::UTILITIES, KSum::TYPE);
}

CClassConstSP const KSum::TYPE = CClass::registerClassLoadMethod(
    "KSum", typeid(KSum), KSum::load);

/******************************/
// for type linking
bool KSumLoad(void){
    return (KSum::TYPE != 0);
}

DRLIB_END_NAMESPACE

