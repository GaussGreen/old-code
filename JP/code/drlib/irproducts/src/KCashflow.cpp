//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KCashflow.cpp
//
//   Description : KCashflow component
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/KCashflow.hpp"

DRLIB_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
//  KCashflowTree
//  FDModel's pricing product/calc 
//  for the KCashflow component
//-----------------------------------------------------------------------------

class KCashflowTree : public FDProduct {

    /********************* methods ********************/
public:
    KCashflowTree(const KCashflowConstSP &inst, FDModel* model);

    virtual string getOutputName() const { return inst->outputName; }
    virtual void init(Control*) const;
    virtual void initProd(void);
    virtual void update(int &step, FDProduct::UpdateType updateType);
    virtual const TreeSlice& getValue(int step, DateTime eventDate) const;
    virtual void recordOutput(Control* ctrl, YieldCurveConstSP, Results* results);

    /******************** variables ********************/
private:
    KCashflowConstSP  inst;
    FDProductSP fxProd;
    TreeSliceSP mainSlice;
};


KCashflowTree::KCashflowTree(const KCashflowConstSP &inst, FDModel* model) 
: FDProduct(model), inst(inst) {
    try {
        // create implied FX product from IndexSpecs
        fxProd = model->createProduct(inst->fx);
        fxProd->addModelResetDates(inst->dates);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KCashflowTree::init(Control*) const{
    try {
        model->addCritDates(inst->dates);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KCashflowTree::initProd(void) {
    try {
        string discCurve = inst->getActualDiscount()->getName();

        mainSlice = model->createSlice(discCurve);
        mainSlice->name = inst->outputName + "_mainSlice";
        // setting mainSlice to zero and starting DEV means the mainSlice product
        // runs in "fast scalar DEV" mode, which means the DEV doesn't apply any
        // DEV calculation but updates the slice counter and is a nice way to 
        // initialize slice
        *mainSlice = 0.0;
        startDEV(mainSlice);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


const TreeSlice& KCashflowTree::getValue(int step, DateTime eventDate) const { 
    try {
        DateTime currentDate = model->getDate(step);

        if (eventDate != currentDate) {
            throw ModelException("Cannot be valued at a date != currentDate");
        }
        return *mainSlice;  // in domestic currency
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KCashflowTree::update(int &step, FDProduct::UpdateType /*updateType*/) {
    try {
        DateTime currDate = model->getDate(step);
        DateTime valueDate = model->getValueDate();
        
        if (currDate <= valueDate)
            return;

        // cashflows may pay on the same date
        // only cashflows that pay after value date are considered
        for (int i = 0; i < inst->dates.size(); i++) {
            if (currDate == inst->dates[i]) {
                const TreeSlice& fx = fxProd->getValue(step, currDate);
                *mainSlice += inst->amounts[i] * fx;
            }
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


void KCashflowTree::recordOutput(Control* ctrl, YieldCurveConstSP, Results* results) {
    try {
        const TreeSlice &price = getValue(0,model->getDate(0));
        const string ccy = inst->getActualDiscount()->getCcy();
        
        recordSliceToOutputName(
            ctrl, results, model, inst->isTopComponent(), inst->outputName, 
            ccy, price);

        inst->recordExtraOutput(ctrl, results);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }

}

//-----------------------------------------------------------------------------
//  KCashflow instrument
//  Multiple cashflow instrument
//-----------------------------------------------------------------------------

FDProductSP KCashflow::createProduct(FDModel * model) const {
    try {
        if (!setupCalled) 
            throw ModelException("steup() has not been called by parent component");
        return FDProductSP(new KCashflowTree(KCashflowConstSP(this), model));
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


// called via KComponent in the getMarket function
void KCashflow::setup(const IModel* model, const MarketData* market) {
    try {
        KComponent::setup(model, market);

        if (dates.size() != amounts.size())
            throw ModelException("\"dates\" and \"amounts\" should have same size");

        fx.reset(new IndexSpecFX(discount->getName(), getActualDiscount()->getName()));
        fx->addResetDates(dates);
        fx->setup(model, market);
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}


DateTime KCashflow::getLastDate(void) const {
    return (dates.empty() ? DateTime() : dates.back());
}


void KCashflow::reportCashFlows(CashflowInfoArray &cashflowInfos, bool amountsNeeded ) const {
    try {
        cashflowInfos.clear();

        // optional: gives a hint on how many cells will be needed in the array
        cashflowInfos.reserve(dates.size());

        // populate "cashFlows" and "cashflowInfos" arrays
        for (int i=0; i<dates.size(); ++i) {
            CashflowInfoSP cfi(new CashflowInfo);

            cfi->date = dates[i];
            cfi->amount = amounts[i];
            cfi->componentName = outputName;
            cfi->cfType = cfType;
            cfi->amountType = CashflowInfo::AmountType::KNOWN;

            cashflowInfos.push_back(cfi);
        }
    } 
    catch (exception& e) {
        throw makeException(e, __FUNCTION__);     
    }
}

void KCashflow::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(KCashflow, clazz);
    SUPERCLASS(KComponent);
    IMPLEMENTS(FDModel::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultConstructor);

    FIELD(dates,   "Array of cashflow payment dates");
    FIELD(amounts, "Cashflow payment amounts in the payment currency");
    FIELD(cfType,  "Type of cashflow generated to report in cashflow events - no effect on pricing");
    FIELD_MAKE_OPTIONAL(cfType);
    FIELD(fx,      "Internal indexSpecFX to convert from payment to pricing currency");
    FIELD_MAKE_TRANSIENT(fx);
 }


CClassConstSP const KCashflow::TYPE = CClass::registerClassLoadMethod(
    "KCashflow", typeid(KCashflow), KCashflow::load);

DEFINE_TEMPLATE_TYPE(KCashflowArray);

bool KCashflowLoad() { return KCashflow::TYPE != 0; }

DRLIB_END_NAMESPACE
